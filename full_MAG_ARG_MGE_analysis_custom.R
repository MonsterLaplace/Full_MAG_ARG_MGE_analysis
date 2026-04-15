###############################
## Full customized R workflow v4.2
## For representative MAG + ARG + MGE analysis
## Dataset: Yak vs White-lipped deer
###############################

rm(list = ls())
options(stringsAsFactors = FALSE)

###############################
## 0. Packages
###############################

pkg_needed <- c(
  "tidyverse", "data.table", "vegan", "ggrepel", "pheatmap",
  "ComplexHeatmap", "circlize", "randomForest", "pROC",
  "igraph", "ggraph", "ape"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(data.table)
library(vegan)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(randomForest)
library(pROC)
library(igraph)
library(ggraph)
library(ape)

if (requireNamespace("Maaslin2", quietly = TRUE)) library(Maaslin2)
if (requireNamespace("ggtree", quietly = TRUE)) library(ggtree)

###############################
## 1. Output directories
###############################

dir.create("R_out", showWarnings = FALSE)
dir.create("R_out/fig", showWarnings = FALSE, recursive = TRUE)
dir.create("R_out/tab", showWarnings = FALSE, recursive = TRUE)
dir.create("R_out/rds", showWarnings = FALSE, recursive = TRUE)

###############################
## 2. Helper functions
###############################

check_required_cols <- function(df, req, df_name = "dataframe"){
  miss <- setdiff(req, colnames(df))
  if(length(miss) > 0){
    stop(paste0("Missing columns in ", df_name, ": ", paste(miss, collapse = ", ")))
  }
}

clr_transform <- function(x, pseudocount = 1e-6){
  x <- as.matrix(x)
  x <- x + pseudocount
  gm <- exp(rowMeans(log(x)))
  log(x / gm)
}

scale01 <- function(x){
  if(all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)){
    return(rep(0, length(x)))
  } else {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
}

extract_taxon <- function(x, prefix = "g__"){
  out <- stringr::str_extract(x, paste0(prefix, "[^;]+"))
  out <- gsub(prefix, "", out)
  out[is.na(out)] <- "Unclassified"
  out
}

plot_box_jitter <- function(df, xvar, yvar, fillvar = NULL, outfile = NULL, width = 4, height = 4){
  fill_use <- ifelse(is.null(fillvar), xvar, fillvar)

  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[fill_use]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = .data[[fill_use]]), width = 0.15, alpha = 0.7, size = 2, show.legend = FALSE) +
    theme_bw() +
    labs(x = NULL, y = yvar)

  if(!is.null(outfile)) ggsave(outfile, p, width = width, height = height)
  return(p)
}

has_exactly_two_groups <- function(group_vec){
  g <- unique(as.character(group_vec[!is.na(group_vec)]))
  length(g) == 2
}

group_counts <- function(group_vec){
  table(as.character(group_vec), useNA = "no")
}

safe_wilcox_by_group <- function(df, value_col, group_col, out_txt = NULL, context = "", verbose = FALSE){
  tmp <- df[, c(value_col, group_col), drop = FALSE]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]

  if(nrow(tmp) == 0){
    msg <- paste0("Skipping ", context, ": no complete cases.")
    if(verbose) message(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  g <- as.character(tmp[[group_col]])
  v <- tmp[[value_col]]

  if(length(unique(g)) != 2){
    msg <- paste0("Skipping ", context, ": ", group_col, " does not contain exactly two groups.")
    if(verbose) message(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  if(length(unique(v)) <= 1 || stats::sd(v, na.rm = TRUE) == 0){
    msg <- paste0("Skipping ", context, ": ", value_col, " has zero variance.")
    if(verbose) message(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  wt <- tryCatch(
    wilcox.test(v ~ g),
    error = function(e) e
  )

  if(inherits(wt, "error")){
    msg <- paste0("Skipping ", context, ": wilcox.test failed for ", value_col, " ~ ", group_col, "; ", wt$message)
    if(verbose) message(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  if(!is.null(out_txt)) capture.output(wt, file = out_txt)
  return(wt)
}

safe_group_compare_table <- function(df, value_col, group_col){
  tmp <- df[, c(value_col, group_col), drop = FALSE]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]

  tmp %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n = n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

filter_nonzero_sd_cols <- function(df){
  keep <- sapply(df, function(x){
    if(!is.numeric(x)) return(FALSE)
    stats::sd(x, na.rm = TRUE) > 0
  })
  df[, keep, drop = FALSE]
}

run_rf_classifier <- function(feature_mat, meta, outcome = "host", prefix = "RF_model", ntree = 1000, seed = 123){
  set.seed(seed)

  dat <- data.frame(sample = rownames(feature_mat), feature_mat, check.names = FALSE) %>%
    left_join(meta, by = "sample")

  dat[[outcome]] <- factor(dat[[outcome]])

  if(length(unique(dat[[outcome]])) != 2){
    writeLines(paste0("Skipping ", prefix, ": outcome does not have exactly two classes."),
               paste0("R_out/tab/", prefix, "_skipped.txt"))
    return(NULL)
  }

  train_idx <- sample(seq_len(nrow(dat)), size = floor(0.7 * nrow(dat)))
  train_dat <- dat[train_idx, ]
  test_dat  <- dat[-train_idx, ]

  x_train <- train_dat %>% select(-sample, -host, -season, -sex, -age, -site)
  y_train <- train_dat[[outcome]]

  x_test <- test_dat %>% select(-sample, -host, -season, -sex, -age, -site)
  y_test <- test_dat[[outcome]]

  if(ncol(x_train) < 2){
    writeLines(paste0("Skipping ", prefix, ": fewer than 2 features after filtering."),
               paste0("R_out/tab/", prefix, "_skipped.txt"))
    return(NULL)
  }

  rf_fit <- randomForest(
    x = x_train,
    y = y_train,
    ntree = ntree,
    importance = TRUE
  )

  pred_prob <- predict(rf_fit, x_test, type = "prob")[, levels(y_train)[2]]
  pred_cls  <- predict(rf_fit, x_test, type = "response")

  roc_obj <- pROC::roc(
    response = y_test,
    predictor = pred_prob,
    levels = levels(y_train),
    direction = "<"
  )
  auc_val <- pROC::auc(roc_obj)

  imp <- importance(rf_fit, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    arrange(desc(MeanDecreaseAccuracy))

  write.csv(imp, paste0("R_out/tab/", prefix, "_importance.csv"), row.names = FALSE)

  pdf(paste0("R_out/fig/", prefix, "_ROC.pdf"), width = 5, height = 5)
  plot(roc_obj, main = paste0(prefix, " AUC = ", round(auc_val, 3)))
  dev.off()

  list(
    model = rf_fit,
    auc = auc_val,
    importance = imp,
    pred = data.frame(sample = test_dat$sample, truth = y_test, pred = pred_cls, prob = pred_prob)
  )
}

###############################
## 3. Read metadata
###############################

meta <- fread("00_meta/samples.tsv") %>% as.data.frame()
check_required_cols(meta, c("sample", "host", "site", "season", "sex", "age"), "samples.tsv")

meta <- meta %>%
  mutate(
    host = case_when(
      host == "BCL" ~ "deer",
      host == "Yak" ~ "yak",
      TRUE ~ host
    ),
    host = factor(host, levels = c("deer", "yak")),
    season = factor(season, levels = c("lowgrass", "highgrass")),
    site = factor(site),
    sex = factor(sex),
    age = factor(age)
  )

write.csv(meta, "R_out/tab/metadata_cleaned.csv", row.names = FALSE)

###############################
## 4. Read abundance tables
###############################

arg_mat <- fread("13_tables/sample_arg_abundance.tsv") %>% as.data.frame()
check_required_cols(arg_mat, c("sample_id"), "sample_arg_abundance.tsv")
colnames(arg_mat)[1] <- "sample"
rownames(arg_mat) <- arg_mat$sample
arg_mat <- arg_mat[, setdiff(colnames(arg_mat), "sample"), drop = FALSE]

mag_mat <- fread("13_tables/sample_mag_abundance.tsv") %>% as.data.frame()
check_required_cols(mag_mat, c("sample_id"), "sample_mag_abundance.tsv")
colnames(mag_mat)[1] <- "sample"
rownames(mag_mat) <- mag_mat$sample
mag_mat <- mag_mat[, setdiff(colnames(mag_mat), "sample"), drop = FALSE]

arg_mat <- arg_mat[meta$sample, , drop = FALSE]
mag_mat <- mag_mat[meta$sample, , drop = FALSE]

###############################
## 5. Read annotation tables
###############################

mag_arg <- fread("13_tables/mag_arg_catalog.tsv") %>% as.data.frame()
check_required_cols(mag_arg, c("MAG", "ARG", "Class"), "mag_arg_catalog.tsv")

gtdb <- fread("08_gtdbtk/classification_pplacer.tsv") %>% as.data.frame()
check_required_cols(gtdb, c("user_genome", "classification"), "classification_pplacer.tsv")
colnames(gtdb)[colnames(gtdb) == "user_genome"] <- "MAG"

if (!file.exists("13_tables/mag_trait_catalog.tsv")) {
  stop("mag_trait_catalog.tsv is required in this workflow.")
}
mag_trait <- fread("13_tables/mag_trait_catalog.tsv") %>% as.data.frame()
if (!"MAG" %in% colnames(mag_trait)) colnames(mag_trait)[1] <- "MAG"

###############################
## 6. Construct MAG summary table
###############################

mag_arg_summary <- mag_arg %>%
  group_by(MAG) %>%
  summarise(
    ARG_record_count_from_arg = n(),
    ARG_class_n_from_arg = n_distinct(Class),
    ARG_subclass_n_from_arg = ifelse("Subclass" %in% colnames(mag_arg), n_distinct(Subclass), NA),
    ARG_group_n_from_arg = ifelse("ARG_group" %in% colnames(mag_arg), n_distinct(ARG_group), NA),
    high_conf_ARG_count_from_arg = ifelse("High_confidence" %in% colnames(mag_arg), sum(High_confidence == 1, na.rm = TRUE), NA),
    validated_by_amrfinder_count_from_arg = ifelse("Validated_by_AMRFinder" %in% colnames(mag_arg), sum(Validated_by_AMRFinder == 1, na.rm = TRUE), NA),
    validated_by_rgi_count_from_arg = ifelse("Validated_by_RGI" %in% colnames(mag_arg), sum(Validated_by_RGI == 1, na.rm = TRUE), NA),
    validated_by_both_count_from_arg = ifelse("Validated_by_both" %in% colnames(mag_arg), sum(Validated_by_both == 1, na.rm = TRUE), NA),
    mobile_ARG_loose_count_from_arg = ifelse("putative_mobile_ARG_loose" %in% colnames(mag_arg), sum(putative_mobile_ARG_loose == 1, na.rm = TRUE), NA),
    mobile_ARG_strict_count_from_arg = ifelse("putative_mobile_ARG_strict" %in% colnames(mag_arg), sum(putative_mobile_ARG_strict == 1, na.rm = TRUE), NA),
    ARG_list_from_arg = paste(unique(ARG), collapse = ";"),
    Class_list_from_arg = paste(unique(Class), collapse = ";"),
    .groups = "drop"
  )

mag_info <- mag_trait %>%
  left_join(mag_arg_summary, by = "MAG") %>%
  left_join(gtdb[, c("MAG", "classification")], by = "MAG")

num_cols <- colnames(mag_info)[sapply(mag_info, is.numeric)]
for (cc in num_cols) {
  mag_info[[cc]][is.na(mag_info[[cc]])] <- 0
}

if (!"has_ARG" %in% colnames(mag_info)) {
  mag_info$has_ARG <- ifelse(mag_info$ARG_count > 0, 1, 0)
}
if (!"has_mobile_ARG_loose" %in% colnames(mag_info)) {
  mag_info$has_mobile_ARG_loose <- ifelse(mag_info$mobile_ARG_loose_count > 0, 1, 0)
}
if (!"has_mobile_ARG_strict" %in% colnames(mag_info)) {
  mag_info$has_mobile_ARG_strict <- ifelse(mag_info$mobile_ARG_strict_count > 0, 1, 0)
}
if (!"has_vfdb_hit" %in% colnames(mag_info)) {
  mag_info$has_vfdb_hit <- ifelse(mag_info$vfdb_hit_count > 0, 1, 0)
}

mag_info <- mag_info %>%
  mutate(
    is_ARG_host = has_ARG == 1 | ARG_count > 0,
    phylum = extract_taxon(classification, "p__"),
    genus  = extract_taxon(classification, "g__")
  )

write.csv(mag_info, "R_out/tab/mag_info_merged.csv", row.names = FALSE)

###############################
## 7. Transform abundance matrices
###############################

arg_rel <- arg_mat / rowSums(arg_mat)
arg_rel[is.na(arg_rel)] <- 0

mag_rel <- mag_mat / rowSums(mag_mat)
mag_rel[is.na(mag_rel)] <- 0

arg_clr <- clr_transform(arg_mat)
mag_clr <- clr_transform(mag_mat)

###############################
## 8. Resistome alpha diversity
###############################

resistome_alpha <- data.frame(
  sample = rownames(arg_mat),
  ARG_burden = rowSums(arg_mat),
  ARG_richness = rowSums(arg_mat > 0),
  ARG_shannon = diversity(arg_mat, index = "shannon")
) %>%
  left_join(meta, by = "sample")

write.csv(resistome_alpha, "R_out/tab/resistome_alpha.csv", row.names = FALSE)

plot_box_jitter(resistome_alpha, "host", "ARG_burden", "host", "R_out/fig/resistome_ARG_burden_by_host.pdf")
plot_box_jitter(resistome_alpha, "host", "ARG_richness", "host", "R_out/fig/resistome_ARG_richness_by_host.pdf")
plot_box_jitter(resistome_alpha, "host", "ARG_shannon", "host", "R_out/fig/resistome_ARG_shannon_by_host.pdf")
plot_box_jitter(resistome_alpha, "season", "ARG_burden", "season", "R_out/fig/resistome_ARG_burden_by_season.pdf")

###############################
## 9. Resistome beta diversity
###############################

arg_bray <- vegdist(arg_rel, method = "bray")
arg_pcoa <- cmdscale(arg_bray, k = 2, eig = TRUE)

arg_pcoa_df <- data.frame(
  sample = rownames(arg_rel),
  PC1 = arg_pcoa$points[,1],
  PC2 = arg_pcoa$points[,2]
) %>% left_join(meta, by = "sample")

p_arg_pcoa <- ggplot(arg_pcoa_df, aes(PC1, PC2, color = host, shape = season)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    data = arg_pcoa_df,
    aes(PC1, PC2, color = host, group = host),
    inherit.aes = FALSE,
    level = 0.95
  ) +
  theme_bw()

ggsave("R_out/fig/resistome_PCoA_bray.pdf", p_arg_pcoa, width = 6, height = 5)

arg_euc <- dist(arg_clr)
perm_arg_bray <- adonis2(arg_bray ~ host * season, data = meta, permutations = 999)
perm_arg_euc  <- adonis2(arg_euc ~ host * season, data = meta, permutations = 999)

capture.output(perm_arg_bray, file = "R_out/tab/resistome_PERMANOVA_bray.txt")
capture.output(perm_arg_euc, file = "R_out/tab/resistome_PERMANOVA_aitchison.txt")

###############################
## 10. Resistome heatmap
###############################

top_arg <- names(sort(colSums(arg_mat), decreasing = TRUE))[1:min(30, ncol(arg_mat))]
arg_z <- scale(arg_rel[, top_arg, drop = FALSE])
arg_z <- as.matrix(arg_z)

meta_hm <- meta[match(rownames(arg_z), meta$sample), ]

if(nrow(arg_z) != nrow(meta_hm)){
  stop("Heatmap row number does not match metadata row number.")
}

ha_row <- rowAnnotation(
  host = meta_hm$host,
  season = meta_hm$season,
  col = list(
    host = c(deer = "#1f78b4", yak = "#e31a1c"),
    season = c(lowgrass = "#33a02c", highgrass = "#ff7f00")
  )
)

pdf("R_out/fig/resistome_heatmap_top30.pdf", width = 10, height = 8)
Heatmap(
  arg_z,
  name = "Z-score",
  left_annotation = ha_row,
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_gp = grid::gpar(fontsize = 8)
)
dev.off()

###############################
## 11. Differential ARG (MaAsLin2)
###############################

if ("Maaslin2" %in% loadedNamespaces()) {
  arg_prev <- colMeans(arg_mat > 0)
  arg_filt <- arg_mat[, arg_prev >= 0.1, drop = FALSE]

  meta_maas <- meta
  rownames(meta_maas) <- meta_maas$sample

  Maaslin2(
    input_data = as.data.frame(arg_filt),
    input_metadata = meta_maas,
    output = "R_out/maaslin2_arg",
    fixed_effects = c("host", "season"),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM"
  )

  if (file.exists("R_out/maaslin2_arg/significant_results.tsv")) {
    maas_arg <- fread("R_out/maaslin2_arg/significant_results.tsv") %>% as.data.frame()
    write.csv(maas_arg, "R_out/tab/maaslin2_ARG_significant.csv", row.names = FALSE)
  }
}

###############################
## 12. MAG community analysis
###############################

mag_bray <- vegdist(mag_rel, method = "bray")
mag_pcoa <- cmdscale(mag_bray, k = 2, eig = TRUE)

mag_pcoa_df <- data.frame(
  sample = rownames(mag_rel),
  PC1 = mag_pcoa$points[,1],
  PC2 = mag_pcoa$points[,2]
) %>% left_join(meta, by = "sample")

p_mag_pcoa <- ggplot(mag_pcoa_df, aes(PC1, PC2, color = host, shape = season)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    data = mag_pcoa_df,
    aes(PC1, PC2, color = host, group = host),
    inherit.aes = FALSE,
    level = 0.95
  ) +
  theme_bw()

ggsave("R_out/fig/MAG_PCoA_bray.pdf", p_mag_pcoa, width = 6, height = 5)

perm_mag_bray <- adonis2(mag_bray ~ host * season, data = meta, permutations = 999)
capture.output(perm_mag_bray, file = "R_out/tab/MAG_PERMANOVA_bray.txt")

###############################
## 13. Differential MAG (MaAsLin2)
###############################

maas_mag_host <- NULL

if ("Maaslin2" %in% loadedNamespaces()) {
  mag_prev <- colMeans(mag_mat > 0)
  mag_filt <- mag_mat[, mag_prev >= 0.1, drop = FALSE]

  meta_maas <- meta
  rownames(meta_maas) <- meta_maas$sample

  Maaslin2(
    input_data = as.data.frame(mag_filt),
    input_metadata = meta_maas,
    output = "R_out/maaslin2_mag",
    fixed_effects = c("host", "season"),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM"
  )

  if (file.exists("R_out/maaslin2_mag/significant_results.tsv")) {
    maas_mag <- fread("R_out/maaslin2_mag/significant_results.tsv") %>% as.data.frame()
    maas_mag_host <- maas_mag %>% filter(grepl("host", metadata))
    write.csv(maas_mag_host, "R_out/tab/host_enriched_MAGs.csv", row.names = FALSE)
  }
}

###############################
## 14. ARG-host MAG abundance per sample
###############################

arg_host_mags <- intersect(unique(mag_arg$MAG), colnames(mag_mat))

arghost_burden <- data.frame(
  sample = rownames(mag_mat),
  ARGhost_MAG_abundance = if(length(arg_host_mags) > 0) rowSums(mag_mat[, arg_host_mags, drop = FALSE]) else 0
) %>%
  left_join(meta, by = "sample")

write.csv(arghost_burden, "R_out/tab/ARGhost_MAG_abundance_by_sample.csv", row.names = FALSE)
plot_box_jitter(arghost_burden, "host", "ARGhost_MAG_abundance", "host", "R_out/fig/ARGhost_MAG_abundance_by_host.pdf")

###############################
## 15. Shared vs host-specific MAG
###############################

mag_presence <- (mag_mat > 0) * 1
yak_samples <- meta$sample[meta$host == "yak"]
deer_samples <- meta$sample[meta$host == "deer"]

mag_shared <- data.frame(
  MAG = colnames(mag_mat),
  yak_prev = colMeans(mag_presence[yak_samples, , drop = FALSE]),
  deer_prev = colMeans(mag_presence[deer_samples, , drop = FALSE])
) %>%
  mutate(
    occurrence_type = case_when(
      yak_prev > 0 & deer_prev > 0 ~ "shared",
      yak_prev > 0 & deer_prev == 0 ~ "yak_only",
      yak_prev == 0 & deer_prev > 0 ~ "deer_only",
      TRUE ~ "rare"
    )
  )

mag_info2 <- mag_info %>%
  left_join(mag_shared, by = "MAG")

###############################
## 16. Add host enrichment labels
###############################

if (!is.null(maas_mag_host) && nrow(maas_mag_host) > 0) {
  mag_host_effect <- maas_mag_host %>%
    select(feature, coef, qval) %>%
    rename(MAG = feature, host_coef = coef, host_q = qval)

  mag_info2 <- mag_info2 %>%
    left_join(mag_host_effect, by = "MAG") %>%
    mutate(
      host_enrichment = case_when(
        !is.na(host_q) & host_q < 0.05 & host_coef > 0 ~ "yak_enriched",
        !is.na(host_q) & host_q < 0.05 & host_coef < 0 ~ "deer_enriched",
        TRUE ~ "non_sig"
      )
    )
} else {
  mag_info2 <- mag_info2 %>% mutate(host_enrichment = "non_sig")
}

###############################
## 17. Group overview tables
###############################

mag_info2$ARG_host_group <- ifelse(mag_info2$is_ARG_host, "ARG_host_MAG", "non_ARG_host_MAG")

arg_host_group_tab <- data.frame(
  group = names(table(mag_info2$ARG_host_group)),
  n = as.integer(table(mag_info2$ARG_host_group))
)
write.csv(arg_host_group_tab, "R_out/tab/ARG_host_group_overview.csv", row.names = FALSE)

host_enrich_tab <- data.frame(
  group = names(table(mag_info2$host_enrichment)),
  n = as.integer(table(mag_info2$host_enrichment))
)
write.csv(host_enrich_tab, "R_out/tab/host_enrichment_overview.csv", row.names = FALSE)

###############################
## 18. Risk score
###############################

for (cc in c("ARG_count", "mobile_ARG_loose_count", "mobile_ARG_strict_count",
             "plasmid_like_contig_count", "virus_like_contig_count",
             "integron_contig_count", "integrase_contig_count",
             "IS_contig_count", "is_element_total", "transposase_total",
             "vfdb_hit_count", "vfdb_gene_count", "vfdb_vfg_count")) {
  if (!cc %in% colnames(mag_info2)) mag_info2[[cc]] <- 0
  mag_info2[[cc]][is.na(mag_info2[[cc]])] <- 0
}

mag_info2 <- mag_info2 %>%
  mutate(
    shared_score = ifelse(occurrence_type == "shared", 1, 0),
    yak_score = ifelse(host_enrichment == "yak_enriched", 1, 0),
    risk_score =
      0.18 * scale01(ARG_count) +
      0.18 * scale01(mobile_ARG_strict_count) +
      0.08 * scale01(mobile_ARG_loose_count) +
      0.10 * scale01(plasmid_like_contig_count) +
      0.08 * scale01(integrase_contig_count) +
      0.08 * scale01(IS_contig_count) +
      0.08 * scale01(transposase_total) +
      0.10 * scale01(vfdb_vfg_count) +
      0.05 * shared_score +
      0.07 * yak_score
  )

write.csv(mag_info2, "R_out/tab/mag_info_with_risk_score.csv", row.names = FALSE)

highrisk_mag <- mag_info2 %>%
  filter(is_ARG_host) %>%
  arrange(desc(risk_score))

write.csv(highrisk_mag, "R_out/tab/high_risk_ARG_host_MAGs.csv", row.names = FALSE)

top20 <- highrisk_mag %>% slice(1:min(20, n()))
if(nrow(top20) > 0){
  p_toprisk <- ggplot(top20, aes(reorder(MAG, risk_score), risk_score, fill = host_enrichment)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = NULL, y = "Risk score")
  ggsave("R_out/fig/top20_highrisk_MAGs.pdf", p_toprisk, width = 7, height = 6)
}

###############################
## 19. ARG host vs non-ARG host MGE comparison
###############################

mge_test_cols <- intersect(c("mobile_ARG_loose_count", "mobile_ARG_strict_count",
                             "plasmid_like_contig_count", "integrase_contig_count",
                             "IS_contig_count", "transposase_total"),
                           colnames(mag_info2))

write.csv(
  data.frame(group = names(group_counts(mag_info2$ARG_host_group)),
             n = as.integer(group_counts(mag_info2$ARG_host_group))),
  "R_out/tab/ARG_host_group_counts.csv",
  row.names = FALSE
)

for (cc in mge_test_cols) {

  out_file <- paste0("R_out/tab/wilcox_", cc, "_ARGhost_vs_nonARGhost.txt")

  wt <- safe_wilcox_by_group(
    df = mag_info2,
    value_col = cc,
    group_col = "ARG_host_group",
    out_txt = out_file,
    context = cc,
    verbose = FALSE
  )

  if(!is.null(wt)){
    p_tmp <- ggplot(mag_info2, aes(x = ARG_host_group, y = .data[[cc]], fill = ARG_host_group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      theme_bw() +
      labs(x = "ARG host group", y = cc)

    ggsave(paste0("R_out/fig/", cc, "_ARGhost_compare.pdf"), p_tmp, width = 4, height = 4)

    summary_tab <- safe_group_compare_table(mag_info2, cc, "ARG_host_group")
    write.csv(summary_tab,
              paste0("R_out/tab/", cc, "_ARGhost_compare_summary.csv"),
              row.names = FALSE)
  }
}

###############################
## 19A. CAZy functional coupling analysis
###############################

cazy_total_cols <- intersect(c("CAZy_count", "CAZyme_count", "cazy_count"), colnames(mag_info2))
cazy_family_cols <- intersect(c("CAZy_family_count"), colnames(mag_info2))
cazy_sub_cols <- intersect(c("GH_count", "GT_count", "CE_count", "CBM_count", "PL_count", "AA_count"),
                           colnames(mag_info2))
cazy_sub_family_cols <- intersect(c("GH_family_count", "GT_family_count", "CE_family_count",
                                    "CBM_family_count", "PL_family_count", "AA_family_count"),
                                  colnames(mag_info2))

if(length(cazy_total_cols) > 0){

  cazy_total_col <- cazy_total_cols[1]

  ## 9.1 ARG-host vs non-ARG-host
  cazy_compare1 <- mag_info2 %>%
    mutate(ARG_host_group = ifelse(is_ARG_host, "ARG_host_MAG", "non_ARG_host_MAG"))

  write.csv(
    data.frame(group = names(group_counts(cazy_compare1$ARG_host_group)),
               n = as.integer(group_counts(cazy_compare1$ARG_host_group))),
    "R_out/tab/CAZy_ARGhost_group_counts.csv",
    row.names = FALSE
  )

  if(has_exactly_two_groups(cazy_compare1$ARG_host_group)){

    cazy_compare1_summary <- safe_group_compare_table(cazy_compare1, cazy_total_col, "ARG_host_group")
    write.csv(cazy_compare1_summary,
              "R_out/tab/CAZy_ARGhost_vs_nonARGhost_summary.csv",
              row.names = FALSE)

    safe_wilcox_by_group(
      df = cazy_compare1,
      value_col = cazy_total_col,
      group_col = "ARG_host_group",
      out_txt = "R_out/tab/wilcox_CAZy_ARGhost_vs_nonARGhost.txt",
      context = "ARG-host vs non-ARG-host CAZy comparison",
      verbose = FALSE
    )

    p_cazy1 <- ggplot(cazy_compare1,
                      aes(x = ARG_host_group, y = .data[[cazy_total_col]], fill = ARG_host_group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      theme_bw() +
      labs(x = NULL, y = cazy_total_col)

    ggsave("R_out/fig/CAZy_ARGhost_vs_nonARGhost.pdf", p_cazy1, width = 5, height = 4)

    all_cazy_compare_cols1 <- c(cazy_family_cols, cazy_sub_cols, cazy_sub_family_cols)

    if(length(all_cazy_compare_cols1) > 0){
      cazy_sub_res1 <- lapply(all_cazy_compare_cols1, function(cc){

        tmp <- cazy_compare1[, c(cc, "ARG_host_group"), drop = FALSE]
        tmp <- tmp[complete.cases(tmp), , drop = FALSE]

        if(nrow(tmp) == 0 || length(unique(tmp$ARG_host_group)) != 2 || stats::sd(tmp[[cc]], na.rm = TRUE) == 0){
          return(data.frame(
            feature = cc,
            p_value = NA,
            mean_ARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
            mean_nonARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE),
            median_ARGhost = median(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
            median_nonARGhost = median(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE)
          ))
        }

        wt <- tryCatch(wilcox.test(tmp[[cc]] ~ tmp$ARG_host_group), error = function(e) NULL)

        data.frame(
          feature = cc,
          p_value = ifelse(is.null(wt), NA, wt$p.value),
          mean_ARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
          mean_nonARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE),
          median_ARGhost = median(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
          median_nonARGhost = median(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE)
        )
      }) %>% bind_rows()

      cazy_sub_res1$q_value <- p.adjust(cazy_sub_res1$p_value, method = "BH")
      write.csv(cazy_sub_res1,
                "R_out/tab/CAZy_subclass_ARGhost_vs_nonARGhost.csv",
                row.names = FALSE)

      cazy_sub_plot1 <- cazy_sub_res1 %>%
        pivot_longer(cols = c(mean_ARGhost, mean_nonARGhost),
                     names_to = "group", values_to = "mean_value")

      p_cazy_sub1 <- ggplot(cazy_sub_plot1,
                            aes(x = feature, y = mean_value, fill = group)) +
        geom_col(position = "dodge") +
        theme_bw() +
        labs(x = NULL, y = "Mean count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave("R_out/fig/CAZy_subclass_ARGhost_vs_nonARGhost.pdf",
             p_cazy_sub1, width = 8, height = 4)
    }

  } else {
    msg <- "Skipping ARG-host vs non-ARG-host CAZy comparison: only one ARG_host_group present."
    message(msg)
    writeLines(msg, "R_out/tab/CAZy_ARGhost_vs_nonARGhost_skipped.txt")
  }

  ## 9.2 host-enriched ARG-host MAG
  cazy_compare2 <- mag_info2 %>%
    filter(is_ARG_host, host_enrichment %in% c("yak_enriched", "deer_enriched"))

  write.csv(
    data.frame(group = names(group_counts(cazy_compare2$host_enrichment)),
               n = as.integer(group_counts(cazy_compare2$host_enrichment))),
    "R_out/tab/CAZy_host_enriched_ARGhost_group_counts.csv",
    row.names = FALSE
  )

  if(has_exactly_two_groups(cazy_compare2$host_enrichment)){

    cazy_compare2_summary <- safe_group_compare_table(cazy_compare2, cazy_total_col, "host_enrichment")
    write.csv(cazy_compare2_summary,
              "R_out/tab/CAZy_host_enriched_ARGhost_summary.csv",
              row.names = FALSE)

    safe_wilcox_by_group(
      df = cazy_compare2,
      value_col = cazy_total_col,
      group_col = "host_enrichment",
      out_txt = "R_out/tab/wilcox_CAZy_host_enriched_ARGhost.txt",
      context = "host-enriched ARG-host CAZy comparison",
      verbose = FALSE
    )

    p_cazy2 <- ggplot(cazy_compare2,
                      aes(x = host_enrichment, y = .data[[cazy_total_col]], fill = host_enrichment)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      theme_bw() +
      labs(x = NULL, y = cazy_total_col)

    ggsave("R_out/fig/CAZy_host_enriched_ARGhost_compare.pdf", p_cazy2, width = 5, height = 4)

    all_cazy_compare_cols2 <- c(cazy_family_cols, cazy_sub_cols, cazy_sub_family_cols)

    if(length(all_cazy_compare_cols2) > 0){
      cazy_sub_res2 <- lapply(all_cazy_compare_cols2, function(cc){

        tmp <- cazy_compare2[, c(cc, "host_enrichment"), drop = FALSE]
        tmp <- tmp[complete.cases(tmp), , drop = FALSE]

        if(nrow(tmp) == 0 || length(unique(tmp$host_enrichment)) != 2 || stats::sd(tmp[[cc]], na.rm = TRUE) == 0){
          return(data.frame(
            feature = cc,
            p_value = NA,
            mean_yak = mean(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
            mean_deer = mean(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE),
            median_yak = median(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
            median_deer = median(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE)
          ))
        }

        wt <- tryCatch(wilcox.test(tmp[[cc]] ~ tmp$host_enrichment), error = function(e) NULL)

        data.frame(
          feature = cc,
          p_value = ifelse(is.null(wt), NA, wt$p.value),
          mean_yak = mean(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
          mean_deer = mean(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE),
          median_yak = median(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
          median_deer = median(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE)
        )
      }) %>% bind_rows()

      cazy_sub_res2$q_value <- p.adjust(cazy_sub_res2$p_value, method = "BH")
      write.csv(cazy_sub_res2,
                "R_out/tab/CAZy_subclass_host_enriched_ARGhost.csv",
                row.names = FALSE)

      cazy_sub_plot2 <- cazy_sub_res2 %>%
        pivot_longer(cols = c(mean_yak, mean_deer),
                     names_to = "group", values_to = "mean_value")

      p_cazy_sub2 <- ggplot(cazy_sub_plot2,
                            aes(x = feature, y = mean_value, fill = group)) +
        geom_col(position = "dodge") +
        theme_bw() +
        labs(x = NULL, y = "Mean count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave("R_out/fig/CAZy_subclass_host_enriched_ARGhost.pdf",
             p_cazy_sub2, width = 8, height = 4)
    }

  } else {
    msg <- "Skipping host-enriched ARG-host CAZy comparison: not enough yak/deer enriched ARG-host MAGs."
    message(msg)
    writeLines(msg, "R_out/tab/CAZy_host_enriched_ARGhost_skipped.txt")
  }

} else {
  msg <- "No CAZy total column detected in mag_info2. Skip CAZy functional coupling analysis."
  message(msg)
  writeLines(msg, "R_out/tab/CAZy_analysis_skipped.txt")
}

###############################
## 19B. VFDB coupling analysis
###############################

vf_cols_main <- intersect(c("vfdb_hit_count", "vfdb_gene_count", "vfdb_vfg_count"), colnames(mag_info2))

if(length(vf_cols_main) > 0){

  vf_compare1 <- mag_info2 %>%
    mutate(ARG_host_group = ifelse(is_ARG_host, "ARG_host_MAG", "non_ARG_host_MAG"))

  write.csv(
    data.frame(group = names(group_counts(vf_compare1$ARG_host_group)),
               n = as.integer(group_counts(vf_compare1$ARG_host_group))),
    "R_out/tab/VFDB_ARGhost_group_counts.csv",
    row.names = FALSE
  )

  if(has_exactly_two_groups(vf_compare1$ARG_host_group)){
    vf_res1 <- lapply(vf_cols_main, function(cc){
      tmp <- vf_compare1[, c(cc, "ARG_host_group"), drop = FALSE]
      tmp <- tmp[complete.cases(tmp), , drop = FALSE]

      if(nrow(tmp) == 0 || length(unique(tmp$ARG_host_group)) != 2 || stats::sd(tmp[[cc]], na.rm = TRUE) == 0){
        return(data.frame(
          feature = cc,
          p_value = NA,
          mean_ARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
          mean_nonARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE),
          median_ARGhost = median(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
          median_nonARGhost = median(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE)
        ))
      }

      wt <- tryCatch(wilcox.test(tmp[[cc]] ~ tmp$ARG_host_group), error = function(e) NULL)

      data.frame(
        feature = cc,
        p_value = ifelse(is.null(wt), NA, wt$p.value),
        mean_ARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
        mean_nonARGhost = mean(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE),
        median_ARGhost = median(tmp[[cc]][tmp$ARG_host_group == "ARG_host_MAG"], na.rm = TRUE),
        median_nonARGhost = median(tmp[[cc]][tmp$ARG_host_group == "non_ARG_host_MAG"], na.rm = TRUE)
      )
    }) %>% bind_rows()

    vf_res1$q_value <- p.adjust(vf_res1$p_value, method = "BH")
    write.csv(vf_res1, "R_out/tab/VFDB_ARGhost_vs_nonARGhost.csv", row.names = FALSE)
  } else {
    msg <- "Skipping VFDB ARG-host vs non-ARG-host comparison: only one ARG_host_group present."
    message(msg)
    writeLines(msg, "R_out/tab/VFDB_ARGhost_vs_nonARGhost_skipped.txt")
  }

  vf_compare2 <- mag_info2 %>%
    filter(is_ARG_host, host_enrichment %in% c("yak_enriched", "deer_enriched"))

  write.csv(
    data.frame(group = names(group_counts(vf_compare2$host_enrichment)),
               n = as.integer(group_counts(vf_compare2$host_enrichment))),
    "R_out/tab/VFDB_host_enriched_ARGhost_group_counts.csv",
    row.names = FALSE
  )

  if(has_exactly_two_groups(vf_compare2$host_enrichment)){
    vf_res2 <- lapply(vf_cols_main, function(cc){
      tmp <- vf_compare2[, c(cc, "host_enrichment"), drop = FALSE]
      tmp <- tmp[complete.cases(tmp), , drop = FALSE]

      if(nrow(tmp) == 0 || length(unique(tmp$host_enrichment)) != 2 || stats::sd(tmp[[cc]], na.rm = TRUE) == 0){
        return(data.frame(
          feature = cc,
          p_value = NA,
          mean_yak = mean(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
          mean_deer = mean(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE),
          median_yak = median(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
          median_deer = median(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE)
        ))
      }

      wt <- tryCatch(wilcox.test(tmp[[cc]] ~ tmp$host_enrichment), error = function(e) NULL)

      data.frame(
        feature = cc,
        p_value = ifelse(is.null(wt), NA, wt$p.value),
        mean_yak = mean(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
        mean_deer = mean(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE),
        median_yak = median(tmp[[cc]][tmp$host_enrichment == "yak_enriched"], na.rm = TRUE),
        median_deer = median(tmp[[cc]][tmp$host_enrichment == "deer_enriched"], na.rm = TRUE)
      )
    }) %>% bind_rows()

    vf_res2$q_value <- p.adjust(vf_res2$p_value, method = "BH")
    write.csv(vf_res2, "R_out/tab/VFDB_host_enriched_ARGhost.csv", row.names = FALSE)
  } else {
    msg <- "Skipping VFDB host-enriched ARG-host comparison: not enough yak/deer enriched ARG-host MAGs."
    message(msg)
    writeLines(msg, "R_out/tab/VFDB_host_enriched_ARGhost_skipped.txt")
  }
}

###############################
## 20. Random forest
###############################

arg_prev <- colMeans(arg_mat > 0)
arg_rf <- arg_mat[, arg_prev >= 0.1, drop = FALSE]

mag_prev <- colMeans(mag_mat > 0)
mag_rf <- mag_mat[, mag_prev >= 0.1, drop = FALSE]

if(length(arg_host_mags) > 0){
  arghost_rf <- mag_mat[, arg_host_mags, drop = FALSE]
  arghost_prev <- colMeans(arghost_rf > 0)
  arghost_rf <- arghost_rf[, arghost_prev >= 0.1, drop = FALSE]
}

rf_results <- list()
tmp_rf <- run_rf_classifier(arg_rf, meta, prefix = "RF_ARG")
if(!is.null(tmp_rf)) rf_results$ARG <- tmp_rf

tmp_rf <- run_rf_classifier(mag_rf, meta, prefix = "RF_MAG")
if(!is.null(tmp_rf)) rf_results$MAG <- tmp_rf

if(exists("arghost_rf") && ncol(arghost_rf) > 2){
  tmp_rf <- run_rf_classifier(arghost_rf, meta, prefix = "RF_ARGHOSTMAG")
  if(!is.null(tmp_rf)) rf_results$ARGHOST <- tmp_rf
}

rf_auc <- data.frame(
  model = names(rf_results),
  AUC = sapply(rf_results, function(x) as.numeric(x$auc))
)

if(nrow(rf_auc) > 0){
  write.csv(rf_auc, "R_out/tab/RF_model_AUC_comparison.csv", row.names = FALSE)

  p_auc <- ggplot(rf_auc, aes(model, AUC, fill = model)) +
    geom_col(width = 0.7) +
    ylim(0, 1) +
    theme_bw() +
    labs(x = NULL, y = "AUC")

  ggsave("R_out/fig/RF_AUC_comparison.pdf", p_auc, width = 5, height = 4)
}

if("ARGHOST" %in% names(rf_results)){
  top_imp <- rf_results$ARGHOST$importance %>% slice(1:min(20, n()))
  p_imp <- ggplot(top_imp, aes(reorder(feature, MeanDecreaseAccuracy), MeanDecreaseAccuracy)) +
    geom_col(fill = "#d95f02") +
    coord_flip() +
    theme_bw()
  ggsave("R_out/fig/RF_ARGHOSTMAG_top20_importance.pdf", p_imp, width = 6, height = 6)
}

###############################
## 21. Linear models
###############################

sample_metrics <- resistome_alpha %>%
  left_join(arghost_burden %>% select(sample, ARGhost_MAG_abundance), by = "sample")

m1 <- lm(log1p(ARG_burden) ~ host * season, data = sample_metrics)
m2 <- lm(log1p(ARG_richness) ~ host * season, data = sample_metrics)
m3 <- lm(log1p(ARG_shannon) ~ host * season, data = sample_metrics)
m4 <- lm(log1p(ARGhost_MAG_abundance) ~ host * season, data = sample_metrics)

capture.output(summary(m1), file = "R_out/tab/lm_ARG_burden_host_season.txt")
capture.output(summary(m2), file = "R_out/tab/lm_ARG_richness_host_season.txt")
capture.output(summary(m3), file = "R_out/tab/lm_ARG_shannon_host_season.txt")
capture.output(summary(m4), file = "R_out/tab/lm_ARGhost_abundance_host_season.txt")

###############################
## 22. ARG-MAG bipartite network
###############################

mag_arg_edges <- mag_arg %>%
  select(MAG, ARG, Class) %>%
  distinct()

top_mag <- highrisk_mag %>% slice(1:min(50, n())) %>% pull(MAG)
top_arg2 <- mag_arg %>% count(ARG, sort = TRUE) %>% slice(1:min(30, n())) %>% pull(ARG)

net_edges <- mag_arg_edges %>%
  filter(MAG %in% top_mag, ARG %in% top_arg2)

if (nrow(net_edges) > 0) {
  mag_nodes <- data.frame(
    name = unique(net_edges$MAG),
    node_type = "MAG"
  ) %>%
    left_join(mag_info2 %>% select(MAG, host_enrichment, risk_score), by = c("name" = "MAG"))

  arg_nodes <- data.frame(
    name = unique(net_edges$ARG),
    node_type = "ARG"
  )

  nodes <- bind_rows(mag_nodes, arg_nodes)

  g_bi <- graph_from_data_frame(
    d = net_edges %>% select(from = MAG, to = ARG),
    vertices = nodes,
    directed = FALSE
  )

  p_net <- ggraph(g_bi, layout = "fr") +
    geom_edge_link(alpha = 0.3, colour = "grey50") +
    geom_node_point(aes(color = node_type,
                        size = ifelse(node_type == "MAG", risk_score + 1, 2))) +
    geom_node_text(aes(label = ifelse(node_type == "ARG", name, "")),
                   repel = TRUE, size = 3) +
    scale_color_manual(values = c("MAG" = "#1f78b4", "ARG" = "#e31a1c")) +
    theme_void()

  ggsave("R_out/fig/network_ARG_MAG_bipartite.pdf", p_net, width = 10, height = 8)
}

###############################
## 23. ARG-MGE-MAG tripartite network
###############################

if (file.exists("13_tables/arg_mge_mag_edges.tsv")) {
  arg_mge_mag <- fread("13_tables/arg_mge_mag_edges.tsv") %>% as.data.frame()

  if (all(c("MAG", "ARG", "MGE_type") %in% colnames(arg_mge_mag))) {
    mag_nodes3 <- data.frame(name = unique(arg_mge_mag$MAG), node_type = "MAG") %>%
      left_join(mag_info2 %>% select(MAG, host_enrichment, risk_score), by = c("name" = "MAG"))

    arg_nodes3 <- data.frame(name = unique(arg_mge_mag$ARG), node_type = "ARG")
    mge_nodes3 <- data.frame(name = unique(arg_mge_mag$MGE_type), node_type = "MGE")

    nodes3 <- bind_rows(mag_nodes3, arg_nodes3, mge_nodes3)

    edges_mag_arg <- arg_mge_mag %>% select(from = MAG, to = ARG) %>% distinct()
    edges_arg_mge <- arg_mge_mag %>% select(from = ARG, to = MGE_type) %>% distinct()
    edges3 <- bind_rows(edges_mag_arg, edges_arg_mge)

    g_tri <- graph_from_data_frame(edges3, vertices = nodes3, directed = FALSE)

    p_tri <- ggraph(g_tri, layout = "fr") +
      geom_edge_link(alpha = 0.25, colour = "grey60") +
      geom_node_point(aes(color = node_type,
                          size = case_when(
                            node_type == "MAG" ~ risk_score + 1,
                            node_type == "ARG" ~ 3,
                            TRUE ~ 4
                          ))) +
      geom_node_text(aes(label = ifelse(node_type != "MAG", name, "")),
                     repel = TRUE, size = 3) +
      scale_color_manual(values = c("MAG" = "#1f78b4", "ARG" = "#e31a1c", "MGE" = "#33a02c")) +
      theme_void()

    ggsave("R_out/fig/network_ARG_MGE_MAG_tripartite.pdf", p_tri, width = 11, height = 9)
  }
}

###############################
## 24. MAG trait correlation
###############################

trait_cols <- intersect(c("ARG_count", "mobile_ARG_strict_count", "mobile_ARG_loose_count",
                          "plasmid_like_contig_count", "IS_contig_count",
                          "transposase_total", "vfdb_vfg_count", "CAZy_count", "risk_score"),
                        colnames(mag_info2))

if(length(trait_cols) >= 2){
  cor_df <- mag_info2 %>% select(all_of(trait_cols)) %>% na.omit()
  cor_df <- filter_nonzero_sd_cols(cor_df)

  if(ncol(cor_df) >= 2){
    cor_mat <- cor(cor_df, method = "spearman")

    pdf("R_out/fig/MAG_trait_correlation_heatmap.pdf", width = 6, height = 5)
    pheatmap(cor_mat, display_numbers = TRUE)
    dev.off()
  } else {
    writeLines("Skipping MAG trait correlation heatmap: fewer than 2 non-zero-SD numeric columns.",
               "R_out/tab/MAG_trait_correlation_skipped.txt")
  }
}

###############################
## 25. Optional phylogenetic tree
###############################

tree_file1 <- "08_gtdbtk/bac120.decorated.tree"

if ("ggtree" %in% loadedNamespaces() && file.exists(tree_file1)) {
  tree <- ape::read.tree(tree_file1)

  tree_anno <- mag_info2 %>%
    select(MAG, host_enrichment, ARG_count, mobile_ARG_strict_count, risk_score, classification) %>%
    distinct()

  tree_anno2 <- tree_anno %>% filter(MAG %in% tree$tip.label)
  rownames(tree_anno2) <- tree_anno2$MAG
  tree_anno2 <- tree_anno2[tree$tip.label, , drop = FALSE]

  p_tree <- ggtree::ggtree(tree, layout = "rectangular", size = 0.2) %<+% tree_anno2 +
    ggtree::geom_tippoint(aes(color = host_enrichment, size = risk_score), alpha = 0.8) +
    scale_color_manual(values = c(
      deer_enriched = "#1f78b4",
      non_sig = "grey70",
      yak_enriched = "#e31a1c"
    ))

  ggsave("R_out/fig/MAG_phylogeny_host_risk.pdf", p_tree, width = 10, height = 18)
}

###############################
## 26. Save objects
###############################

saveRDS(meta, "R_out/rds/meta.rds")
saveRDS(arg_mat, "R_out/rds/arg_mat.rds")
saveRDS(mag_mat, "R_out/rds/mag_mat.rds")
saveRDS(mag_info2, "R_out/rds/mag_info2.rds")
saveRDS(highrisk_mag, "R_out/rds/highrisk_mag.rds")
saveRDS(rf_results, "R_out/rds/rf_results.rds")

writeLines(capture.output(sessionInfo()), "R_out/tab/sessionInfo.txt")

cat("All customized analyses completed successfully.\n")
