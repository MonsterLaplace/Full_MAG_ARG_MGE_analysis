###############################
## Full customized R workflow v4.3
## For representative MAG + ARG + MGE analysis
## Dataset: Yak vs White-lipped deer
## Beautified/updated version with smart boxplots
###############################

rm(list = ls())
options(stringsAsFactors = FALSE)

###############################
## 0. Packages
###############################

pkg_needed <- c(
  "tidyverse", "data.table", "vegan", "ggrepel", "pheatmap",
  "ComplexHeatmap", "circlize", "randomForest", "pROC",
  "igraph", "ggraph", "ape", "scales"
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

safe_wilcox_by_group <- function(df, value_col, group_col, out_txt = NULL, context = ""){
  tmp <- df[, c(value_col, group_col), drop = FALSE]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]

  if(nrow(tmp) == 0){
    msg <- paste0("Skipping ", context, ": no complete cases.")
    warning(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  g <- as.character(tmp[[group_col]])
  v <- tmp[[value_col]]

  if(length(unique(g)) != 2){
    msg <- paste0("Skipping ", context, ": ", group_col, " does not contain exactly two groups.")
    warning(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  if(length(unique(v)) <= 1 || stats::sd(v, na.rm = TRUE) == 0){
    msg <- paste0("Skipping ", context, ": ", value_col, " has zero variance.")
    warning(msg)
    if(!is.null(out_txt)) writeLines(msg, out_txt)
    return(NULL)
  }

  wt <- tryCatch(
    wilcox.test(v ~ g),
    error = function(e) e
  )

  if(inherits(wt, "error")){
    msg <- paste0("Skipping ", context, ": wilcox.test failed for ", value_col, " ~ ", group_col, "; ", wt$message)
    warning(msg)
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

should_use_log1p <- function(x,
                             require_nonnegative = TRUE,
                             zero_fraction_threshold = 0.1,
                             ratio_threshold = 20,
                             quantile_ratio_threshold = 10,
                             verbose = FALSE) {

  x <- x[is.finite(x)]
  if(length(x) == 0) return(FALSE)

  if(require_nonnegative && any(x < 0, na.rm = TRUE)) {
    if(verbose) message("Auto-scale: negative values detected, use raw scale.")
    return(FALSE)
  }

  if(length(unique(x)) <= 1 || stats::sd(x, na.rm = TRUE) == 0) {
    if(verbose) message("Auto-scale: zero variance, use raw scale.")
    return(FALSE)
  }

  x_pos <- x[x > 0]

  if(length(x_pos) == 0) {
    if(verbose) message("Auto-scale: all zeros, use raw scale.")
    return(FALSE)
  }

  zero_fraction <- mean(x == 0, na.rm = TRUE)
  max_median_ratio <- max(x, na.rm = TRUE) / max(stats::median(x, na.rm = TRUE), 1e-8)
  q99_q50_ratio <- stats::quantile(x, 0.99, na.rm = TRUE) / max(stats::quantile(x, 0.50, na.rm = TRUE), 1e-8)

  use_log <- (
    zero_fraction >= zero_fraction_threshold ||
      max_median_ratio >= ratio_threshold ||
      q99_q50_ratio >= quantile_ratio_threshold
  )

  if(verbose) {
    message(
      "Auto-scale decision for vector: ",
      "zero_fraction=", round(zero_fraction, 3),
      "; max/median=", round(max_median_ratio, 3),
      "; q99/q50=", round(q99_q50_ratio, 3),
      "; use_log1p=", use_log
    )
  }

  return(use_log)
}

smart_boxplot_with_stats <- function(df, value_col, group_col,
                                     outfile = NULL,
                                     title = NULL,
                                     xlab = NULL,
                                     ylab = NULL,
                                     transform_mode = c("auto", "log1p", "raw"),
                                     width = 5,
                                     height = 4.5,
                                     palette = NULL,
                                     verbose = TRUE) {

  transform_mode <- match.arg(transform_mode)

  plot_df <- df[, c(value_col, group_col), drop = FALSE]
  plot_df <- plot_df[complete.cases(plot_df), , drop = FALSE]

  if(nrow(plot_df) == 0){
    warning(paste0("No complete cases for plotting: ", value_col, " ~ ", group_col))
    return(NULL)
  }

  plot_df[[group_col]] <- as.character(plot_df[[group_col]])
  raw_y <- plot_df[[value_col]]

  use_log1p <- FALSE
  if(transform_mode == "log1p"){
    use_log1p <- TRUE
  } else if(transform_mode == "raw"){
    use_log1p <- FALSE
  } else if(transform_mode == "auto"){
    use_log1p <- should_use_log1p(raw_y, verbose = verbose)
  }

  if(use_log1p){
    plot_df$plot_value <- log1p(raw_y)
    ylab_use <- ifelse(is.null(ylab), paste0("log1p(", value_col, ")"), ylab)
    scale_note <- "log1p"
  } else {
    plot_df$plot_value <- raw_y
    ylab_use <- ifelse(is.null(ylab), value_col, ylab)
    scale_note <- "raw"
  }

  xlab_use <- ifelse(is.null(xlab), group_col, xlab)
  title_use <- ifelse(is.null(title), value_col, title)

  group_levels <- unique(plot_df[[group_col]])
  plot_df[[group_col]] <- factor(plot_df[[group_col]], levels = group_levels)

  if(is.null(palette)){
    if(length(group_levels) == 2){
      palette <- c("#E64B35", "#4DBBD5")
      names(palette) <- group_levels
    } else {
      palette <- scales::hue_pal()(length(group_levels))
      names(palette) <- group_levels
    }
  } else {
    if(is.null(names(palette))){
      if(length(palette) >= length(group_levels)){
        names(palette) <- group_levels
      }
    }
    palette <- palette[group_levels]
  }

  wt <- NULL
  p_label <- "Wilcoxon p = NA"
  if(length(unique(plot_df[[group_col]])) == 2 && stats::sd(raw_y, na.rm = TRUE) > 0){
    wt <- tryCatch(
      wilcox.test(raw_y ~ plot_df[[group_col]]),
      error = function(e) NULL
    )
    if(!is.null(wt)){
      p_label <- paste0("Wilcoxon p = ", signif(wt$p.value, 3))
    }
  }

  y_max <- max(plot_df$plot_value, na.rm = TRUE)
  y_min <- min(plot_df$plot_value, na.rm = TRUE)
  y_span <- max(1e-8, y_max - y_min)
  y_annot <- y_max + 0.08 * y_span

  subtitle_use <- paste0("scale: ", scale_note)

  p <- ggplot(plot_df, aes(x = .data[[group_col]], y = plot_value, fill = .data[[group_col]])) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.85,
      color = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      aes(color = .data[[group_col]]),
      width = 0.15,
      alpha = 0.75,
      size = 2,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_color_manual(values = palette, drop = FALSE) +
    annotate("text", x = 1.5, y = y_annot, label = p_label, size = 4) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 15, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5)
    ) +
    labs(
      x = xlab_use,
      y = ylab_use,
      title = title_use,
      subtitle = subtitle_use
    ) +
    coord_cartesian(clip = "off")

  if(!is.null(outfile)){
    ggsave(outfile, p, width = width, height = height)
  }

  if(verbose){
    message("Plot created: ", value_col, " ~ ", group_col, " [", scale_note, "]")
  }

  return(list(
    plot = p,
    wilcox = wt,
    use_log1p = use_log1p,
    scale_mode = scale_note
  ))
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

if (!file.exists("13_tables/all.mag_trait_catalog.with_vfdb_cazy.tsv")) {
  stop("mag_trait_catalog.tsv is required in this workflow.")
}
mag_trait <- fread("13_tables/all.mag_trait_catalog.with_vfdb_cazy.tsv") %>% as.data.frame()
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
  ARG_shannon =vegan::diversity(arg_mat, index = "shannon")
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

arg_bray <- vegan::vegdist(arg_rel, method = "bray")
arg_pcoa <- ape::pcoa(arg_bray)
arg_pcoa_df <- data.frame(
  sample = rownames(arg_rel),
  PC1 = arg_pcoa$vectors[,1],
  PC2 = arg_pcoa$vectors[,2]
) %>% 
  left_join(meta, by = "sample")
pc1_exp <- round(arg_pcoa$values$Relative_eig[1] * 100, 2)
pc2_exp <- round(arg_pcoa$values$Relative_eig[2] * 100, 2)
p_arg_pcoa <- ggplot(arg_pcoa_df, aes(PC1, PC2, color = host, shape = season)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    data = arg_pcoa_df,
    aes(PC1, PC2, color = host, group = host),
    inherit.aes = FALSE,
    level = 0.95
  ) +
  labs(
    x = paste0("PCoA1 (", pc1_exp, "%)"),
    y = paste0("PCoA2 (", pc2_exp, "%)")
  ) +
  theme_bw()

ggsave("R_out/fig/resistome_PCoA_bray.pdf", p_arg_pcoa, width = 6, height = 5)

arg_euc <- dist(arg_clr)
perm_arg_bray <- vegan::adonis2(arg_bray ~ host * season, data = meta, permutations = 999)
perm_arg_euc  <- vegan::adonis2(arg_euc ~ host * season, data = meta, permutations = 999)

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

mag_bray <- vegan::vegdist(mag_rel, method = "bray")
mag_pcoa <- cmdscale(mag_bray, k = 2, eig = TRUE)

eig_prop <- mag_pcoa$eig / sum(mag_pcoa$eig[mag_pcoa$eig > 0])

mag_pcoa_df <- data.frame(
  sample = rownames(mag_rel),
  PC1 = mag_pcoa$points[,1],
  PC2 = mag_pcoa$points[,2]
) %>% 
  left_join(meta, by = "sample")

p_mag_pcoa <- ggplot(mag_pcoa_df, aes(PC1, PC2, color = host, shape = season)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = host),
    level = 0.95
  ) +
  labs(
    x = sprintf("PCoA1 (%.2f%%)", eig_prop[1] * 100),
    y = sprintf("PCoA2 (%.2f%%)", eig_prop[2] * 100)
  ) +
  theme_bw()

ggsave("R_out/fig/MAG_PCoA_bray.pdf", p_mag_pcoa, width = 6, height = 5)

perm_mag_bray <- vegan::adonis2(mag_bray ~ host * season, data = meta, permutations = 999)
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

## 输出结果：保存文件
write.csv(mag_shared, "R_out/tab/mag_shared_host_specific.csv", row.names = FALSE)
write.csv(mag_info2, "R_out/tab/mag_info2_with_occurrence.csv", row.names = FALSE)

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

## 输出结果：保存最终结果
write.csv(mag_info2, "R_out/tab/mag_info2_with_host_enrichment.csv", row.names = FALSE)

###############################
## 17. Risk score
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
## 18. ARG host vs non-ARG host MGE comparison
###############################

mge_test_cols <- intersect(c("mobile_ARG_loose_count", "mobile_ARG_strict_count",
                             "plasmid_like_contig_count", "integrase_contig_count",
                             "IS_contig_count", "transposase_total"),
                           colnames(mag_info2))

mag_info2$ARG_host_group <- ifelse(mag_info2$is_ARG_host, "ARG_host_MAG", "non_ARG_host_MAG")

write.csv(
  data.frame(group = names(group_counts(mag_info2$ARG_host_group)),
             n = as.integer(group_counts(mag_info2$ARG_host_group))),
  "R_out/tab/ARG_host_group_counts.csv",
  row.names = FALSE
)

plot_info_list <- list()

for (cc in mge_test_cols) {

  out_file <- paste0("R_out/tab/wilcox_", cc, "_ARGhost_vs_nonARGhost.txt")

  wt <- safe_wilcox_by_group(
    df = mag_info2,
    value_col = cc,
    group_col = "ARG_host_group",
    out_txt = out_file,
    context = cc
  )

  if(!is.null(wt)){

    res_plot <- smart_boxplot_with_stats(
      df = mag_info2,
      value_col = cc,
      group_col = "ARG_host_group",
      outfile = paste0("R_out/fig/", cc, "_ARGhost_compare.pdf"),
      title = paste0(cc, " comparison"),
      xlab = "ARG host group",
      ylab = NULL,
      transform_mode = "auto",
      width = 5,
      height = 4.5,
      palette = c(
        "ARG_host_MAG" = "#E64B35",
        "non_ARG_host_MAG" = "#4DBBD5"
      ),
      verbose = TRUE
    )

    summary_tab <- safe_group_compare_table(mag_info2, cc, "ARG_host_group")
    write.csv(summary_tab,
              paste0("R_out/tab/", cc, "_ARGhost_compare_summary.csv"),
              row.names = FALSE)

    plot_info_list[[cc]] <- data.frame(
      feature = cc,
      scale_mode = res_plot$scale_mode,
      use_log1p = res_plot$use_log1p,
      stringsAsFactors = FALSE
    )
  }
}

if(length(plot_info_list) > 0){
  plot_info_df <- do.call(rbind, plot_info_list)
  write.csv(plot_info_df,
            "R_out/tab/ARGhost_compare_plot_scale_modes.csv",
            row.names = FALSE)
}

###############################
## 18A. CAZy functional coupling analysis
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
      context = "ARG-host vs non-ARG-host CAZy comparison"
    )

    p_cazy1 <- smart_boxplot_with_stats(
      df = cazy_compare1,
      value_col = cazy_total_col,
      group_col = "ARG_host_group",
      outfile = "R_out/fig/CAZy_ARGhost_vs_nonARGhost.pdf",
      title = "CAZy comparison: ARG-host vs non-ARG-host",
      xlab = "ARG host group",
      ylab = NULL,
      transform_mode = "auto",
      width = 5,
      height = 4.5,
      palette = c(
        "ARG_host_MAG" = "#E64B35",
        "non_ARG_host_MAG" = "#4DBBD5"
      ),
      verbose = TRUE
    )

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
    warning(msg)
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
      context = "host-enriched ARG-host CAZy comparison"
    )

    p_cazy2 <- smart_boxplot_with_stats(
      df = cazy_compare2,
      value_col = cazy_total_col,
      group_col = "host_enrichment",
      outfile = "R_out/fig/CAZy_host_enriched_ARGhost_compare.pdf",
      title = "CAZy comparison in host-enriched ARG-host MAGs",
      xlab = "Host enrichment",
      ylab = NULL,
      transform_mode = "auto",
      width = 5,
      height = 4.5,
      palette = c(
        "yak_enriched" = "#00A087",
        "deer_enriched" = "#3C5488"
      ),
      verbose = TRUE
    )

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
    warning(msg)
    writeLines(msg, "R_out/tab/CAZy_host_enriched_ARGhost_skipped.txt")
  }

} else {
  msg <- "No CAZy total column detected in mag_info2. Skip CAZy functional coupling analysis."
  warning(msg)
  writeLines(msg, "R_out/tab/CAZy_analysis_skipped.txt")
}

###############################
## 18B. VFDB coupling analysis
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

    # VFDB plots: ARG-host vs non-ARG-host
    for(cc in vf_cols_main){
      tmp_plot <- vf_compare1[, c(cc, "ARG_host_group"), drop = FALSE]
      tmp_plot <- tmp_plot[complete.cases(tmp_plot), , drop = FALSE]

      if(nrow(tmp_plot) > 0 && length(unique(tmp_plot$ARG_host_group)) == 2){
        smart_boxplot_with_stats(
          df = tmp_plot,
          value_col = cc,
          group_col = "ARG_host_group",
          outfile = paste0("R_out/fig/", cc, "_VFDB_ARGhost_compare.pdf"),
          title = paste0(cc, ": ARG-host vs non-ARG-host"),
          xlab = "ARG host group",
          ylab = NULL,
          transform_mode = "auto",
          width = 5,
          height = 4.5,
          palette = c(
            "ARG_host_MAG" = "#E64B35",
            "non_ARG_host_MAG" = "#4DBBD5"
          ),
          verbose = TRUE
        )
      }
    }

  } else {
    writeLines("Skipping VFDB ARG-host vs non-ARG-host comparison: only one ARG_host_group present.",
               "R_out/tab/VFDB_ARGhost_vs_nonARGhost_skipped.txt")
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

    # VFDB plots: yak-enriched vs deer-enriched ARG-host MAGs
    for(cc in vf_cols_main){
      tmp_plot <- vf_compare2[, c(cc, "host_enrichment"), drop = FALSE]
      tmp_plot <- tmp_plot[complete.cases(tmp_plot), , drop = FALSE]

      if(nrow(tmp_plot) > 0 && length(unique(tmp_plot$host_enrichment)) == 2){
        smart_boxplot_with_stats(
          df = tmp_plot,
          value_col = cc,
          group_col = "host_enrichment",
          outfile = paste0("R_out/fig/", cc, "_VFDB_host_enriched_ARGhost_compare.pdf"),
          title = paste0(cc, ": yak-enriched vs deer-enriched ARG-host MAGs"),
          xlab = "Host enrichment",
          ylab = NULL,
          transform_mode = "auto",
          width = 5,
          height = 4.5,
          palette = c(
            "yak_enriched" = "#00A087",
            "deer_enriched" = "#3C5488"
          ),
          verbose = TRUE
        )
      }
    }

  } else {
    writeLines("Skipping VFDB host-enriched ARG-host comparison: not enough yak/deer enriched ARG-host MAGs.",
               "R_out/tab/VFDB_host_enriched_ARGhost_skipped.txt")
  }

} else {
  writeLines("Skipping VFDB analysis: no VFDB columns found in mag_info2.",
             "R_out/tab/VFDB_analysis_skipped.txt")
}

###############################
## 19. Random forest
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
## 20. Linear models
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
## 21. ARG-MAG bipartite network
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
## 22. ARG-MGE-MAG tripartite network
###############################

if (file.exists("13_tables/arg_mge_mag_edges.tsv")) {

  arg_mge_mag <- fread("13_tables/arg_mge_mag_edges.tsv") %>% as.data.frame()

  # 若没有 MGE_type，则从现有MGE证据列自动构建
  if (all(c("MAG", "ARG", "MGE_type") %in% colnames(arg_mge_mag))) {

    edge_long <- arg_mge_mag %>%
      dplyr::select(MAG, ARG, MGE_type) %>%
      dplyr::filter(!is.na(MAG), !is.na(ARG), !is.na(MGE_type), MGE_type != "")

  } else if (all(c("MAG", "ARG") %in% colnames(arg_mge_mag))) {

    edge_list <- list()

    if ("plasmid_like" %in% colnames(arg_mge_mag)) {
      edge_list[["plasmid"]] <- arg_mge_mag %>%
        dplyr::filter(plasmid_like == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "Plasmid")
    }

    if ("virus_like" %in% colnames(arg_mge_mag)) {
      edge_list[["virus"]] <- arg_mge_mag %>%
        dplyr::filter(virus_like == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "Virus")
    }

    if ("integron_present" %in% colnames(arg_mge_mag)) {
      edge_list[["integron"]] <- arg_mge_mag %>%
        dplyr::filter(integron_present == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "Integron")
    }

    if ("integrase_present" %in% colnames(arg_mge_mag)) {
      edge_list[["integrase"]] <- arg_mge_mag %>%
        dplyr::filter(integrase_present == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "Integrase")
    }

    if ("IS_present" %in% colnames(arg_mge_mag)) {
      edge_list[["IS"]] <- arg_mge_mag %>%
        dplyr::filter(IS_present == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "IS element")
    }

    if ("transposase_count" %in% colnames(arg_mge_mag)) {
      edge_list[["transposase"]] <- arg_mge_mag %>%
        dplyr::filter(transposase_count > 0) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "Transposase")
    }

    if ("TIR_present" %in% colnames(arg_mge_mag)) {
      edge_list[["TIR"]] <- arg_mge_mag %>%
        dplyr::filter(TIR_present == 1) %>%
        dplyr::transmute(MAG, ARG, MGE_type = "TIR")
    }

    edge_long <- dplyr::bind_rows(edge_list) %>%
      dplyr::distinct()

  } else {
    edge_long <- NULL
  }

  if (!is.null(edge_long) && nrow(edge_long) > 0) {

    # 若有 ARG_std，优先用简化名称
    if ("ARG_std" %in% colnames(arg_mge_mag)) {
      arg_label_map <- arg_mge_mag %>%
        dplyr::select(ARG, ARG_std) %>%
        dplyr::filter(!is.na(ARG), !is.na(ARG_std), ARG_std != "") %>%
        dplyr::distinct()

      edge_long <- edge_long %>%
        left_join(arg_label_map, by = "ARG") %>%
        dplyr::mutate(ARG_label = ifelse(!is.na(ARG_std) & ARG_std != "", ARG_std, ARG))
    } else {
      edge_long$ARG_label <- edge_long$ARG
    }

    # 节点表
    mag_nodes3 <- data.frame(name = unique(edge_long$MAG), node_type = "MAG") %>%
      left_join(
        mag_info2 %>% dplyr::select(MAG, host_enrichment, risk_score),
        by = c("name" = "MAG")
      )

    arg_nodes3 <- data.frame(
      name = unique(edge_long$ARG_label),
      node_type = "ARG"
    )

    mge_nodes3 <- data.frame(
      name = unique(edge_long$MGE_type),
      node_type = "MGE"
    )

    nodes3 <- bind_rows(mag_nodes3, arg_nodes3, mge_nodes3)

    # 边表
    edges_mag_arg <- edge_long %>%
      dplyr::select(from = MAG, to = ARG_label) %>%
      dplyr::distinct()

    edges_arg_mge <- edge_long %>%
      dplyr::select(from = ARG_label, to = MGE_type) %>%
      dplyr::distinct()

    edges3 <- bind_rows(edges_mag_arg, edges_arg_mge) %>%
      dplyr::distinct()

    if (nrow(edges3) > 0) {

      g_tri <- graph_from_data_frame(edges3, vertices = nodes3, directed = FALSE)

      V(g_tri)$degree <- igraph::degree(g_tri)

      node_df <- igraph::as_data_frame(g_tri, what = "vertices") %>%
        dplyr::mutate(
          risk_score = ifelse(is.na(risk_score), 0, risk_score)
        )

      # 将 MAG 风险值分级，用于 size 图例
      mag_risk_vals <- node_df$risk_score[node_df$node_type == "MAG"]

      if (length(mag_risk_vals[is.finite(mag_risk_vals)]) >= 3 &&
          length(unique(mag_risk_vals[is.finite(mag_risk_vals)])) >= 3) {

        qs <- quantile(mag_risk_vals, probs = c(0.33, 0.67), na.rm = TRUE)

        node_df <- node_df %>%
          dplyr::mutate(
            size_group = dplyr::case_when(
              node_type == "MAG" & risk_score <= qs[1] ~ "MAG: low risk",
              node_type == "MAG" & risk_score > qs[1] & risk_score <= qs[2] ~ "MAG: medium risk",
              node_type == "MAG" & risk_score > qs[2] ~ "MAG: high risk",
              node_type == "ARG" ~ "ARG",
              node_type == "MGE" ~ "MGE",
              TRUE ~ "Other"
            )
          )
      } else {
        node_df <- node_df %>%
          dplyr::mutate(
            size_group = dplyr::case_when(
              node_type == "MAG" ~ "MAG",
              node_type == "ARG" ~ "ARG",
              node_type == "MGE" ~ "MGE",
              TRUE ~ "Other"
            )
          )
      }

      # 标签策略：MGE全部标；ARG只标高连接度；MAG不标
      arg_deg <- node_df$degree[node_df$node_type == "ARG"]
      arg_degree_cutoff <- if (length(arg_deg) > 0) quantile(arg_deg, 0.75, na.rm = TRUE) else Inf
      if (!is.finite(arg_degree_cutoff)) arg_degree_cutoff <- 2

      node_df <- node_df %>%
        dplyr::mutate(
          label = dplyr::case_when(
            node_type == "MGE" ~ name,
            node_type == "ARG" & degree >= arg_degree_cutoff ~ name,
            TRUE ~ NA_character_
          )
        )

      g_tri2 <- graph_from_data_frame(edges3, vertices = node_df, directed = FALSE)

      # size 图例定义
      if (any(grepl("^MAG: ", node_df$size_group))) {
        size_values <- c(
          "MAG: low risk" = 4,
          "MAG: medium risk" = 6,
          "MAG: high risk" = 8,
          "ARG" = 4.3,
          "MGE" = 5.3
        )
        size_breaks <- c("MAG: low risk", "MAG: medium risk", "MAG: high risk", "ARG", "MGE")
      } else {
        size_values <- c(
          "MAG" = 6,
          "ARG" = 4.3,
          "MGE" = 5.3
        )
        size_breaks <- c("MAG", "ARG", "MGE")
      }

      p_tri <- ggraph(g_tri2, layout = "stress") +
        geom_edge_link(
          colour = "grey75",
          alpha = 0.35,
          linewidth = 0.45
        ) +
        geom_node_point(
          aes(fill = node_type, shape = node_type, size = size_group),
          color = "black",
          alpha = 0.95,
          stroke = 0.45
        ) +
        geom_node_text(
          aes(label = label),
          repel = TRUE,
          size = 3.2,
          color = "black",
          max.overlaps = Inf
        ) +
        scale_fill_manual(
          name = "Node type",
          values = c(
            "MAG" = "#4C78A8",
            "ARG" = "#E45756",
            "MGE" = "#54A24B"
          )
        ) +
        scale_shape_manual(
          name = "Node type",
          values = c(
            "MAG" = 21,
            "ARG" = 22,
            "MGE" = 24
          )
        ) +
        scale_size_manual(
          name = "Node importance",
          values = size_values,
          breaks = size_breaks
        ) +
        guides(
          fill = guide_legend(
            order = 1,
            override.aes = list(
              shape = 21,
              size = 5,
              color = "black",
              alpha = 1
            )
          ),
          shape = guide_legend(
            order = 2,
            override.aes = list(
              fill = "grey70",
              size = 5,
              color = "black",
              alpha = 1
            )
          ),
          size = guide_legend(
            order = 3,
            override.aes = list(
              shape = 21,
              fill = "grey70",
              color = "black",
              alpha = 1
            )
          )
        ) +
        labs(
          title = "ARG–MGE–MAG tripartite network",
          subtitle = "Node fill and shape indicate category; MAG node size reflects risk level",
          caption = "ARG labels are simplified using ARG_std when available. Only high-degree ARGs are labeled."
        ) +
        theme_void(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
          plot.caption = element_text(size = 9, color = "grey40"),
          legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.2, "cm")
        )

      ggsave(
        "R_out/fig/network_ARG_MGE_MAG_tripartite.pdf",
        p_tri,
        width = 12,
        height = 9
      )

      write.table(
        edge_long,
        "R_out/tab/ARG_MGE_MAG_tripartite_edges_used.tsv",
        sep = "\t", quote = FALSE, row.names = FALSE
      )

    } else {
      writeLines(
        "Tripartite network skipped: no valid edges after processing.",
        "R_out/tab/ARG_MGE_MAG_tripartite_skipped.txt"
      )
    }

  } else {
    writeLines(
      "Tripartite network skipped: no MGE_type column and no positive MGE evidence columns found.",
      "R_out/tab/ARG_MGE_MAG_tripartite_skipped.txt"
    )
  }

} else {
  writeLines(
    "Tripartite network skipped: file 13_tables/arg_mge_mag_edges.tsv not found.",
    "R_out/tab/ARG_MGE_MAG_tripartite_skipped.txt"
  )
}

###############################
## 23. MAG trait correlation
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
## 24. Optional phylogenetic tree
###############################

tree_file1 <- "08_gtdbtk/phylophlan.tree"

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
## 25. Save objects
###############################

saveRDS(meta, "R_out/rds/meta.rds")
saveRDS(arg_mat, "R_out/rds/arg_mat.rds")
saveRDS(mag_mat, "R_out/rds/mag_mat.rds")
saveRDS(mag_info2, "R_out/rds/mag_info2.rds")
saveRDS(highrisk_mag, "R_out/rds/highrisk_mag.rds")
saveRDS(rf_results, "R_out/rds/rf_results.rds")

writeLines(capture.output(sessionInfo()), "R_out/tab/sessionInfo.txt")

cat("All customized analyses completed successfully.\n")
