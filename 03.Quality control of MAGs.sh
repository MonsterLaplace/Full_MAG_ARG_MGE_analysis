# This process is the downstream analysis of the MAGs generated from assembly and binning
#1. 按照以下的目录整理结果
#########################################################################################
#    project/
#    ├── 00_meta/
#    │   └── samples.tsv
#    ├── 01_raw/
#    ├── 02_qc/
#    ├── 03_host_removed/
#    ├── 04_assembly/         
#    ├── 05_bins_raw/         
#    ├── 06_mag_qc/
#    ├── 07_drep/
#    ├── 08_gtdbtk/
#    ├── 09_mapping/
#    ├── 10_annotation/
#    │   ├── bakta/
#    │   ├── proteins/
#    │   ├── arg/
#    │   ├── cazy/
#    │   ├── vfdb/
#    ├── 11_mge/
#    │   ├── genomad/
#    │   ├── integron/
#    │   ├── isescan/
#    ├── 12_contig_level/
#    ├── 13_tables/
#    └── 14_fig/
#    └── 15_R_out/
#    │   ├── fig/
#    │   ├── rds/
#    │   ├── tab/
#########################################################################################
#样本列表：samples.tsv
#列：sample_id host site season sex age
#宿主：yak / deer

#2. MAG 目录升级：质量评估、去污染、去冗余、分类
#已经有 MAGs 组装结果，以“所有 bins”作为输入，重建一个高质量 non-redundant catalog。
#假设：每个 bin 一个 fasta：05_bins_raw/*.fa
#2.1 CheckM2
conda activate checkm2
mkdir -p 06_mag_qc/checkm2
checkm2 predict \
  --input 05_bins_raw/ \
  --output-directory 06_mag_qc/checkm2 \
  --threads 32
#输出核心表通常在：quality_report.tsv  
#2.2 GUNC 污染检查
mkdir -p 06_mag_qc/gunc
gunc run \
  --input_dir 05_bins_raw \
  --out_dir 06_mag_qc/gunc \
  --threads 32 \
  -r /path/to/gnucdb/gunc_db_progenomes2.1.dmnd
  mkdir -p 06_mag_qc/selected

awk -F'\t' 'NR==1 || ($2>=70 && $3<=5)' 06_mag_qc/checkm2/quality_report.tsv > 06_mag_qc/high_quality.tsv
cut -f1 06_mag_qc/high_quality.tsv | tail -n +2 > 06_mag_qc/high_quality_ids.txt
#2.3 筛选高质量MAG：
#completeness ≥ 70
#contamination ≤ 5
#GUNC 通过

mkdir -p 06_mag_qc/selected

awk -F'\t' 'NR==1 || ($2>=70 && $3<=5)' 06_mag_qc/checkm2/quality_report.tsv > 06_mag_qc/high_quality.tsv
cut -f1 06_mag_qc/high_quality.tsv | tail -n +2 > 06_mag_qc/high_quality_ids.txt

while read id; do
  cp "05_bins_raw/${id}.fa" 06_mag_qc/selected/
done < 06_mag_qc/high_quality_ids.txt
#如果 GUNC 有可疑污染 bins，再排掉。
#2.4 dRep 去冗余形成非冗余 MAG catalog
mkdir -p 07_drep
dRep dereplicate 07_drep/out \
  -g 06_mag_qc/selected/*.fa \
  -p 32 \
  -comp 70 -con 5 \
  -sa 0.95
