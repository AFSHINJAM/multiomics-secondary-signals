###############################################################################
####  MULTI-OMICS DATA PREPROCESSING                                       ####
####                                                                       ####
####  Preprocesses six omics layers from the Aracena et al. 2024 dataset  ####
####  (Nature Genetics) for multi-omics integration analysis.             ####
####                                                                       ####
####  Input:  Raw RNA-seq counts, ATAC-seq and ChIP-seq peak files        ####
####  Output: Normalized feature matrices (samples x 1000 features each)  ####
####          saved in GAUDI-compatible format + workspace.RData           ####
####                                                                       ####
####  Omics layers: RNA-seq, ATAC-seq, H3K27ac, H3K4me3, H3K4me1,        ####
####                H3K27me3                                               ####
###############################################################################

rm(list = ls())
options(warn = 1)

# ==============================================================================
# CONFIGURATION — update these paths before running
# ==============================================================================

# Root directory containing raw data
DATA_RAW_DIR <- "/path/to/raw/data"          # e.g. folder with RNAseq_counts.dat

# Directory containing ATAC-seq and ChIP-seq peak files
PEAK_BASE_DIR <- "/path/to/peak/files"       # subdirs: ATACseq_peaks, H3K27ac_peaks, etc.

# Directory containing pre-computed merged peak BED files
MERGED_PEAKS_DIR <- "/path/to/merged/peaks"  # merged_peaks.bed, merged_H3K27ac_peaks.bed, etc.

# Output directory for processed data and results
OUTPUT_DIR <- "/path/to/output"

# Directory containing GAUDI/runfactorization.R script
SCRIPTS_DIR <- "/path/to/scripts"

# Number of top variable features to retain per omics layer
N_FEATURES <- 1000

# ==============================================================================
# LIBRARIES
# ==============================================================================

library(edgeR)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(ggsci)
library(ggforce)

source(file.path(SCRIPTS_DIR, "runfactorization.R"))

# ==============================================================================
# STEP 1: RNA-SEQ — Load, normalize (TMM + log-CPM), select top variable genes
# ==============================================================================

cat("\n=== Step 1: RNA-seq ===\n")

rna_data <- read.table(file.path(DATA_RAW_DIR, "RNAseq_counts.dat"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Keep only NI and Flu samples (exclude Mock)
samples <- colnames(rna_data)[7:ncol(rna_data)]
samples <- samples[!grepl("Mock$", samples)]

split_samples <- strsplit(samples, "_")
first_part    <- sapply(split_samples, `[`, 1)
cluster       <- substr(first_part, 1, 2)
subject       <- substr(first_part, 3, nchar(first_part))
time_point    <- sapply(split_samples, function(x) tail(x, 1))

metadata_rna <- data.frame(
  sample     = samples,
  cluster    = cluster,
  subject    = subject,
  time_point = time_point,
  stringsAsFactors = FALSE
)
metadata_rna <- metadata_rna[metadata_rna$time_point %in% c("NI", "Flu"), ]

rna_data_filtered <- rna_data[, c(1, 2, which(colnames(rna_data) %in% metadata_rna$sample))]
rna_data_numeric  <- rna_data_filtered[, 3:ncol(rna_data_filtered)]
rna_data_t        <- t(rna_data_numeric)
colnames(rna_data_t) <- rna_data_filtered$Symbol

# Remove zero-variance genes
nonzero_var   <- apply(rna_data_t, 2, var) > 0
rna_data_t    <- rna_data_t[, nonzero_var]

# TMM normalization + log-CPM
dge           <- DGEList(counts = t(rna_data_t))
dge           <- calcNormFactors(dge, method = "TMM")
log_cpm       <- cpm(dge, log = TRUE, prior.count = 1)
rna_normalized <- t(log_cpm)

# Select top N_FEATURES most variable genes
gene_variances      <- apply(rna_normalized, 2, var)
top_variable_genes  <- names(sort(gene_variances, decreasing = TRUE)[1:N_FEATURES])
rna_normalized_var  <- rna_normalized[, top_variable_genes]

cat("RNA-seq prepared:", nrow(rna_normalized_var), "samples,",
    ncol(rna_normalized_var), "features\n")

# ==============================================================================
# STEP 2: ATAC-SEQ — Build count matrix from peak overlaps, normalize, select
# ==============================================================================

cat("\n=== Step 2: ATAC-seq ===\n")

atac_peak_dir <- file.path(PEAK_BASE_DIR, "ATACseq_peaks")
atac_files    <- list.files(atac_peak_dir, pattern = "\\.narrowPeak\\.bed$", full.names = TRUE)
sample_info   <- stringr::str_match(basename(atac_files),
                                    "^([A-Za-z]{2})(\\d+)_(NI|Flu|Mock)_peaks\\.narrowPeak\\.bed$")

metadata_atac <- data.frame(
  file       = atac_files,
  cluster    = sample_info[, 2],
  subject    = paste0(sample_info[, 2], sample_info[, 3]),
  time_point = sample_info[, 4],
  stringsAsFactors = FALSE
) %>%
  filter(time_point %in% c("NI", "Flu")) %>%
  mutate(
    cluster = ifelse(cluster %in% c("EU", "Eu"), "EU",
              ifelse(cluster %in% c("AF", "Af"), "AF", cluster))
  )

# Helper: read a BED/narrowPeak file and return a data.frame
read_peaks <- function(file) {
  peaks <- import.bed(file)
  as.data.frame(peaks) %>%
    mutate(peak_id = paste(seqnames, start, end, sep = "_"),
           width   = end - start) %>%
    dplyr::select(peak_id, seqnames, start, end, width, name, score)
}

peak_list <- lapply(metadata_atac$file, read_peaks)
names(peak_list) <- paste0(metadata_atac$subject, "_", metadata_atac$time_point)

# Load consensus (merged) peaks
merged_peaks    <- import.bed(file.path(MERGED_PEAKS_DIR, "merged_peaks.bed"))
merged_peak_ids <- paste0(seqnames(merged_peaks), "_",
                          start(merged_peaks), "_", end(merged_peaks))

# Build sample x peak count matrix using overlap scores
count_matrix <- matrix(0,
                       nrow = length(merged_peaks),
                       ncol = nrow(metadata_atac),
                       dimnames = list(merged_peak_ids, names(peak_list)))

for (i in seq_len(nrow(metadata_atac))) {
  sp  <- peak_list[[i]]
  sgr <- GRanges(seqnames = sp$seqnames,
                 ranges   = IRanges(start = sp$start, end = sp$end))
  ov  <- findOverlaps(merged_peaks, sgr)
  if (length(ov) > 0) {
    scores <- tapply(sp$score[subjectHits(ov)], queryHits(ov), max)
    count_matrix[as.numeric(names(scores)), i] <- scores
  }
}

# Filter: keep peaks present in at least 10% of samples
keep_peaks     <- rowSums(count_matrix > 0) >= ncol(count_matrix) * 0.1
count_filtered <- count_matrix[keep_peaks, ]

# TMM normalization + log-CPM
dge_atac  <- DGEList(counts = count_filtered)
dge_atac  <- calcNormFactors(dge_atac, method = "TMM")
tmm       <- cpm(dge_atac, log = TRUE, prior.count = 1)
top_peaks <- names(sort(apply(tmm, 1, var), decreasing = TRUE)[1:N_FEATURES])
atac_normalized_var <- t(tmm[top_peaks, ])

cat("ATAC-seq prepared:", nrow(atac_normalized_var), "samples,",
    ncol(atac_normalized_var), "features\n")

# ==============================================================================
# STEP 2B: HISTONE MARKS — same pipeline as ATAC-seq, applied to each mark
# ==============================================================================

cat("\n=== Step 2B: Histone marks ===\n")

process_histone_mark <- function(mark_name, file_pattern, merged_bed_file, top_n = 1000) {

  cat(sprintf("\nProcessing %s...\n", mark_name))

  histone_files <- list.files(file.path(PEAK_BASE_DIR, paste0(mark_name, "_peaks")),
                              pattern = file_pattern, full.names = TRUE)

  sample_info <- stringr::str_match(basename(histone_files),
                                    "^([A-Za-z]{2})(\\d+)_(NI|Flu|Mock)_")

  metadata_histone <- data.frame(
    file       = histone_files,
    cluster    = sample_info[, 2],
    subject    = paste0(sample_info[, 2], sample_info[, 3]),
    time_point = sample_info[, 4],
    stringsAsFactors = FALSE
  ) %>%
    filter(time_point %in% c("NI", "Flu")) %>%
    mutate(cluster = ifelse(cluster %in% c("EU", "Eu"), "EU",
                     ifelse(cluster %in% c("AF", "Af"), "AF", cluster)))

  peak_list <- lapply(metadata_histone$file, read_peaks)
  names(peak_list) <- paste0(metadata_histone$subject, "_", metadata_histone$time_point)

  merged_peaks    <- import.bed(merged_bed_file)
  merged_peak_ids <- paste0(seqnames(merged_peaks), "_",
                            start(merged_peaks), "_", end(merged_peaks))

  count_matrix <- matrix(0,
                         nrow = length(merged_peaks),
                         ncol = nrow(metadata_histone),
                         dimnames = list(merged_peak_ids, names(peak_list)))

  for (i in seq_len(nrow(metadata_histone))) {
    sp  <- peak_list[[i]]
    sgr <- GRanges(seqnames = sp$seqnames,
                   ranges   = IRanges(start = sp$start, end = sp$end))
    ov  <- findOverlaps(merged_peaks, sgr)
    if (length(ov) > 0) {
      scores <- tapply(sp$score[subjectHits(ov)], queryHits(ov), max)
      count_matrix[as.numeric(names(scores)), i] <- scores
    }
  }

  keep_peaks     <- rowSums(count_matrix > 0) >= ncol(count_matrix) * 0.1
  count_filtered <- count_matrix[keep_peaks, ]

  dge  <- DGEList(counts = count_filtered)
  dge  <- calcNormFactors(dge, method = "TMM")
  tmm  <- cpm(dge, log = TRUE, prior.count = 1)
  top  <- names(sort(apply(tmm, 1, var), decreasing = TRUE)[1:top_n])
  normalized_var <- t(tmm[top, ])

  cat(sprintf("%s prepared: %d samples, %d features\n",
              mark_name, nrow(normalized_var), ncol(normalized_var)))
  return(normalized_var)
}

h3k27ac_final <- process_histone_mark(
  mark_name      = "H3K27ac",
  file_pattern   = "_H3K27ac_peaks\\.narrowPeak\\.bed$",
  merged_bed_file = file.path(MERGED_PEAKS_DIR, "merged_H3K27ac_peaks.bed"),
  top_n = N_FEATURES
)

h3k4me3_final <- process_histone_mark(
  mark_name      = "H3K4me3",
  file_pattern   = "_H3K4me3_peaks\\.narrowPeak\\.bed$",
  merged_bed_file = file.path(MERGED_PEAKS_DIR, "merged_H3K4me3_peaks.bed"),
  top_n = N_FEATURES
)

h3k4me1_final <- process_histone_mark(
  mark_name      = "H3K4me1",
  file_pattern   = "_H3K4me1_peaks\\.narrowPeak\\.bed$",
  merged_bed_file = file.path(MERGED_PEAKS_DIR, "merged_H3K4me1_peaks.bed"),
  top_n = N_FEATURES
)

h3k27me3_final <- process_histone_mark(
  mark_name      = "H3K27me3",
  file_pattern   = "_H3K27me3_peaks\\.broadPeak\\.bed$",
  merged_bed_file = file.path(MERGED_PEAKS_DIR, "merged_H3K27me3_peaks.bed"),
  top_n = N_FEATURES
)

# ==============================================================================
# STEP 3: ALIGN SAMPLES — keep only samples present in ALL six omics layers
# ==============================================================================

cat("\n=== Step 3: Aligning samples across omics ===\n")

common_samples <- Reduce(intersect, list(
  rownames(rna_normalized_var),
  rownames(atac_normalized_var),
  rownames(h3k27ac_final),
  rownames(h3k4me3_final),
  rownames(h3k4me1_final),
  rownames(h3k27me3_final)
))
cat("Common samples across all 6 omics:", length(common_samples), "\n")

rna_final    <- rna_normalized_var[common_samples, ]
atac_final   <- atac_normalized_var[common_samples, ]
h3k27ac_final <- h3k27ac_final[common_samples, ]
h3k4me3_final <- h3k4me3_final[common_samples, ]
h3k4me1_final <- h3k4me1_final[common_samples, ]
h3k27me3_final <- h3k27me3_final[common_samples, ]

# Build metadata table
sample_names    <- strsplit(common_samples, "_")
sample_metadata <- data.frame(
  sample     = common_samples,
  subject    = sapply(sample_names, function(x) x[1]),
  time_point = sapply(sample_names, function(x) x[2]),
  cluster    = substr(sapply(sample_names, function(x) x[1]), 1, 2),
  stringsAsFactors = FALSE
)
sample_metadata$group <- paste0(sample_metadata$cluster, "_", sample_metadata$time_point)

cat("\nSample distribution:\n")
print(table(sample_metadata$time_point))
print(table(sample_metadata$cluster))
print(table(sample_metadata$group))

# ==============================================================================
# STEP 4: SAVE — write all matrices in GAUDI-compatible format
# ==============================================================================

cat("\n=== Step 4: Saving data ===\n")

data_folder <- file.path(OUTPUT_DIR, "data")
dir.create(data_folder, showWarnings = FALSE, recursive = TRUE)

# Helper: write matrix in GAUDI tab-delimited format (features as rows)
write_gaudi_format <- function(data, file, sep = "\t") {
  write.table(cbind(probe = colnames(data), t(data)), file,
              sep = sep, append = FALSE, quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

write_gaudi_format(rna_final,    file.path(data_folder, "omics1_RNAseq.txt"))
write_gaudi_format(atac_final,   file.path(data_folder, "omics2_ATACseq.txt"))
write_gaudi_format(h3k27ac_final, file.path(data_folder, "omics3_H3K27ac.txt"))
write_gaudi_format(h3k4me3_final, file.path(data_folder, "omics4_H3K4me3.txt"))
write_gaudi_format(h3k4me1_final, file.path(data_folder, "omics5_H3K4me1.txt"))
write_gaudi_format(h3k27me3_final, file.path(data_folder, "omics6_H3K27me3.txt"))

write.table(sample_metadata[, c("sample", "time_point")],
            file.path(data_folder, "clusters_timepoint.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE,
            col.names = c("sample.id", "cluster.id"))

write.table(sample_metadata[, c("sample", "cluster")],
            file.path(data_folder, "clusters_ancestry.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE,
            col.names = c("sample.id", "cluster.id"))

write.table(sample_metadata[, c("sample", "group")],
            file.path(data_folder, "clusters_combined.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE,
            col.names = c("sample.id", "cluster.id"))

cat("Data saved to:", data_folder, "\n")

# ==============================================================================
# STEP 5: RUN FACTORIZATION METHODS
# ==============================================================================

cat("\n=== Step 5: Running factorization methods ===\n")

results_folder <- file.path(OUTPUT_DIR, "results")
dir.create(results_folder, showWarnings = FALSE, recursive = TRUE)

factorization_methods <- c("RGCCA", "MCIA", "MOFA", "intNMF", "JIVE", "tICA", "UBMI")

cl_data          <- read.table(file.path(data_folder, "clusters_combined.txt"),
                               sep = "\t", header = TRUE, stringsAsFactors = FALSE)
min_cluster_size <- min(table(cl_data$cluster.id))

out <- runfactorization(
  data_folder,
  c("omics1_RNAseq.txt", "omics2_ATACseq.txt", "omics3_H3K27ac.txt",
    "omics4_H3K4me3.txt", "omics5_H3K4me1.txt", "omics6_H3K27me3.txt"),
  num.factors      = 2,
  n_components_conc = 2,
  compute_features = TRUE,
  min_pts          = max(5, min_cluster_size),
  sep              = "\t",
  filtering        = "none",
  methods          = factorization_methods
)

all_results <- list(Combined = out)

# Save workspace for downstream analyses (simulation, SHAP)
save.image(file.path(data_folder, "workspace.RData"))
cat("\nWorkspace saved to:", file.path(data_folder, "workspace.RData"), "\n")
cat("=== Preprocessing complete ===\n")
