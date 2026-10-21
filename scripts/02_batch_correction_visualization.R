###############################################################################
####    VISUALIZATION: Clusters Before vs After Batch Correction           ####
####                                                                        ####
####    Left panel:  Time point (Flu/NI) clusters — before ComBat          ####
####    Right panel: Ancestry (AF/EU) clusters   — after ComBat            ####
####                                                                        ####
####    Also generates manuscript summary figures and clustering metrics    ####
####    table comparing all integration methods.                            ####
###############################################################################

# Clear workspace
rm(list = ls())
set.seed(42)

# ==============================================================================
# CONFIGURATION — set these paths before running
# ==============================================================================

# Directory containing workspace.RData (output of 01_data_preprocessing.R)
WORKSPACE_PATH <- "/path/to/data/workspace.RData"

# Path to runfactorization.R (GAUDI/UBMI integration script)
GAUDI_SCRIPT   <- "/path/to/scripts/runfactorization.R"

# Directory where all figures and tables will be saved
OUTPUT_DIR     <- "/path/to/output/figures"

# Number of parallel cores for ComBat batch correction
NUM_CORES <- 8

# ==============================================================================
# PACKAGES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)
library(RGCCA)
library(omicade4)
library(IntNMF)
library(r.jive)
library(reticulate)
library(MOFA2)
library(sva)
library(readr)
library(future.apply)
library(mclust)
library(tidyr)

plan(multicore, workers = NUM_CORES)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")
load(WORKSPACE_PATH)

setwd(dirname(GAUDI_SCRIPT))
source(basename(GAUDI_SCRIPT))
setwd(dirname(WORKSPACE_PATH))

# Prepare list of omics matrices
omics_data_original <- list(
  "RNA-seq"   = rna_final,
  "ATAC-seq"  = atac_final,
  "H3K27ac"   = h3k27ac_final,
  "H3K4me3"   = h3k4me3_final,
  "H3K4me1"   = h3k4me1_final,
  "H3K27me3"  = h3k27me3_final
)

# Standardise ancestry labels
if (!"cluster" %in% colnames(sample_metadata)) {
  sample_metadata$cluster <- substr(sample_metadata$subject, 1, 2)
}
sample_metadata$ancestry <- ifelse(
  sample_metadata$cluster %in% c("EU", "Eu"), "EU",
  ifelse(sample_metadata$cluster %in% c("AF", "Af"), "AF",
         sample_metadata$cluster)
)

cat("Samples:", nrow(sample_metadata), "\n")
cat("Time points:", paste(names(table(sample_metadata$time_point)),
                          table(sample_metadata$time_point),
                          sep = "=", collapse = ", "), "\n")
cat("Ancestry:", paste(names(table(sample_metadata$ancestry)),
                       table(sample_metadata$ancestry),
                       sep = "=", collapse = ", "), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Apply ComBat batch correction to a single omics matrix
apply_combat <- function(omics_data, batch_labels) {
  omics_clean <- omics_data
  omics_clean[is.na(omics_clean)] <- 0
  omics_clean[is.infinite(omics_clean)] <- 0

  corrected <- tryCatch({
    t(ComBat(dat = t(omics_clean), batch = batch_labels,
             mod = NULL, par.prior = TRUE, prior.plots = FALSE))
  }, error = function(e) {
    cat("ComBat error:", e$message, "\n")
    omics_clean
  })

  corrected[is.na(corrected)] <- 0
  corrected[is.infinite(corrected)] <- 0
  corrected
}

# Write a matrix in the tab-separated format expected by runfactorization()
write_gaudi_format <- function(data, file) {
  df_out <- as.data.frame(cbind(probe = colnames(data), t(data)))
  readr::write_tsv(df_out, file)
}

# Run MOFA2 on a list of omics matrices and return a factors matrix
run_mofa_func <- function(data_list, num_factors = 2) {
  temp_hdf5 <- tempfile(pattern = "mofa_", tmpdir = tempdir(), fileext = ".hdf5")

  result <- tryCatch({
    mofa_data <- lapply(data_list, function(x) {
      mat <- t(as.matrix(x))
      mat[is.na(mat)] <- 0
      mat[is.infinite(mat)] <- 0
      mat
    })
    names(mofa_data) <- paste0("view_", seq_along(mofa_data))

    mofa_obj <- create_mofa(mofa_data)

    model_opts <- get_default_model_options(mofa_obj)
    model_opts$num_factors <- num_factors

    train_opts <- get_default_training_options(mofa_obj)
    train_opts$convergence_mode <- "fast"
    train_opts$verbose <- FALSE
    train_opts$seed <- 42

    mofa_obj <- prepare_mofa(mofa_obj,
                             model_options  = model_opts,
                             training_options = train_opts)
    mofa_obj <- MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE, outfile = temp_hdf5)

    get_factors(mofa_obj)[[1]]

  }, error = function(e) {
    cat("MOFA error:", e$message, "\n")
    NULL
  })

  if (file.exists(temp_hdf5)) file.remove(temp_hdf5)
  result
}

# Run all integration methods via runfactorization() and return a named list
# of factor matrices (samples × factors)
run_integration_methods <- function(data_list, metadata, methods_to_run,
                                    temp_base, label) {
  temp_folder <- file.path(temp_base,
                           paste0("temp_", label, "_", sample(10000:99999, 1)))
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)

  omics_files <- c()
  for (i in seq_along(data_list)) {
    filename <- sprintf("omics%d.txt", i)
    write_gaudi_format(data_list[[i]], file.path(temp_folder, filename))
    omics_files <- c(omics_files, filename)
  }

  write.table(metadata[, c("sample", "ancestry")],
              file.path(temp_folder, "clusters_ancestry.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE,
              col.names = c("sample.id", "cluster.id"))

  min_cluster_size <- min(table(metadata$ancestry))

  integration_result <- tryCatch({
    runfactorization(temp_folder, omics_files,
                     num.factors         = 2,
                     n_components_conc   = 2,
                     compute_features    = FALSE,
                     min_pts             = max(5, min_cluster_size),
                     sep                 = "\t",
                     filtering           = "none",
                     methods             = methods_to_run)
  }, error = function(e) {
    cat("GAUDI error:", e$message, "\n")
    NULL
  })

  results <- list()

  if (!is.null(integration_result)) {
    for (i in seq_along(integration_result$method)) {
      method_name <- integration_result$method[i]
      factors     <- integration_result$factorizations[[i]]

      # Unwrap nested lists (different methods return different structures)
      if (is.list(factors) && !is.data.frame(factors)) {
        for (nm in c("factors", "score", "composite", "global",
                     "global.factors", "joint")) {
          if (nm %in% names(factors)) { factors <- factors[[nm]]; break }
        }
        if (is.list(factors) && length(factors) > 0) factors <- factors[[1]]
      }

      if (is.data.frame(factors)) factors <- as.matrix(factors)

      if (is.matrix(factors) && nrow(factors) >= 2 && ncol(factors) >= 1) {
        results[[method_name]] <- factors
      }
    }
  }

  unlink(temp_folder, recursive = TRUE)
  results
}

# Compute ARI, silhouette, and cluster purity from a factor matrix
calculate_metrics <- function(factors, true_labels) {

  set.seed(42)
  k        <- length(unique(true_labels))
  clusters <- kmeans(factors[, 1:2], centers = k, nstart = 20)$cluster

  # Cluster purity
  cluster_purity <- function(clusters, true_labels) {
    n <- length(clusters)
    purity_sum <- 0
    for (k in unique(clusters)) {
      mask <- clusters == k
      if (sum(mask) > 0)
        purity_sum <- purity_sum + max(table(true_labels[mask]))
    }
    purity_sum / n
  }

  purity        <- cluster_purity(clusters, true_labels)
  ari           <- adjustedRandIndex(clusters, true_labels)
  dist_matrix   <- dist(factors[, 1:2])
  sil           <- silhouette(clusters, dist_matrix)
  silhouette_sc <- mean(sil[, 3])

  list(purity = purity, ari = ari, silhouette = silhouette_sc)
}

# ==============================================================================
# PART 1 — BEFORE / AFTER BATCH CORRECTION CLUSTER PLOTS
# ==============================================================================

cat("Running integration methods...\n\n")

methods_to_run <- c("RGCCA", "MCIA", "intNMF", "MOFA", "JIVE", "tICA", "UBMI")

# ── BEFORE: load pre-computed factors from workspace ──────────────────────────
cat("Loading pre-computed factors from workspace (BEFORE batch correction)...\n")

factors_before <- list()

for (i in seq_along(out$method)) {
  method_name <- out$method[i]
  factors     <- out$factorizations[[i]]

  if (is.list(factors) && !is.data.frame(factors)) {
    for (nm in c("factors", "score", "composite", "global", "global.factors", "joint")) {
      if (nm %in% names(factors)) { factors <- factors[[nm]]; break }
    }
    if (is.list(factors) && length(factors) > 0) factors <- factors[[1]]
  }

  if (is.data.frame(factors)) factors <- as.matrix(factors)

  if (is.matrix(factors) && nrow(factors) >= 2 && ncol(factors) >= 1) {
    factors_before[[method_name]] <- factors
    cat(sprintf("  ✓ %s: %d samples × %d factors\n",
                method_name, nrow(factors), ncol(factors)))
  }
}

cat(sprintf("\nLoaded %d methods from workspace\n\n", length(factors_before)))

# ── AFTER: apply ComBat to remove time-point (infection) signal ───────────────
cat("Applying ComBat batch correction (batch = time_point)...\n")
corrected_list <- future_lapply(omics_data_original, function(x) {
  apply_combat(x, sample_metadata$time_point)
}, future.seed = TRUE)

cat("Running integration methods AFTER batch correction...\n")
factors_after <- run_integration_methods(corrected_list, sample_metadata,
                                         methods_to_run, OUTPUT_DIR, "after")

cat("Running MOFA AFTER batch correction...\n")
mofa_factors_after <- run_mofa_func(corrected_list, num_factors = 2)
if (!is.null(mofa_factors_after)) {
  factors_after[["MOFA"]] <- mofa_factors_after
}

# ── Generate side-by-side PDF (one page per method) ───────────────────────────
cat("\nGenerating cluster plots PDF...\n")

timepoint_colors <- c("NI" = "#4DBBD5FF", "Flu" = "#E64B35FF")
ancestry_colors  <- c("AF" = "#00A087FF", "EU"  = "#F39B7FFF")
ancestry_shapes  <- c("AF" = 15,          "EU"  = 18)

pdf_file <- file.path(OUTPUT_DIR, "cluster_comparison_final.pdf")
pdf(pdf_file, width = 14, height = 6)

for (method_name in names(factors_before)) {

  if (is.null(factors_before[[method_name]]) ||
      is.null(factors_after[[method_name]])) next

  factors_b <- factors_before[[method_name]]
  factors_a <- factors_after[[method_name]]

  metrics_before <- calculate_metrics(factors_b, sample_metadata$time_point)
  metrics_after  <- calculate_metrics(factors_a, sample_metadata$ancestry)

  df_before <- data.frame(
    Factor1   = factors_b[, 1],
    Factor2   = factors_b[, 2],
    TimePoint = sample_metadata$time_point,
    Ancestry  = sample_metadata$ancestry
  )

  df_after <- data.frame(
    Factor1  = factors_a[, 1],
    Factor2  = factors_a[, 2],
    Ancestry = sample_metadata$ancestry
  )

  subtitle_before <- sprintf(
    "Color: TimePoint (batch) | Shape: Ancestry\nPurity: %.3f | ARI: %.3f | Silhouette: %.3f",
    metrics_before$purity, metrics_before$ari, metrics_before$silhouette
  )

  p1 <- ggplot(df_before,
               aes(x = Factor1, y = Factor2, color = TimePoint, shape = Ancestry)) +
    geom_point(size = 3.5, alpha = 0.8) +
    scale_color_manual(values = timepoint_colors) +
    scale_shape_manual(values = ancestry_shapes) +
    labs(title    = paste(method_name, "- BEFORE Batch Correction"),
         subtitle = subtitle_before,
         x = "Factor 1", y = "Factor 2",
         color = "Time Point", shape = "Ancestry") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom",
          plot.title       = element_text(face = "bold", size = 16),
          plot.subtitle    = element_text(size = 11, color = "gray30"))

  subtitle_after <- sprintf(
    "Ancestry Clusters (AF/EU)\nPurity: %.3f | ARI: %.3f | Silhouette: %.3f",
    metrics_after$purity, metrics_after$ari, metrics_after$silhouette
  )

  p2 <- ggplot(df_after,
               aes(x = Factor1, y = Factor2, color = Ancestry, shape = Ancestry)) +
    geom_point(size = 3.5, alpha = 0.8) +
    scale_color_manual(values = ancestry_colors) +
    scale_shape_manual(values = ancestry_shapes) +
    labs(title    = paste(method_name, "- AFTER Batch Correction"),
         subtitle = subtitle_after,
         x = "Factor 1", y = "Factor 2",
         color = "Ancestry", shape = "Ancestry") +
    theme_bw(base_size = 14) +
    theme(legend.position  = "bottom",
          plot.title       = element_text(face = "bold", size = 16),
          plot.subtitle    = element_text(size = 11, color = "gray30"))

  grid.arrange(p1, p2, ncol = 2)
}

dev.off()
cat("✓ Cluster plots saved:", pdf_file, "\n\n")

# ==============================================================================
# PART 2 — MANUSCRIPT SUMMARY FIGURES AND METRICS TABLE
# ==============================================================================

cat("Collecting metrics for all methods...\n")

all_methods  <- intersect(names(factors_before), names(factors_after))
metrics_rows <- list()

for (method_name in all_methods) {
  fb <- factors_before[[method_name]]
  fa <- factors_after[[method_name]]
  if (is.null(fb) || is.null(fa)) next

  mb <- calculate_metrics(fb, sample_metadata$time_point)
  ma <- calculate_metrics(fa, sample_metadata$ancestry)

  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    Method    = method_name, Condition = "Before",
    ARI       = mb$ari, Silhouette = mb$silhouette, Purity = mb$purity,
    stringsAsFactors = FALSE
  )
  metrics_rows[[length(metrics_rows) + 1]] <- data.frame(
    Method    = method_name, Condition = "After",
    ARI       = ma$ari, Silhouette = ma$silhouette, Purity = ma$purity,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-8s  BEFORE: ARI=%.3f Sil=%.3f Pur=%.3f\n",
              method_name, mb$ari, mb$silhouette, mb$purity))
  cat(sprintf("           AFTER:  ARI=%.3f Sil=%.3f Pur=%.3f\n",
              ma$ari, ma$silhouette, ma$purity))
}

metrics_df <- do.call(rbind, metrics_rows)
metrics_df$Method    <- factor(metrics_df$Method,    levels = all_methods)
metrics_df$Condition <- factor(metrics_df$Condition, levels = c("Before", "After"))
cat("\n✓ Metrics collected\n\n")

# ── Shared ggplot theme ───────────────────────────────────────────────────────
ms_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "#f0f0f0", colour = NA),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 9.5),
    legend.position    = "top",
    legend.key.size    = unit(0.45, "cm"),
    panel.spacing      = unit(0.9, "lines")
  )

# ── Figure 1: Before & After bars for all methods and all 3 metrics ──────────
long_df        <- pivot_longer(metrics_df,
                               cols      = c(ARI, Silhouette, Purity),
                               names_to  = "Metric",
                               values_to = "Value")
long_df$Metric <- factor(long_df$Metric, levels = c("ARI", "Silhouette", "Purity"))

fig1 <- ggplot(long_df, aes(x = Method, y = Value, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
  facet_wrap(~ Metric, ncol = 1, scales = "free_y",
             labeller = labeller(Metric = c(
               ARI        = "Adjusted Rand Index (ARI)",
               Silhouette = "Silhouette Score",
               Purity     = "Cluster Purity"
             ))) +
  scale_fill_manual(
    values = c(Before = "#7FB3D3", After = "#E07B54"),
    name   = "Batch correction:"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.06))) +
  labs(x = NULL, y = NULL) +
  ms_theme

ggsave(file.path(OUTPUT_DIR, "Fig_summary_metrics_before_after.pdf"),
       fig1, width = 7.5, height = 9, device = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "Fig_summary_metrics_before_after.png"),
       fig1, width = 7.5, height = 9, dpi = 300)
cat("✓ Figure 1 (Fig_summary_metrics_before_after) saved\n")

# ── Figure 2: Delta lollipop — change per method ─────────────────────────────
delta_df <- metrics_df %>%
  pivot_wider(names_from = Condition, values_from = c(ARI, Silhouette, Purity)) %>%
  mutate(
    dARI        = ARI_After        - ARI_Before,
    dSilhouette = Silhouette_After - Silhouette_Before,
    dPurity     = Purity_After     - Purity_Before
  ) %>%
  dplyr::select(Method, dARI, dSilhouette, dPurity) %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Delta") %>%
  mutate(
    Metric = gsub("^d", "", Metric),
    Metric = factor(Metric, levels = c("ARI", "Silhouette", "Purity"))
  )

fig2 <- ggplot(delta_df, aes(x = reorder(Method, Delta), y = Delta)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey50", linetype = "dashed") +
  geom_segment(aes(xend = reorder(Method, Delta), yend = 0),
               colour = "#5B9DC4", linewidth = 0.9, alpha = 0.7) +
  geom_point(colour = "#5B9DC4", size = 3.8) +
  facet_wrap(~ Metric, ncol = 3, scales = "free_x",
             labeller = labeller(Metric = c(
               ARI        = "Δ ARI",
               Silhouette = "Δ Silhouette",
               Purity     = "Δ Purity"
             ))) +
  coord_flip() +
  labs(x = NULL, y = "Δ metric value") +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "#f0f0f0", colour = NA),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "Fig_delta_lollipop.pdf"),
       fig2, width = 10, height = 4.5, device = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "Fig_delta_lollipop.png"),
       fig2, width = 10, height = 4.5, dpi = 300)
cat("✓ Figure 2 (Fig_delta_lollipop) saved\n")

# ── Table: Full metrics CSV ───────────────────────────────────────────────────
table_wide <- metrics_df %>%
  pivot_wider(names_from  = Condition,
              values_from = c(ARI, Silhouette, Purity)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(Method,
                ARI_Before, ARI_After,
                Silhouette_Before, Silhouette_After,
                Purity_Before, Purity_After)

write.csv(table_wide,
          file.path(OUTPUT_DIR, "Table_clustering_metrics.csv"),
          row.names = FALSE)

cat("\n=== CLUSTERING METRICS TABLE ===\n")
cat(sprintf("%-8s | %14s | %22s | %18s\n",
            "Method", "ARI (B→A)", "Silhouette (B→A)", "Purity (B→A)"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(table_wide))) {
  cat(sprintf("%-8s | %5.3f → %5.3f | %5.3f → %5.3f          | %5.3f → %5.3f\n",
              table_wide$Method[i],
              table_wide$ARI_Before[i],        table_wide$ARI_After[i],
              table_wide$Silhouette_Before[i], table_wide$Silhouette_After[i],
              table_wide$Purity_Before[i],     table_wide$Purity_After[i]))
}

# ==============================================================================
# DONE
# ==============================================================================
cat("\n=== COMPLETE — all files saved to:", OUTPUT_DIR, "===\n")
cat("  Supplementary: cluster_comparison_final.pdf\n")
cat("  Figure 1:      Fig_summary_metrics_before_after.pdf/.png\n")
cat("  Figure 2:      Fig_delta_lollipop.pdf/.png\n")
cat("  Table:         Table_clustering_metrics.csv\n")
