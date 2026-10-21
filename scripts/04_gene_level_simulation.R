###############################################################################
####   GENE-LEVEL SIMULATION - COMPARING GENE SETS                         ####
####   Testing 238 genes (≥3 omics) vs 54 genes (≥4 omics)                ####
####   Metrics: ARI, Silhouette, Purity                                    ####
####   Fig3 generated automatically at the end                             ####
###############################################################################

rm(list = ls())
set.seed(42)

# ==============================================================================
# CONFIGURATION — set these paths before running
# ==============================================================================

RUN_ID   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
RUN_NAME <- "gene_level_comparison"

# Path to gene selection results RData (gene_selection_results object)
GENE_SELECTION_PATH   <- "/path/to/data/gene_selection_feng_method.RData"

# Path to simulation target genes RData (simulation_target_genes object)
SIMULATION_GENES_PATH <- "/path/to/data/simulation_target_genes.RData"

# Path to runfactorization.R (GAUDI/UBMI integration script)
GAUDI_SCRIPT          <- "/path/to/scripts/runfactorization.R"

# Path to results_raw.csv from simulation_random_feature_injection.R
# (used to build Panel A of Figure 3 — set to NULL to skip Panel A)
RANDOM_INJECTION_CSV  <- "/path/to/results/simulations/results_raw.csv"

# Directory where per-run results are saved
OUTPUT_DIR <- file.path("/path/to/results/simulations",
                        paste0(RUN_ID, "_", RUN_NAME))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Directory where the final Fig3 is saved
FIGURES_DIR <- "/path/to/results/figures"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

config <- list(
  run_id   = RUN_ID,
  run_name = RUN_NAME,

  # Simulation parameters
  n_simulations = 10,
  prop_affected = c(0.2, 0.5, 1.0),
  effect_sizes  = c(0.4, 0.8),
  gene_sets     = c("genes_3plus", "genes_4plus"),

  # Methods passed to runfactorization()
  methods = c("RGCCA", "MCIA", "intNMF", "JIVE", "tICA", "UBMI"),
  run_mofa_separately = TRUE,

  # Paths (stored for reference in log)
  output_dir            = OUTPUT_DIR,
  gene_selection_path   = GENE_SELECTION_PATH,
  simulation_genes_path = SIMULATION_GENES_PATH,
  gaudi_script          = GAUDI_SCRIPT
)

cat("========================================\n")
cat("SIMULATION RUN:", RUN_ID, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("========================================\n\n")

# ==============================================================================
# PACKAGES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(sva)
library(cluster)
library(mclust)
library(fpc)
library(pheatmap)
library(RGCCA)
library(omicade4)
library(IntNMF)
library(r.jive)
library(reticulate)
library(MOFA2)

# ==============================================================================
# LOGGING
# ==============================================================================

LOG_FILE <- file.path(OUTPUT_DIR, "LOG.txt")

write_log <- function(message, also_print = TRUE) {
  timestamp   <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_message <- paste(timestamp, message)
  cat(log_message, "\n", file = LOG_FILE, append = TRUE)
  if (also_print) cat(log_message, "\n")
}

write_log("Starting gene-level simulation comparing gene sets...")

# ==============================================================================
# LOAD DATA
# ==============================================================================

write_log("Loading gene selection results...")
load(config$gene_selection_path)

write_log("Loading simulation target genes...")
load(config$simulation_genes_path)

write_log("Sourcing GAUDI functions...")
setwd(dirname(GAUDI_SCRIPT))
source(basename(GAUDI_SCRIPT))
setwd(dirname(GENE_SELECTION_PATH))

# Build omics list from gene-level matrices
omics_data_original <- list(
  "RNA-seq"   = gene_selection_results$rna_genes_selected,
  "ATAC-seq"  = gene_selection_results$atac_genes_selected,
  "H3K27ac"   = gene_selection_results$h3k27ac_genes_selected,
  "H3K4me3"   = gene_selection_results$h3k4me3_genes_selected,
  "H3K4me1"   = gene_selection_results$h3k4me1_genes_selected,
  "H3K27me3"  = gene_selection_results$h3k27me3_genes_selected
)

write_log("Omics data (gene-level) loaded:")
for (name in names(omics_data_original)) {
  write_log(sprintf("  %s: %d samples x %d genes", name,
                    nrow(omics_data_original[[name]]),
                    ncol(omics_data_original[[name]])))
}

# Standardise ancestry labels
sample_metadata <- gene_selection_results$sample_metadata

if (!"cluster" %in% colnames(sample_metadata)) {
  sample_metadata$cluster <- substr(sample_metadata$subject, 1, 2)
}
sample_metadata$ancestry <- ifelse(
  sample_metadata$cluster %in% c("EU", "Eu"), "EU",
  ifelse(sample_metadata$cluster %in% c("AF", "Af"), "AF",
         sample_metadata$cluster)
)

write_log(paste("Ancestry:", paste(names(table(sample_metadata$ancestry)),
                                   table(sample_metadata$ancestry),
                                   sep = "=", collapse = ", ")))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

get_ancestry_encoding <- function(metadata) {
  ifelse(metadata$ancestry == "AF", 1, 0)
}

# Inject ancestry signal into a subset of target genes within one omics layer.
# For each affected gene g in sample i:
#
#   x_tilde_ig = x_ig + ES * (a_i - a_bar) * sigma_g + epsilon_ig
#
# where ES = effect size (equals Cohen's d because signal is scaled by sigma_g),
#       a_i in {0,1} is ancestry indicator (1=AF, 0=EU),
#       a_bar = mean(a_i)  →  centered_ancestry ∈ {+0.5, -0.5},
#       sigma_g = SD of gene g across samples,
#       epsilon_ig ~ N(0, 0.1 * sigma_g)  [10% noise relative to gene SD]
inject_ancestry_signal_genes <- function(omics_data, ancestry_encoding,
                                         target_genes_all, prop_affected,
                                         effect_size) {
  available_genes   <- colnames(omics_data)
  available_targets <- intersect(target_genes_all, available_genes)

  if (length(available_targets) == 0)
    return(list(data = omics_data, affected_genes = character(0)))

  n_samples      <- nrow(omics_data)
  n_affected     <- max(1, round(length(available_targets) * prop_affected))
  affected_genes <- sample(available_targets, n_affected, replace = FALSE)

  perturbed_data    <- omics_data
  centered_ancestry <- ancestry_encoding - mean(ancestry_encoding, na.rm = TRUE)

  for (gene_name in affected_genes) {
    gene_values <- omics_data[, gene_name]
    gene_sd     <- sd(gene_values, na.rm = TRUE)
    if (is.na(gene_sd) || gene_sd == 0) gene_sd <- 1
    noise        <- rnorm(n_samples, mean = 0, sd = 0.1 * gene_sd)
    perturbation <- effect_size * centered_ancestry * gene_sd + noise
    perturbed_data[, gene_name] <- gene_values + perturbation
  }

  list(data = perturbed_data, affected_genes = affected_genes)
}

apply_combat <- function(omics_data, batch_labels) {
  omics_clean <- omics_data
  omics_clean[is.na(omics_clean)]      <- 0
  omics_clean[is.infinite(omics_clean)] <- 0

  corrected <- tryCatch({
    t(ComBat(dat = t(omics_clean), batch = batch_labels,
             mod = NULL, par.prior = TRUE, prior.plots = FALSE))
  }, error = function(e) {
    write_log(paste("ComBat error:", e$message))
    omics_clean
  })

  corrected[is.na(corrected)]      <- 0
  corrected[is.infinite(corrected)] <- 0
  corrected
}

# Compute ARI, Silhouette, and Purity from a factor matrix and true labels
calculate_metrics <- function(factors, true_labels) {
  if (is.null(factors))
    return(list(ARI = NA, Silhouette = NA, Purity = NA))

  if (is.list(factors) && !is.data.frame(factors)) {
    for (nm in c("factors", "score", "composite", "global",
                 "global.factors", "joint")) {
      if (nm %in% names(factors)) { factors <- factors[[nm]]; break }
    }
    if (is.list(factors) && length(factors) > 0) factors <- factors[[1]]
  }

  if (is.data.frame(factors)) factors <- as.matrix(factors)
  if (!is.matrix(factors) || nrow(factors) < 2 || ncol(factors) < 1)
    return(list(ARI = NA, Silhouette = NA, Purity = NA))

  n_factors   <- min(5, ncol(factors))
  factor_data <- factors[, 1:n_factors, drop = FALSE]
  factor_data[is.na(factor_data)]      <- 0
  factor_data[is.infinite(factor_data)] <- 0

  valid_idx <- !is.na(true_labels)
  if (sum(valid_idx) < 10)
    return(list(ARI = NA, Silhouette = NA, Purity = NA))

  factor_data <- factor_data[valid_idx, , drop = FALSE]
  true_labels <- true_labels[valid_idx]

  k <- length(unique(true_labels))
  if (k < 2) k <- 2

  set.seed(123)
  km <- tryCatch(
    kmeans(factor_data, centers = k, nstart = 25, iter.max = 1000),
    error = function(e) NULL
  )
  if (is.null(km))
    return(list(ARI = NA, Silhouette = NA, Purity = NA))

  pred    <- km$cluster
  metrics <- list(ARI = NA, Silhouette = NA, Purity = NA)

  tryCatch({ metrics$ARI <- adjustedRandIndex(pred, true_labels) },
           error = function(e) {})
  tryCatch({
    metrics$Silhouette <- mean(
      silhouette(pred, dist(factor_data))[, 3], na.rm = TRUE)
  }, error = function(e) {})
  tryCatch({
    cm <- table(pred, true_labels)
    metrics$Purity <- sum(apply(cm, 1, max)) / length(pred)
  }, error = function(e) {})

  metrics
}

write_gaudi_format <- function(data, file) {
  write.table(cbind(probe = colnames(data), t(data)), file,
              sep = "\t", append = FALSE, quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

run_mofa_func <- function(data_list, num_factors = 2) {
  temp_hdf5 <- tempfile(pattern = "mofa_", tmpdir = tempdir(), fileext = ".hdf5")

  result <- tryCatch({
    mofa_data <- lapply(data_list, function(x) {
      mat <- t(as.matrix(x))
      mat[is.na(mat)]      <- 0
      mat[is.infinite(mat)] <- 0
      mat
    })
    names(mofa_data) <- paste0("view_", seq_along(mofa_data))

    mofa_obj   <- create_mofa(mofa_data)
    model_opts <- get_default_model_options(mofa_obj)
    model_opts$num_factors <- num_factors

    train_opts <- get_default_training_options(mofa_obj)
    train_opts$convergence_mode <- "fast"
    train_opts$verbose          <- FALSE
    train_opts$seed             <- 42

    mofa_obj <- prepare_mofa(mofa_obj,
                             model_options    = model_opts,
                             training_options = train_opts)
    mofa_obj <- MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE,
                                outfile = temp_hdf5)
    get_factors(mofa_obj)[[1]]

  }, error = function(e) {
    write_log(paste("MOFA error:", e$message))
    NULL
  })

  if (file.exists(temp_hdf5)) file.remove(temp_hdf5)
  result
}

run_integration_methods <- function(data_list, metadata, config,
                                    temp_base, label) {
  temp_folder <- file.path(temp_base,
                           paste0("temp_", label, "_", sample(10000:99999, 1)))
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)

  omics_files <- c()
  for (i in seq_along(data_list)) {
    filename    <- sprintf("omics%d.txt", i)
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
                     num.factors       = 2,
                     n_components_conc = 2,
                     compute_features  = FALSE,
                     min_pts           = max(5, min_cluster_size),
                     sep               = "\t",
                     filtering         = "none",
                     methods           = config$methods)
  }, error = function(e) {
    write_log(paste("GAUDI error:", e$message))
    NULL
  })

  results <- data.frame()

  if (!is.null(integration_result)) {
    for (i in seq_along(integration_result$method)) {
      method_name <- integration_result$method[i]
      factors     <- integration_result$factorizations[[i]]
      if (is.list(factors) && length(factors) > 0) factors <- factors[[1]]

      metrics <- calculate_metrics(factors, metadata$ancestry)

      results <- rbind(results, data.frame(
        method     = method_name,
        ARI        = metrics$ARI,
        Silhouette = metrics$Silhouette,
        Purity     = metrics$Purity,
        stringsAsFactors = FALSE
      ))
    }
  }

  unlink(temp_folder, recursive = TRUE)

  if (config$run_mofa_separately) {
    mofa_factors <- run_mofa_func(data_list, num_factors = 2)
    mofa_metrics <- calculate_metrics(mofa_factors, metadata$ancestry)

    results <- rbind(results, data.frame(
      method     = "MOFA",
      ARI        = mofa_metrics$ARI,
      Silhouette = mofa_metrics$Silhouette,
      Purity     = mofa_metrics$Purity,
      stringsAsFactors = FALSE
    ))
  }

  results
}

# ==============================================================================
# RUN SCENARIO
# ==============================================================================

run_scenario <- function(params, omics_original, metadata, config,
                         simulation_target_genes) {

  set.seed(42 + params$row_idx * 1000 + params$sim_id)

  target_genes <- simulation_target_genes[[params$gene_set]]
  n_genes      <- length(target_genes)

  write_log(sprintf("  Gene set: %s (%d genes)", params$gene_set, n_genes))

  ancestry_encoding  <- get_ancestry_encoding(metadata)
  perturbed_list     <- list()
  all_affected_genes <- list()

  # Inject ancestry signal into every omics layer using target genes
  for (omics_name in names(omics_original)) {
    result <- inject_ancestry_signal_genes(
      omics_original[[omics_name]], ancestry_encoding,
      target_genes, params$prop_affected, params$effect_size
    )
    perturbed_list[[omics_name]]     <- result$data
    all_affected_genes[[omics_name]] <- result$affected_genes
    if (length(result$affected_genes) > 0)
      write_log(sprintf("    %s: %d genes affected",
                        omics_name, length(result$affected_genes)))
  }

  # Remove dominant time-point signal via ComBat
  write_log("  Applying ComBat...")
  corrected_list <- lapply(perturbed_list, function(x) {
    apply_combat(x, metadata$time_point)
  })

  # Run all integration methods and collect metrics
  write_log("  Running integration methods...")
  results <- run_integration_methods(corrected_list, metadata, config,
                                     config$output_dir,
                                     paste0("sim_", params$row_idx))

  results$prop_affected <- params$prop_affected
  results$effect_size   <- params$effect_size
  results$gene_set      <- params$gene_set
  results$n_genes       <- n_genes
  results$sim_id        <- params$sim_id

  write_log(sprintf("    Mean ARI: %.3f", mean(results$ARI, na.rm = TRUE)))

  results
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================

write_log("")
write_log("========================================")
write_log("STARTING SIMULATION")
write_log("========================================")

param_grid <- expand.grid(
  prop_affected = config$prop_affected,
  effect_size   = config$effect_sizes,
  gene_set      = config$gene_sets,
  sim_id        = 1:config$n_simulations,
  stringsAsFactors = FALSE
)
param_grid$row_idx <- 1:nrow(param_grid)

n_total <- nrow(param_grid)
write_log(sprintf("Total scenarios: %d", n_total))
write_log("")

start_time  <- Sys.time()
all_results <- data.frame()

for (i in 1:n_total) {
  params <- param_grid[i, ]

  msg <- sprintf("[%d/%d] Effect=%.2f, Prop=%.0f%%, GeneSet=%s, Sim=%d",
                 i, n_total, params$effect_size, params$prop_affected * 100,
                 params$gene_set, params$sim_id)
  write_log(msg)

  result <- tryCatch({
    run_scenario(params, omics_data_original, sample_metadata,
                 config, simulation_target_genes)
  }, error = function(e) {
    write_log(sprintf("  ERROR: %s", e$message))
    data.frame()
  })

  all_results <- rbind(all_results, result)

  # Save checkpoint every 5 scenarios
  if (i %% 5 == 0) {
    elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    remaining <- (elapsed / i) * (n_total - i)
    write_log(sprintf("  Progress: %.1f min elapsed, ~%.1f min remaining",
                      elapsed, remaining))
    write.csv(all_results,
              file.path(OUTPUT_DIR, "results_intermediate.csv"),
              row.names = FALSE)
  }
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write_log("")
write_log("SAVING RESULTS...")

results_file <- file.path(OUTPUT_DIR, "results_raw.csv")
write.csv(all_results, results_file, row.names = FALSE)
write_log(sprintf("Raw results: %s", results_file))

summary_results <- all_results %>%
  group_by(method, effect_size, prop_affected, gene_set, n_genes) %>%
  summarise(
    mean_ARI        = mean(ARI,        na.rm = TRUE),
    sd_ARI          = sd(ARI,          na.rm = TRUE),
    mean_Silhouette = mean(Silhouette, na.rm = TRUE),
    mean_Purity     = mean(Purity,     na.rm = TRUE),
    n               = n(),
    .groups = "drop"
  )

summary_file <- file.path(OUTPUT_DIR, "results_summary.csv")
write.csv(summary_results, summary_file, row.names = FALSE)
write_log(sprintf("Summary: %s", summary_file))

# ==============================================================================
# SUPPLEMENTARY PLOTS — per gene set (ARI, Silhouette, Purity)
# ==============================================================================

write_log("GENERATING SUPPLEMENTARY PLOTS...")

method_colors <- c(
  "RGCCA"  = "#E64B35FF", "MCIA"   = "#4DBBD5FF", "intNMF" = "#3C5488FF",
  "JIVE"   = "#F39B7FFF", "tICA"   = "#8491B4FF", "UBMI"   = "#91D1C2FF",
  "MOFA"   = "#00A087FF"
)

metric_formulas <- list(
  ARI        = "ARI = (RI - Expected_RI) / (max(RI) - Expected_RI)",
  Silhouette = "Silhouette = (b - a) / max(a, b)",
  Purity     = "Purity = (1/N) * Σ max_j |cluster_i ∩ class_j|"
)

create_metric_boxplot <- function(data, metric_col, metric_name,
                                  formula, colors, n_genes) {
  medians <- data %>%
    group_by(method, effect_size, prop_affected) %>%
    summarise(median_val = median(.data[[metric_col]], na.rm = TRUE),
              .groups = "drop")

  data$prop_label    <- factor(sprintf("Prop=%.0f%%", data$prop_affected * 100),
                               levels = c("Prop=20%", "Prop=50%", "Prop=100%"))
  medians$prop_label <- factor(sprintf("Prop=%.0f%%", medians$prop_affected * 100),
                               levels = c("Prop=20%", "Prop=50%", "Prop=100%"))

  ggplot(data, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    geom_text(data = medians,
              aes(x = method, y = median_val,
                  label = sprintf("%.2f", median_val)),
              vjust = -0.8, size = 3, fontface = "bold", color = "black") +
    facet_grid(effect_size ~ prop_label,
               labeller = labeller(
                 effect_size = function(x) paste0("Effect = ", x))) +
    scale_fill_manual(values = colors) +
    labs(title    = sprintf("%s (%d genes)", metric_name, n_genes),
         subtitle = formula,
         x = "", y = metric_name) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      legend.position  = "none",
      strip.text.x     = element_text(size = 11, face = "bold"),
      strip.text.y     = element_text(size = 11, face = "bold", angle = 0),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.subtitle    = element_text(size = 9, color = "gray30", face = "italic"),
      panel.spacing    = unit(1, "lines")
    )
}

for (gs in config$gene_sets) {
  subset_data   <- all_results %>% filter(gene_set == gs)
  n_genes_label <- unique(subset_data$n_genes)

  pdf_file_gs <- file.path(OUTPUT_DIR,
                           sprintf("simulation_feng_%s_%dgenes.pdf",
                                   gs, n_genes_label))
  pdf(pdf_file_gs, width = 18, height = 10)

  for (metric in c("ARI", "Silhouette", "Purity")) {
    print(create_metric_boxplot(subset_data, metric,
                                metric, metric_formulas[[metric]],
                                method_colors, n_genes_label))
  }

  # Summary heatmap for this gene set
  heatmap_long <- subset_data %>%
    group_by(method) %>%
    summarise(
      ARI        = mean(ARI,        na.rm = TRUE),
      Silhouette = mean(Silhouette, na.rm = TRUE),
      Purity     = mean(Purity,     na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(ARI, Silhouette, Purity),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("ARI", "Silhouette", "Purity")))

  p_heat <- ggplot(heatmap_long,
                   aes(x = Metric, y = reorder(method, Value), fill = Value)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              fontface = "bold", size = 4) +
    scale_fill_gradient2(low = "#E64B35FF", mid = "white", high = "#4DBBD5FF",
                         midpoint = 0.5, limits = c(-0.1, 1), name = "Score") +
    labs(title    = sprintf("Performance Heatmap (%d genes)", n_genes_label),
         subtitle = "Average across all effect sizes and proportions",
         x = "Metric", y = "Method") +
    theme_bw() +
    theme(axis.text.x     = element_text(angle = 0, hjust = 0.5, size = 10),
          axis.text.y     = element_text(size = 10),
          plot.title      = element_text(size = 14, face = "bold"),
          legend.position = "right")
  print(p_heat)

  dev.off()
  write_log(sprintf("Supplementary plots saved: %s", pdf_file_gs))
}

# ==============================================================================
# FIGURE 3 — MAIN PAPER FIGURE
# ==============================================================================

write_log("")
write_log("GENERATING FIGURE 3...")

fig3_colors <- c(
  RGCCA  = "#E64B35FF", MCIA   = "#4DBBD5FF", intNMF = "#3C5488FF",
  JIVE   = "#F39B7FFF", tICA   = "#8491B4FF", UBMI   = "#91D1C2FF",
  MOFA   = "#00A087FF"
)
fig3_shapes <- c(MCIA = 16, intNMF = 17, JIVE = 15, UBMI = 18,
                 RGCCA = 8, tICA = 10, MOFA = 4)

# ── Panel A: heatmap from random-injection results ────────────────────────────
# (requires RANDOM_INJECTION_CSV to point to results_raw.csv from
#  simulation_random_feature_injection.R; set to NULL to skip this panel)
pA <- NULL

if (!is.null(RANDOM_INJECTION_CSV) && file.exists(RANDOM_INJECTION_CSV)) {
  pA <- tryCatch({
    raw <- read.csv(RANDOM_INJECTION_CSV, stringsAsFactors = FALSE)

    heat_data <- raw %>%
      group_by(method) %>%
      summarise(
        ARI        = mean(ARI,        na.rm = TRUE),
        Silhouette = mean(Silhouette, na.rm = TRUE),
        Purity     = mean(Purity,     na.rm = TRUE),
        .groups = "drop"
      )

    method_order <<- heat_data %>% arrange(desc(ARI)) %>% pull(method)

    heat_long <- heat_data %>%
      pivot_longer(-method, names_to = "Metric", values_to = "Score") %>%
      mutate(
        Metric = factor(Metric, levels = c("ARI", "Silhouette", "Purity")),
        method = factor(method, levels = rev(method_order))
      )

    ggplot(heat_long, aes(x = Metric, y = method, fill = Score)) +
      geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.2f", Score)),
                size = 4.2, fontface = "bold") +
      scale_fill_gradient2(
        low = "#E64B35FF", mid = "white", high = "#4DBBD5FF",
        midpoint = 0.5, limits = c(-0.05, 1), name = "Score") +
      labs(
        title    = "A   Average simulation performance",
        subtitle = "Random feature injection\n(n = 10 replicates per condition)",
        x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(colour = "grey40", size = 9),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right",
        panel.grid      = element_blank()
      )
  }, error = function(e) {
    write_log(paste("WARNING: Could not build Panel A —", e$message))
    NULL
  })
} else {
  write_log("Panel A skipped (RANDOM_INJECTION_CSV not set or file not found)")
}

# ── Panel B: gene-level ARI from this run ────────────────────────────────────

# Use method order from Panel A if available, otherwise rank by mean ARI
if (!exists("method_order")) {
  method_order <- all_results %>%
    group_by(method) %>%
    summarise(ARI = mean(ARI, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ARI)) %>%
    pull(method)
}

feng238_avg <- all_results %>%
  filter(gene_set == "genes_3plus") %>%
  group_by(method) %>%
  summarise(ARI = mean(ARI, na.rm = TRUE), .groups = "drop") %>%
  mutate(gene_set = "Large gene set\n(238 genes, 3+ omics)")

feng54_avg <- all_results %>%
  filter(gene_set == "genes_4plus") %>%
  group_by(method) %>%
  summarise(ARI = mean(ARI, na.rm = TRUE), .groups = "drop") %>%
  mutate(gene_set = "Small gene set\n(54 genes, 4+ omics)")

gene_data <- bind_rows(feng238_avg, feng54_avg) %>%
  mutate(
    method   = factor(method, levels = method_order),
    gene_set = factor(gene_set,
                      levels = c("Large gene set\n(238 genes, 3+ omics)",
                                 "Small gene set\n(54 genes, 4+ omics)"))
  )

y_min <- max(floor(min(gene_data$ARI, na.rm = TRUE) * 10) / 10 - 0.05, -0.15)
y_max <- min(ceiling(max(gene_data$ARI, na.rm = TRUE) * 10) / 10 + 0.05,  0.85)

pB <- ggplot(gene_data,
             aes(x = method, y = ARI, colour = method, shape = method)) +
  geom_hline(yintercept = 0, colour = "grey40",
             linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 4.5, stroke = 1.2) +
  geom_text(aes(label = sprintf("%.2f", ARI),
                vjust = ifelse(ARI >= 0, -0.9, 1.6)),
            size = 3.4, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ gene_set, ncol = 2) +
  scale_colour_manual(values = fig3_colors, guide = "none") +
  scale_shape_manual(values  = fig3_shapes,  guide = "none") +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = c(-0.1, 0, 0.2, 0.4, 0.6, 0.8)) +
  labs(
    title    = "B   Gene-level signal detection (candidate genes)",
    subtitle = "Mean ARI across effect sizes and feature proportions",
    x = NULL,
    y = "Adjusted Rand Index (ARI)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(colour = "grey40", size = 9),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.x        = element_text(size = 10, face = "bold",
                                      angle = 30, hjust = 1),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none"
  )

# ── Save Figure 3 ─────────────────────────────────────────────────────────────

out_pdf <- file.path(FIGURES_DIR, "Fig3_simulation_summary.pdf")
out_png <- file.path(FIGURES_DIR, "Fig3_simulation_summary.png")

if (!is.null(pA)) {
  pdf(out_pdf, width = 13, height = 5.5)
  grid.arrange(pA, pB, ncol = 2, widths = c(1, 1.5))
  dev.off()

  png(out_png, width = 13, height = 5.5, units = "in", res = 300)
  grid.arrange(pA, pB, ncol = 2, widths = c(1, 1.5))
  dev.off()
} else {
  # Save Panel B alone when Panel A data is unavailable
  pdf(out_pdf, width = 8, height = 5.5)
  print(pB)
  dev.off()
  png(out_png, width = 8, height = 5.5, units = "in", res = 300)
  print(pB)
  dev.off()
  write_log("Note: Fig3 saved with Panel B only (Panel A data unavailable)")
}

write_log(sprintf("Fig3 saved:\n  %s\n  %s", out_pdf, out_png))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

write_log("")
write_log("========================================")
write_log("COMPLETE")
write_log("========================================")
write_log(sprintf("Runtime: %.1f minutes", total_time))
write_log("")
write_log("METHOD RANKING BY GENE SET (by mean ARI):")

for (gs in config$gene_sets) {
  write_log(sprintf("\n%s:", gs))
  ranking <- all_results %>%
    filter(gene_set == gs) %>%
    group_by(method) %>%
    summarise(
      mean_ARI = mean(ARI, na.rm = TRUE),
      sd_ARI   = sd(ARI,   na.rm = TRUE),
      n_genes  = unique(n_genes),
      .groups  = "drop"
    ) %>%
    arrange(desc(mean_ARI))

  for (j in 1:nrow(ranking)) {
    write_log(sprintf("  %d. %s: ARI=%.4f (SD=%.4f) [%d genes]",
                      j, ranking$method[j], ranking$mean_ARI[j],
                      ranking$sd_ARI[j], ranking$n_genes[j]))
  }
}

cat("\n========================================\n")
cat("SIMULATION COMPLETE\n")
cat(sprintf("Output:  %s\n", OUTPUT_DIR))
cat(sprintf("Fig3:    %s\n", out_png))
cat(sprintf("Runtime: %.1f minutes\n", total_time))
cat("========================================\n")
