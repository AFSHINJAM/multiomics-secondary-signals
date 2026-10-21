###############################################################################
####   MULTI-OMICS INTEGRATION SIMULATION STUDY                             ####
####   After Batch Correction Only                                           ####
####   Metrics: ARI, Silhouette, Purity                                     ####
###############################################################################

rm(list = ls())
set.seed(42)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(sva)
library(cluster)
library(mclust)
library(fpc)
library(pheatmap)

# GAUDI method dependencies
library(RGCCA)
library(omicade4)    # for MCIA
library(IntNMF)
library(r.jive)      # for JIVE

# MOFA - use basilisk (no manual Python setup needed)
library(reticulate)
library(MOFA2)

cat("MOFA will use basilisk - no manual Python setup needed\n\n")

# ==============================================================================
# CONFIGURATION — set these paths before running
# ==============================================================================

RUN_ID   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
RUN_NAME <- "ancestry_simulation_after_only"

# Path to workspace.RData produced by 01_data_preprocessing.R
WORKSPACE_PATH <- "/path/to/data/workspace.RData"

# Path to runfactorization.R (GAUDI/UBMI integration script)
GAUDI_SCRIPT   <- "/path/to/scripts/runfactorization.R"

# Directory where all results and plots will be saved
OUTPUT_DIR <- file.path("/path/to/results/simulations",
                        paste0(RUN_ID, "_", RUN_NAME))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Simulation parameters
config <- list(
  run_id   = RUN_ID,
  run_name = RUN_NAME,

  # Simulation parameters
  n_simulations    = 10,
  prop_affected    = c(0.05, 0.10),
  effect_sizes     = c(0.4, 0.8),
  n_omics_affected = c(2, 6),

  # Methods passed to runfactorization()
  methods = c("RGCCA", "MCIA", "intNMF", "JIVE", "tICA", "UBMI"),
  run_mofa_separately = TRUE,

  # Paths (stored for reference in log)
  output_dir     = OUTPUT_DIR,
  workspace_path = WORKSPACE_PATH,
  gaudi_script   = GAUDI_SCRIPT
)

cat("========================================\n")
cat("SIMULATION RUN:", RUN_ID, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("========================================\n\n")

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

write_log("Starting simulation (After Batch Correction Only)...")

# ==============================================================================
# LOAD DATA
# ==============================================================================

write_log("Loading workspace...")
load(config$workspace_path)

write_log("Sourcing GAUDI functions...")
setwd(dirname(GAUDI_SCRIPT))
source(basename(GAUDI_SCRIPT))
setwd(dirname(WORKSPACE_PATH))

# Prepare omics data
omics_data_original <- list(
  "RNA-seq"   = rna_final,
  "ATAC-seq"  = atac_final,
  "H3K27ac"   = h3k27ac_final,
  "H3K4me3"   = h3k4me3_final,
  "H3K4me1"   = h3k4me1_final,
  "H3K27me3"  = h3k27me3_final
)

write_log("Omics data loaded:")
for (name in names(omics_data_original)) {
  write_log(sprintf("  %s: %d x %d", name,
                    nrow(omics_data_original[[name]]),
                    ncol(omics_data_original[[name]])))
}

# Standardise ancestry labels
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

# Signal injection using mean-centered ancestry encoding scaled by feature SD.
# For each affected feature j in sample i:
#
#   x_tilde_ij = x_ij + ES * (a_i - a_bar) * sigma_j + epsilon_ij
#
# where ES = effect size (equals Cohen's d because signal is scaled by sigma_j),
#       a_i in {0,1} is ancestry indicator (1=AF, 0=EU),
#       a_bar = mean(a_i)  →  centered_ancestry ∈ {+0.5, -0.5},
#       sigma_j = SD of feature j,
#       epsilon_ij ~ N(0, 0.1 * sigma_j)  [10% noise relative to feature SD]
inject_ancestry_signal <- function(omics_data, ancestry_encoding,
                                   prop_affected, effect_size) {
  n_samples  <- nrow(omics_data)
  n_features <- ncol(omics_data)
  n_affected <- round(n_features * prop_affected)

  affected_features <- sample(1:n_features, n_affected, replace = FALSE)
  perturbed_data    <- omics_data
  feature_sds       <- apply(omics_data, 2, sd)
  centered_ancestry <- ancestry_encoding - mean(ancestry_encoding, na.rm = TRUE)

  for (feat_idx in affected_features) {
    feat_sd <- feature_sds[feat_idx]
    if (is.na(feat_sd) || feat_sd == 0) feat_sd <- 1
    noise       <- rnorm(n_samples, mean = 0, sd = 0.1 * feat_sd)
    perturbation <- effect_size * centered_ancestry * feat_sd + noise
    perturbed_data[, feat_idx] <- omics_data[, feat_idx] + perturbation
  }

  list(data = perturbed_data, affected_features = affected_features)
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
  if (is.null(factors)) return(list(ARI = NA, Silhouette = NA, Purity = NA))

  if (is.list(factors) && !is.data.frame(factors)) {
    for (nm in c("factors", "score", "composite", "global", "global.factors", "joint")) {
      if (nm %in% names(factors)) { factors <- factors[[nm]]; break }
    }
    if (is.list(factors) && length(factors) > 0) factors <- factors[[1]]
  }

  if (is.data.frame(factors)) factors <- as.matrix(factors)
  if (!is.matrix(factors) || nrow(factors) < 2 || ncol(factors) < 1) {
    return(list(ARI = NA, Silhouette = NA, Purity = NA))
  }

  n_factors   <- min(5, ncol(factors))
  factor_data <- factors[, 1:n_factors, drop = FALSE]
  factor_data[is.na(factor_data)]      <- 0
  factor_data[is.infinite(factor_data)] <- 0

  valid_idx <- !is.na(true_labels)
  if (sum(valid_idx) < 10) return(list(ARI = NA, Silhouette = NA, Purity = NA))

  factor_data <- factor_data[valid_idx, , drop = FALSE]
  true_labels <- true_labels[valid_idx]

  k <- length(unique(true_labels))
  if (k < 2) k <- 2

  set.seed(123)
  km <- tryCatch(
    kmeans(factor_data, centers = k, nstart = 25, iter.max = 1000),
    error = function(e) NULL
  )
  if (is.null(km)) return(list(ARI = NA, Silhouette = NA, Purity = NA))

  pred    <- km$cluster
  metrics <- list(ARI = NA, Silhouette = NA, Purity = NA)

  tryCatch({ metrics$ARI <- adjustedRandIndex(pred, true_labels) }, error = function(e) {})
  tryCatch({
    metrics$Silhouette <- mean(silhouette(pred, dist(factor_data))[, 3], na.rm = TRUE)
  }, error = function(e) {})
  tryCatch({
    cm             <- table(pred, true_labels)
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

    mofa_obj <- create_mofa(mofa_data)

    model_opts               <- get_default_model_options(mofa_obj)
    model_opts$num_factors   <- num_factors

    train_opts                    <- get_default_training_options(mofa_obj)
    train_opts$convergence_mode   <- "fast"
    train_opts$verbose            <- FALSE
    train_opts$seed               <- 42

    mofa_obj <- prepare_mofa(mofa_obj,
                             model_options    = model_opts,
                             training_options = train_opts)
    mofa_obj <- MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE, outfile = temp_hdf5)

    get_factors(mofa_obj)[[1]]

  }, error = function(e) {
    write_log(paste("MOFA error:", e$message))
    NULL
  })

  if (file.exists(temp_hdf5)) file.remove(temp_hdf5)
  result
}

run_integration_methods <- function(data_list, metadata, config, temp_base, label) {

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

run_scenario <- function(params, omics_original, metadata, config) {

  set.seed(42 + params$row_idx * 1000 + params$sim_id)

  ancestry_encoding <- get_ancestry_encoding(metadata)

  # Select which omics layers to perturb
  if (params$n_omics_affected == 6) {
    affected_omics <- names(omics_original)
  } else {
    affected_omics <- sample(names(omics_original),
                             params$n_omics_affected, replace = FALSE)
  }

  write_log(sprintf("  Affecting omics: %s", paste(affected_omics, collapse = ", ")))

  # Inject ancestry signal into selected omics
  perturbed_list <- list()
  for (omics_name in names(omics_original)) {
    if (omics_name %in% affected_omics) {
      result <- inject_ancestry_signal(
        omics_original[[omics_name]], ancestry_encoding,
        params$prop_affected, params$effect_size
      )
      perturbed_list[[omics_name]] <- result$data
      write_log(sprintf("    %s: %d features affected",
                        omics_name, length(result$affected_features)))
    } else {
      perturbed_list[[omics_name]] <- omics_original[[omics_name]]
    }
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

  results$prop_affected    <- params$prop_affected
  results$effect_size      <- params$effect_size
  results$n_omics_affected <- params$n_omics_affected
  results$sim_id           <- params$sim_id

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
  prop_affected    = config$prop_affected,
  effect_size      = config$effect_sizes,
  n_omics_affected = config$n_omics_affected,
  sim_id           = 1:config$n_simulations,
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

  msg <- sprintf("[%d/%d] Effect=%.2f, Prop=%.0f%%, N_omics=%d, Sim=%d",
                 i, n_total, params$effect_size, params$prop_affected * 100,
                 params$n_omics_affected, params$sim_id)
  write_log(msg)

  result <- tryCatch({
    run_scenario(params, omics_data_original, sample_metadata, config)
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
  group_by(method, effect_size, prop_affected, n_omics_affected) %>%
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
# PLOTS — ARI, Silhouette, Purity
# ==============================================================================

write_log("GENERATING PLOTS...")

method_colors <- c(
  "RGCCA"  = "#E64B35FF",
  "MCIA"   = "#4DBBD5FF",
  "intNMF" = "#3C5488FF",
  "JIVE"   = "#F39B7FFF",
  "tICA"   = "#8491B4FF",
  "UBMI"   = "#91D1C2FF",
  "MOFA"   = "#00A087FF"
)

pdf_file <- file.path(OUTPUT_DIR, "plots.pdf")
pdf(pdf_file, width = 16, height = 10)

# Helper: boxplot with median labels
create_metric_boxplot <- function(data, metric_col, metric_name, colors) {

  medians <- data %>%
    group_by(method, effect_size, prop_affected, n_omics_affected) %>%
    summarise(median_val = median(.data[[metric_col]], na.rm = TRUE), .groups = "drop")

  data$param_label    <- sprintf("Effect=%.1f\nProp=%.0f%%\nOmics=%d",
                                  data$effect_size,
                                  data$prop_affected * 100,
                                  data$n_omics_affected)
  medians$param_label <- sprintf("Effect=%.1f\nProp=%.0f%%\nOmics=%d",
                                  medians$effect_size,
                                  medians$prop_affected * 100,
                                  medians$n_omics_affected)

  ggplot(data, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    geom_text(data = medians,
              aes(x = method, y = median_val,
                  label = sprintf("%.2f", median_val)),
              vjust = -0.8, size = 3, fontface = "bold", color = "black") +
    facet_wrap(~ param_label, ncol = 4) +
    scale_fill_manual(values = colors) +
    labs(title    = metric_name,
         subtitle = "After Batch Correction - Parameter Combinations",
         x = "", y = metric_name) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          legend.position  = "none",
          strip.text       = element_text(size = 11, face = "bold"),
          plot.title       = element_text(size = 14, face = "bold"),
          plot.subtitle    = element_text(size = 10, color = "gray40"))
}

print(create_metric_boxplot(all_results, "ARI",        "ARI (Adjusted Rand Index)", method_colors))
print(create_metric_boxplot(all_results, "Silhouette", "Silhouette Score",          method_colors))
print(create_metric_boxplot(all_results, "Purity",     "Purity Score",              method_colors))

# Combined heatmap
heatmap_long <- all_results %>%
  group_by(method) %>%
  summarise(
    ARI        = mean(ARI,        na.rm = TRUE),
    Silhouette = mean(Silhouette, na.rm = TRUE),
    Purity     = mean(Purity,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(ARI, Silhouette, Purity),
               names_to = "Metric", values_to = "Value")

heatmap_long$Metric <- factor(heatmap_long$Metric,
                               levels = c("ARI", "Silhouette", "Purity"))

p_heatmap <- ggplot(heatmap_long,
                    aes(x = Metric, y = reorder(method, Value), fill = Value)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Value)), fontface = "bold", size = 4) +
  scale_fill_gradient2(low = "#E64B35FF", mid = "white", high = "#4DBBD5FF",
                       midpoint = 0.5, limits = c(-0.1, 1), name = "Score") +
  labs(title    = "Performance Heatmap After Batch Correction",
       subtitle = "ARI, Silhouette, and Purity across integration methods",
       x = "Metric", y = "Method") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title  = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "right")
print(p_heatmap)

dev.off()
write_log(sprintf("Plots: %s", pdf_file))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

write_log("")
write_log("========================================")
write_log("COMPLETE")
write_log("========================================")
write_log(sprintf("Runtime: %.1f minutes", total_time))

ranking <- all_results %>%
  group_by(method) %>%
  summarise(
    mean_ARI = mean(ARI, na.rm = TRUE),
    sd_ARI   = sd(ARI,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_ARI))

write_log("")
write_log("METHOD RANKING (by mean ARI):")
for (i in 1:nrow(ranking)) {
  write_log(sprintf("  %d. %s: ARI=%.4f (SD=%.4f)",
                    i, ranking$method[i], ranking$mean_ARI[i], ranking$sd_ARI[i]))
}

write_log(sprintf("\n  BEST METHOD: %s (ARI=%.4f)",
                  ranking$method[1], ranking$mean_ARI[1]))

cat("\n========================================\n")
cat("SIMULATION COMPLETE\n")
cat(sprintf("Output: %s\n", OUTPUT_DIR))
cat(sprintf("Runtime: %.1f minutes\n", total_time))
cat("========================================\n")
