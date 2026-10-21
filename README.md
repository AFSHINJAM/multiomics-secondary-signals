# multiomics-secondary-signals

Benchmarking multi-omics integration methods (GAUDI/UBMI, MCIA, JIVE, RGCCA, intNMF, MOFA, tICA) for detecting secondary signals (ancestry) after removing a dominant infection response.

Data: Aracena et al. 2024, *Nature Genetics* — RNA-seq, ATAC-seq and ChIP-seq (H3K27ac, H3K4me3, H3K4me1, H3K27me3) from non-infected (NI) and influenza-infected (Flu) samples of African (AF) and European (EU) ancestry.

## Scripts

Run in order. Each script has a **CONFIGURATION** block at the top — set the paths there before running.

| Script | Description |
|---|---|
| `scripts/01_data_preprocessing.R` | Normalizes the six omics layers (TMM + log-CPM), keeps the top 1000 variable features per layer, writes GAUDI-format matrices, runs the factorization methods and saves `workspace.RData`. |
| `scripts/02_batch_correction_visualization.R` | Removes the time-point (infection) signal with ComBat and compares clusters before (time point) and after (ancestry) correction; produces summary figures and a clustering metrics table. |
| `scripts/03_random_feature_injection_simulation.R` | Simulation injecting ancestry signal into randomly selected features across varying effect sizes, feature proportions and numbers of affected omics. |
| `scripts/04_gene_level_simulation.R` | Simulation injecting ancestry signal into candidate gene sets (238 genes in ≥3 omics, 54 genes in ≥4 omics); generates Figure 3. |

Clustering performance is evaluated with ARI, silhouette score and cluster purity.

## Requirements

R with: `edgeR`, `sva`, `rtracklayer`, `GenomicRanges`, `IRanges`, `tidyverse`, `ggplot2`, `gridExtra`, `ggsci`, `ggforce`, `cluster`, `mclust`, `fpc`, `pheatmap`, `future.apply`, `RGCCA`, `omicade4`, `IntNMF`, `r.jive`, `MOFA2`, `reticulate`.

The factorization methods are run through `runfactorization.R` from GAUDI (set its location via `GAUDI_SCRIPT` / `SCRIPTS_DIR`).

## License

MIT
