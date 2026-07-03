# =============================================================================
# STEP_03_GSE_S2_interpretability.R
#
# Purpose:
#   Establish an interpretability bridge between GSE embedding dimensions
#   and Sentinel-2 spectral bands by computing pairwise Spearman correlations
#   between the top-ranked GSE dimensions and Sentinel-2 band-percentile
#   combinations, stratified by seagrass canopy morphology class.
#
#   Outputs correspond to Figure 6 (GSE-Sentinel-2 linkage heatmaps) in
#   the manuscript.
#
# Requires:
#   - Output CSVs from STEP_02: GSE_sensitivity_*.csv
#
# Study:
#   Interpretable Estimation of Seagrass Percent Cover Across Extensive
#   Coastal Environments in Indonesia Using Google Satellite Embeddings
#   International Journal of Digital Earth (2026)
#   Manuscript ID: TJDE-2026-0306
#
# Authors: Muhammad Hafizt, Stuart Phinn, Pramaditya Wicaksono,
#          Udhi Eko Hernawan, Huwaida Nur Salsabila,
#          Mitchell Lyons, Kathryn McMahon, Chris Roelfsema
#
# NOTE:
#   Field training data are not included in this repository.
#   Update the USER CONFIGURATION section before running.
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(tidyr)
library(stringr)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input dataset (field SPC + Sentinel-2 + GSE values)
data_path <- "path/to/your/input_data.csv"

# Paths to GSE sensitivity results from STEP_02
gse_mono_path        <- "path/to/your/output_directory/GSE_sensitivity_mono_species.csv"
gse_mixed_short_path <- "path/to/your/output_directory/GSE_sensitivity_mixed_short.csv"
gse_mixed_long_path  <- "path/to/your/output_directory/GSE_sensitivity_mixed_long.csv"

# Output directory
out_dir <- "path/to/your/output_directory"

# Number of top GSE dimensions to include in linkage analysis
TOP_K <- 5

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df              <- read_csv(data_path, show_col_types = FALSE)
gse_mono        <- read_csv(gse_mono_path, show_col_types = FALSE)
gse_mixed_short <- read_csv(gse_mixed_short_path, show_col_types = FALSE)
gse_mixed_long  <- read_csv(gse_mixed_long_path, show_col_types = FALSE)

cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. DEFINE VARIABLES
# =============================================================================

# Sentinel-2 band-percentile combinations
s2_bands <- c(
  "B2_p0",  "B2_p20",  "B2_p40",  "B2_p60",  "B2_p80",  "B2_p100",
  "B3_p0",  "B3_p20",  "B3_p40",  "B3_p60",  "B3_p80",  "B3_p100",
  "B4_p0",  "B4_p20",  "B4_p40",  "B4_p60",  "B4_p80",  "B4_p100",
  "B8_p0",  "B8_p20",  "B8_p40",  "B8_p60",  "B8_p80",  "B8_p100"
)

# =============================================================================
# 3. GSE-S2 LINKAGE FUNCTION
# =============================================================================

run_gse_s2_link <- function(df_sub, gse_rank, morph_label) {

  # Select top-K GSE dimensions by combined sensitivity rank
  top_gse <- gse_rank %>%
    arrange(rank_mean) %>%
    slice(1:TOP_K) %>%
    pull(band)

  cat("\nMorphology:", morph_label, "\n")
  cat("Top GSE dimensions:", paste(top_gse, collapse = ", "), "\n")

  # Compute pairwise Spearman correlations: GSE x S2
  cor_mat <- map_dfr(
    top_gse,
    function(gse_band) {
      map_dfr(
        s2_bands,
        function(s2_band) {
          ct <- cor.test(
            df_sub[[gse_band]],
            df_sub[[s2_band]],
            method = "spearman",
            exact  = FALSE
          )
          tibble(
            morphology = morph_label,
            GSE_band   = gse_band,
            S2_band    = s2_band,
            rho        = ct$estimate
          )
        }
      )
    }
  )

  # Heatmap
  p <- ggplot(cor_mat, aes(x = S2_band, y = GSE_band, fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      name     = "Spearman rho"
    ) +
    labs(
      title = paste("GSE to Sentinel-2 linkage:", morph_label),
      x     = "Sentinel-2 bands",
      y     = "GSE embedding dimensions"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)

  # Save correlation table
  write_csv(
    cor_mat,
    file.path(out_dir, paste0("GSE_S2_link_", morph_label, ".csv"))
  )

  return(cor_mat)
}

# =============================================================================
# 4. RUN PER MORPHOLOGY CLASS
# =============================================================================

# Mono-species
df_mono <- df %>% filter(sg_morpho == "mono")

gse_mono_rank <- gse_mono %>%
  group_by(band) %>%
  summarise(rank_mean = mean(rank_mean, na.rm = TRUE), .groups = "drop")

link_mono <- run_gse_s2_link(df_mono, gse_mono_rank, "mono")

# Mixed short-leaved
df_mixed_short <- df %>% filter(sg_morpho == "mixed_short")

gse_mixed_short_rank <- gse_mixed_short %>%
  select(band, rank_mean)

link_mixed_short <- run_gse_s2_link(df_mixed_short, gse_mixed_short_rank, "mixed_short")

# Mixed long-leaved
df_mixed_long <- df %>% filter(sg_morpho == "mixed_long")

gse_mixed_long_rank <- gse_mixed_long %>%
  select(band, rank_mean)

link_mixed_long <- run_gse_s2_link(df_mixed_long, gse_mixed_long_rank, "mixed_long")

cat("STEP-03 completed: GSE to Sentinel-2 interpretability heatmaps saved\n")
