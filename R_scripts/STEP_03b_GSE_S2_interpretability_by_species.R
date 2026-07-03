# =============================================================================
# STEP_03b_GSE_S2_interpretability_by_species.R
#
# Purpose:
#   Extend the GSE-Sentinel-2 interpretability analysis to the species level.
#   Pairwise Spearman correlations between top-ranked GSE dimensions and
#   Sentinel-2 band-percentile combinations are computed for selected
#   mono-species under high-cover conditions (SPC >= 80%).
#
#   Only species with sufficient sample size (n >= 30) are included.
#   In the study, only Thalassia hemprichii met this threshold.
#
#   Outputs correspond to Figure 7 (species-level GSE-S2 linkage) in
#   the manuscript.
#
# Requires:
#   - Output CSV from STEP_02: GSE_sensitivity_mono_species.csv
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
library(stringr)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input dataset
data_path <- "path/to/your/input_data.csv"

# Path to GSE sensitivity results from STEP_02 (mono-species)
gse_rank_path <- "path/to/your/output_directory/GSE_sensitivity_mono_species.csv"

# Output directory
out_dir <- "path/to/your/output_directory"

# Minimum sample size for species-level analysis
MIN_SAMPLES <- 30

# Number of top GSE dimensions to include
TOP_K <- 5

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df       <- read_csv(data_path, show_col_types = FALSE)
gse_rank <- read_csv(gse_rank_path, show_col_types = FALSE)

cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. DEFINE VARIABLES
# =============================================================================

# Target species for species-level analysis
# Restricted to species with sufficient high-cover mono-species samples
target_species <- c("Ea_SPC", "Th_SPC", "Cr_SPC")

species_labels <- c(
  Ea_SPC = "Enhalus acoroides",
  Th_SPC = "Thalassia hemprichii",
  Cr_SPC = "Cymodocea rotundata"
)

# Sentinel-2 band-percentile combinations
s2_bands <- c(
  "B2_p0",  "B2_p20",  "B2_p40",  "B2_p60",  "B2_p80",  "B2_p100",
  "B3_p0",  "B3_p20",  "B3_p40",  "B3_p60",  "B3_p80",  "B3_p100",
  "B4_p0",  "B4_p20",  "B4_p40",  "B4_p60",  "B4_p80",  "B4_p100",
  "B8_p0",  "B8_p20",  "B8_p40",  "B8_p60",  "B8_p80",  "B8_p100"
)

# =============================================================================
# 3. SPECIES-LEVEL LINKAGE FUNCTION
# =============================================================================

run_species_link <- function(sp) {

  cat("\n============================================\n")
  cat("Species:", species_labels[sp], "\n")

  # Filter to high-cover mono-species quadrats for this species
  df_sp <- df %>%
    filter(
      sg_morpho == "mono",
      .data[[sp]] >= 80
    )

  cat("Samples (SPC >= 80%):", nrow(df_sp), "\n")

  if (nrow(df_sp) < MIN_SAMPLES) {
    message(paste("Skipped: insufficient samples (n <", MIN_SAMPLES, ")"))
    return(NULL)
  }

  # Top GSE dimensions for this species from STEP_02 sensitivity results
  top_gse <- gse_rank %>%
    filter(response == sp) %>%
    arrange(rank_mean) %>%
    slice(1:TOP_K) %>%
    pull(band)

  cat("Top GSE dimensions:", paste(top_gse, collapse = ", "), "\n")

  # Pairwise Spearman correlations: GSE x S2
  cor_df <- map_dfr(
    top_gse,
    function(gse_band) {
      map_dfr(
        s2_bands,
        function(s2_band) {
          ct <- cor.test(
            df_sp[[gse_band]],
            df_sp[[s2_band]],
            method = "spearman",
            exact  = FALSE
          )
          tibble(
            species  = species_labels[sp],
            GSE_band = gse_band,
            S2_band  = s2_band,
            rho      = ct$estimate
          )
        }
      )
    }
  )

  # Heatmap
  p <- ggplot(cor_df, aes(x = S2_band, y = GSE_band, fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      name     = "Spearman rho"
    ) +
    labs(
      title    = paste("GSE to Sentinel-2 linkage:", species_labels[sp]),
      subtitle = "Mono-species, SPC >= 80%",
      x        = "Sentinel-2 bands",
      y        = "GSE embedding dimensions"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p)

  write_csv(
    cor_df,
    file.path(out_dir, paste0("GSE_S2_link_", sp, "_high_cover.csv"))
  )

  return(cor_df)
}

# =============================================================================
# 4. RUN FOR TARGET SPECIES
# =============================================================================

species_links <- map(target_species, run_species_link) %>%
  compact() %>%
  bind_rows()

cat("STEP-03b completed: species-level GSE to Sentinel-2 linkage saved\n")
