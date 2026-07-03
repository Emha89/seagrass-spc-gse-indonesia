# =============================================================================
# STEP_01_S2_spectral_sensitivity.R
#
# Purpose:
#   Assess Sentinel-2 spectral band sensitivity to seagrass percent cover (SPC)
#   across three canopy morphology classes:
#     (1) Mono-species (per species SPC)
#     (2) Mixed short-leaved (total SPC)
#     (3) Mixed long-leaved (total SPC)
#
#   Sensitivity is quantified using:
#     - Spearman rank correlation (Spearman 1904)
#     - Generalised Additive Models (GAM; Wood 2017)
#     - Combined ranking metric (mean of Spearman and GAM ranks)
#
#   Outputs correspond to Figure 4 (Sentinel-2 spectral sensitivity) and
#   Supplementary Figure S1 (GAM smooth plots) in the manuscript.
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
#   Data are available from the corresponding author upon reasonable request,
#   subject to BRIN data governance requirements (see Data Availability Statement).
#   Update the USER CONFIGURATION section before running.
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(mgcv)
library(ggplot2)
library(forcats)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input dataset (field SPC + Sentinel-2 band values)
# Column structure: loc, year, sg_morpho, total_SPC, species SPC columns,
#                  S2 band percentiles (e.g. B2_p0 to B8_p100)
data_path <- "path/to/your/input_data.csv"

# Output directory for CSV and figure files
out_dir <- "path/to/your/output_directory"

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df <- read_csv(data_path, show_col_types = FALSE)
cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. DEFINE VARIABLES
# =============================================================================

# Species SPC columns (used for mono-species analysis)
species_spc <- c(
  "Ea_SPC", "Th_SPC", "Cr_SPC", "Cs_SPC", "Si_SPC",
  "Hu_SPC", "Ho_SPC", "Hp_SPC", "Tc_SPC",
  "Hm_SPC", "Hs_SPC", "Hd_SPC"
)

# Sentinel-2 band-percentile combinations
s2_bands <- c(
  "B2_p0",  "B2_p20",  "B2_p40",  "B2_p60",  "B2_p80",  "B2_p100",
  "B3_p0",  "B3_p20",  "B3_p40",  "B3_p60",  "B3_p80",  "B3_p100",
  "B4_p0",  "B4_p20",  "B4_p40",  "B4_p60",  "B4_p80",  "B4_p100",
  "B8_p0",  "B8_p20",  "B8_p40",  "B8_p60",  "B8_p80",  "B8_p100"
)

# Logical band ordering for figures (band x percentile)
band_order <- expand.grid(
  band = c("B2", "B3", "B4", "B8"),
  perc = c("p0", "p20", "p40", "p60", "p80", "p100")
) %>%
  mutate(band_full = paste0(band, "_", perc)) %>%
  pull(band_full)

# Top N bands to retain per morphology class for figure
TOP_N <- 10

# =============================================================================
# 3. SENSITIVITY ANALYSIS FUNCTION
# =============================================================================

run_sensitivity <- function(data, response_var, label) {

  cat("\n============================================\n")
  cat("Analysis:", label, "| Response:", response_var, "\n")
  cat("Number of samples:", nrow(data), "\n")
  cat("============================================\n")

  if (nrow(data) < 30) {
    message("Skipped: insufficient samples (n < 30)")
    return(NULL)
  }

  # Spearman rank correlation
  spearman <- map_dfr(
    s2_bands,
    function(band) {
      ct <- cor.test(
        data[[band]],
        data[[response_var]],
        method = "spearman",
        exact = FALSE
      )
      tibble(
        group        = label,
        response     = response_var,
        band         = band,
        spearman_rho = ct$estimate,
        spearman_p   = ct$p.value
      )
    }
  ) %>%
    arrange(desc(abs(spearman_rho)))

  cat("\nTop 5 bands by |Spearman rho|:\n")
  print(head(spearman, 5))

  # GAM: deviance explained
  gam_res <- map_dfr(
    s2_bands,
    function(band) {
      m <- gam(
        as.formula(paste(response_var, "~ s(", band, ", k = 5)", sep = "")),
        data   = data,
        method = "REML"
      )
      s <- summary(m)
      tibble(
        group               = label,
        response            = response_var,
        band                = band,
        deviance_explained  = s$dev.expl,
        edf                 = s$s.table[1, "edf"],
        gam_p               = s$s.table[1, "p-value"]
      )
    }
  ) %>%
    arrange(desc(deviance_explained))

  cat("\nTop 5 bands by GAM deviance explained:\n")
  print(head(gam_res, 5))

  # Combined ranking (mean of Spearman rank and GAM rank)
  final_res <- spearman %>%
    left_join(gam_res, by = c("group", "response", "band")) %>%
    mutate(
      rank_spearman = rank(-abs(spearman_rho)),
      rank_gam      = rank(-deviance_explained),
      rank_mean     = (rank_spearman + rank_gam) / 2
    ) %>%
    arrange(rank_mean)

  cat("\nTop 5 combined ranking:\n")
  print(head(final_res, 5))

  return(final_res)
}

# =============================================================================
# 4. MONO-SPECIES ANALYSIS
# =============================================================================

df_mono <- df %>% filter(sg_morpho == "mono")

mono_results <- map(
  species_spc,
  function(sp) {
    df_sub <- df_mono %>%
      filter(.data[[sp]] > 0, .data[[sp]] <= 100)
    run_sensitivity(data = df_sub, response_var = sp, label = "mono_species")
  }
) %>%
  compact() %>%
  bind_rows()

# =============================================================================
# 5. MIXED SHORT-LEAVED ANALYSIS
# =============================================================================

df_mixed_short <- df %>%
  filter(sg_morpho == "mixed_short", total_SPC > 0, total_SPC <= 100)

mixed_short_results <- run_sensitivity(
  data         = df_mixed_short,
  response_var = "total_SPC",
  label        = "mixed_short"
)

# =============================================================================
# 6. MIXED LONG-LEAVED ANALYSIS
# =============================================================================

df_mixed_long <- df %>%
  filter(sg_morpho == "mixed_long", total_SPC > 0, total_SPC <= 100)

mixed_long_results <- run_sensitivity(
  data         = df_mixed_long,
  response_var = "total_SPC",
  label        = "mixed_long"
)

# =============================================================================
# 7. SAVE SENSITIVITY RESULTS
# =============================================================================

write_csv(mono_results,        file.path(out_dir, "S2_sensitivity_mono_species.csv"))
write_csv(mixed_short_results, file.path(out_dir, "S2_sensitivity_mixed_short.csv"))
write_csv(mixed_long_results,  file.path(out_dir, "S2_sensitivity_mixed_long.csv"))

cat("STEP-01 completed: S2 sensitivity results saved\n")

# =============================================================================
# 8. FIGURE 4 — Sentinel-2 spectral sensitivity to seagrass percent cover
# =============================================================================

# Aggregate mono-species results across species
fig_mono <- mono_results %>%
  group_by(band) %>%
  summarise(rank_mean = mean(rank_mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = "Mono-species")

fig_mixed_short <- mixed_short_results %>%
  select(band, rank_mean) %>%
  mutate(group = "Mixed short")

fig_mixed_long <- mixed_long_results %>%
  select(band, rank_mean) %>%
  mutate(group = "Mixed long")

fig_data <- bind_rows(fig_mono, fig_mixed_short, fig_mixed_long) %>%
  group_by(group) %>%
  slice_min(rank_mean, n = TOP_N) %>%
  ungroup() %>%
  mutate(band = factor(band, levels = band_order))

fig_plot <- ggplot(
  fig_data,
  aes(x = rank_mean, y = fct_reorder(band, rank_mean))
) +
  geom_point(size = 2.5, colour = "black") +
  facet_wrap(~ group, scales = "free_y") +
  scale_x_continuous(limits = c(0, 25)) +
  labs(
    x = "Combined sensitivity rank (lower = more influential)",
    y = "Sentinel-2 band"
  ) +
  theme_minimal() +
  theme(
    strip.text    = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x  = element_text(margin = margin(t = 10))
  )

print(fig_plot)

ggsave(
  filename = file.path(out_dir, "Figure_4_S2_spectral_sensitivity.png"),
  plot     = fig_plot,
  width    = 12, height = 5, dpi = 300
)

# Save top-10 bands per morphology to CSV
fig_export <- fig_data %>%
  mutate(band = as.character(band)) %>%
  arrange(group, rank_mean)

write_csv(fig_export, file.path(out_dir, "S2_top10_spectral_bands_by_morphology.csv"))

cat("Figure 4 saved\n")
