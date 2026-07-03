# =============================================================================
# STEP_04_S2_spectral_profiles.R
#
# Purpose:
#   Generate Sentinel-2 multispectral reflectance profiles for mono-species
#   seagrass under high percent cover conditions (SPC >= 80%), stratified
#   by water depth class. Profiles summarise the median and interquartile
#   range of Sentinel-2 p60 surface reflectance per species and depth class.
#
#   Outputs correspond to Figure 3 (spectral profiles) in the manuscript.
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
library(tidyr)
library(ggplot2)
library(stringr)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input dataset (field SPC + Sentinel-2 + depth values)
# Required columns: sg_morpho, species SPC columns, S2 band columns, depth
data_path <- "path/to/your/input_data.csv"

# Output directory
out_dir <- "path/to/your/output_directory"

# Sentinel-2 percentile to use for spectral profiles
PERC_TARGET <- 60

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df <- read_csv(data_path, show_col_types = FALSE)
cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. DEFINE VARIABLES
# =============================================================================

# Validate percentile
stopifnot(PERC_TARGET %in% c(0, 20, 40, 60, 80, 100))

# Species SPC columns (mono-species only)
species_spc <- c(
  "Ea_SPC", "Th_SPC", "Cr_SPC", "Cs_SPC", "Si_SPC",
  "Hu_SPC", "Ho_SPC", "Hp_SPC", "Tc_SPC",
  "Hm_SPC", "Hs_SPC", "Hd_SPC"
)

# Sentinel-2 bands at the target percentile
s2_bands <- paste0(c("B2", "B3", "B4", "B8"), "_p", PERC_TARGET)

# Display labels for band axis
band_labels <- setNames(
  c("Blue (B2)", "Green (B3)", "Red (B4)", "NIR (B8)"),
  s2_bands
)

# =============================================================================
# 3. FILTER: MONO-SPECIES HIGH COVER
# =============================================================================

df_mono <- df %>% filter(sg_morpho == "mono")

# =============================================================================
# 4. DEPTH STRATIFICATION
# =============================================================================

# Depth derived from Allen Coral Atlas bathymetric product
# Reported RMSE: up to 1.9 m (Li et al. 2021)
df_mono <- df_mono %>%
  mutate(
    depth_m = abs(depth),
    depth_class = case_when(
      depth_m <= 2  ~ "0-2 m",
      depth_m <= 5  ~ "2-5 m",
      depth_m <= 10 ~ "5-10 m",
      TRUE          ~ ">10 m"
    ),
    depth_class = factor(depth_class, levels = c("0-2 m", "2-5 m", "5-10 m", ">10 m"))
  )

# =============================================================================
# 5. BUILD SPECTRAL PROFILES PER SPECIES AND DEPTH CLASS
# =============================================================================

profile_list <- lapply(species_spc, function(sp) {

  df_sp <- df_mono %>%
    filter(
      .data[[sp]] >= 80,
      .data[[sp]] <= 100
    )

  if (nrow(df_sp) < 3) {
    message("Skipped ", sp, ": insufficient samples")
    return(NULL)
  }

  df_long <- df_sp %>%
    select(all_of(s2_bands), depth_class) %>%
    pivot_longer(
      cols      = all_of(s2_bands),
      names_to  = "band",
      values_to = "value"
    )

  df_long %>%
    group_by(depth_class, band) %>%
    summarise(
      species = sp,
      n       = n(),
      median  = median(value, na.rm = TRUE),
      q25     = quantile(value, 0.25, na.rm = TRUE),
      q75     = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
})

profile_df <- bind_rows(profile_list)
cat("Spectral profile rows:", nrow(profile_df), "\n")

# =============================================================================
# 6. FIGURE 3 — Sentinel-2 spectral profiles (depth-stratified)
# =============================================================================

profile_df$band <- factor(
  profile_df$band,
  levels = s2_bands,
  labels = band_labels[s2_bands]
)

p <- ggplot(profile_df,
            aes(x = band, y = median, group = species, color = species)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.15,
    alpha = 0.4
  ) +
  facet_wrap(~ depth_class, nrow = 1) +
  labs(
    x     = "Sentinel-2 band",
    y     = paste0("Surface reflectance (dimensionless DN/10000, p", PERC_TARGET, ")"),
    color = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 20, hjust = 1),
    legend.position  = "bottom"
  )

print(p)

ggsave(
  filename = file.path(out_dir, paste0("Figure_3_S2_spectral_profiles_p", PERC_TARGET, ".png")),
  plot     = p,
  width    = 14, height = 6, dpi = 300
)

# =============================================================================
# 7. SAVE OUTPUTS
# =============================================================================

write_csv(
  profile_df %>% mutate(band = as.character(band)),
  file.path(out_dir, paste0("S2_spectral_profiles_mono_high_cover_p", PERC_TARGET, ".csv"))
)

cat("STEP-04 completed: spectral profiles saved\n")
