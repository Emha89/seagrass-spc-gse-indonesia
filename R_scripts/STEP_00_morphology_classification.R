# =============================================================================
# STEP_00_morphology_classification.R
#
# Purpose:
#   Classify seagrass canopy morphology into three classes based on
#   species-level percent cover (SPC) per quadrat observation.
#   The resulting morphology column (sg_morpho) and three-class
#   reclassification (morph3) are used as inputs to all subsequent
#   analysis scripts (STEP_01 to STEP_05).
#
#   Three morphology classes (sg_morpho):
#     mono        - single species per quadrat
#     mixed_short - multi-species assemblage without Enhalus acoroides
#     mixed_long  - multi-species assemblage including Enhalus acoroides
#
#   Three-class reclassification for modelling (morph3):
#     mono_Ea                     - mono Enhalus acoroides meadow
#     mixed_short_plus_mono_short - short-leaved canopy (mono non-Ea or mixed_short)
#     mixed_long                  - mixed canopy including Enhalus acoroides
#
#   Classification is based on:
#     - Number of species with SPC > 0 per quadrat
#     - Presence of Enhalus acoroides (Ea_SPC > 0)
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
#   Data are available from the corresponding author upon reasonable
#   request, subject to BRIN data governance requirements.
#   Update the USER CONFIGURATION section before running.
# =============================================================================

library(dplyr)
library(readr)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input field dataset (per-species SPC per quadrat observation)
# Required columns: gee_id, Ea_SPC, Th_SPC, Cr_SPC, Cs_SPC, Si_SPC,
#                   Hu_SPC, Ho_SPC, Hp_SPC, Tc_SPC, Hm_SPC, Hs_SPC, Hd_SPC
data_path <- "path/to/your/input_field_data.csv"

# Output path for classified dataset
out_path <- "path/to/your/output_directory/field_data_with_morphology.csv"

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df <- read_csv(data_path, show_col_types = FALSE)
cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. DEFINE SPECIES COLUMNS
# =============================================================================

species_cols <- c(
  "Ea_SPC", "Th_SPC", "Cr_SPC", "Cs_SPC", "Si_SPC",
  "Hu_SPC", "Ho_SPC", "Hp_SPC", "Tc_SPC",
  "Hm_SPC", "Hs_SPC", "Hd_SPC"
)

# Replace NA in species columns with 0 (absent)
df[species_cols] <- df[species_cols] %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# =============================================================================
# 3. CLASSIFY sg_morpho (PRIMARY 3-CLASS MORPHOLOGY)
# =============================================================================
#
# Classification rules:
#   mono        - only one species has SPC > 0
#   mixed_long  - two or more species with SPC > 0, including Ea_SPC > 0
#   mixed_short - two or more species with SPC > 0, without Ea_SPC
#
# =============================================================================

df <- df %>%
  mutate(
    n_species_present = rowSums(across(all_of(species_cols)) > 0, na.rm = TRUE),

    sg_morpho = case_when(
      n_species_present == 1                        ~ "mono",
      n_species_present > 1 & Ea_SPC > 0           ~ "mixed_long",
      n_species_present > 1 & Ea_SPC == 0          ~ "mixed_short",
      TRUE                                           ~ NA_character_
    )
  )

cat("\nsg_morpho distribution:\n")
print(df %>% count(sg_morpho))

# =============================================================================
# 4. RECLASSIFY TO morph3 (THREE-CLASS FOR RF MODELLING)
# =============================================================================
#
# morph3 reclassification (consistent with STEP_05 R script and
# GEE deployment script):
#
#   mono_Ea                      - mono quadrat with Ea_SPC > 0
#                                  (Enhalus acoroides dominant meadow)
#   mixed_short_plus_mono_short  - mixed_short OR mono without Ea_SPC
#                                  (short-leaved canopy)
#   mixed_long                   - mixed assemblage including Ea_SPC
#                                  (long-leaved, structurally complex canopy)
#
# =============================================================================

morph_levels <- c(
  "mixed_short_plus_mono_short",
  "mixed_long",
  "mono_Ea"
)

df <- df %>%
  mutate(
    morph3 = case_when(
      sg_morpho == "mixed_long"                          ~ "mixed_long",
      sg_morpho == "mono" & Ea_SPC > 0                  ~ "mono_Ea",
      sg_morpho %in% c("mixed_short", "mono")           ~ "mixed_short_plus_mono_short",
      TRUE                                               ~ NA_character_
    ),
    morph3 = factor(morph3, levels = morph_levels)
  )

cat("\nmorph3 distribution:\n")
print(df %>% count(morph3) %>% mutate(prop = round(n / sum(n), 3)))

# =============================================================================
# 5. CALCULATE TOTAL SPC
# =============================================================================

# total_SPC = maximum species cover per quadrat
# This represents the best estimate of total canopy cover when
# species covers are recorded independently and may overlap
df <- df %>%
  mutate(
    total_SPC = pmax(!!!syms(species_cols), na.rm = TRUE)
  )

cat("\ntotal_SPC summary:\n")
print(summary(df$total_SPC))

# =============================================================================
# 6. SAVE OUTPUT
# =============================================================================

write_csv(df, out_path)
cat("\nSTEP-00 completed: morphology classification saved to", out_path, "\n")
cat("Columns added: sg_morpho, morph3, n_species_present, total_SPC\n")
