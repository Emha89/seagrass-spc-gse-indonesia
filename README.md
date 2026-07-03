# Seagrass Percent Cover Estimation Using Google Satellite Embeddings

R analysis scripts for:

> Hafizt, M., Phinn, S., Wicaksono, P., Hernawan, U.E., Salsabila, H.N., Lyons, M., McMahon, K., & Roelfsema, C. (2026). **Interpretable Estimation of Seagrass Percent Cover Across Extensive Coastal Environments in Indonesia Using Google Satellite Embeddings**. *International Journal of Digital Earth*. Manuscript ID: TJDE-2026-0306.

---

## Contents

```
R_scripts/
├── STEP_01_S2_spectral_sensitivity.R       # Sentinel-2 spectral sensitivity to SPC (Figure 4)
├── STEP_02_GSE_spectral_sensitivity.R      # GSE embedding sensitivity to SPC (Figure 5)
├── STEP_03_GSE_S2_interpretability.R       # GSE–Sentinel-2 linkage per morphology (Figure 6)
├── STEP_03b_GSE_S2_interpretability_by_species.R  # Species-level GSE–S2 linkage (Figure 7)
├── STEP_04_S2_spectral_profiles.R          # Sentinel-2 spectral profiles per depth class (Figure 3)
└── STEP_05_two_stage_RF_LOLOCV.R          # Two-stage RF model + LOLO-CV (Tables 2–3, Figure 8)

GEE_scripts/
├── 01_S2_percentile_composite.js           # Sentinel-2 annual percentile composite generation
├── 02_GSE_feature_extraction.js            # GSE annual layer extraction per field survey year
└── 03_SPC_spatial_deployment.js            # Spatial SPC prediction deployment using trained model
```

---

## Methods Overview

The analysis follows a five-step framework:

1. **STEP 01** — Spearman correlation and GAM-based sensitivity analysis of Sentinel-2 spectral bands to seagrass percent cover (SPC), stratified by canopy morphology class.
2. **STEP 02** — Equivalent sensitivity analysis applied to 64 GSE embedding dimensions.
3. **STEP 03** — Interpretability bridge: pairwise Spearman correlations between top-ranked GSE dimensions and Sentinel-2 band-percentile combinations, per morphology class and at the species level.
4. **STEP 04** — Sentinel-2 multispectral reflectance profiles for mono-species seagrass under high cover conditions, stratified by water depth class.
5. **STEP 05** — Two-stage Random Forest modelling framework under leave-one-location-out (LOLO) cross-validation:
   - **Stage A**: morphology classification (3 classes) using GSE dimensions
   - **Stage B**: SPC regression using GSE dimensions + morphology class probabilities

---

## Requirements

### R packages

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "purrr",
  "mgcv", "randomForest", "Metrics", "caret",
  "ggplot2", "forcats", "stringr", "tibble"
))
```

R version >= 4.0 recommended.

### Google Earth Engine

A Google Earth Engine account is required to run the GEE scripts. Access the GEE Code Editor at: https://code.earthengine.google.com

---

## Data

Field training data are **not included** in this repository. Data are available from the corresponding author (Muhammad Hafizt, m.hafizt@uq.edu.au) upon reasonable request, subject to data governance requirements of the National Research and Innovation Agency of Indonesia (BRIN).

### Expected input data structure

The input CSV file (`input_data.csv`) should contain the following columns:

| Column | Description |
|--------|-------------|
| `loc` | Study site identifier (Ayau, Bintan, Karimunjawa, Komodo, Rote) |
| `year` | Survey year (2018, 2019, 2021, 2022) |
| `gee_id` | Unique observation ID |
| `total_SPC` | Total seagrass percent cover (0–100%) |
| `sg_morpho` | Canopy morphology class (mono, mixed_short, mixed_long) |
| `Ea_SPC` | *Enhalus acoroides* SPC (0–100%) |
| `Th_SPC`, `Cr_SPC`, ... | Per-species SPC columns |
| `B2_p0` to `B8_p100` | Sentinel-2 band-percentile values (DN scaled × 10000) |
| `GSE_A00` to `GSE_A63` | GSE embedding dimensions (64 total) |
| `depth` | Water depth (m, negative values from Allen Coral Atlas) |

---

## Usage

Update the `USER CONFIGURATION` section at the top of each script with your local file paths before running. Scripts are designed to be run sequentially: STEP 01 → 02 → 03 → 03b → 04 → 05A.

---

## Spatial Prediction

A complete time series of spatial SPC prediction results (2017–2024) for all five study locations is available at:
https://muhammadhafizt.users.earthengine.app/view/seagrasspercentcoverindonesia

---

## License

MIT License. See LICENSE file for details.

---

## Contact

Muhammad Hafizt  
School of the Environment, The University of Queensland, Brisbane, Australia  
National Research and Innovation Agency of Indonesia (BRIN)  
m.hafizt@uq.edu.au
