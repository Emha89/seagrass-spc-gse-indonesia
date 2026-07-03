# =============================================================================
# STEP_05A_two_stage_RF_LOLOCV.R
#
# Purpose:
#   Implement the two-stage Random Forest (RF) modelling framework for
#   seagrass percent cover (SPC) estimation using Google Satellite Embeddings
#   (GSE), evaluated under leave-one-location-out (LOLO) cross-validation.
#
#   Stage A: RF classification of seagrass canopy morphology (3 classes)
#            using GSE dimensions as predictors. Morphology class membership
#            probabilities are extracted from the training set (OOF votes)
#            and from the test set predictions.
#
#   Stage B: RF regression of total SPC using GSE dimensions combined with
#            morphology class probabilities from Stage A as additional
#            predictors.
#
#   Outputs correspond to:
#     - Table 2 (morphology classification accuracy, LOLO)
#     - Table 3 (SPC regression performance, LOLO)
#     - Figure 8 (spatial deployment and LOLO performance)
#     - Supplementary Table S2 (SPC performance by cover interval)
#     - Supplementary Figure S3 (predictor importance stability)
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
library(randomForest)
library(Metrics)
library(caret)
library(purrr)
library(tidyr)
library(ggplot2)
library(tibble)

set.seed(42)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to input dataset
# Required columns: loc, year, gee_id, total_SPC, sg_morpho, Ea_SPC,
#                   GSE_A00 to GSE_A63
data_path <- "path/to/your/input_data.csv"

# Output directory
out_dir <- "path/to/your/output_directory"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df <- read_csv(data_path, show_col_types = FALSE)
cat("Data loaded:", nrow(df), "rows\n")

# =============================================================================
# 2. CONFIGURATION
# =============================================================================

gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))

required_cols <- c(
  "loc", "year", "gee_id", "total_SPC",
  "sg_morpho", "Ea_SPC", gse_vars
)

# RF hyperparameters
ntree_cls <- 700
ntree_reg <- 700

# Minimum test samples per LOLO fold
min_samples_test <- 5

# Grid search values for mtry (classification)
mtry_cls_grid <- unique(pmax(
  2, c(8, 12, 16, floor(sqrt(length(gse_vars))), floor(length(gse_vars) / 3))
))

nodesize_cls <- 3

# Grid search values for mtry and nodesize (regression)
mtry_reg_grid <- unique(pmax(
  2, c(8, 12, 16,
       floor(sqrt(length(gse_vars) + 3)),
       floor((length(gse_vars) + 3) / 3))
))
nodesize_reg_grid <- c(3, 5)

# =============================================================================
# 3. CLEAN AND FILTER DATA
# =============================================================================

df <- df %>%
  select(any_of(required_cols)) %>%
  mutate(
    year      = as.integer(year),
    gee_id    = as.character(gee_id),
    loc       = as.character(loc),
    sg_morpho = as.character(sg_morpho)
  ) %>%
  filter(
    !is.na(loc),
    !is.na(total_SPC),
    total_SPC >= 0,
    total_SPC <= 100,
    !is.na(sg_morpho)
  ) %>%
  filter(if_all(all_of(gse_vars), ~ !is.na(.x)))

cat("After filtering:", nrow(df), "rows\n")

# =============================================================================
# 4. MORPHOLOGY RECLASSIFICATION (3 CLASSES)
# =============================================================================

# Three canopy morphology classes:
#   mixed_short_plus_mono_short: short-leaved assemblages (no Enhalus acoroides)
#   mixed_long:                  assemblages including Enhalus acoroides
#   mono_Ea:                     pure Enhalus acoroides mono-species meadows

morph_levels <- c(
  "mixed_short_plus_mono_short",
  "mixed_long",
  "mono_Ea"
)

df <- df %>%
  mutate(
    morph3 = case_when(
      sg_morpho == "mixed_long"                                ~ "mixed_long",
      sg_morpho == "mono" & !is.na(Ea_SPC) & Ea_SPC > 0      ~ "mono_Ea",
      sg_morpho %in% c("mixed_short", "mono")                 ~ "mixed_short_plus_mono_short",
      TRUE ~ NA_character_
    ),
    morph3 = factor(morph3, levels = morph_levels)
  ) %>%
  filter(!is.na(morph3))

cat("\nMorphology class distribution:\n")
print(df %>% count(morph3) %>% mutate(prop = n / sum(n)))

# =============================================================================
# 5. VALID LOCATIONS FOR LOLO
# =============================================================================

valid_locs <- df %>%
  count(loc) %>%
  filter(n >= min_samples_test) %>%
  pull(loc)

cat("\nValid LOLO locations:\n")
print(valid_locs)

if (length(valid_locs) == 0) stop("No valid locations for LOLO.")

# =============================================================================
# 6. HELPER FUNCTION
# =============================================================================

# Ensure all morphology probability columns are present (fill with 0 if absent)
ensure_prob_cols <- function(prob_df, labels) {
  for (cc in paste0("P_", labels)) {
    if (!cc %in% names(prob_df)) prob_df[[cc]] <- 0
  }
  prob_df[, paste0("P_", labels)]
}

# =============================================================================
# 7. STORAGE
# =============================================================================

cls_metrics_list <- list()
reg_metrics_list <- list()
confusion_list   <- list()
oof_prob_list    <- list()
reg_imp_list     <- list()
pred_spc_test_df <- list()

# =============================================================================
# 8. LOLO CROSS-VALIDATION LOOP
# =============================================================================

for (test_loc in valid_locs) {

  cat("\n--------------------------------------------------\n")
  cat("LOLO iteration | Test location:", test_loc, "\n")
  cat("--------------------------------------------------\n")

  train_df <- df %>% filter(loc != test_loc)
  test_df  <- df %>% filter(loc == test_loc)

  # -----------------------------------------------------------------
  # STAGE A: MORPHOLOGY CLASSIFICATION
  # -----------------------------------------------------------------

  # Grid search over mtry to minimise OOB error
  best_cls <- NULL
  best_oob <- Inf

  for (mtry_val in mtry_cls_grid) {
    rf_try <- randomForest(
      x        = train_df[, gse_vars],
      y        = train_df$morph3,
      ntree    = ntree_cls,
      mtry     = min(mtry_val, length(gse_vars)),
      nodesize = nodesize_cls
    )
    oob <- tail(rf_try$err.rate[, "OOB"], 1)
    if (oob < best_oob) {
      best_oob <- oob
      best_cls <- rf_try
    }
  }

  cat("Best OOB error (morphology RF):", round(best_oob, 4), "\n")

  # Extract OOF class probabilities for training set (from RF votes)
  train_probs <- ensure_prob_cols(
    as.data.frame(best_cls$votes) %>% setNames(paste0("P_", colnames(.))),
    morph_levels
  )

  # Predict class probabilities for test set
  test_probs <- ensure_prob_cols(
    as.data.frame(predict(best_cls, test_df[, gse_vars], type = "prob")) %>%
      setNames(paste0("P_", colnames(.))),
    morph_levels
  )

  # Evaluate classification performance
  pred_class <- predict(best_cls, test_df[, gse_vars])
  cm         <- confusionMatrix(pred_class, test_df$morph3)

  cat("Class-wise accuracy (test set):\n")
  print(round(prop.table(cm$table, margin = 2), 2))

  cm_df <- as.data.frame(cm$table)
  names(cm_df)[1:3] <- c("Predicted", "Reference", "Count")

  confusion_list[[test_loc]] <- cm_df %>% mutate(TestLocation = test_loc)

  cls_metrics_list[[test_loc]] <- tibble(
    TestLocation = test_loc,
    Accuracy     = as.numeric(cm$overall["Accuracy"])
  )

  oof_prob_list[[test_loc]] <- bind_cols(
    test_df %>% select(gee_id, year, loc, total_SPC, morph3),
    test_probs
  )

  # -----------------------------------------------------------------
  # STAGE B: SPC REGRESSION (GSE + morphology probabilities)
  # -----------------------------------------------------------------

  # Combine GSE dimensions and morphology class probabilities as predictors
  train_x <- bind_cols(train_df[, gse_vars], train_probs)
  test_x  <- bind_cols(test_df[, gse_vars], test_probs)

  # Grid search over mtry and nodesize to minimise OOB MSE
  best_reg <- NULL
  best_mse <- Inf

  for (mtry_val in mtry_reg_grid) {
    for (nodesize_val in nodesize_reg_grid) {
      rf_try <- randomForest(
        x        = train_x,
        y        = train_df$total_SPC,
        ntree    = ntree_reg,
        mtry     = min(mtry_val, ncol(train_x)),
        nodesize = nodesize_val,
        importance = TRUE
      )
      mse <- tail(rf_try$mse, 1)
      if (mse < best_mse) {
        best_mse <- mse
        best_reg <- rf_try
      }
    }
  }

  cat("Best OOB MSE (SPC regression):", round(best_mse, 2), "\n")

  pred_spc <- predict(best_reg, test_x)

  pred_spc_test_df[[test_loc]] <- tibble(
    TestLocation = test_loc,
    SPC_obs      = test_df$total_SPC,
    SPC_pred     = pred_spc
  )

  reg_metrics_list[[test_loc]] <- tibble(
    TestLocation = test_loc,
    RMSE         = rmse(test_df$total_SPC, pred_spc),
    MAE          = mae(test_df$total_SPC, pred_spc),
    R2           = caret::R2(pred_spc, test_df$total_SPC)
  )

  print(reg_metrics_list[[test_loc]])

  reg_imp_list[[test_loc]] <- importance(best_reg, type = 2) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    rename(importance = IncNodePurity) %>%
    mutate(TestLocation = test_loc)
}

# =============================================================================
# 9. COMPILE RESULTS
# =============================================================================

cls_metrics_df    <- bind_rows(cls_metrics_list)
reg_metrics_df    <- bind_rows(reg_metrics_list)
confusion_df      <- bind_rows(confusion_list)
oof_probs_df      <- bind_rows(oof_prob_list)
reg_importance_df <- bind_rows(reg_imp_list)
spc_pred_all      <- bind_rows(pred_spc_test_df)

# =============================================================================
# 10. COMPUTE MACRO-F1 FROM CONFUSION MATRICES (Table 2)
# =============================================================================

macro_f1_df <- bind_rows(confusion_list) %>%
  group_by(TestLocation) %>%
  group_modify(~ {
    cm      <- xtabs(Count ~ Predicted + Reference, data = .x)
    classes <- intersect(rownames(cm), colnames(cm))
    f1_per_class <- sapply(classes, function(cl) {
      tp    <- cm[cl, cl]
      fp    <- sum(cm[cl, ]) - tp
      fn    <- sum(cm[, cl]) - tp
      denom <- 2 * tp + fp + fn
      if (denom == 0) return(NA_real_)
      2 * tp / denom
    })
    tibble(MacroF1 = mean(f1_per_class, na.rm = TRUE))
  })

cls_metrics_df <- cls_metrics_df %>%
  left_join(macro_f1_df, by = "TestLocation") %>%
  left_join(
    confusion_df %>%
      group_by(TestLocation) %>%
      summarise(n_test = sum(Count), .groups = "drop"),
    by = "TestLocation"
  )

cat("\nTable 2 — Morphology classification (LOLO):\n")
print(cls_metrics_df)

cat("\nTable 3 — SPC regression (LOLO):\n")
print(reg_metrics_df)

# =============================================================================
# 11. SAVE OUTPUTS
# =============================================================================

write_csv(cls_metrics_df,    file.path(out_dir, "Table2_LOLO_morphology_classification.csv"))
write_csv(reg_metrics_df,    file.path(out_dir, "Table3_LOLO_SPC_regression.csv"))
write_csv(confusion_df,      file.path(out_dir, "LOLO_morphology_confusion_matrices.csv"))
write_csv(oof_probs_df,      file.path(out_dir, "LOLO_morphology_OOF_probabilities.csv"))
write_csv(spc_pred_all,      file.path(out_dir, "LOLO_SPC_test_predictions.csv"))
write_csv(reg_importance_df, file.path(out_dir, "LOLO_SPC_regression_importance.csv"))

# =============================================================================
# 12. PREDICTOR IMPORTANCE STABILITY (Supplementary Figure S3)
# =============================================================================

importance_stability <- reg_importance_df %>%
  group_by(variable) %>%
  summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance),
    .groups         = "drop"
  ) %>%
  arrange(desc(mean_importance))

write_csv(
  importance_stability,
  file.path(out_dir, "FigS3_importance_stability.csv")
)

p_imp <- ggplot(
  importance_stability %>% slice_head(n = 20),
  aes(x = reorder(variable, mean_importance), y = mean_importance)
) +
  geom_col(fill = "grey70") +
  geom_errorbar(
    aes(ymin = mean_importance - sd_importance,
        ymax = mean_importance + sd_importance),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    x = "Predictor",
    y = "Importance (IncNodePurity)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "FigS3_importance_stability.png"),
  p_imp, width = 7, height = 6, dpi = 300
)

# =============================================================================
# 13. SPC PERFORMANCE BY COVER INTERVAL (Supplementary Table S2)
# =============================================================================

spc_interval_metrics <- spc_pred_all %>%
  mutate(
    SPC_bin = case_when(
      SPC_obs < 20 ~ "0-20",
      SPC_obs < 40 ~ "20-40",
      SPC_obs < 60 ~ "40-60",
      SPC_obs < 80 ~ "60-80",
      TRUE         ~ "80-100"
    )
  ) %>%
  group_by(SPC_bin) %>%
  summarise(
    n    = n(),
    RMSE = rmse(SPC_obs, SPC_pred),
    R2   = caret::R2(SPC_pred, SPC_obs),
    .groups = "drop"
  )

cat("\nSupplementary Table S2 — SPC performance by cover interval:\n")
print(spc_interval_metrics)

write_csv(
  spc_interval_metrics,
  file.path(out_dir, "TableS2_SPC_performance_by_interval.csv")
)

# =============================================================================
# 14. SAMPLE COMPOSITION SUMMARY (Supplementary Table S1)
# =============================================================================

n_total_dataset <- nrow(df)

location_summary <- df %>%
  group_by(loc) %>%
  summarise(
    n_total_location = n(),
    years_sampled    = paste(sort(unique(year)), collapse = ", "),
    .groups          = "drop"
  ) %>%
  filter(loc %in% valid_locs) %>%
  mutate(
    n_train_model = n_total_dataset - n_total_location,
    n_test        = n_total_location
  ) %>%
  select(loc, n_total_location, n_train_model, n_test, years_sampled) %>%
  arrange(desc(n_total_location))

cat("\nSupplementary Table S1 — Sample composition per location:\n")
print(location_summary)

write_csv(
  location_summary,
  file.path(out_dir, "TableS1_sample_composition_per_location.csv")
)

cat("\nSTEP-05A completed\n")
