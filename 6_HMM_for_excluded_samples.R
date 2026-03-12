# ================================
# STEP 0: Load required libraries
# ================================
library(lubridate)
library(dplyr)
library(reshape2)
library(ggplot2)
library(caret)
library(depmixS4)
library(pROC)
library(sva)

rm(list = ls())

set.seed(1234)
# ================================
# STEP 1: Load and clean phenotype data
# ================================
pheno_all <- read.csv("pheno_data_345.csv", header = TRUE, fill = TRUE)

## ---- Batch1: Training ----
batch1_names <- read.csv("batch1.csv", stringsAsFactors = FALSE)
batch1_samples <- subset(pheno_all, samples %in% batch1_names$sample_id)
# Clean time strings — remove empty or invalid ones
batch1_samples <- batch1_samples %>%
  mutate(time_window_clean = trimws(time_window)) %>%
  mutate(time_window_clean = ifelse(grepl("^\\d{1,4}:\\d{2}:\\d{2}$", time_window_clean), time_window_clean, NA))
# Parse time safely
batch1_samples <- batch1_samples %>%  mutate(time_parsed = lubridate::hms(time_window_clean))
pheno_test <- batch1_samples[!(batch1_samples$samples %in% pheno_train$samples), ]
pheno_test <- pheno_test[!(pheno_test$final_diagnosis == ""), ]
pheno_test <- subset(pheno_test, final_diagnosis != "Hemorrhage")#comment this for CAT1
row.names(pheno_test) <- pheno_test$samples

# Subset to only genes in glist first
expr_test <- cpm_train[, colnames(cpm_train) %in% genes_to_use, drop = FALSE]

# Keep only genes common to ALL three datasets
expr_test <- expr_test[, common_genes, drop = FALSE]

# ================================
# Prepare training and validation matrices
# ================================
test_df <- expr_test
test_df <- as.data.frame(test_df)
test_df <- test_df[rownames(test_df) %in% rownames(pheno_test), ]
test_df$Cat1 <- as.numeric(pheno_test$final_diagnosis == "Stroke")  # change to Hemorrhage for CAT1

test_scaled <- predict(scaler, test_df[, genes_to_use])
test_scaled$Cat1 <- test_df$Cat1

test_model <- depmix(response = formula_obj, data = test_scaled, nstates = best_k, family = binomial())
test_model <- setpars(test_model, getpars(fit_train_model))

test_post <- posterior(test_model, type="viterbi")
#val_states <- val_post$state
#val_target_state <- which.max(tapply(val_scaled$Cat1, val_states, mean))
test_probs <- test_post[, paste0("S", train_target_state)]

roc_test<- roc(test_scaled$Cat1, test_probs)
test_pred <- ifelse(test_probs >= best_threshold, 1, 0)
test_confmat <- confusionMatrix(factor(test_pred, levels=c(0,1)), factor(test_scaled$Cat1, levels=c(0,1)), positive="1")
test_confmat
