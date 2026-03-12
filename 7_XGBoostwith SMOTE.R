###############################################################################
### 1️⃣  Load packages
###############################################################################
library(DMwR)         # For SMOTE oversampling
library(xgboost)
library(lubridate)
library(dplyr)
library(reshape2)
library(ggplot2)
library(caret)
library(depmixS4)
library(pROC)
library(sva)
library(hms)
library(readr)
library(corrplot)
library(glmnet)
library(stringr)
library(plyr)
library(skimr)
library(caretEnsemble)
#library(caretSDM)
library(rminer)
library(MLeval)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)

rm(list = ls())

set.seed(1234)
# ==============================================================================
# -------------------Load and clean phenotype data------------------------------
# ==============================================================================
pheno_all <- read.csv("pheno_data_345.csv", header = TRUE, fill = TRUE)

## ---- Batch1: Training ----
batch1_names <- read.csv("batch1.csv", stringsAsFactors = FALSE)
batch1_samples <- subset(pheno_all, samples %in% batch1_names$sample_id)
## Clean time strings — remove empty or invalid ones
batch1_samples <- batch1_samples %>%
  mutate(time_window_clean = trimws(time_window)) %>%
  mutate(time_window_clean = ifelse(grepl("^\\d{1,4}:\\d{2}:\\d{2}$", time_window_clean), time_window_clean, NA))
## Parse time safely
batch1_samples <- batch1_samples %>%  mutate(time_parsed = lubridate::hms(time_window_clean))
## filter
batch1_samples_clean <- subset(batch1_samples, clinical_diagnosis != "TIA")
batch1_samples_clean <- subset(batch1_samples_clean, final_diagnosis != "" & !is.na(final_diagnosis)) # (S0408)
batch1_samples_clean <- subset(batch1_samples_clean, treatment != "MT")
batch1_samples_clean <- subset(batch1_samples_clean, samples != "S0455") #unknown diagnosis
batch1_samples_clean <- subset(batch1_samples_clean, samples != "S0085") #confusion TPA/Mimic
batch1_samples_clean <- subset(batch1_samples_clean, !(samples %in% c("S0145", "S0227", "S0247", "S0580", "S0613"))) # HC conversion, second or multiple stroke
threshold <- dhours(4) + dminutes(30)  # or dhours(4.5)
batch1_samples_clean <- subset(batch1_samples_clean, !(treatment == "NONE" & final_diagnosis == "Stroke" & time_parsed > threshold))
batch1_samples_clean$severity_score <- as.numeric(batch1_samples_clean$severity_score)
batch1_samples_clean <- subset(batch1_samples_clean, !(treatment == "NONE" & final_diagnosis == "Stroke" & severity_score < 5))
batch1_samples_clean <- subset(batch1_samples_clean, !(samples %in% c("S0001", "S0092", "S0178", "S0213", "S0331", "S0629", "S0657"))) # low read depth (1M)
row.names(batch1_samples_clean) <- batch1_samples_clean$samples
pheno_train <- batch1_samples_clean

## ---- Batch2 + Batch3: Validation ----
batch2_names <- read.csv("batch2.csv", stringsAsFactors = FALSE)
batch3_names <- read.csv("batch3.csv", stringsAsFactors = FALSE)

pheno_val <- subset(pheno_all, samples %in% c(batch2_names$sample_id, batch3_names$sample_id)) %>%
  filter(final_diagnosis != "" & !is.na(final_diagnosis))

## Excluded samples (Reasons for exclusion are detailed in pheno_data_345.csv)
exclude_samples <- c("S0675","S0699","S0702","S0811","S0802","S0807","S0866","S0873","S0964", "S0791",
                     "S0815","S0826","S0838","S0839","S0844","S0849","S0852","S0854","S0865","S0896","S0912","S0948","S0777","S0198")

pheno_val <- pheno_val[!pheno_val$samples %in% exclude_samples, ]
row.names(pheno_val) <- pheno_val$samples

# ==============================================================================
# ------------------------ Load expression data --------------------------------
# ==============================================================================
df_train <- read.table("4_transcript_count_matrix_batch1.csv", sep = ",", header = TRUE)
df_val2 <- read.table("5_transcript_count_matrix_batch2.csv", sep = ",", header = TRUE)
df_val3 <- read.table("6_transcript_count_matrix_batch3.csv", sep = ",", header = TRUE)

## Clean gene names
clean_genes <- function(df){
  df$transcript_id <- gsub("[.|\\-|]", "_", df$transcript_id)
  rownames(df) <- df$transcript_id
  df$transcript_id <- NULL
  names(df) <- substr(names(df), 1, 5)
  return(df)
}

df_train <- clean_genes(df_train)
df_val2 <- clean_genes(df_val2)
df_val3 <- clean_genes(df_val3)

## Filter lowly expressed genes 
filter_low_counts <- function(df, min_count = 10, min_samples_prop = 0.1) {
  min_samples <- ceiling(ncol(df) * min_samples_prop)
  keep <- rowSums(df >= min_count) >= min_samples
  df_filtered <- df[keep, ]
  message("Genes retained: ", sum(keep), " / ", nrow(df))
  return(df_filtered)
}

df_train <- filter_low_counts(df_train, min_count = 10, min_samples_prop = 0.1)
df_val2  <- filter_low_counts(df_val2,  min_count = 10, min_samples_prop = 0.1)
df_val3  <- filter_low_counts(df_val3,  min_count = 10, min_samples_prop = 0.1)

## perform CPM calculation
set.seed(1234)
to_cpm <- function(df) {
  counts <- as.matrix(df)
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  return(cpm)
}

## Filtered CPM matrices
cpm_train <- to_cpm(df_train)
cpm_val2  <- to_cpm(df_val2)
cpm_val3  <- to_cpm(df_val3)

cpm_train <- t(cpm_train)
cpm_val2  <- t(cpm_val2)
cpm_val3  <- t(cpm_val3)

# ==============================================================================
# ------------Subset genes of interest from differential expression-------------
# ==============================================================================
genes1 <- read.csv("results/1_HS_vs_TPA_50_t.csv")
genes2 <- read.csv("results/2_TPA_vs_SM_50_t.csv")
genes3 <- read.csv("results/3_HS_vs_SM_50_t.csv")

g <- setdiff(genes1$X, genes2$X)
glist <- setdiff(g, genes3$X)

glist  <- gsub("[.|-]", "_", glist)
genes_to_use <- glist

## Subset CPM to selected genes
expr_train <- cpm_train[, colnames(cpm_train) %in% genes_to_use, drop = FALSE]
expr_val2  <- cpm_val2[,  colnames(cpm_val2)  %in% genes_to_use, drop = FALSE]
expr_val3  <- cpm_val3[,  colnames(cpm_val3)  %in% genes_to_use, drop = FALSE]

## Keep only common genes
common_genes <- Reduce(intersect, list(colnames(expr_train),
                                       colnames(expr_val2),
                                       colnames(expr_val3)))
expr_train <- expr_train[, common_genes, drop = FALSE]
expr_val2  <- expr_val2[,  common_genes, drop = FALSE]
expr_val3  <- expr_val3[,  common_genes, drop = FALSE]

# ==============================================================================
# ---------------Merge validation batches and batch correction -----------------
# ==============================================================================
expr_val <- rbind(expr_val2, expr_val3)
common_val_samples <- intersect(rownames(expr_val), rownames(pheno_val))
expr_val <- expr_val[common_val_samples, ]

pheno_val$batch <- ifelse(
  pheno_val$samples %in% rownames(expr_val2), "batch2",
  ifelse(pheno_val$samples %in% rownames(expr_val3), "batch3", NA)
)

mod_val <- model.matrix(~ final_diagnosis, data = pheno_val)                    #Includes final_diagnosis in the model preserves biological differences
expr_val_corrected <- ComBat(dat = t(expr_val), batch = pheno_val$batch, mod = mod_val,
                             par.prior = TRUE, prior.plots = FALSE)             # batch1 or training set as reference cause over-correction wrongly preserve biological signal
expr_val_corrected <- t(expr_val_corrected)

# ==============================================================================
# --------------- Prepare training and validation data sets --------------------
# ==============================================================================
train_df <- as.data.frame(expr_train)
train_df <- train_df[rownames(train_df) %in% rownames(pheno_train), ]
train_df$Cat1 <- factor(ifelse(pheno_train$final_diagnosis == "Hemorrhage", "Y", "N"), levels = c("N", "Y"))

val_df <- as.data.frame(expr_val_corrected)
val_df <- val_df[rownames(val_df) %in% rownames(pheno_val), ]
val_df$Cat1 <- factor(ifelse(pheno_val$final_diagnosis == "Hemorrhage", "Y", "N"), levels = c("N", "Y"))

# ==============================================================================
# ---------------------------- Scale data sets ---------------------------------
# ==============================================================================
scaler <- preProcess(train_df[, common_genes], method = c("center","scale"))
train_scaled <- predict(scaler, train_df[, common_genes])
train_scaled$Cat1 <- train_df$Cat1

val_scaled <- predict(scaler, val_df[, common_genes])
val_scaled$Cat1 <- val_df$Cat1

# ==============================================================================
# --------------- Apply SMOTE to balanced training data set --------------------
# ==============================================================================
set.seed(1234)
train_smote <- SMOTE(Cat1 ~ ., data = train_scaled, perc.over = 600, perc.under = 150)
table(train_smote$Cat1)  # Confirm balance


## Train XGBoost model using caret
fitControl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

xgb_grid <- expand.grid(
  nrounds = 200,
  max_depth = c(3, 5, 7),
  eta = c(0.05, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

xgb_model <- train(
  Cat1 ~ .,
  data = train_smote,
  method = "xgbTree",
  metric = "ROC",
  trControl = fitControl,
  tuneGrid = xgb_grid
)

print(xgb_model)
best_model <- xgb_model$finalModel

## Predict on validation data
val_prob <- predict(xgb_model, newdata = val_scaled, type = "prob")[, "Y"]
roc_obj <- roc(val_scaled$Cat1, val_prob, levels = c("N","Y"), direction = "<")
auc_val <- auc(roc_obj)
cat("Validation AUC:", auc_val, "\n")

##  Optimize threshold using Youden’s J
best_thresh <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
best_thresh <- as.numeric(best_thresh)
cat("Optimal threshold (Youden’s J):", best_thresh, "\n")

## Predict classes and evaluate
val_pred <- ifelse(val_prob > best_thresh, "Y", "N")
val_pred <- factor(val_pred, levels = c("N", "Y"))
val_Cat1_factor <- factor(val_scaled$Cat1, levels = c("N", "Y"))

cm <- confusionMatrix(val_pred, val_Cat1_factor, positive = "Y")
print(cm)

# ROC plot
## Extract CV predictions for the best model
train_preds <- xgb_model$pred

## Keep only best-tuned hyperparameters
best_tune <- xgb_model$bestTune

train_preds <- train_preds[
  train_preds$nrounds == best_tune$nrounds &
    train_preds$max_depth == best_tune$max_depth &
    train_preds$eta == best_tune$eta &
    train_preds$colsample_bytree == best_tune$colsample_bytree &
    train_preds$min_child_weight == best_tune$min_child_weight &
    train_preds$subsample == best_tune$subsample,
]

## Training ROC
roc_train <- roc(
  response = train_preds$obs,
  predictor = train_preds$Y,
  levels = c("N", "Y"),
  direction = "<"
)

auc_train <- pROC::auc(roc_train)

roc_val <- roc(
  response = val_scaled$Cat1,
  predictor = val_prob,
  levels = c("N", "Y"),
  direction = "<"
)

auc_val <- pROC::auc(roc_val)

tiff(filename = "plots/Fig_4C.tiff", width = 6, height = 6, units = "in", res = 300)
##Plot Training vs Validation ROC curves
plot(
  roc_train,
  col = "blue",
  lwd = 2,
  main = "ROC Curves: Training vs Validation",
  legacy.axes = TRUE
)

plot(
  roc_val,
  col = "red",
  lwd = 2, lty = 2, 
  add = TRUE
)

legend(
  "bottomright",
  legend = c(
    paste0("Training AUC = ", round(auc_train, 3)),
    paste0("Validation AUC = ", round(auc_val, 3))
  ),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)

dev.off()