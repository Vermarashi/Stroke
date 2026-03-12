###############################################################################
###############################################################################
# Title: Stroke Transcriptomics Analysis: Caret based classifier using SMOTE
# Author: Rashi Verma
# Date: 2025-11-13
# Description: 
#   - Load phenotype and RNA-seq count data
#   - Filter lowly-expressed genes
#   - CPM normalization
#   - Subset genes to DETs
#   - Batch correction for validation data
#   - Apply SMOTE to training set to balance the samples
#   - Train Caret Classifiers on train using differential expressed transcripts
#   - Evaluate on independent validation
#   - Plot expression, ROC and summary statistics
###############################################################################
###############################################################################
# Load packages
library(DMwR)         
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

# ==============================================================================
# ------------- Perform Caret modelling on SMOTE training data set -------------
# ==============================================================================
## Define caret training control
fitControl <- trainControl(
  method = 'repeatedcv',
  repeats = 10,
  number = 10,
  savePredictions = 'final',
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

## Train multiple models on SMOTE data
set.seed(1234)
model_mars <- train(Cat1 ~ ., data=train_smote, method='earth', tuneLength=20, trControl=fitControl, metric="ROC")
model_rf   <- train(Cat1 ~ ., data=train_smote, method='rf', tuneLength=20, trControl=fitControl, importance=TRUE, metric="ROC")
model_svm  <- train(Cat1 ~ ., data=train_smote, method='svmRadial', tuneLength=20, trControl=fitControl, metric="ROC")
model_gbm  <- train(Cat1 ~ ., data=train_smote, method='gbm', tuneLength=10, trControl=fitControl, metric="ROC")
model_nnet <- train(Cat1 ~ ., data=train_smote, method='nnet', tuneLength=20, trControl=fitControl, metric="ROC")
model_pls  <- train(Cat1 ~ ., data=train_smote, method='pls', tuneLength=10, trControl=fitControl, metric="ROC")
model_rpart<- train(Cat1 ~ ., data=train_smote, method='rpart', tuneLength=10, trControl=fitControl, metric="ROC")
model_ranger<- train(Cat1 ~ ., data=train_smote, method='ranger', tuneLength=10, trControl=fitControl, metric="ROC")
model_glmnet<- train(Cat1 ~ ., data=train_smote, method='glmnet', tuneLength=10, trControl=fitControl, metric="ROC")


## Compare models and evaluate performance
models_compare <- resamples(list(RF=model_rf,MARS=model_mars,GLMNET=model_glmnet, 
                                 GBM=model_gbm,SVM=model_svm,RANGER=model_ranger, 
                                 RPART=model_rpart, PLS=model_pls, NNET=model_nnet))

## save session (optional)
save.image(file = "caret_classifiers.RData")

Caret_AUC_ROC <- evalm(list(RF=model_rf,MARS=model_mars,GLMNET=model_glmnet, 
                    GBM=model_gbm,SVM=model_svm,RANFER=model_ranger, 
                    RPART=model_rpart, PLS=model_pls, NNET=model_nnet),
               gnames=c('rf', 'mars', 'glmnet', 'GBM', 'svm', 'ranger', 'rpart', 'pls','nnet'),
               rlinethick=0.8,fsize=8,
               plots='r')

## Plot Fig.4A
tiff(file="plots/Fig_4A.tiff", unit= "in", res = 600, width = 4, height = 4)
Caret_AUC_ROC$roc
dev.off()

# ==============================================================================
# ------------- Evaluate models on training and validation data set ------------
# ==============================================================================
pred_models <- list(MARS=model_mars, RF=model_rf, GBM=model_gbm, SVM=model_svm,
                    RANGER=model_ranger, RPART=model_rpart, PLS=model_pls, NNET=model_nnet,
                    GLMNET=model_glmnet)

## Training Dataset
for (mod_name in names(pred_models)) {
  cat("===== Model:", mod_name, "=====\n")
  predicted_train <- predict(pred_models[[mod_name]], newdata=train_scaled)
  print(confusionMatrix(reference=train_scaled$Cat1, data=predicted_train, mode='everything', positive='Y'))
}

## Validation Dataset
for (mod_name in names(pred_models)) {
  cat("===== Model:", mod_name, "=====\n")
  predicted_val <- predict(pred_models[[mod_name]], newdata=val_scaled)
  print(confusionMatrix(reference=val_scaled$Cat1, data=predicted_val, mode='everything', positive='Y'))
}

# ==============================================================================
# ----- Refit and evaluate best models on training and validation data set -----
# ==============================================================================
## Refit the best (GBM) model on the training data set
model_gbm<- train(Cat1 ~ ., data=train_smote, method='gbm', tuneLength=10, trControl=fitControl, metric="ROC")

### train data set
cat("gbm model on train\n")
predicted_train <- predict(model_gbm, newdata=train_scaled)
print(confusionMatrix(reference=train_scaled$Cat1, data=predicted_train, mode='everything', positive='Y'))

### validation data set
cat("gbm model on validation\n")
predicted_val <- predict(model_gbm, newdata=val_scaled)
print(confusionMatrix(reference = val_scaled$Cat1, data = predicted_val, mode='everything', positive='Y'))

## plot Fig.4B
## AUC_ROC plot of best model (GBM) on training and validation data set
### Get predicted probabilities for GBM
train_probs <- predict(pred_models[["GBM"]], newdata=train_scaled, type="prob")[, "Y"]
val_probs   <- predict(pred_models[["GBM"]], newdata=val_scaled, type="prob")[, "Y"]

### Compute ROC objects
roc_train <- roc(response = train_scaled$Cat1, predictor = train_probs, levels=c("N","Y"), direction="<")
roc_val   <- roc(response = val_scaled$Cat1, predictor = val_probs, levels=c("N","Y"), direction="<")

### Convert to data frames for ggplot
roc_train_df <- data.frame(
  FPR = 1 - roc_train$specificities,
  TPR = roc_train$sensitivities,
  Dataset = "Train"
)

roc_val_df <- data.frame(
  FPR = 1 - roc_val$specificities,
  TPR = roc_val$sensitivities,
  Dataset = "Validation"
)

roc_df <- bind_rows(roc_train_df, roc_val_df)

tiff(file="plots/Fig_4B.tiff", unit= "in", res = 600, width = 4, height = 4)
p <-ggplot(roc_df, aes(x = FPR, y = TPR, color = Dataset, linetype = Dataset)) +
  geom_path(size=1.2) +
  scale_color_manual(values = c("Train" = "blue", "Validation" = "red")) +
  scale_linetype_manual(values = c("Train" = "solid", "Validation" = "dashed")) +
  labs(#title = "ROC Curve for GBM: Train vs Validation",
    x = "False Positive Rate",
    y = "True Positive Rate") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(color = "black", face="bold", hjust=0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "NONE"
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "gray")

dev.off()

####################################################################################################################
# Supplementary Figure 2

####################################################################################################################
####################################################################################################################
####################################################################################################################

# ==============================================================================
# -------------------------- Optional: Statistics ------------------------------
# ==============================================================================
# Compute ROC and AUC for train and validation data set
gbm_probs_t <- predict(pred_models[["GBM"]], newdata = train_scaled, type = "prob")[, "Y"]
roc_obj_t <- roc(response = train_scaled$Cat1, predictor = gbm_probs_t, levels = c("N","Y"), direction = "<")
auc_train <- auc(roc_obj_t)
cat("GBM AUC:", auc_train, "\n")

gbm_probs_v <- predict(pred_models[["GBM"]], newdata = val_scaled, type = "prob")[, "Y"]
roc_obj_v <- roc(response = val_scaled$Cat1, predictor = gbm_probs_v, levels = c("N","Y"), direction = "<")
auc_val <- auc(roc_obj_v)
cat("GBM AUC:", auc_val, "\n")

# Find best threshold (Youden's J)
best_coords <- coords(roc_obj_t, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
best_thresh <- as.numeric(best_coords["threshold"])
best_sens   <- as.numeric(best_coords["sensitivity"])
best_spec   <- as.numeric(best_coords["specificity"])

cat("Best threshold for GBM (Youden's J):", best_thresh, "\n")
cat("Sensitivity:", best_sens, "Specificity:", best_spec, "\n")

# predict classes using best threshold
val_scaled$Cat1 <- factor(val_scaled$Cat1, levels=c("N","Y"))
predicted_class <- ifelse(gbm_probs_v >= best_thresh, "Y", "N")
predicted_class <- factor(predicted_class, levels=c("N","Y"))
conf_mat <- confusionMatrix(predicted_class, reference=val_scaled$Cat1, positive="Y")
print(conf_mat)

# ==============================================================================
# ------------------ GBM with threshold optimization -------------------------
# ==============================================================================
## No change observed after optimizing the threshold

## Train RPART model on SMOTE-balanced data
set.seed(1234)
model_gbm <- train(
  Cat1 ~ .,
  data = train_smote,
  method = "gbm",
  tuneLength = 10,
  trControl = fitControl,
  metric = "ROC"
)

cat("\n===== TRAIN CONFUSION MATRIX (default threshold 0.5) =====\n")
predicted_train_class <- predict(model_gbm, newdata = train_scaled)
print(confusionMatrix(
  reference = train_scaled$Cat1,
  data = predicted_train_class,
  mode = "everything",
  positive = "Y"
))

## Optimize THRESHOLD on TRAIN using probabilities
train_prob <- predict(model_gbm, newdata = train_scaled, type = "prob")[, "Y"]

thresholds <- seq(0, 1, 0.01)

youden_J <- sapply(thresholds, function(th) {
  pred <- ifelse(train_prob > th, "Y", "N")
  cm <- confusionMatrix(
    factor(pred, levels = c("N","Y")),
    train_scaled$Cat1,
    positive = "Y"
  )
  J <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  return(J)
})

opt_threshold <- thresholds[which.max(youden_J)]

cat("\n===== OPTIMAL THRESHOLD (TRAIN) =====\n")
print(opt_threshold)

## Confusion Matrix on TRAIN using optimal threshold

train_pred_opt <- ifelse(train_prob > opt_threshold, "Y", "N")

cat("\n===== TRAIN CONFUSION MATRIX (optimal threshold) =====\n")
print(confusionMatrix(
  factor(train_pred_opt, levels = c("N","Y")),
  train_scaled$Cat1,
  mode = "everything",
  positive = "Y"
))

## Apply same threshold to VALIDATION
val_prob <- predict(model_gbm, newdata = val_scaled, type = "prob")[, "Y"]
val_pred_opt <- ifelse(val_prob > opt_threshold, "Y", "N")

cat("\n===== VALIDATION CONFUSION MATRIX (optimal threshold) =====\n")
print(confusionMatrix(
  factor(val_pred_opt, levels = c("N","Y")),
  val_scaled$Cat1,
  mode = "everything",
  positive = "Y"
))


# ==============================================================================
# ---------------- All Caret model with threshold optimization -----------------
# ==============================================================================
pred_models <- list(
  MARS = model_mars,
  RF = model_rf,
  GBM = model_gbm,
  SVM = model_svm,
  RANGER = model_ranger,
  RPART = model_rpart,
  PLS = model_pls,
  NNET = model_nnet,
  GLMNET = model_glmnet
)

# Sequence of thresholds we will test
thresholds <- seq(0, 1, 0.01)


for (mod_name in names(pred_models)) {
  
  cat("\n=====================\n")
  cat("MODEL:", mod_name, "\n")
  cat("=====================\n")
  
  model <- pred_models[[mod_name]]
  
  # ---------------------------------------------------------
  # 1. TRAIN Probabilities
  # ---------------------------------------------------------
  train_prob <- predict(model, newdata = train_scaled, type = "prob")[, "Y"]
  
  # -------- Optimize threshold using Youden’s J ----------
  J_values <- sapply(thresholds, function(th) {
    pred_class <- ifelse(train_prob > th, "Y", "N")
    cm <- confusionMatrix(
      factor(pred_class, levels=c("N","Y")),
      train_scaled$Cat1,
      positive="Y"
    )
    cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  })
  
  opt_th <- thresholds[which.max(J_values)]
  cat("Optimal threshold =", opt_th, "\n")
  
  # ---------------------------------------------------------
  # 2. Confusion matrix on TRAIN using optimal threshold
  # ---------------------------------------------------------
  train_pred_opt <- ifelse(train_prob > opt_th, "Y", "N")
  cat("\nTRAIN confusion matrix (Optimized threshold)\n")
  print(confusionMatrix(
    factor(train_pred_opt, levels=c("N","Y")),
    train_scaled$Cat1,
    positive="Y"
  ))
  
  # ---------------------------------------------------------
  # 3. Validation using same threshold
  # ---------------------------------------------------------
  val_prob <- predict(model, newdata = val_scaled, type = "prob")[, "Y"]
  val_pred_opt <- ifelse(val_prob > opt_th, "Y", "N")
  
  cat("\nVALIDATION confusion matrix (Optimized threshold)\n")
  print(confusionMatrix(
    factor(val_pred_opt, levels=c("N","Y")),
    val_scaled$Cat1,
    positive="Y"
  ))
}
