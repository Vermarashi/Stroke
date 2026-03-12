###############################################################################
###############################################################################
# Title: Stroke Transcriptomics Analysis: HMM-based Gene Selection & Validation
# Ischemic Stroke vs. Stroke Mimic
# Author: Rashi Verma
# Date: 2025-11-13
# Description: 
#   - Load phenotype and RNA-seq count data
#   - Filter lowly-expressed genes
#   - CPM normalization
#   - Subset genes to DETs
#   - Batch correction for validation data
#   - Train Hidden Markov Models (HMMs) per gene and multi-gene
#   - Evaluate on independent validation
#   - Plot expression, ROC, heatmaps, and summary statistics
###############################################################################
###############################################################################

# Load packages
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
batch1_samples_clean <- subset(batch1_samples_clean, final_diagnosis != "Hemorrhage")
row.names(batch1_samples_clean) <- batch1_samples_clean$samples
pheno_train <- batch1_samples_clean

## ---- Batch2 + Batch3: Validation ----
batch2_names <- read.csv("batch2.csv", stringsAsFactors = FALSE)
batch3_names <- read.csv("batch3.csv", stringsAsFactors = FALSE)

pheno_val <- subset(pheno_all, samples %in% c(batch2_names$sample_id, batch3_names$sample_id)) %>%
  filter(final_diagnosis != "" & !is.na(final_diagnosis))

## Excluded samples (Reasons for exclusion are detailed in pheno_data_345.csv)
exclude_samples <- c("S0675","S0699","S0702","S0807","S0811","S0866","S0964","S0791","S0802", "S0873",
                     "S0815","S0826","S0838","S0839","S0844","S0849","S0852","S0854","S0865","S0896","S0912","S0948","S0777","S0198")

pheno_val <- pheno_val[!pheno_val$samples %in% exclude_samples, ]
pheno_val <- subset(pheno_val, final_diagnosis != "Hemorrhage")
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
filter_low_counts <- function(df, min_count = 10, min_samples_prop = 0.10) {
  min_samples <- ceiling(ncol(df) * min_samples_prop)
  keep <- rowSums(df >= min_count) >= min_samples
  df_filtered <- df[keep, ]
  message("Genes retained: ", sum(keep), " / ", nrow(df))
  return(df_filtered)
}

df_train <- filter_low_counts(df_train, min_count = 10, min_samples_prop = 0.10)
df_val2  <- filter_low_counts(df_val2,  min_count = 10, min_samples_prop = 0.10)
df_val3  <- filter_low_counts(df_val3,  min_count = 10, min_samples_prop = 0.10)


## CPM transformation
to_cpm <- function(df) {
  counts <- as.matrix(df)
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  return(cpm)
}

cpm_train <- to_cpm(df_train)
cpm_val2 <- to_cpm(df_val2)
cpm_val3 <- to_cpm(df_val3)

cpm_train <- t(cpm_train)
cpm_val2  <- t(cpm_val2)
cpm_val3  <- t(cpm_val3)

# ==============================================================================
# ------------Subset genes of interest from differential expression-------------
# ==============================================================================
genes1=read.csv("results/1_HS_vs_TPA_50_t.csv")
genes1[1:20,1:5]; dim(genes1) # 100x7

genes2=read.csv("results/2_TPA_vs_SM_50_t.csv")
genes2[1:20,1:5]; dim(genes2) # 100x7

genes3=read.csv("results/3_HS_vs_SM_50_t.csv")
genes3[1:20,1:5]; dim(genes3) # 100x7

## Combine all genes into one vector
g<- setdiff(genes2$X, genes3$X)
glist <-setdiff(g, genes1$X)
length(glist) # 80 gene IDs. # 82 gene IDs

glist  <- gsub("[.]","_", glist)
glist  <- gsub("[|]","_", glist)
glist  <- gsub("[-]","_", glist)

genes_to_use <- glist

## Subset to only genes in glist first
expr_train <- cpm_train[, colnames(cpm_train) %in% genes_to_use, drop = FALSE]
expr_val2  <- cpm_val2[,  colnames(cpm_val2)  %in% genes_to_use, drop = FALSE]
expr_val3  <- cpm_val3[,  colnames(cpm_val3)  %in% genes_to_use, drop = FALSE]

## Keep only genes common to ALL three datasets
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

mod_val <- model.matrix(~ final_diagnosis, data = pheno_val)
expr_val_corrected <- ComBat(dat = t(expr_val), batch = pheno_val$batch, mod = mod_val, par.prior = TRUE, prior.plots = FALSE)
expr_val_corrected <- t(expr_val_corrected)

# ==============================================================================
# --------------- Prepare training and validation data sets --------------------
# ==============================================================================
train_df <- expr_train
train_df <- as.data.frame(train_df)
train_df <- train_df[rownames(train_df) %in% rownames(pheno_train), ]
train_df$Cat1 <- as.numeric(pheno_train$final_diagnosis == "Stroke")


val_df <- expr_val_corrected
val_df <- as.data.frame(val_df)
val_df <- val_df[rownames(val_df) %in% rownames(pheno_val), ]
val_df$Cat1 <- as.numeric(pheno_val$final_diagnosis == "Stroke")

# ==============================================================================
# ---------- Sanity check between pheno data and expression data ---------------
# ==============================================================================
identical(rownames(train_df), pheno_train$samples)
identical(rownames(val_df), pheno_val$samples)

# ==============================================================================
# ---------------------------- Scale data sets ---------------------------------
# ==============================================================================
genes_to_use <- common_genes
scaler <- preProcess(train_df[, genes_to_use], method = c("center","scale"))
train_scaled <- predict(scaler, train_df[, genes_to_use])
train_scaled$Cat1 <- train_df$Cat1

val_scaled <- predict(scaler, val_df[, genes_to_use])
val_scaled$Cat1 <- val_df$Cat1

# ==============================================================================
# ----------------- Per-gene HMM analysis for feature selection ----------------
# ==============================================================================
set.seed(1234)

gene_cols <- genes_to_use
formula_obj <- as.formula(paste("Cat1 ~", paste(gene_cols, collapse = " + ")))

## Choose best number of states (2-4)
candidate_states <- 2:4
model_info <- data.frame(States = integer(), AIC = numeric(), BIC = numeric())

for (nst in candidate_states){
  mod <- depmix(response = formula_obj, data = train_scaled, nstates = nst, family = binomial())
  fit_mod <- try(fit(mod, verbose = FALSE), silent = TRUE)
  if(inherits(fit_mod, "try-error")) next
  model_info <- rbind(model_info,
                      data.frame(States = nst, AIC = AIC(fit_mod), BIC = BIC(fit_mod)))
}

best_k <- model_info$States[which.min(model_info$BIC)]
cat("Optimal number of states (BIC):", best_k, "\n")

results_list <- list()

for (g in gene_cols) {
  cat("Processing:", g, "\n")
  
  tryCatch({
    formula_obj <- as.formula(paste("Cat1 ~", g))
    
    model <- depmix(response = formula_obj, data = train_scaled, nstates = best_k, family = binomial())
    fit_model <- fit(model, verbose = FALSE)
    
    post_train <- posterior(fit_model, type = "viterbi")
    states <- post_train$state
    actual <- train_scaled$Cat1
    
    if (length(unique(states)) < 2) stop("Only one state detected.")
    
    state_mean <- tapply(actual, states, mean)
    target_state <- which.max(state_mean)
    probs <- post_train[, paste0("S", target_state)]
    
    roc_train <- roc(actual, probs, quiet = TRUE)
    best_coords <- suppressWarnings(coords(roc_train, "best", best.method = "youden", transpose = FALSE))
    best_threshold <- as.numeric(best_coords["threshold"])
    if (is.na(best_threshold)) best_threshold <- 0.5
    
    pred_train <- ifelse(probs >= best_threshold, 1, 0)
    if (length(unique(pred_train)) < 2) stop("Single-class predictions.")
    
    cm_train <- confusionMatrix(factor(pred_train, levels = c(0,1)), factor(actual, levels = c(0,1)), positive = "1")
    
    train_sens <- cm_train$byClass["Sensitivity"]
    train_spec <- cm_train$byClass["Specificity"]
    train_J <- as.numeric(train_sens + train_spec - 1)
    
    results_list[[g]] <- data.frame(
      Gene = g,
      Train_Accuracy = cm_train$overall["Accuracy"],
      Train_AUC = as.numeric(auc(roc_train)),
      Train_Sensitivity = train_sens,
      Train_Specificity = train_spec,
      Train_YoudenJ = train_J,
      Threshold_Used = best_threshold,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    cat("⚠️ Skipping", g, "due to error:", e$message, "\n")
  })
}

results_df <- do.call(rbind, results_list)
results_df <- results_df[order(-results_df$Train_YoudenJ), ]

## Select top features with perfect accuracy
perfect_genes <- subset(results_df, Train_Specificity == 1 & Train_Sensitivity == 1)

# ==============================================================================
# --------- Systematic Screening Across Combinations (min-2; max 4) ------------
# ==============================================================================
## loop provide an idea of promising transcripts combinations, but the performance estimates 
## are not fully reliable because the loop can occasionally reuse scaling parameters 
## from previous iterations. Therefore, each selected combination must be rerun 
## independently with fresh pre-processing to obtain accurate and unbiased results.

perfect_genes <- subset(results_df, Train_Specificity == 1 & Train_Sensitivity == 1)$Gene

## Initialize results data frame
perfect_df <- data.frame(
  GeneCombo = character(),
  Train_Accuracy = numeric(),
  Train_Sensitivity = numeric(),
  Train_Specificity = numeric(),
  Train_YoudenJ = numeric(),
  Val_Accuracy = numeric(),
  Val_Sensitivity = numeric(),
  Val_Specificity = numeric(),
  Val_YoudenJ = numeric(),
  stringsAsFactors = FALSE
)

best_k <- 2  # chosen number of hidden states

## Generate combinations
gene_combos <- unlist(
  lapply(2:5, function(k) {
    combn(perfect_genes, k, simplify = FALSE)
  }),
  recursive = FALSE
)

cat("Generating combinations of 2 to 5 genes →",
    length(gene_combos), "total\n")


## Loop through each combo
for (glist in gene_combos) {
  glist <- unlist(glist)  # ensure character vector
  
  tryCatch({
    
    # Skip if any gene has zero variance
    if(any(apply(train_df[, glist], 2, var, na.rm = TRUE) == 0)) next
    
    ## ---------------- Scaling ----------------
    scaler <- preProcess(train_df[, glist], method = c("center", "scale"))
    train_scaled <- predict(scaler, train_df[, glist])
    train_scaled$Cat1 <- train_df$Cat1
    
    val_scaled <- predict(scaler, val_df[, glist])
    val_scaled$Cat1 <- val_df$Cat1
    
    ## ---------------- Model Formula ----------------
    formula_obj <- as.formula(paste("Cat1 ~", paste(glist, collapse = " + ")))
    
    ## ---------------- Train HMM ----------------
    fit_train_model <- tryCatch({
      train_model <- depmix(response = formula_obj, data = train_scaled, 
                            nstates = best_k, family = binomial())
      fit(train_model, verbose = FALSE)
    }, error = function(e) {
      message("⚠️ Skipping combo (train fit error): ", paste(glist, collapse = ","))
      return(NULL)
    })
    
    if(is.null(fit_train_model)) next
    
    ## ---------------- Posterior on Train ----------------
    train_post <- posterior(fit_train_model, type = "viterbi")
    if(any(is.na(train_post))) next
    
    train_states <- train_post$state
    ## Skip degenerate states
    if(any(table(train_states, train_scaled$Cat1) == 0)) next
    
    train_target_state <- which.max(tapply(train_scaled$Cat1, train_states, mean))
    train_probs <- train_post[, paste0("S", train_target_state)]
    
    roc_train <- roc(train_scaled$Cat1, train_probs)
    best_threshold <- coords(roc_train, "best", best.method = "youden", transpose = FALSE)$threshold
    if(is.na(best_threshold)) best_threshold <- 0.5
    
    train_pred <- ifelse(train_probs >= best_threshold, 1, 0)
    cm_train <- confusionMatrix(
      factor(train_pred, levels=c(0,1)),
      factor(train_scaled$Cat1, levels=c(0,1)),
      positive="1"
    )
    youden_train <- cm_train$byClass["Sensitivity"] + cm_train$byClass["Specificity"] - 1
    
    ## ---------------- Validation ----------------
    val_model <- depmix(response = formula_obj, data=val_scaled, nstates=best_k, family=binomial())
    val_model <- setpars(val_model, getpars(fit_train_model))
    
    val_post <- tryCatch({
      posterior(val_model, type="viterbi")
    }, error = function(e) {
      message("⚠️ Skipping combo (validation error): ", paste(glist, collapse = ","))
      return(NULL)
    })
    if(is.null(val_post) || any(is.na(val_post))) next
    
    val_states <- val_post$state
    if(any(table(val_states, val_scaled$Cat1) == 0)) next  # skip degenerate
    
    val_target_state <- which.max(tapply(val_scaled$Cat1, val_states, mean))
    val_probs <- val_post[, paste0("S", train_target_state)]
    val_pred <- ifelse(val_probs >= best_threshold, 1, 0)
    
    cm_val <- confusionMatrix(
      factor(val_pred, levels=c(0,1)),
      factor(val_scaled$Cat1, levels=c(0,1)),
      positive="1"
    )
    youden_val <- cm_val$byClass["Sensitivity"] + cm_val$byClass["Specificity"] - 1
    
    ## ---------------- Save results ----------------
    perfect_df <- rbind(perfect_df, data.frame(
      GeneCombo = paste(glist, collapse=","),
      Train_Accuracy = cm_train$overall["Accuracy"],
      Train_Sensitivity = cm_train$byClass["Sensitivity"],
      Train_Specificity = cm_train$byClass["Specificity"],
      Train_YoudenJ = youden_train,
      Val_Accuracy = cm_val$overall["Accuracy"],
      Val_Sensitivity = cm_val$byClass["Sensitivity"],
      Val_Specificity = cm_val$byClass["Specificity"],
      Val_YoudenJ = youden_val,
      stringsAsFactors = FALSE
    ))
    
    cat("Completed combo:", paste(glist, collapse=","), "\n")
    
  }, error=function(e){
    message("⚠️ Skipping combo (other error): ", paste(glist, collapse=","), " -> ", e$message)
  })
  
} # end loop

## Save results
write.csv(perfect_df, "plots/panel_CAT2.csv", row.names = FALSE)

# ==============================================================================
# Run Multi-transcript (Best combination) HMM on training and validation data set 
# ==============================================================================
## To ensure that the predictive accuracy of each gene set is evaluated without 
## any influence from prior iterations, it is necessary to re-run each promising 
## combination independently. Conducting a fresh scaling step and model initialization 
## for every selected combination guarantees that the resulting classification metrics 
## reflect the intrinsic behavior of that specific gene set. This independent evaluation 
## step removes any residual computational artifacts and ensures that the final reported 
## performance is accurate, reproducible, and methodologically sound.

## We selected 4-gene combinations because they provide a reliable, actionable 
## signature that can be quickly quantified in an emergency setting. Models with 
## 1–3 genes lack sufficient signal integration, leading to unstable thresholds, 
## while using more genes would increase complexity without substantial benefit, 
## whereas four genes are sufficient to classify patients as ischemic or 
## hemorrhagic stroke, enabling timely administration of tPA
best_k <- 2
perfect_genes <- c("MSTRG_7936_41", "MSTRG_10294_21", "ENST00000397762", "ENST00000418132")
glist_top40 <- perfect_genes
glist_top40

# ==============================================================================
# ---------------------------- Scale data sets ---------------------------------
# ==============================================================================
scaler <- preProcess(train_df[, glist_top40], method = c("center","scale"))
train_scaled <- predict(scaler, train_df[, glist_top40])
train_scaled$Cat1 <- train_df$Cat1

val_scaled <- predict(scaler, val_df[, glist_top40])
val_scaled$Cat1 <- val_df$Cat1

# ==============================================================================
# -------- Perform Multi-transcripts HMM modelling on training data set --------
# ==============================================================================
set.seed(1234)
formula_obj <- as.formula(paste("Cat1 ~", paste(glist_top40, collapse = " + ")))

train_model <- depmix(response = formula_obj, data = train_scaled, nstates = best_k, family = binomial())
fit_train_model <- fit(train_model, verbose = FALSE)

# ==============================================================================
# ------------- Evaluate models on training and validation data set ------------
# ==============================================================================
## Training data Set
## Posterior on training
train_post <- posterior(fit_train_model, type = "viterbi")
train_states <- train_post$state
train_target_state <- which.max(tapply(train_scaled$Cat1, train_states, mean))
train_probs <- train_post[, paste0("S", train_target_state)]

roc_train <- roc(train_scaled$Cat1, train_probs)
best_threshold <- coords(roc_train, "best", best.method = "youden", transpose = FALSE)$threshold
train_pred <- ifelse(train_probs >= best_threshold, 1, 0)

train_confmat <- confusionMatrix(factor(train_pred, levels=c(0,1)), factor(train_scaled$Cat1, levels=c(0,1)), positive="1")
train_confmat
cat("Best threshold:", best_threshold, "\n")

## Validation data Set
val_model <- depmix(response = formula_obj, data = val_scaled, nstates = best_k, family = binomial())
val_model <- setpars(val_model, getpars(fit_train_model))

val_post <- posterior(val_model, type="viterbi")
val_probs <- val_post[, paste0("S", train_target_state)]

roc_val <- roc(val_scaled$Cat1, val_probs)
val_pred <- ifelse(val_probs >= best_threshold, 1, 0)
val_confmat <- confusionMatrix(factor(val_pred, levels=c(0,1)), factor(val_scaled$Cat1, levels=c(0,1)), positive="1")
val_confmat

## Plot Fig.5A
my_clean_theme <- theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Helper function to make ROC data frame
clean_roc_df <- function(roc_obj, dataset_name) {
  df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Dataset = dataset_name
  ) %>%
    arrange(FPR, TPR) %>%
    distinct(FPR, .keep_all = TRUE)  # ensure unique FPR
  return(df)
}

## Compute ROC objects
roc_train_plot <- roc(train_scaled$Cat1, train_probs)
roc_val_plot   <- roc(val_scaled$Cat1, val_probs)

## Convert to data frames for ggplot
roc_train_df <- clean_roc_df(roc_train_plot, "Train")
roc_val_df   <- clean_roc_df(roc_val_plot, "Validation")

roc_df <- bind_rows(roc_train_df, roc_val_df)
roc_df$Dataset <- factor(roc_df$Dataset, levels = c("Train", "Validation"))

tiff(file="plots/Fig_5A.tiff", unit= "in", res = 600, width = 2.307, height = 2.239)
ggplot(roc_df, aes(x = FPR, y = TPR, color = Dataset, linetype = Dataset)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "solid") +  # center diagonal
  scale_color_manual(values = c("Train" = "blue", "Validation" = "red")) +
  scale_linetype_manual(values = c("Train" = "solid", "Validation" = "dashed")) +
  labs(#title = "B) ROC Curves",
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_minimal() +
  my_clean_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    #legend.position = "none"
  )

dev.off()

## Compute AUC
auc_train <- auc(roc_train_plot)
auc_val   <- auc(roc_val_plot)
cat("AUC (Train):", round(auc_train, 3), "\n")
cat("AUC (Validation):", round(auc_val, 3), "\n")

## Plot Fig.5B
val_actual <- ifelse(val_scaled$Cat1 == 1, 1, 0)
val_states <- val_post$state
val_pred <- ifelse(val_probs >= best_threshold, 1, 0)

val_plot_df <- data.frame(
  actual = val_actual,
  predicted = val_pred,
  state = factor(val_states)
)

tiff(file="plots/Fig_5B.tiff", unit= "in", res = 600, width = 2.3, height = 2.3)
set.seed(1234)
ggplot(val_plot_df, aes(x = actual, y = predicted)) +
  geom_jitter(aes(color = factor(actual)), width = 0.07, height = 0.07, size = 2, alpha = 0.8) +
  geom_vline(xintercept = best_threshold, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = best_threshold, linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = best_threshold + 0.05, y = 1.05, 
           label = paste0("FP threshold: ", round(best_threshold, 3)),
           size = 3, color = "black") +
  annotate("text", x = 1.05, y = best_threshold - 0.05, 
           label = paste0("FN threshold: ", round(best_threshold, 3)),
           size = 3, color = "black", angle = 90) +
  labs(
    title = "G) Actual vs Predicted Cat1 (Binary) at Best Threshold",
    subtitle = paste("Threshold =", round(best_threshold, 3)),
    x = "Actual",
    y = "Predicted"
  ) +
  scale_color_manual(
    values = c("0" = "olivedrab", "1" = "#1f77b4"),
    name = "Stroke Type",
    labels = c("0" = "Stroke Mimic", "1" = "Stroke") 
  ) +
  theme_minimal() + 
  my_clean_theme +
  theme(legend.position = "none")

dev.off()

