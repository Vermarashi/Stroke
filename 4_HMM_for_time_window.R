###############################################################################
###############################################################################
# Title: Stroke Transcriptomics Analysis: HMM-based Gene Selection & Validation
# Ischemic Stroke - With in Time Window vs. Out of time Window
# Author: Rashi Verma
# Date: 2025-11-14
# Description: 
#   - Load phenotype and RNA-seq count data
#   - Filter lowly-expressed genes
#   - CPM normalization
#   - Subset highly correlated transcripts to time window (continuous)
#   - Batch correction for validation data
#   - Train Hidden Markov Models (HMMs) per gene and multi-gene
#   - Evaluate on independent validation
#   - Plot expression, ROC, heatmaps, and summary statistics
###############################################################################
###############################################################################
# Load packages
library(lubridate)
library(dplyr)
library(hms)
library(readr)
library(corrplot)
library(caret)
library(glmnet)
library(stringr)
library(plyr)
library(dplyr)
library(skimr)
library(caretEnsemble)
library(reshape2)
library(pROC)
#library(caretSDM)
library(rminer)
library(MLeval)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(depmixS4)
library(inflection)
library(gridExtra)
library(colorRamp2)
library(RColorBrewer)
library(statmod)
library(sva)
library(MatchIt)
library(tidyr)
library(ggplot2)
library(limma)
library(edgeR)
library(MLmetrics)
library(ggfortify)
library(ComplexHeatmap)
library(grid)
library(patchwork)  
library(boot)
library(ggpubr)
library(Metrics)  
library(verification)
library(biomaRt)


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
# filter
batch1_samples_clean <- subset(batch1_samples, clinical_diagnosis == "AIS")
batch1_samples_clean <- subset(batch1_samples_clean, time_window != "" & !is.na(time_window)) # (S0408)
batch1_samples_clean <- subset(batch1_samples_clean, !(samples %in% c("S0145", "S0227", "S0247", "S0580", "S0613", "S0216", "S0355", "S0423"))) # HC conversion, second or multiple stroke
#batch1_samples_clean$severity_score <- as.numeric(batch1_samples_clean$severity_score) 216,355,370
#batch1_samples_clean <- subset(batch1_samples_clean, !(treatment == "NONE" & final_diagnosis == "Stroke" & severity_score < 5))
batch1_samples_clean <- subset(batch1_samples_clean, !(samples %in% c("S0001", "S0092", "S0178", "S0213", "S0331", "S0629", "S0657"))) # low read depth (1M)
row.names(batch1_samples_clean) <- batch1_samples_clean$samples
pheno_train <- batch1_samples_clean

## ---- Batch2 + Batch3: Validation ----
batch2_names <- read.csv("batch2.csv", stringsAsFactors = FALSE)
batch3_names <- read.csv("batch3.csv", stringsAsFactors = FALSE)

pheno_val <- subset(pheno_all, samples %in% c(batch2_names$sample_id, batch3_names$sample_id)) %>%
  filter(time_window != "" & !is.na(time_window))

pheno_val <- subset(pheno_val, clinical_diagnosis == "AIS")

pheno_val <- pheno_val %>%
  mutate(time_window_clean = trimws(time_window)) %>%
  mutate(time_window_clean = ifelse(grepl("^\\d{1,4}:\\d{2}:\\d{2}$", time_window_clean), time_window_clean, NA))
# Parse time safely
pheno_val <- pheno_val %>%  mutate(time_parsed = lubridate::hms(time_window_clean),
                                   time_flag = if_else(as.numeric(time_parsed, "hours") >= 3.5, 1, 0))

# Remove known problematic samples
exclude_samples <- c("S0675", "S0699","S0702","S0811","S0866","S0964","S0791","S0802", "S0873",
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

## CPM transformation 
to_cpm <- function(df) {
  counts <- as.matrix(df)
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  return(cpm)
}

cpm_train <- to_cpm(df_train)
cpm_val2  <- to_cpm(df_val2)
cpm_val3  <- to_cpm(df_val3)

cpm_train <- t(cpm_train)
cpm_val2  <- t(cpm_val2)
cpm_val3  <- t(cpm_val3)

# ==============================================================================
# --------------- Subset correlated transcripts to time window -----------------
# ==============================================================================
trainData_counts<- cpm_train
pheno_train$time_hours <- as.numeric(pheno_train$time_parsed, units = "hours")

common_samples <- intersect(rownames(cpm_train), pheno_train$samples)
## Subset counts matrix
cpm_train <- cpm_train[common_samples, ]

## Reorder phenotype data to match (optional)
pheno_train <- pheno_train[match(common_samples, pheno_train$samples), ]
stopifnot(all(rownames(cpm_train) == pheno_train$samples))

## Perform correlation
correlation_results <- apply(t(cpm_train), 1, function(x) {
  cor.test(x, pheno_train$time_hours, method = "pearson")
})

## Extract the correlation coefficients and p-values
correlations <- sapply(correlation_results, function(res) res$estimate)
p_values <- sapply(correlation_results, function(res) res$p.value)

adjusted_p_values_bh <- p.adjust(p_values, method = "BH")

results <- data.frame(
  correlation = correlations,
  p_value = p_values,
  p_value_adj_bh = adjusted_p_values_bh
)

head(results)
results
summary(results$correlation)

## subset to absolute correlation > 0.6
results[1:4,]
corgenes<-results[abs(results$correlation)>0.6 ,]  # 254
corgenes<-corgenes[abs(corgenes$p_value_adj_bh)<0.05 ,]  # 44
corgenes2<-corgenes[complete.cases(corgenes),] # removes rows with missing data
corgenes2$correlation
corgenes2$genes <- sub("\\.cor$", "", rownames(corgenes2)) # remove .cor from row names and add a new column as gene
corgenes2 
dim(corgenes2)

correlation<-data.frame(cor=corgenes2$correlation,genes=corgenes2$genes)
glist<-correlation$genes

## subset correlated transcripts to only genes in glist first
genes_to_use <- glist

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

mod_val <- model.matrix(~ time_flag, data = pheno_val) #Including final_diagnosis in the model preserves biological differences
expr_val_corrected <- ComBat(dat = t(expr_val), batch = pheno_val$batch, mod = mod_val, par.prior = TRUE, prior.plots = FALSE) # batch1 or training set as reference cause over-correction wrongly preserve biological signal
expr_val_corrected <- t(expr_val_corrected)

# ==============================================================================
# --------------- Prepare training and validation data sets --------------------
# ==============================================================================
train_df <- expr_train
train_df <- as.data.frame(train_df)
train_df <- train_df[rownames(train_df) %in% rownames(pheno_train), ]
threshold <- dhours(3) + dminutes(30)  # or dhours(3.5)
train_df$Cat1 <- as.numeric(pheno_train$time_parsed <= threshold)


val_df <- expr_val_corrected
val_df <- as.data.frame(val_df)
val_df <- val_df[rownames(val_df) %in% rownames(pheno_val), ]
threshold <- dhours(3) + dminutes(30)  # or dhours(3.5)
val_df$Cat1 <- as.numeric(pheno_val$time_parsed <= threshold)

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

## run loop over each gene train model and test it on train data set
rm(confusionMatrix)
confmat_fn <- caret::confusionMatrix

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
    
    roc_train <- pROC::roc(actual, probs, quiet = TRUE)
    best_coords <- suppressWarnings(pROC::coords(roc_train, "best", best.method = "youden", transpose = FALSE))
    best_threshold <- as.numeric(best_coords["threshold"])
    if (is.na(best_threshold)) best_threshold <- 0.5
    
    pred_train <- ifelse(probs >= best_threshold, 1, 0)
    if (length(unique(pred_train)) < 2) stop("Single-class predictions.")
    
    # Explicitly named arguments
    cm_train <- confmat_fn(
      data = factor(pred_train, levels = c(0,1)),
      reference = factor(actual, levels = c(0,1)),
      positive = "1"
    )
    
    train_sens <- cm_train$byClass["Sensitivity"]
    train_spec <- cm_train$byClass["Specificity"]
    train_J <- as.numeric(train_sens + train_spec - 1)
    
    results_list[[g]] <- data.frame(
      Gene = g,
      Train_Accuracy = cm_train$overall["Accuracy"],
      Train_AUC = as.numeric(pROC::auc(roc_train)),
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

# ==============================================================================
# --------- Systematic Screening Across Combinations (min-2; max 4) ------------
# ==============================================================================
## loop provide an idea of promising transcripts combinations, but the performance estimates 
## are not fully reliable because the loop can occasionally reuse scaling parameters 
## from previous iterations. Therefore, each selected combination must be rerun 
## independently with fresh pre-processing to obtain accurate and unbiased results.
# Select top features with perfect accuracy
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

#filter <- subset(perfect_df, Train_Accuracy > 0.95 & Val_Accuracy > 0.95)
## Save results
write.csv(perfect_df, "plots/panel_TW.csv", row.names = FALSE)

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
## whereas four genes are sufficient to classify ischemic patients as with in window 
## and out of window, enabling timely administration of tPA
best_k <- 2
perfect_genes <- c("ENST00000531913", "MSTRG_11018_26", "MSTRG_19969_3", "MSTRG_54881_4")
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

## Plot Fig.6A
top_genes <- correlation[correlation$genes %in% common_genes, ]

tiff(file="plots/Fig_6A.tiff", unit= "in", res = 600, width = 2.3, height = 2.3)
p <- ggplot(top_genes, aes(x = reorder(genes, cor), y = cor, color = cor > 0)) +
  geom_point(size = 2) +
  coord_flip() +
  scale_color_manual(values = c("firebrick", "steelblue")) +
  labs(
    title = "Top 20 Genes Correlated with Stroke Onset Time",
    x = "Transcript",
    y = "Pearson Correlation Coefficient"
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )
dev.off()
print(p)
# To save the plot
#ggsave("plots/Fig_6A.png", plot = p, width = 3.5, height = 2.7, dpi = 600)

## Plot Fig. 6B
## GO enrichment of correlated genes
##save unique transcript list of CAT1, 2 and 3 and retrieve ensemble gene name for ENST Ids
glist <- top_genes #common_genes

## Filter out MSTRG IDs
glist_clean <- glist[!grepl("^MSTR", glist)]

##Keep only valid gene symbols
glist_symbols <- unique(glist_clean)
#glist_symbols <- glist_symbols[glist_symbols != ""]  # remove empty symbols
glist_symbols <- glist_symbols[!is.na(glist_symbols) & glist_symbols != ""]

## Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Map transcripts to gene info
tx2gene <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_transcript_id",
  values = glist_symbols,
  mart = mart
)

head(tx2gene)

## Run GO enrichment
go_enrichment <- enrichGO(
  gene          = tx2gene$hgnc_symbol,
  keyType       = "SYMBOL",
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(go_enrichment@result)

go_df <- as.data.frame(go_enrichment)
dim(go_df)

## Calculate fold enrichment using GeneRatio and BgRatio
go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) eval(parse(text=x)))  # Convert GeneRatio to numeric
go_df$BgRatio_num <- sapply(go_df$BgRatio, function(x) eval(parse(text=x)))  # Convert BgRatio to numeric
go_df$FoldEnrichment <- go_df$GeneRatio_num / go_df$BgRatio_num  # Calculate fold enrichment
go_df$GO_type<-go_df$ONTOLOGY
go_df_sorted <- go_df[order(go_df$FoldEnrichment),  ]  # Sort in descending order

## Extract the top 18 GO pathways based on FoldEnrichment
top_20_go_pathways <- head(go_df_sorted, 50)

## Print the top 18 pathways
print(top_20_go_pathways)

## Convert p.adjust to -log10 scale for better visualization
top_20_go_pathways$log10_padj <- -log10(top_20_go_pathways$p.adjust)

## Ensure GO type is a factor (optional, useful if you later want to facet)
top_20_go_pathways$GO_type <- factor(top_20_go_pathways$GO_type, levels = c("BP", "CC", "MF"))

## top_20_bp_pathways <- subset(top_20_go_pathways, GO_type == "BP")
top_20_go_pathways$Description_wrapped <- str_wrap(top_20_go_pathways$Description, width = 80)

## Create the dotplot without splitting
tiff(file="plots/Fig_6B.tiff", unit= "in", res = 600, width = 2.3, height = 2.3)
p <- ggplot(top_20_go_pathways, 
            aes(x = FoldEnrichment, 
                y = reorder(Description_wrapped, FoldEnrichment), 
                size = Count, 
                color = log10_padj)) +
  geom_point() +
  scale_color_gradientn(colors = c("#D94F4F", "#B080A0", "#3C68C8")) +
  scale_size(range = c(1.5, 3.5)) +  # Adjust point sizes
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),#face = "bold"
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  
  labs(
    x = "Fold Enrichment",
    y = "GO Term",
    color = "-log10(p.adjust)",
    size = "Gene Count"
  )

ggsave("plots/Fig_6B.tiff", plot = p, width = 5.0, height = 7, dpi = 600)
print(p)
#dev.off()

## Plot Fig.6C
train_state_prob_df <- data.frame(
  Prob = train_probs,
  TimeGroup = factor(
    train_scaled$Cat1,
    levels = c(0, 1),
    labels = c("> 3.5 h", "≤ 3.5 h")
  )
)

p<-ggplot(train_state_prob_df, aes(x = TimeGroup, y = Prob, fill = TimeGroup)) +
  
  ## Violin plot
  geom_violin(
    alpha = 0.7,
    trim = FALSE,
    color = NA
  ) +
  
  ## Boxplot
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  
  ## Jittered points
  geom_jitter(
    aes(color = TimeGroup),
    width = 0.12,
    alpha = 0.4,
    size = 1
  ) +
  
  ## Classification threshold
  geom_hline(
    yintercept = best_threshold,
    linetype = "dashed",
    linewidth = 0.9,
    color = "black"
  ) +
  
  ## Threshold annotation
  annotate(
    "text",
    x = 1.5,
    y = best_threshold + 0.05,
    label = paste0("Threshold = ", round(best_threshold, 2)),
    size = 4,
    color = "black"
  ) +
  
  ## Fixed x-axis labels with sample sizes
  scale_x_discrete(
    labels = c(
      "> 3.5 h\n(n = 15)",
      "≤ 3.5 h\n(n = 15)"
    )
  ) +
  
  ## Force y-axis from 0 to 1
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  
  scale_fill_manual(
    values = c("> 3.5 h" = "#56B4E9", "≤ 3.5 h" = "#E69F00")
  ) +
  scale_color_manual(
    values = c("> 3.5 h" = "#56B4E9", "≤ 3.5 h" = "#E69F00")
  ) +
  
  labs(
    x = "Time to Event",
    y = paste0("Posterior P(State ", train_target_state, ")"),
    #title = "Posterior Probability of High-Risk Hidden State"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 11, color = "black",
                                margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, color = "black",
                                margin = margin(r = 8)),
    axis.text = element_text(color = "black")
  )
ggsave("plots/Fig_6D.tiff", plot = p, width = 3.232, height = 3.159, dpi = 600)


## Plot Fig.6D
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

# Plot ROC curves
tiff(file="plots/Fig_6D", unit= "in", res = 600, width = 2.307, height = 2.239)
p<-ggplot(roc_df, aes(x = FPR, y = TPR, color = Dataset, linetype = Dataset)) +
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
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "none"
  )

dev.off()
ggsave("plots/Fig_6D.tiff", plot = p, width = 3.195, height = 3.095, dpi = 600)

## Compute AUC
auc_train <- pROC::auc(roc_train_plot)
auc_val   <- pROC::auc(roc_val_plot)

cat("AUC (Train):", round(auc_train, 3), "\n")
cat("AUC (Validation):", round(auc_val, 3), "\n")

#################################################################################
## Additional Confusion matrix Plot
val_actual <- ifelse(val_scaled$Cat1 == 1, 1, 0)
val_states <- val_post$state
val_pred <- ifelse(val_probs >= best_threshold, 1, 0)

val_plot_df <- data.frame(
  actual = val_actual,
  predicted = val_pred,
  state = factor(val_states)
)

tiff(file="plots/TimeWindow_HMM_CM.tiff", unit= "in", res = 600, width = 2.3, height = 2.3)
set.seed(1234)
ggplot(val_plot_df, aes(x = actual, y = predicted)) +
  geom_jitter(aes(color = factor(actual)), width = 0.07, height = 0.07, size = 2, alpha = 0.8) +
  geom_vline(xintercept = best_threshold, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = best_threshold, linetype = "dashed", color = "red", size = 0.8) +
  #annotate("text", x = best_threshold + 0.05, y = 1.05, 
  #label = paste0("FP threshold: ", round(best_threshold, 3)),
  #size = 3, color = "black") +
  #annotate("text", x = 1.05, y = best_threshold - 0.05, 
  #label = paste0("FN threshold: ", round(best_threshold, 3)),
  #size = 3, color = "black", angle = 90) +
  labs(
    #title = "G) Actual vs Predicted Cat1 (Binary) at Best Threshold",
    #subtitle = paste("Threshold =", round(best_threshold, 3)),
    x = "Actual",
    y = "Predicted"
  ) +
  scale_color_manual(
    values = c("0" = "red", "1" = "blue"),
    name = "Stroke Type",
    labels = c("0" = "Stroke Mimic", "1" = "Stroke") # <---- this line changes legend labels
  ) +
  theme_minimal() + 
  my_clean_theme +
  theme(legend.position = "none")

dev.off()