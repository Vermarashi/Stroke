# ==============================================================================
# Loop over predefined gene combinations and evaluate depmix models
#
# For each gene set listed in panel_CAT1.csv, this loop:
#   1) Checks that all genes are present in the training data and have non-zero variance
#   2) Scales training and validation data using training-set parameters
#   3) Fits a binomial depmixS4 model with a fixed number of hidden states (best_k)
#   4) Identifies the target hidden state most associated with Cat1 = 1
#   5) Determines the optimal classification threshold on the training set via Youden’s index
#   6) Evaluates model performance on both training and validation data
#   7) Appends confusion matrices, thresholds, and gene lists to a single output file
#
# Failed model fits are logged without interrupting the loop.
# ==============================================================================
# Read Excel combos
library(readxl)
combo_df <- read.csv("plot/panel_CAT1.csv")
combo_list <- strsplit(combo_df[[1]], ",")

# Define ONE output file
out_file <- "plot/panel4_CAT1.txt"
if (file.exists(out_file)) file.remove(out_file)

# Loop and APPEND results
library(caret)
library(depmixS4)
library(pROC)

for (i in seq_along(combo_list)) {
  
  glist_top40 <- combo_list[[i]]
  cat("\nRunning combo", i, "\n")
  
  if (!all(glist_top40 %in% colnames(train_df))) next
  if (any(apply(train_df[, glist_top40], 2, var, na.rm = TRUE) == 0)) next
  
  tryCatch({
    
    # ==============================================================================
    # Scaling
    # ==============================================================================
    scaler <- preProcess(train_df[, glist_top40], method = c("center","scale"))
    train_scaled <- predict(scaler, train_df[, glist_top40])
    train_scaled$Cat1 <- train_df$Cat1
    
    val_scaled <- predict(scaler, val_df[, glist_top40])
    val_scaled$Cat1 <- val_df$Cat1
    
    # ==============================================================================
    # Model
    # ==============================================================================
    set.seed(1234)
    formula_obj <- as.formula(paste("Cat1 ~", paste(glist_top40, collapse = " + ")))
    
    train_model <- depmix(response = formula_obj,
                          data = train_scaled,
                          nstates = best_k,
                          family = binomial())
    
    fit_train_model <- fit(train_model, verbose = FALSE)
    
    # ==============================================================================
    # Training evaluation
    # ==============================================================================
    train_post <- posterior(fit_train_model, type = "viterbi")
    train_states <- train_post$state
    train_target_state <- which.max(tapply(train_scaled$Cat1, train_states, mean))
    train_probs <- train_post[, paste0("S", train_target_state)]
    
    roc_train <- roc(train_scaled$Cat1, train_probs)
    best_threshold <- coords(
      roc_train, "best", best.method = "youden", transpose = FALSE
    )$threshold
    
    train_pred <- ifelse(train_probs >= best_threshold, 1, 0)
    train_confmat <- confusionMatrix(
      factor(train_pred, levels = c(0,1)),
      factor(train_scaled$Cat1, levels = c(0,1)),
      positive = "1"
    )
    
    # ==============================================================================
    # Validation evaluation
    # ==============================================================================
    val_model <- depmix(response = formula_obj,
                        data = val_scaled,
                        nstates = best_k,
                        family = binomial())
    
    val_model <- setpars(val_model, getpars(fit_train_model))
    val_post <- posterior(val_model, type = "viterbi")
    val_probs <- val_post[, paste0("S", train_target_state)]
    
    val_pred <- ifelse(val_probs >= best_threshold, 1, 0)
    val_confmat <- confusionMatrix(
      factor(val_pred, levels = c(0,1)),
      factor(val_scaled$Cat1, levels = c(0,1)),
      positive = "1"
    )
    
    # ==============================================================================
    # APPEND results to ONE TXT
    # ==============================================================================
    sink(out_file, append = TRUE)
    
    cat("\n============================================================\n")
    cat("Combo:", i, "\n")
    cat("Genes:\n", paste(glist_top40, collapse = ", "), "\n\n")
    
    cat("Best threshold:", best_threshold, "\n\n")
    
    cat("Training Confusion Matrix:\n")
    print(train_confmat)
    
    cat("\nValidation Confusion Matrix:\n")
    print(val_confmat)
    
    sink()
    
  }, error = function(e) {
    sink(out_file, append = TRUE)
    cat("\nCombo", i, "FAILED:", e$message, "\n")
    sink()
  })
}

