###############################################################################
###############################################################################
# Title: Differential Expression Analysis for CAT1: Hemorrhagic Stroke
# Author: Rashi Verma
# Date: 2025-11-13
# Description: 
#   - Load phenotype and RNA-seq count data
#   - Clean and filter phenotype data
#   - Time-window and clinical exclusions
#   - Match hemorrhage vs ischemic stroke samples on age and severity
#   - Filter/normalize counts (edgeR TMM)
#   - Model with voom-limma
#   - identify significantly differentially expressed transcripts with multiple-testing correction and fold-change thresholds.
###############################################################################
###############################################################################

library("statmod")
library("edgeR")
library("limma")
library("dplyr")
library("tibble")
library("ggplot2")
library("cobalt")
library("MatchIt")
library("matrixStats")
library("pheatmap")
library("ggrepel")
library("viridis")
library("lubridate")
library("hms")

# If RESULTS and PLOTS folders don't exist, create them.
if (!dir.exists("plots")) {
  dir.create("plots")
}  
if (!dir.exists("results")) {
  dir.create("results")
}


set.seed(99)

# Train Dataset
pheno_1 <- read.csv("pheno_data_345.csv", header = TRUE, fill = TRUE); typeofdata="Phenotype"
batch1_names <- read.csv("batch1.csv", header = TRUE, stringsAsFactors = FALSE)
batch1_ids <- batch1_names$sample_id
batch1_samples <- subset(pheno_1, samples %in% batch1_ids)
# Clean time strings — remove empty or invalid ones
batch1_samples <- batch1_samples %>%
  mutate(time_window_clean = trimws(time_window)) %>%
  mutate(time_window_clean = ifelse(grepl("^\\d{1,4}:\\d{2}:\\d{2}$", time_window_clean), time_window_clean, NA))
# Parse time safely
batch1_samples <- batch1_samples %>%  mutate(time_parsed = lubridate::hms(time_window_clean))
# filter
batch1_samples_clean <- subset(batch1_samples, clinical_diagnosis != "TIA")
batch1_samples_clean <- subset(batch1_samples_clean, final_diagnosis != "" & !is.na(final_diagnosis)) # (S0408)
batch1_samples_clean <- subset(batch1_samples_clean, treatment != "MT")
batch1_samples_clean <- subset(batch1_samples_clean, samples != "S0455") #unknown diagnosis
batch1_samples_clean <- subset(batch1_samples_clean, samples != "S0085") #confusion TPA/Mimic
batch1_samples_clean <- subset(batch1_samples_clean, !(samples %in% c("S0145", "S0227", "S0247", "S0580", "S0613", "S0423"))) # HC conversion, second or multiple stroke
threshold <- dhours(4) + dminutes(30)  # or dhours(4.5)
batch1_samples_clean <- subset(batch1_samples_clean, !(treatment == "NONE" & final_diagnosis == "Stroke" & time_parsed > threshold))
batch1_samples_clean$severity_score <- as.numeric(batch1_samples_clean$severity_score)
pheno_1_clean <- subset(batch1_samples_clean, !(treatment == "NONE" & final_diagnosis == "Stroke" & severity_score < 5))
targets <- subset(pheno_1_clean, !(samples %in% c("S0001", "S0092", "S0178", "S0213", "S0331", "S0629", "S0657"))) # low read depth (1M)
targets <- targets %>% filter(!is.na(age) & !is.na(severity_score))

dim(targets)
filtered_targets_clean <- targets %>% filter(treatment %in% c("TPA", "TPA+MT") | final_diagnosis == "Hemorrhage")
dim(filtered_targets_clean)

countsTable <- read.csv(file='transcript_count_matrix_batch1.csv')
dim(countsTable)

filtered_gene_counts <- countsTable %>% select(transcript_id, all_of(filtered_targets_clean$samples)) %>% column_to_rownames(var = "transcript_id")
head(filtered_gene_counts)
dim(filtered_gene_counts)

all(colnames(filtered_gene_counts) %in% filtered_targets_clean$samples)
all(colnames(filtered_gene_counts) == filtered_targets_clean$samples)

exp_counts <- filtered_gene_counts[(rowSums(filtered_gene_counts) > 37),]
dim(exp_counts)

filtered_targets_clean$final_diagnosis_binary <- ifelse(filtered_targets_clean$final_diagnosis == "Hemorrhage", 1, 0)
dim(filtered_targets_clean) 

exp_counts_t<-data.frame(t(exp_counts))

exp_counts_t$AGE<-filtered_targets_clean$age
exp_counts_t$SCORE<-filtered_targets_clean$severity_score
exp_counts_t$DIAGNOSIS<-filtered_targets_clean$final_diagnosis_binary
dim(exp_counts_t) 

m.out <- matchit(as.factor(exp_counts_t$DIAGNOSIS) ~ AGE + SCORE, data = exp_counts_t, 
                 method = "full", distance = "mahalanobis")
summary(m.out)
matched_samples <- match.data(m.out)
dim(matched_samples) 

exp_Matched <- matched_samples %>% select(-AGE, -SCORE, -DIAGNOSIS, -weights, -subclass) %>% 
  select(where(is.numeric))
colnames(exp_Matched) <- rownames(exp_counts)
dim(exp_Matched)

count_final <- (t(exp_Matched))

# Filtering: Only once
dge <- DGEList(counts = count_final)
min_samples <- ceiling(ncol(dge) / 3)
keep <- rowSums(cpm(dge) > 10) >= min_samples
dge <- dge[keep, , keep.lib.size = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Match and align metadata
final_targets <- filtered_targets_clean %>% filter(samples %in% colnames(dge))
final_targets <- final_targets[match(colnames(dge), final_targets$samples), ]

group <- as.factor(final_targets$final_diagnosis)
sex <- final_targets$sex
age <- final_targets$age # Assuming age is numeric
race <- final_targets$race
score <- final_targets$severity_score
treatment <- final_targets$treatment

design <- model.matrix(~0 + group + age + sex, data = final_targets)
colnames(design) <- gsub("group", "", colnames(design))

# Voom
v <- voom(dge, design, weights = final_targets$weights, plot = TRUE)

# Linear modeling
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(ISvsHS = Hemorrhage - Stroke, levels = design)
fit <-lmFit(v, design)
# fit to contrast matix, to identify contrast variables
cfit <- contrasts.fit(fit, contrasts=contrast.matrix)
# finally apply Bayesian correction
efit <- eBayes(cfit)
# Mason paper used 1.5 fold changes 
tfit <- treat(efit, lfc=(log2(1.2)))
dt <- decideTests(tfit)
summary(dt)

# DE table
top.table<- topTable(efit, sort.by = "P", n = Inf, adjust.method = "BH")

length(which(top.table$adj.P.Val < 0.05))
length(which(top.table$P.Value < 0.05))

filtered_genes <- top.table[top.table$P.Value < 0.05 & abs(top.table$logFC) > log2(1.5), ]
#write.csv(filtered_genes, "results/1_HS_vs_TPA_all_t.csv", row.names = TRUE)

new_column_name <- "ID"
filtered_genes[, new_column_name] <- rownames(filtered_genes)
str(filtered_genes)

# Define upregulated and downregulated genes
upregulated_genes <- filtered_genes[filtered_genes$logFC > log2(1.5) & filtered_genes$P.Value < 0.05, ]
downregulated_genes <- filtered_genes[filtered_genes$logFC < -log2(1.5) & filtered_genes$P.Value < 0.05, ]

# Sort the upregulated and downregulated genes by P.Value
sorted_upregulated_genes <- upregulated_genes[order(upregulated_genes$P.Value), ]
sorted_downregulated_genes <- downregulated_genes[order(downregulated_genes$P.Value), ]

# Select the top 25 upregulated and 25 downregulated genes
top_upregulated_genes <- head(sorted_upregulated_genes, 50)
top_downregulated_genes <- head(sorted_downregulated_genes, 50)

# Combine the top upregulated and downregulated genes
sorted_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
gene_names <- sorted_genes$ID
#write.csv(sorted_genes, "results/1_HS_vs_TPA_50_t.csv", row.names = TRUE)

# Select the gene expression data for the sorted genes
sorted_gene_expression <- cpm(dge)[sorted_genes$ID, ]

# Reorder the samples (columns) so that "Stroke" samples come first
# Create a vector with the reordered column names
group_vector <- ifelse(final_targets$final_diagnosis == "Stroke", "Stroke", "Hemorrhage")
group_df <- data.frame(Groups = group_vector)
rownames(group_df) <- colnames(sorted_gene_expression)

# Define custom colors for annotation
ann_colors = list(Groups = c("Stroke" = "#0072B2", "Hemorrhage" = "darkorange"))

# Reorder the samples (columns) so that "Stroke" samples come first
# Create a vector with the reordered column names
ordered_samples <- order(group_df$Groups)  # "Stroke" samples come first based on group vector

# Reorder the columns in the gene expression matrix
sorted_gene_expression_ordered <- sorted_gene_expression[, ordered_samples]
group_df_ordered <- group_df[ordered_samples, ]

# Check if the column names of the gene expression matrix match the row names of the annotation
# Ensure column names in sorted_gene_expression_ordered match row names in group_df_ordered
if (!all(colnames(sorted_gene_expression_ordered) == rownames(group_df_ordered))) {
  stop("Column names in the gene expression matrix do not match row names in the annotation data.")
}

# Reorder group_df (annotation) according to the column order in the expression matrix
group_df_ordered <- group_df[match(colnames(sorted_gene_expression_ordered), rownames(group_df)), , drop = FALSE]

tiff(file="plots/1_heatmap_CAT1.tiff", unit= "in", res = 600, width = 2.6, height = 2.5)
# Create the heatmap
pheatmap(
  sorted_gene_expression_ordered,     # Use the sorted gene expression data
  cluster_cols = FALSE,               # Don't cluster columns (samples)
  cluster_rows = FALSE,                # Cluster rows (genes)
  scale = 'row',                      # Scale rows (genes) across samples
  annotation_col = group_df_ordered,   # Add sample group annotations
  annotation_colors = ann_colors,      # Apply custom colors to annotations
  show_colnames = FALSE,               # Show column names (sample IDs)
  show_rownames = FALSE,               # Show row names (gene IDs)
  border_color = NA,             # Set border color to white
  cex.axis = 0.8,                     # Adjust axis text size
  cex.names = 0.8,                    # Adjust row/column names size
  #main = "Hemorrhage vs. Stroke",  # Add title
  fontsize = 6,                      # Adjust overall font size
  fontsize_row = 4,                   # Adjust row name font size
  fontsize_col= 6,                   # Adjust column name font size
  colorRampPalette(c("purple", "white","yellow"))(500),  # More color steps for fine details #FFD700", "white", "#0072B2
  breaks = seq(-1, 1, length.out = 501),  # More granular breaks
  #legend = TRUE,                      # Display the color legend
  #annotation_legend = TRUE,           # Show the annotation legend
  drop_levels = TRUE                  # Drop unused factors from annotations
)

# Perform PCA
pca_data <- prcomp(t(sorted_gene_expression), center = TRUE, scale. = TRUE)

# Calculate variance explained
pca_var <- pca_data$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)  # Convert to percentage

# Extract PCA components
pca_df <- data.frame(
  PC1 = pca_data$x[,1], 
  PC2 = pca_data$x[,2], 
  Sample = colnames(sorted_gene_expression)
)

# Add diagnosis information (replace 'sample_diagnosis' with your actual diagnosis vector)
pca_df$Diagnosis <- final_targets$final_diagnosis  # Vector with "Hemorrhage" or "Stroke" corresponding to samples

# Define shape mapping: Circle (16) for hemorrhage, Diamond (18) for stroke
shape_mapping <- c("Hemorrhage" = 24, "Stroke" = 21)
# Add diagnosis information (replace 'sample_diagnosis' with your actual diagnosis vector)
pca_df$Diagnosis <- final_targets$final_diagnosis  # Vector with "Hemorrhage" or "Stroke" corresponding to samples

# Modify color mapping: Hemorrhage (red), Stroke (blue)
fill_mapping <- c("Hemorrhage" = "darkorange", "Stroke" = "#0072B2")

#tiff(file="plots/1_PCA_CAT1.tiff", unit= "in", res = 600, width = 2.8, height = 2.5)
# Plot PCA with variance percentage in axis labels
ggplot(pca_df, aes(x = PC1, y = PC2, shape = Diagnosis, fill = Diagnosis)) +
  geom_point(size = 2, color = "black", stroke = 0.5) +
  scale_shape_manual(values = shape_mapping) +
  scale_fill_manual(values = fill_mapping) +
  #ggtitle("Hemorrhage vs. Stroke") +
  xlab(paste0("PC1 (", pca_var_perc[1], "% variance)")) +
  ylab(paste0("PC2 (", pca_var_perc[2], "% variance)")) +
  theme_minimal(base_size = 14) +
  theme(
    #legend.position = "NONE",            # remove legend
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.border = element_blank(),
    axis.line = element_line(color = "black"), # keep axis lines
    plot.title = element_text(hjust = 0.5, color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# Volcano no legend
# Define thresholds
fc_threshold <- 1.5
pval_threshold <- 0.05

# Classify significance
top.table <- top.table %>%
  mutate(Significance = case_when(
    P.Value < pval_threshold & abs(logFC) >= fc_threshold ~ "Significant: FC & p-value",
    P.Value < pval_threshold ~ "Significant: p-value",
    abs(logFC) >= fc_threshold ~ "Significant: FC",
    TRUE ~ "Not Significant"
  ))

# Custom colors
custom_colors <- c(
  "Not Significant" = "grey60",
  "Significant: FC" = "darkgreen",
  "Significant: p-value" = "blue4",
  "Significant: FC & p-value" = "red4"
)

# Volcano plot without legend
volcano_plot <- ggplot(top.table, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(size = 0.5, shape = 16) +
  scale_color_manual(values = custom_colors) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Hemorrhage vs. Stroke",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none"  # <- Remove the legend
  )
ggsave("plots/1_volcano_CAT1.tiff", plot = volcano_plot, width = 2.4, height = 2.4, units = "in", dpi = 600)

################################################################################################################
# Supplementary Figure 1
## love plot to prove matchit results
## Covariates
covariates <- c("AGE", "SCORE")

## Std. Mean Differences
smd_before <- c(0.0285, 0.1001)
smd_after  <- c(0.2606, 0.0182)

## Create data frame
df <- data.frame(
  Covariate = rep(covariates, 2),
  SMD = c(smd_before, smd_after),
  Status = rep(c("Before Matching", "After Matching"), each = 2)
)

## Love plot with only axis lines
love_plot <- ggplot(df, aes(x = Covariate, y = SMD, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray") +
  labs(#title = "Covariate Balance Before and After Matching",
    y = "Standardized Mean Difference (SMD)",
    x = "Covariate") +
  scale_fill_manual(values = c("Before Matching" = "red", "After Matching" = "blue")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
ggsave("plots/Suppl_Fig1_A.tiff", plot = love_plot, width = 2.4, height = 2.4, units = "in", dpi = 300)

## Assume full matching object m.out exists
## Add weights to original data
exp_counts_t$weights <- m.out$weights
exp_counts_t$group <- ifelse(exp_counts_t$DIAGNOSIS == 1, "Treated", "Control")

## Reshape data to long format for faceting
df_long <- exp_counts_t %>%
  select(AGE, SCORE, group, weights) %>%
  pivot_longer(cols = c(AGE, SCORE), names_to = "Covariate", values_to = "Value")

## Plot weighted density with facets
weight_plot <-ggplot(df_long, aes(x = Value, color = group, fill = group, weight = weights)) +
  geom_density(alpha = 0.3, size = 1) +
  facet_wrap(~Covariate, scales = "free") +
  scale_color_manual(values = c("Control" = "#0072B2", "Treated" = "darkorange")) +
  scale_fill_manual(values = c("Control" = "#0072B2", "Treated" = "darkorange")) +
  labs(#title = "Weighted Density Plots of Matched Covariates",
    x = "Value",
    y = "Weighted Density") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
ggsave("plots/Suppl_Fig1_B.tiff", plot = weight_plot, width = 4, height = 4, units = "in", dpi = 300)