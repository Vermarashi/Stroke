###############################################################################
###############################################################################
# Title: Gene Ontology Enrichment Analysis for Category-Specific Transcripts
# Author: Rashi Verma
# Date: 2025-11-13
# Description:
#   - Load differential transcript expression results for three pairwise
#     comparisons (HS vs TPA, TPA vs SM, HS vs SM)
#   - Identify category-specific (non-overlapping) top 50 transcripts for:
#       • CAT1: HS-specific transcripts
#       • CAT2: TPA-specific transcripts
#   - Filter Ensembl transcript IDs and remove non-ENST entries
#   - Map Ensembl transcript IDs to HGNC gene symbols using Ensembl BioMart
#   - Perform Gene Ontology (GO) enrichment analysis (BP, CC, MF) using
#     clusterProfiler with multiple-testing correction (BH)
#   - Calculate fold enrichment from GeneRatio and BgRatio
#   - Rank GO terms by fold enrichment
#   - Visualize top enriched GO terms using publication-quality dot plots
#     (fold enrichment vs GO term, colored by –log10 adjusted p-value)
###############################################################################
###############################################################################

library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(stringr)

rm(list=ls())

set.seed(1234)
# Save unique transcript list of CAT1, 2 and 3 and retrieve ensemble gene name for ENST Ids
## CAT1
genes1=read.csv("results/1_HS_vs_TPA_50_t.csv")
genes1[1:20,1:5]; dim(genes1) # 100x7

genes2=read.csv("results/2_TPA_vs_SM_50_t.csv")
genes2[1:20,1:5]; dim(genes2) # 100x7

genes3=read.csv("results/3_HS_vs_SM_50_t.csv")
genes3[1:20,1:5]; dim(genes3) # 100x7

# Combine all genes into one vector
# Subset non-overlapping genes.
g<- setdiff(genes1$X, genes2$X)
glist <-setdiff(g, genes3$X)
length(glist) # 80 gene IDs. # 82 gene IDs

# Filter out MSTRG IDs
glist_clean <- glist[grepl("^ENST", glist)]

# Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl transcript IDs to gene symbols
transcript2symbol <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = glist_clean,
  mart = mart
)
#write.csv(transcript2symbol, "results_DTE/Cat1_50_GO.csv", row.names = FALSE)

# Keep only valid gene symbols
glist_symbols <- unique(transcript2symbol$hgnc_symbol)
glist_symbols <- glist_symbols[!is.na(glist_symbols) & glist_symbols != ""]

# Run GO enrichment
go_enrichment <- enrichGO(
  gene          = glist_symbols,
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

# Calculate fold enrichment using GeneRatio and BgRatio
go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) eval(parse(text=x)))  # Convert GeneRatio to numeric
go_df$BgRatio_num <- sapply(go_df$BgRatio, function(x) eval(parse(text=x)))  # Convert BgRatio to numeric
go_df$FoldEnrichment <- go_df$GeneRatio_num / go_df$BgRatio_num  # Calculate fold enrichment
go_df$GO_type<-go_df$ONTOLOGY
go_df_sorted <- go_df[order(go_df$FoldEnrichment,  decreasing = TRUE), ]  # Sort in descending order

# Extract the top 18 GO pathways based on FoldEnrichment
top_go_pathways <- head(go_df_sorted, 50)
print(top_go_pathways)

# Convert p.adjust to -log10 scale for better visualization
top_go_pathways$log10_padj <- -log10(top_go_pathways$p.adjust)

# Ensure GO type is a factor (optional, useful if you later want to facet)
top_go_pathways$GO_type <- factor(top_go_pathways$GO_type, levels = c("BP", "CC", "MF"))

top_go_pathways$Description_wrapped <- str_wrap(top_go_pathways$Description, width = 50)

# Create the dotplot without splitting
tiff(file="plots/Fig_6B.tiff", unit= "in", res = 600, width = 7, height = 5)
ggplot(top_go_pathways, 
            aes(x = FoldEnrichment, 
                y = reorder(Description_wrapped, FoldEnrichment), 
                size = Count, 
                color = log10_padj)) +
  geom_point() +
  scale_color_gradientn(colors = c("#D94F4F", "#B080A0", "#3C68C8")) +
  scale_size(range = c(1.5, 3.5)) +  # Adjust point sizes
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Fold Enrichment",
    y = "GO Term",
    color = "-log10(p.adjust)",
    size = "Gene Count"
  )

dev.off()

#################################################################################
## CAT2
genes1=read.csv("results/1_HS_vs_TPA_50_t.csv")
genes1[1:20,1:5]; dim(genes1) # 100x7

genes2=read.csv("results/2_TPA_vs_SM_50_t.csv")
genes2[1:20,1:5]; dim(genes2) # 100x7

genes3=read.csv("results/3_HS_vs_SM_50_t.csv")
genes3[1:20,1:5]; dim(genes3) # 100x7

# Combine all genes into one vector
# Subset non-overlapping genes.
g<- setdiff(genes2$X, genes1$X)
glist <-setdiff(g, genes3$X)
length(glist) # 80 gene IDs. # 82 gene IDs

# Filter out MSTRG IDs
glist_clean <- glist[grepl("^ENST", glist)]

# Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl transcript IDs to gene symbols
transcript2symbol <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = glist_clean,
  mart = mart
)
#write.csv(transcript2symbol, "results_DTE/Cat1_50_GO.csv", row.names = FALSE)

# Keep only valid gene symbols
glist_symbols <- unique(transcript2symbol$hgnc_symbol)
glist_symbols <- glist_symbols[!is.na(glist_symbols) & glist_symbols != ""]

# Run GO enrichment
go_enrichment <- enrichGO(
  gene          = glist_symbols,
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

# Calculate fold enrichment using GeneRatio and BgRatio
go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) eval(parse(text=x)))  # Convert GeneRatio to numeric
go_df$BgRatio_num <- sapply(go_df$BgRatio, function(x) eval(parse(text=x)))  # Convert BgRatio to numeric
go_df$FoldEnrichment <- go_df$GeneRatio_num / go_df$BgRatio_num  # Calculate fold enrichment
go_df$GO_type<-go_df$ONTOLOGY
go_df_sorted <- go_df[order(go_df$FoldEnrichment,  decreasing = TRUE), ]  # Sort in descending order

# Extract the top 18 GO pathways based on FoldEnrichment
top_go_pathways <- head(go_df_sorted, 50)
print(top_go_pathways)

# Convert p.adjust to -log10 scale for better visualization
top_go_pathways$log10_padj <- -log10(top_go_pathways$p.adjust)

# Ensure GO type is a factor (optional, useful if you later want to facet)
top_go_pathways$GO_type <- factor(top_go_pathways$GO_type, levels = c("BP", "CC", "MF"))

top_go_pathways$Description_wrapped <- str_wrap(top_go_pathways$Description, width = 50)

# Create the dotplot without splitting
tiff(file="plots/GO_Cat2_50.tiff", unit= "in", res = 600, width = 7, height = 5)
ggplot(top_go_pathways, 
       aes(x = FoldEnrichment, 
           y = reorder(Description_wrapped, FoldEnrichment), 
           size = Count, 
           color = log10_padj)) +
  geom_point() +
  scale_color_gradientn(colors = c("#D94F4F", "#B080A0", "#3C68C8")) +
  scale_size(range = c(1.5, 3.5)) +  # Adjust point sizes
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Fold Enrichment",
    y = "GO Term",
    color = "-log10(p.adjust)",
    size = "Gene Count"
  )

dev.off()
