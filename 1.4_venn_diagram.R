###############################################################################
###############################################################################
# Title: Venn Diagram Analysis of Differentially Expressed Transcripts
# Author: Rashi Verma
# Date: 2025-11-13
# Description:
#   - Load differential expression result tables from pairwise comparisons
#     (HS vs TPA, TPA vs SM, HS vs SM)
#   - Extract unique transcript/gene identifiers from each comparison
#   - Generate high-resolution Venn diagrams for:
#       • All significant transcripts
#       • Top 50 up- and down-regulated transcripts
###############################################################################
###############################################################################

library(VennDiagram)
library(readr)

# All transcripts
# Read CSV files (replace with your actual filenames)
data1 <- read.csv("results/1_HS_vs_TPA_all_t.csv")
data2 <- read.csv("results/2_TPA_vs_SM_all_t.csv")
data3 <- read.csv("results/3_HS_vs_SM_all_t.csv")

# Extract gene lists (replace 'Gene' with the actual column name)
genes1 <- unique(data1$X)
genes2 <- unique(data2$X)
genes3 <- unique(data3$X)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(
    Set1 = genes1,
    Set2 = genes2,
    Set3 = genes3
  ),
  filename = "plots/venn_diagram_all.tiff",
  output = TRUE,
  imagetype = "tiff",
  height = 3000,
  width = 3000,
  resolution = 600,
  col = "transparent",
  fill = c("#FFA000", "#283593", "#41C9F8"), # Customize colors
  alpha = 0.6,
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.fontface = "bold",
  cat.pos = 0
)

# witout label set1 , 2 and 3
venn.plot <- venn.diagram(
  x = list(
    HS_vs_TPA = genes1,
    TPA_vs_SM = genes2,
    HS_vs_SM = genes3
  ),
  filename = "plots/venn_diagram_publication.tiff",
  output = TRUE,
  imagetype = "tiff",
  units = "in",
  height = 4,
  width = 4,
  resolution = 600,
  col = "transparent",
  fill = c("#FFA000", "#283593", "#41C9F8"),
  alpha = 0.6,
  cex = 2,
  fontface = "bold",
  cat.cex = 0,        # Hide set labels
  cat.pos = 0
)

#######################################################################
# Top 50 up and down
# Read CSV files (replace with your actual filenames)
data1 <- read.csv("results/1_HS_vs_TPA_50_t.csv")
data2 <- read.csv("results/2_TPA_vs_SM_50_t.csv")
data3 <- read.csv("results/3_HS_vs_SM_50_t.csv")

# Extract gene lists (replace 'Gene' with the actual column name)
genes1 <- unique(data1$X)
genes2 <- unique(data2$X)
genes3 <- unique(data3$X)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(
    Set1 = genes1,
    Set2 = genes2,
    Set3 = genes3
  ),
  filename = "plots/venn_diagram_unique_50.tiff",
  output = TRUE,
  imagetype = "tiff",
  height = 3000,
  width = 3000,
  resolution = 600,
  col = "transparent",
  fill = c("#FFA000", "#283593", "#41C9F8"), # Customize colors
  alpha = 0.6,
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.fontface = "bold",
  cat.pos = 0
)

# witout label set1 , 2 and 3
venn.plot <- venn.diagram(
  x = list(
    HS_vs_TPA = genes1,
    TPA_vs_SM = genes2,
    HS_vs_SM = genes3
  ),
  filename = "plots/venn_diagram_unique_50.tiff",
  output = TRUE,
  imagetype = "tiff",
  units = "in",
  height = 4,
  width = 4,
  resolution = 600,
  col = "transparent",
  fill = c("#FFA000", "#283593", "#41C9F8"),
  alpha = 0.6,
  cex = 2,
  fontface = "bold",
  cat.cex = 0,        # Hide set labels
  cat.pos = 0
)
##################################################################################
