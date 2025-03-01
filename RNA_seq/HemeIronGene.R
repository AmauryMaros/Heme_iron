# Load required libraries
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readxl)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(grid)
library(extrafont)
# font_import()
loadfonts()
loadfonts(device = "pdf")
par(family = "Arial", font = 2)

police <- "Arial"

# Load data----

SampleInfo <- read.csv("HRG1_K562_SampleInfo.csv", row.names = 1, stringsAsFactors = TRUE)
HRG1.K562 <- readRDS("HRG1-K562.hg38_12samples.rds")

# Rename colnames of HRG1.K562
colnames(HRG1.K562) <- sapply(strsplit(colnames(HRG1.K562), "\\."), `[`, 1)

# Ensure colnames of count table match rownames of sample metadata
metadata <- SampleInfo %>%
  mutate(Treatment = recode(Treatment, "WT_U" = "Ctl", "KO_U" = "Ctl", "WT_H" = "SB", "KO_H" = "SB")) %>%
  column_to_rownames('SampleID')

metadata$Genotype_treatment <- paste(metadata$Genotype, metadata$Treatment, sep = "_")
metadata$Genotype_treatment <- factor(metadata$Genotype_treatment)

# Check if colnames of count data match rownames of metadata
stopifnot(all(colnames(HRG1.K562) %in% rownames(metadata)))
stopifnot(all(colnames(HRG1.K562) == rownames(metadata)))

# KO_U vs WT_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_Ctl") # set reference 
dds <- DESeqDataSetFromMatrix(countData = assay(HRG1.K562),  
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_U <- results(dds, contrast = c("Genotype_treatment", "KO_Ctl", "WT_Ctl"), pAdjustMethod="BH")

# lfcShrink ashr helps stabilize the fold changes for genes with low counts
res_U <- lfcShrink(dds, coef = "Genotype_treatment_KO_Ctl_vs_WT_Ctl", type = "ashr")

# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("HemeIronGene_plot/MAplot_KO_Ctl_vs_WT_Ctl.pdf", width = 8, height = 6)
par(family = police, font=2)
plotMA(res_U, main = "WT vs KO in Control (Ctl)", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("HemeIronGene_plot/volcano_plot_KO_Ctl_vs_WT_Ctl.pdf", width = 8, height = 6, family = police)
EnhancedVolcano(res_U,
                lab = rownames(res_U),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "KO_Ctl vs WT_Ctl",
                subtitle = "Ref. WT_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
dev.off()


# KO_SB vs WT_SB----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_SB")
dds <- DESeqDataSetFromMatrix(countData = assay(HRG1.K562),  
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)
res_SB <- results(dds, contrast = c("Genotype_treatment", "KO_SB", "WT_SB"), pAdjustMethod="BH")
res_SB <- lfcShrink(dds, coef = "Genotype_treatment_KO_SB_vs_WT_SB", type = "ashr")

pdf("HemeIronGene_plot/MAplot_KO_SB_vs_WT_SB.pdf", width = 8, height = 6)
par(family = police, font=2)
plotMA(res_SB, main = "WT vs KO in Treated (SB) Condition", ylim = c(-5, 5))
dev.off()

# Volcano plot
pdf("HemeIronGene_plot/volcano_plot_KO_SB_vs_WT_SB.pdf", width = 8, height = 6, , family = police)
EnhancedVolcano(res_SB,
                lab = rownames(res_SB),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "KO_SB vs WT_SB",
                subtitle = "Ref. WT_SB",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
dev.off()


# KO_SB vs KO_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "KO_Ctl") # set reference 

dds <- DESeqDataSetFromMatrix(countData = assay(HRG1.K562), #change count matrix if needed 
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_koSB_koU <- results(dds, contrast = c("Genotype_treatment", "KO_SB", "KO_Ctl"), pAdjustMethod="BH")


# lfcShrink ashr helps stabilize the fold changes for genes with low counts

res_koSB_koU <- lfcShrink(dds, coef = "Genotype_treatment_KO_SB_vs_KO_Ctl", type = "ashr")

# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("HemeIronGene_plot/MAplot_KO_SB_vs_KO_Ctl.pdf", width = 8, height = 6)
par(family = police, font=2)
plotMA(res_koSB_koU, main = "Treated vs Untreated in KO Condition", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("HemeIronGene_plot/volcano_plot_KO_SB_vs_KO_Ctl.pdf", width = 8, height = 6, family = police)
EnhancedVolcano(res_koSB_koU,
                lab = rownames(res_koSB_koU),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "KO_SB vs KO_Ctl",
                subtitle = "Ref. KO_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
dev.off()


# Save FC tables
KO_SB_vs_KO_U_all <- as.data.frame(res_koSB_koU)

KO_SB_vs_KO_U_filter <- KO_SB_vs_KO_U_all[KO_SB_vs_KO_U_all$log2FoldChange > 3 | KO_SB_vs_KO_U_all$log2FoldChange < -3 | (-log10(KO_SB_vs_KO_U_all$padj) > 10), ]

KO_SB_vs_KO_U_filter = KO_SB_vs_KO_U_filter[!is.na(KO_SB_vs_KO_U_filter),]
# check some genes to see if the filtering is ok
# KO_U_vs_WT_U_filter[rownames(KO_U_vs_WT_U_filter ) == 'FZD3',]

write.csv(KO_SB_vs_KO_U_all, "FC table/KO_SB_vs_KO_U_all.csv", row.names = TRUE)
write.csv(KO_SB_vs_KO_U_filter, "FC table/KO_SB_vs_KO_U_filter.csv", row.names = TRUE)

# WT_SB vs WT_U----

# Construct a DESeqDataSet object
metadata$Genotype_treatment <- relevel(factor(metadata$Genotype_treatment), ref = "WT_Ctl") # set reference 

dds <- DESeqDataSetFromMatrix(countData = assay(HRG1.K562), #change count matrix if needed 
                              colData = metadata,
                              design = ~ Genotype_treatment)
#  Filter Out Low-Count Genes
dds<-dds[rowSums(counts(dds))>1,]

# Run DESeq2
dds <- DESeq(dds)

# Extract Differential Expression Results
res_wtSB_wtU <- results(dds, contrast = c("Genotype_treatment", "WT_SB", "WT_Ctl"), pAdjustMethod="BH")


# lfcShrink ashr helps stabilize the fold changes for genes with low counts

res_wtSB_wtU <- lfcShrink(dds, coef = "Genotype_treatment_WT_SB_vs_WT_Ctl", type = "ashr")

# MA plots display log fold changes (LFC) vs. mean expression levels
pdf("HemeIronGene_plot/MAplot_WT_SB_vs_WT_Ctl.pdf", width = 8, height = 6)
par(family = police, font=2)
plotMA(res_wtSB_wtU, main = "Treated vs Untreated in WT Condition", ylim = c(-5, 5))
dev.off()


# Volcano plot
pdf("HemeIronGene_plot/volcano_plot_WT_SB_vs_WT_Ctl.pdf", width = 8, height = 6, family = police)
EnhancedVolcano(res_wtSB_wtU,
                lab = rownames(res_wtSB_wtU),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                title = "WT_SB vs WT_Ctl",
                subtitle = "Ref. WT_Ctl",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                legendLabSize = 10)
dev.off()


# Save FC tables
WT_SB_vs_WT_U_all <- as.data.frame(res_wtSB_wtU)
WT_SB_vs_WT_U_filter <- WT_SB_vs_WT_U_all[WT_SB_vs_WT_U_all$log2FoldChange > 3 | WT_SB_vs_WT_U_all$log2FoldChange < -3 | (-log10(WT_SB_vs_WT_U_all$padj) > 10), ]
WT_SB_vs_WT_U_filter = WT_SB_vs_WT_U_filter[!is.na(WT_SB_vs_WT_U_filter),]
# check some genes to see if the filtering is ok
# KO_U_vs_WT_U_filter[rownames(KO_U_vs_WT_U_filter ) == 'FZD3',]
write.csv(WT_SB_vs_WT_U_all, "FC table/WT_SB_vs_WT_U_all.csv", row.names = TRUE)
write.csv(WT_SB_vs_WT_U_filter, "FC table/WT_SB_vs_WT_U_filter.csv", row.names = TRUE)

# Plot KO_SB/KO_U vs WT_SB/ST_U ----

data.to.plot <- dplyr::full_join(KO_SB_vs_KO_U_all %>% rownames_to_column("Gene"), 
                                WT_SB_vs_WT_U_all %>% rownames_to_column("Gene"), 
                                by = "Gene", 
                                suffix = c("_KO_SB_vs_KO_Ctl", "_WT_SB_vs_WT_Ctl"))

up_down_treshlod <- 1

# Create subdataframe for each corner of the future plots
data_up_up <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
data_up_down <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl > up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
data_down_up <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl > up_down_treshlod)
data_down_down <- subset(data.to.plot,log2FoldChange_KO_SB_vs_KO_Ctl < -up_down_treshlod & log2FoldChange_WT_SB_vs_WT_Ctl < -up_down_treshlod)
data_middle <- data.to.plot[!(data.to.plot$Gene %in% c(data_up_up$Gene, data_up_down$Gene, data_down_up$Gene, data_down_down$Gene)), ]


# Create a matrix with the 4 values (counts) in a 2x2 table
table_data <- matrix(
  c(nrow(data_down_up), nrow(data_down_down), nrow(data_up_up), nrow(data_up_down)),
  nrow = 2, ncol = 2,
  dimnames = list(NULL, NULL)  # No row and column names
)
# Color for text in table (based on groups)
text_colors <- c(
  'orangered', # for data_up_up
  'turquoise', # for data_down_down
  'turquoise', # for data_down_up
  'orangered'  # for data_up_down
)

# Create a table theme for styling (no rownames or colnames)
table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 30, fontface = "bold", col = text_colors),  # Set colors for text
    bg_params = list(fill = "white"),
    lwd = 2,
    lty = 2,
    col = "black"  # Border color
  ),
  colhead = list(
    fg_params = list(fontsize = 30, fontface = "bold", col = "black"),
    bg_params = list(fill = "white")
  ),
  rowhead = list(
    fg_params = list(fontsize = 30, fontface = "bold", col = "black"),
    bg_params = list(fill = "white")
  )
)


pdf("HemeIronGene_plot/WT_SB_WT_Ctl_vs_KO_SB_KO_Ctl.pdf", width = 8, height = 6, family = police)
# Create the plot
ggplot(data_middle, aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl)) +
  geom_point(color="black", alpha=0.5, size = 5) +
  geom_point(data = data_up_up,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl),
             color = 'turquoise', alpha=0.5, size=5) +
  geom_point(data = data_down_down,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'turquoise', alpha=0.5, size=5) +
  geom_point(data = data_down_up,
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  geom_point(data = data_up_down, 
             aes(x = log2FoldChange_WT_SB_vs_WT_Ctl, y = log2FoldChange_KO_SB_vs_KO_Ctl), 
             color = 'orangered', alpha=0.5, size=5) +
  
  geom_hline(yintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -up_down_treshlod, linetype = "dashed", color = "black") +
  geom_vline(xintercept = up_down_treshlod, linetype = "dashed", color = "black") +
  
  labs(x = bquote(Log[2] * "FC (KO SB vs KO Ctl)"), 
       y = bquote(Log[2] * "WT SB vs WT Ctl")) +
  
  scale_x_continuous() +
  scale_y_continuous() +
  
  theme_minimal() +
  theme(
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank(),
    axis.text = element_text(size = 30, color = "black", family = police, face = "bold"),
    axis.title = element_text(size = 30, color = "black", family = police, face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 1)  # Adds tick marks
  ) +
  
  # Add the 2x2 table to the plot, placed in the top-left corner
  annotation_custom(
    grob = tableGrob(
      table_data, 
      theme = table_theme, 
      rows = NULL, # No rownames
      cols = NULL  # No colnames
    ), 
    xmin = -Inf, xmax = -1, ymin = Inf, ymax = 4
  )

dev.off()



# PCA----

# Variance Stabilizing Transformation
vsd <- vst(dds,blind=FALSE)
vsd$SampleID <- factor(as.character(vsd$Genotype_treatment), levels=c("WT_Ctl", "KO_Ctl", "WT_SB", "KO_SB"))

# Regularized Log Transformation
rld<-rlog(dds,blind=FALSE)
rld$SampleID <- factor(as.character(rld$Genotype_treatment), levels=c("WT_Ctl", "KO_Ctl", "WT_SB", "KO_SB"))

# Use rld instead of vsd because relatively small dataset = less than 10–15 samples per condition

# Personalization of PCA plot
fontsize <- 20
cbPalette <- c( "darkgrey", "black","pink", "red")

# Extract PCA coordinates
pca_data <- as.data.frame(plotPCA(rld, intgroup = c("Genotype_treatment"), returnData = TRUE)) %>%
  arrange(Genotype_treatment) %>%
  mutate(Genotype_treatment = str_replace_all(Genotype_treatment, "_", " "),  # Replace underscores with spaces
         Genotype_treatment = str_replace_all(Genotype_treatment, "U", "Ctl"))

percentVar <- round(100 * attr(pca_data, "percentVar"))  # Round to whole numbers

# Define custom offsets for each point
x_nudge <- 10
pca_data$nudge_x <- c(x_nudge,x_nudge,x_nudge,
                      -x_nudge,-x_nudge,-x_nudge,
                      x_nudge,x_nudge,x_nudge,
                      -x_nudge,-x_nudge,-x_nudge)
y_nudge <- 0
pca_data$nudge_y <- c(y_nudge,y_nudge,y_nudge+1,# WT Ctl
                      y_nudge,y_nudge,y_nudge, # KO Ctl
                      y_nudge,y_nudge,y_nudge+2, # WT SB
                      y_nudge,y_nudge,y_nudge) # KO SB

# Create PCA plot with variance in axis labels
pdf("HemeIronGene_plot/PCA_plot.pdf", width = 6, height = 4, family = police)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  # Plot PCA points with specified aesthetics
  geom_point(aes(fill = Genotype_treatment), 
             size = 4, shape = 21, color = "black", stroke = 1) +
  # Set the colors for Genotype_treatment
  scale_fill_manual(values = cbPalette) +
  # Add text labels for Genotype_treatment with adjusted position
  geom_text(aes(label = Genotype_treatment, color = Genotype_treatment,
                x = PC1 + nudge_x, y = PC2 + nudge_y), 
            size = 6, fontface = "bold") +
  # Manually set the text color for labels to match the fill color
  scale_color_manual(values = cbPalette) +
  # Add dashed lines at x = 0 and y = 0 for visual reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Add % variance explained to x and y axis labels
  labs(x = paste0("PC1 (", percentVar[1], "%)"), 
       y = paste0("PC2 (", percentVar[2], "%)")) +
  # Apply a clean theme for the plot
  theme_bw() +
  theme(
    # Remove the legend
    legend.position = "none",  
    # Remove grid lines (both major and minor)
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),#                           Add a thick border around the plot
    axis.line = element_blank(), #                                                                Remove the default x and y axis lines
    axis.text = element_text(size = fontsize, color = "black", family = police, face = "bold"),#  Format axis text
    axis.title = element_text(size = fontsize, color = "black", family = police, face = "bold")#  Format axis titles
  )
dev.off()

# Heatmap Top 35 ----

# Use rld data to plot heatmap of most variable genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 35)
topVarData <- as.data.frame(assay(rld)[topVarGenes, ])

annotation_col <- metadata[, c('Genotype_treatment'), drop = FALSE]
rownames(annotation_col) <- rownames(metadata)

annotation_colors <- list(
  Genotype_treatment = c(
    "KO_Ctl" = "pink",
    "KO_SB" = "red",
    "WT_Ctl" = "grey",
    "WT_SB" = "black"
  )
)


pdf("HemeIronGene_plot/heatmap_top35_genes.pdf", width = 8, height = 6, family = police)
pheatmap(
  topVarData,
  scale = "row",             # Scale data by rows (z-score normalization)
  cluster_rows = TRUE,       # Cluster genes by similarity
  cluster_cols = TRUE,       # Cluster samples by similarity
  show_rownames = TRUE,      # Show gene names
  show_colnames = TRUE,      # Show sample names
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # Color palette
  main = "Top 35 Variable Genes",
  annotation_col = annotation_col,  # Add column annotations
  annotation_colors = annotation_colors  # Apply custom colors to annotations
)
dev.off()


# Heatmap of AB genes list----

# Use list of genes from AB and plot heatmap
AB_genes_list <- read_excel("list_of_genes.xlsx")
AB_genes_list <- AB_genes_list$both
AB_genes_list <- unique(AB_genes_list)
AB_genes_list <- unlist(strsplit(AB_genes_list, "/"))

all_genes_capitalized <- toupper(rownames(assay(rld)))
AB_genes_capitalized <- toupper(AB_genes_list)

my_genes <- c()
for(i in AB_genes_capitalized) {
  if (i %in% all_genes_capitalized) {
    my_genes <- c(my_genes, i)
  }
}

matching_indices <- which(rownames(assay(rld)) %in% my_genes)
data_for_heatmap <- as.data.frame(assay(rld)[matching_indices, ])
data_for_heatmap <- data_for_heatmap[rowSums(data_for_heatmap != 0) > 0, ]

# Plot heatmap
pdf("HemeIronGene_plot/heatmap_AB_genes.pdf", width = 8, height = 6, family = police)
pheatmap(
  data_for_heatmap,
  scale = "row",          # Scale data by rows (z-score normalization)
  cluster_rows = T,       # Cluster genes by similarity
  cluster_cols = T,       # Cluster samples by similarity
  show_rownames = T,      # Show gene names
  show_colnames = T,      # Show sample names
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # Color palette
  main = "Genes from AB list",
  annotation_col = annotation_col,  # Add column annotations
  annotation_colors = annotation_colors,  # Apply custom colors to annotations
  fontsize_row = 5
)
dev.off()

# Check if all the genes from rownames(assay(rld))[matching_indices] exist in res_SB
matching_genes <- rownames(assay(rld))[matching_indices]
valid_genes <- matching_genes[matching_genes %in% rownames(res_SB)]
invalid_genes <- matching_genes[!matching_genes %in% rownames(res_SB)]
print(invalid_genes)

res_SB_filtered <- res_SB[matching_genes, ]

# Volcano plot
pdf("HemeIronGene_plot/volcano_plot_AB_genes.pdf", width = 8, height = 6, family = police)
EnhancedVolcano(res_SB_filtered,
                lab = rownames(res_SB_filtered),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "KO_SB vs WT_SB",
                subtitle = "Ref. WT_SB",
                caption = "Differential Expression",
                colAlpha = 0.7,
                legendPosition = 'top',
                labSize = 7,
                labFace = "italic",
                labCol = 'black',
                drawConnectors = TRUE,
                # widthConnectors = 0.5,
                encircle = FALSE,
                xlab = expression(Log[2] ~ "FC (KO SB vs WT SB)"),
                ylim = c(0, 18),  # Adjust y-axis max (change 10 to desired value)
                xlim = c(-3.5, 4.5) 
) + 
  theme_classic() + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),  # Add border
    axis.title.x = element_text(size = 20, face = "bold", colour = 'black'),  # Increase x-axis title size and bold
    axis.title.y = element_text(size = 20, face = "bold", colour = 'black'),
    axis.text.x = element_text(size = 20, face = "bold", colour = 'black'),   # Increase x-axis tick labels size and bold
    axis.text.y = element_text(size = 20, face = "bold", colour = 'black')    # Increase y-axis tick labels size and bold
  )
dev.off()



# HBE1 plot ----

data_for_plot_HBE1 <- as.data.frame(assay(rld))
data_for_plot_HBE1 <- data_for_plot_HBE1["HBE1", , drop = FALSE]
data_for_plot_HBE1 <- as.data.frame(t(data_for_plot_HBE1))
SampleName <- rownames(data_for_plot_HBE1)
data_for_plot_HBE1$SampleName <- SampleName
data_for_plot_HBE1 <- merge(data_for_plot_HBE1, metadata)

# Define custom colors
genotype_colors <- c(
  "KO_Ctl" = "pink",
  "KO_SB" = "red",
  "WT_Ctl" = "grey",
  "WT_SB" = "black"
)


pdf("HemeIronGene_plot/HBE1.pdf", width = 8, height = 6)
ggplot(data_for_plot_HBE1, aes(x = Genotype_treatment, y = HBE1, color = Genotype_treatment)) +
  geom_point(size = 4)+
  scale_color_manual(values = genotype_colors) +
  theme_minimal() +
  labs(x = "Group", y = "HBE1 (rld)", color = "Group") +  # Rename legend
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
dev.off()


## TEST DU TRUC BIZARRE

# dds <- DESeqDataSetFromMatrix(countData = assay(HRG1.K562),
#                               colData = metadata,
#                               design = ~ Genotype_treatment)
# dds <- DESeq(dds)
# res.U <- results(dds, contrast = c('Genotype_treatment',"WT_U", "KO_U"), pAdjustMethod="BH")
# res.U <- lfcShrink(dds,contrast = c('Genotype_treatment',"WT_U", "KO_U"), res=res.U, type = "ashr")
# 
# 
# res.SB <- results(dds,contrast = c('Genotype_treatment',"WT_SB", "KO_SB"), pAdjustMethod="BH")
# res.SB <- lfcShrink(dds, contrast = c('Genotype_treatment',"WT_SB", "KO_SB"), res=res.SB, type = "ashr")
# 
# 
# df_U <- as.data.frame(res.U)
# colnames(df_U) <- paste0(colnames(df_U), "_U")
# 
# df_SB <- as.data.frame(res.SB)
# colnames(df_SB) <- paste0(colnames(df_SB), "_SB")
# 
# 
# merged_df <- cbind(df_U, df_SB)
# 
# 
# UU <- subset(merged_df,log2FoldChange_U > 1 & log2FoldChange_SB > 1)
# UD <- subset(merged_df,log2FoldChange_U > 1 & log2FoldChange_SB < -1)
# DU <- subset(merged_df,log2FoldChange_U < -1 & log2FoldChange_SB > 1)
# DD <- subset(merged_df,log2FoldChange_U < -1 & log2FoldChange_SB < -1)
# 
# U_NS <- subset(merged_df,between(log2FoldChange_U, -1, 1))
# SB_NS <- subset(merged_df,between(log2FoldChange_SB, -1, 1))
# 
# # Subgroup for labeling
# genename <- row.names(merged_df)
# HemeIronGene <- subset(merged_df, rownames(merged_df) %in% rownames(rld)[matching_indices])#make a subset of gene
# HemeIronGeneName <- row.names(HemeIronGene)#get the gene names of subset for annotation

# 
# pdf("HemeIronGene_plot/test.pdf", width = 8, height = 6)
# # Make a basic ggplot2 object
# vol <- ggplot(merged_df, aes(x = log2FoldChange_U, y = log2FoldChange_SB))# Make a basic ggplot2 object
# My_theme <- theme(legend.position = "right", axis.text=element_text(size=28, family="sans", face="bold"),
#                   axis.title=element_text(size=36,family="sans", face="bold"))
# 
# 
# #Label all genes
# vol +
#   ggtitle(label = "Fold changes", subtitle = "Colored by fold-change direction") +
#   geom_point(data = U_NS, size = 5, colour = "black", alpha = 0.5, na.rm = T) +
#   geom_point(data = SB_NS, size = 5, colour = "black", alpha = 0.5, na.rm = T) +
#   geom_point(data = UU, size = 5, colour = "turquoise", alpha = 0.5, na.rm = T) +
#   geom_point(data = DD, size = 5, colour = "turquoise", alpha = 0.5, na.rm = T) +
#   geom_point(data = UD, size = 5, colour = "orangered", alpha = 0.5, na.rm = T) +
#   geom_point(data = DU, size = 5, colour = "orangered", alpha = 0.5, na.rm = T) +
#   geom_text(aes(label = genename), size = 6, vjust = -1) +
#   #geom_cor(method = "pearson", ypos = 1e5) +
#   theme_classic() +
#   My_theme +
#   xlab(expression(log[2]("WT_U"/"KO_U"))) +
#   ylab(expression(log[2]("WT_SB"/"KO_SB"))) +
#   geom_hline(yintercept = 1, colour="#990000", linetype="dashed") +
#   geom_hline(yintercept = -1, colour="#990000", linetype="dashed") +
#   geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
#   geom_vline(xintercept = -1, colour="#990000", linetype="dashed")
# 
# #Only label HRGs
# vol +
#   ggtitle(label = "Fold changes", subtitle = "Colored by fold-change direction") +
#   geom_point(data = U_NS, size = 5, colour = "black", alpha = 0.5, na.rm = T) +
#   geom_point(data = SB_NS, size = 5, colour = "black", alpha = 0.5, na.rm = T) +
#   geom_point(data = UU, size = 5, colour = "turquoise", alpha = 0.5, na.rm = T) +
#   geom_point(data = DD, size = 5, colour = "turquoise", alpha = 0.5, na.rm = T) +
#   geom_point(data = UD, size = 5, colour = "orangered", alpha = 0.5, na.rm = T) +
#   geom_point(data = DU, size = 5, colour = "orangered", alpha = 0.5, na.rm = T) +
#   geom_text(data = HemeIronGene, aes(label = HemeIronGeneName), size = 10, vjust = -1) +
#   #geom_cor(method = "pearson", ypos = 1e5) +
#   theme_classic() +
#   My_theme +
#   xlab(expression(log[2]("WT_U"/"KO_U"))) +
#   ylab(expression(log[2]("WT_SB"/"KO_SB"))) +
#   geom_hline(yintercept = 1, colour="#990000", linetype="dashed") +
#   geom_hline(yintercept = -1, colour="#990000", linetype="dashed") +
#   geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
#   geom_vline(xintercept = -1, colour="#990000", linetype="dashed")
# 
# dev.off()
# 
# 
# 
# 
