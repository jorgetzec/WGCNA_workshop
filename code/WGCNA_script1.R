###########--------------------------------------------------------------------###########
#                                 WGCNA Analysis (Salt Stress- Session 1)
###########--------------------------------------------------------------------###########

# Prefix for all output files
condition_prefix <- "salt"
# Current date for file naming
date_export <- format(Sys.Date(), "%Y%m%d")

###########----------------------- Package Installation and Loading ---------------

# List of CRAN packages
cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "readr", "reshape2",
  "ggplot2", "ggrepel", "ggpubr", "patchwork", "gridExtra", "pheatmap",
  "here"  
)

# List of Bioconductor packages
bioc_pkgs <- c(
  "DESeq2", "WGCNA", "igraph" # For network visualization if used. Consider include "RCy3"
)

# Install CRAN packages if missing
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if missing
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load libraries
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(DESeq2)
library(WGCNA)
library(here)  
#library(igraph)          

# WGCNA settings
allowWGCNAThreads() # Allow multi-threading

###########-----------------------  Work Directory  ---------------

# Set working directory to the location of this script
# setwd(" ")
# dr_here()
# here()

###########----------------------- Expression Data Loading  ---------------

countData_raw <- read.delim(here("data", "rawData.csv"), header = TRUE, sep = ",", row.names = 1)
head(countData_raw); dim(countData_raw)

###########-----------------------  coldata (PhenoData) Preparation ---------------

coldata<- read.delim(here("data", "coldata.csv"), header = TRUE, sep = ",", row.names = 1)

coldata$time <- as.numeric(coldata$time)
coldata$replicates <- as.integer(coldata$replicates)
coldata$replicates <- as.factor(coldata$replicates)
coldata$condition <- as.factor(coldata$condition)
coldata$treatment <- as.factor(coldata$treatment)

head(coldata); dim(coldata)

###########----------------------- Data Preparation: Filtering and Normalization/Transformtion ---------------
# Analize the filteting options one by one and decide which one should be used.

# OPTION 1. Filtering non-expressed or low-count genes for WGCNA using a common threshold (e.g., mean count > 10)
countData_filtered <- countData_raw %>%
  filter(rowMeans(.) > 10)
dim(countData_filtered) # Check dimensions after filtering


# OPTION 2. Trimming tails to remove outliers (e.g. percentile 99).
gene_means <- rowMeans(countData_raw) # mean calculatarion per gene.
p99 <- quantile(gene_means, probs = 0.95) # percentile calculation
countData_filtered <- countData_raw[gene_means <= p99, ] # filtering genes above percentile 99
dim(countData_filtered) # Check dimensions after filtering

# VST transformation with DESeq2
# For WGCNA, often blind=TRUE is used, or design = ~1 for general variance stabilization.
dds <- DESeqDataSetFromMatrix(countData = countData_filtered,
                              colData = coldata,
                              design = ~ 1) # Using a simple design for VST

    # Normalization steps are not needed here as VST will handle it.
    #dds <- estimateSizeFactors(dds)
    #normalized_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = TRUE)

countData_vst <- assay(vsd)
head(countData_vst, 3); dim(countData_vst)

###########----------------------- Outlier Detection and Removal (Samples) ---------------

# Hierarchical clustering before outlier removal (using VST data for better distance metric)
htree_before_removal <- hclust(dist(t(countData_vst)), method = "average")

plot(htree_before_removal, main = "Sample Clustering (VST data, Before Outlier Removal)", xlab = "", sub = "")
abline(h=30, col="red")

  # pdf(paste0(date_export, "_", condition_prefix, "_dendro_samples_BEFORE_removal.pdf"), width = 7, height = 6)
  # plot again to export
  # dev.off()

# Samples to be removed
samples_to_remove <- c("salt2h_1", "salt72h_1") #Specify sample IDs to remove based on clustering

# Remove specified outlier samples from VST data and coldata
countData_vst_filtered <- countData_vst[, !colnames(countData_vst) %in% samples_to_remove]
coldata_filtered <- coldata[!rownames(coldata) %in% samples_to_remove, ]

dim(countData_vst_filtered)
dim(coldata_filtered)

# Ensure consistency
colnames(countData_vst_filtered) == rownames(coldata_filtered)

# Hierarchical clustering after outlier removal
htree_before_removal <- hclust(dist(t(countData_vst_filtered)), method = "average")
plot(htree_before_removal, main = "Sample Clustering (VST data, Before Outlier Removal)", xlab = "", sub = "")
abline(h=30, col="red")


###########----------------------- Data Properties and Distribution Visualization ---------------

# Properties of data at different stages
message("Dimensions of raw count data: ", paste(dim(countData_raw), collapse = " x "))
message("Dimensions of VST transformed data (before sample removal): ", paste(dim(countData_vst), collapse = " x "))
message("Dimensions of VST transformed data (after sample removal, for WGCNA): ", paste(dim(countData_vst_filtered), collapse = " x "))
message("Dimensions of transposed expression data for WGCNA (datExpr): ", paste(dim(datExpr), collapse = " x "))


# Data distribution plots
# Density: raw counts (original full set)
plot_density_raw <- reshape2::melt(countData_raw, 
                                   variable.name = "treatment", 
                                   value.name = "reads") %>%
  ggplot(aes(x = reads, group = treatment)) +
  geom_density(aes(fill = treatment, 
                   color = treatment), 
               alpha = 0.01) +
  xlab("Raw Reads") + 
  ylab("Density") + 
  theme_classic()+
  theme(legend.position = "none") 


# Density: VST transformed, filtered counts (used for WGCNA)
plot_density_vst_filtered <- reshape2::melt(as.data.frame(countData_vst_filtered), 
                                            variable.name = "treatment", 
                                            value.name = "ExprVST") %>%
  ggplot(aes(x = ExprVST, group = treatment)) +
  geom_density(aes(fill = treatment, 
                   color = treatment), 
               alpha = 0.01) +
  xlab("VST Transformed Counts (Filtered)") + ylab("Density") + 
  theme_classic() +
  theme(legend.position = "none") 


plot_power_conect <- grid.arrange(plot_density_raw, plot_density_vst_filtered, nrow = 1)
print(plot_power_conect)

# ggsave(filename = paste0(date_export, "_", condition_prefix, "_density_vst_filtered_dist.pdf"), plot = plot_power_conect, width = 10, height = 5)

# After checking the density plots, can you see what filtering option is required?
# Go back if necessary and apply the filtering option you prefer.

###########----------------------- Transpose Data for WGCNA ---------------

# Transpose data for WGCNA functions (samples in rows, genes in columns)
datExpr <- t(countData_vst_filtered)
head(datExpr[,1:5]); dim(datExpr)

###########----------------------- 7. Exploratory Data Analysis (Heatmap) ---------------

# Colors definition for plots
colores_time <- c("0_h"= "#E41A1C","2_h"="#377EB8", "4_h"= "#4DAF4A", "8_h"= "#984EA3",
                  "12_h"="#0000EE", "24_h"="#FF7F00", "48_h"="#00EE00", "72_h"="black")


# Calculate Pearson's correlation matrix between samples (columns)
cor_matrix <- cor(countData_vst_filtered, method = "pearson")

# Extracting treatment information for annotation
annotation_df <- data.frame(time = coldata$treatment)
rownames(annotation_df) <- rownames(coldata)

levels_time <- unique(annotation_df$time) #Defining unique levels
ann_colors <- list(time = colores_time) # Defining list for annotation_colors

hm_corSample <- pheatmap(cor_matrix,
                         annotation_col = annotation_df,
                         annotation_colors = ann_colors,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete", #check other methods if needed
                         display_numbers = FALSE,
                         number_format = "%.2f",
                         fontsize_number = 8,
                         #main = "Correlation between samples (Pearson)",
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         border_color = NA)


###########----------------------- 8. WGCNA: Soft Threshold (Power) Definition ---------------
# Reminder: datExpr is t(countData_vst_filtered)
# Check different types of networks before to select one for all the analysis: "signed", "signed hybrid", or "unsigned"


# Choose a set of soft-threshold powers
power_options <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr,
                         powerVector = power_options,
                         networkType = "unsigned", # "signed hybrid" "signed" or "unsigned"
                         verbose = 5)

# Review sft content before decide the beta value (power)
sft_data <- sft$fitIndices

# Plot scale-free topology fit index vs. power
plot_sft_rsq <- ggplot(sft_data, 
                       aes(Power, 
                           SFT.R.sq, 
                           label = Power)) +
  geom_point() +
  geom_text_repel(size = 4, 
                  segment.color = "black", 
                  segment.size = 0.5) +
  geom_hline(yintercept = 0.8, color = 'red') + # Common threshold for R^2
  labs(x = 'Power', y = expression('Scale free topology model fit, signed ' * R^2)) 

# Plot mean connectivity vs. power
plot_sft_meank <- ggplot(sft_data, 
                         aes(Power, 
                             mean.k., 
                             label = Power)) +
            geom_point() +
            geom_text_repel(size = 4, 
                            segment.color = "black", 
                            segment.size = 0.5) +
            labs(x = 'Power', y = 'Mean Connectivity') 

# Arrange plots
plot_sft_combined <- grid.arrange(plot_sft_rsq, plot_sft_meank, nrow = 1)

# Select soft power (e.g., first power to reach R^2 > 0.8)
soft_power_val <- sft_data %>%
                  filter(SFT.R.sq >= 0.8) %>% # Adjust threshold if needed
                  pull(Power) %>%
                  head(1)

# soft_power_val <-6 # If necessary set manually after reviewing the plots

###########----------------------- 10. WGCNA: Network and Module Construction ---------------
# Ensure WGCNA::cor is used if conflicts exist
temp_cor_holder <- cor
cor <- WGCNA::cor

# Adjust maxBlockSize if needed based on number of genes in datExpr
# dim(datExpr) # Check number of genes (ncol(datExpr))
# maxBlockSize should be >= number of genes if you want a single block. 
max_block_size <- ncol(datExpr) + 100

# Let's review the documentation before proceed. Add the parameters you need.
?blockwiseModules # Check documentation for blockwiseModules

bwnet <- blockwiseModules(datExpr,
                          maxBlockSize = max_block_size, # Max number of genes in a block
                          saveTOMs = TRUE, # Save TOMs for later use
                          # Requires definition
                          )

cor <- temp_cor_holder # Restore original cor function

# Module colors assigned to genes
module_colors_assigned <- bwnet$colors
table(module_colors_assigned)

# Plot gene counts per module
module_gene_counts <- as.data.frame(table(module_colors_assigned)) %>%
  setNames(c("Module", "Genes")) %>%
  filter(Module != "grey") %>% # Exclude unassigned (grey) module
  arrange(desc(Genes))

# Use Module name itself for fill if actual colors
module_gene_counts$ColorFill <- ifelse(
  as.character(module_gene_counts$Module) %in% standardColors(100),
  as.character(module_gene_counts$Module),
  "grey"
)

plot_module_counts_bar <- ggplot(module_gene_counts, 
                                 aes(x = reorder(Module, -Genes), 
                                     y = Genes, 
                                     fill = ColorFill)) +
        geom_bar(stat = "identity", colour = "black", linewidth = 0.5) +
        geom_text(aes(label = Genes), vjust = -0.5, size = 3.5) +
        scale_fill_identity(guide = "none") + # Use actual color names for fill
        labs(x = "Module", y = "Gene Count") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

plot_module_counts_bar


# Plot dendrogram and module colors
# pdf(paste0(date_export, "_", condition_prefix, "_dendro_modules.pdf"), width = 8, height = 3) # Adjusted size
plotDendroAndColors(bwnet$dendrograms[[1]],
                    #cbind(bwnet$unmergedColors, bwnet$colors),
                    #c("unmerged", "merged"),
                    bwnet$colors, # Merged colors
                    "Merged Modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05,
                    main = NULL) #main = "Gene Dendrogram and Module Colors"
# dev.off()

###########----------------------- 12. Module Profile Plots (Example Highlighting a Gene) ---------------
# This section shows how to plot all gene profiles within a module and highlight one or more specific genes.
# It is reusable by changing `module_to_plot` and `gene_to_highlight`.
# The x-axis will show individual sample IDs.

# Genes of interest: Identifying MYB transcription factors (TFs) in modules

# Cre03.g197100 
# Cre02.g108350 
# Cre03.g197350 
# Cre07.g345350 
# Cre14.g632176 
# Cre02.g103450 

# Define module and gene for highlighting (ensure these are valid for your data)
module_to_plot <- "blue"
gene_to_highlight <- c("Cre14.g632176")
gene_highlight_color <- "blue"     # Color for the highlighted gene line

# --- 1. Data Preparation ---
# Get genes in the selected module
gene_module_key <- tibble::enframe(bwnet$colors, name = "Genes", value = "module")

genes_in_module_list <- gene_module_key %>%
  filter(module == module_to_plot) %>%
  pull(Genes)

# Subset expression data for these genes
# drop=FALSE ensures it remains a matrix even if only one gene is selected or one sample remains
module_expression_data <- countData_vst_filtered[rownames(countData_vst_filtered) %in% genes_in_module_list, , drop = FALSE]

# Calculate mean expression profile for the module (per sample)
mean_profile_for_module <- colMeans(module_expression_data)

# Combine with original module data for melting, adding the mean profile with a unique Gene ID
module_expression_data_with_mean <- rbind(module_expression_data,
                                          "MeanProfileInternalID" = mean_profile_for_module)

# Reshape to long format
module_expr_long_df <- as.data.frame(module_expression_data_with_mean) %>%
  rownames_to_column("Gene") %>%
  melt(id.vars = "Gene", variable.name = "id", value.name = "Expression")

# Define plotting groups
module_expr_long_df <- module_expr_long_df %>%
  mutate(
    PlotGroup = dplyr::case_when(
      Gene == "MeanProfileInternalID" ~ "Mean Profile",
      Gene %in% gene_to_highlight     ~ "Highlighted Gene",
      TRUE                            ~ "Other Genes"
    ),
    PlotGroup = factor(PlotGroup, levels = c("Other Genes", "Highlighted Gene", "Mean Profile"))
  )

sample_order_for_plot <- colnames(countData_vst_filtered)
module_expr_long_df$id <- factor(module_expr_long_df$id, levels = sample_order_for_plot)

# Define aesthetic mappings
color_palette_module_plot <- c(
  "Other Genes" = "grey80",
  "Highlighted Gene" = gene_highlight_color,
  "Mean Profile" = "black"
)

linewidth_palette_module_plot <- c(
  "Other Genes" = 0.3,
  "Highlighted Gene" = 1.0,
  "Mean Profile" = 1.0
)

linetype_palette_module_plot <- c(
  "Other Genes" = "solid",
  "Highlighted Gene" = "solid",
  "Mean Profile" = "dashed"
)

legend_labels_module_plot <- c(
  "Other Genes" = paste0("Other Genes (", module_to_plot, ")"),
  "Highlighted Gene" = paste0("Highlighted Genes\n (", paste(gene_to_highlight, collapse = ", "), ")"),
  "Mean Profile" = paste0("Mean Profile (", module_to_plot, ")")
)

# --- Plotting with Explicit Layering ---

plot_module_profile_samples_viz <- ggplot(module_expr_long_df,
                                          aes(x = id, y = Expression, group = Gene)) + # Basic aes for all layers
  
  # Layer 1: "Other Genes" - drawn first (at the bottom)
  geom_line(data = . %>% filter(PlotGroup == "Other Genes"),
            aes(color = PlotGroup, linewidth = PlotGroup, , linetype = PlotGroup)) +
  
  # Layer 2: "Highlighted Gene" - drawn on top of "Other Genes"
  geom_line(data = . %>% filter(PlotGroup == "Highlighted Gene"),
            aes(color = PlotGroup, linewidth = PlotGroup, linetype = PlotGroup)) +
  
  # Layer 3: "Mean Profile" - drawn on top of "Highlighted Gene" (and "Other Genes")
  geom_line(data = . %>% filter(PlotGroup == "Mean Profile"),
            aes(color = PlotGroup, linewidth = PlotGroup, linetype = PlotGroup)) +
  
  scale_color_manual(
    values = color_palette_module_plot,
    labels = legend_labels_module_plot,
    name = "Line Type",
    breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")
  ) +
  scale_linewidth_manual(
    values = linewidth_palette_module_plot,
    labels = legend_labels_module_plot,
    name = "Line Type",
    breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")
  ) +
  scale_linetype_manual(
    values = linetype_palette_module_plot,
    labels = legend_labels_module_plot,
    name = "Line Type",
    breaks = c("Other Genes", "Highlighted Gene", "Mean Profile")
  ) +
  guides(linewidth = "none") + 
  labs(x = "Time treatment under 200 mM NaCl (h)", y = "Normalized expression (VST)",
       #title = paste("Expression Profile of Module:", module_to_plot)
       ) +theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

print(plot_module_profile_samples_viz)


###########----------------------- 13. Module-Trait Relationships ---------------
# Module Eigengenes (MEs)
module_eigengenes <- bwnet$MEs # Modules in columns, samples in rows

# Create a numeric trait matrix. For time-course, you might use time points as traits.
# Example: Binary traits for each time point (1 if sample is at that time, 0 otherwise)
module_eigengenes <- bwnet$MEs %>%
  select(-MEgrey)

nSamples <- nrow(datExpr)

sample_headers <- coldata_filtered$id
traits <- as.data.frame(diag(1, nrow = length(sample_headers), ncol = length(sample_headers)))
rownames(traits) <- sample_headers
colnames(traits) <- sample_headers

# 'module_eigengenes' contains the eigengenes of the modules
# and 'traits' is a data frame with treatment variables
module_trait_cor <- cor(module_eigengenes, traits, use = "p")

# 'nSamples' is the number of samples
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples)

# Create a text matrix with correlations and p-values.
textMatrix <- paste(signif(module_trait_cor, 2), "\n(",
                    signif(module_trait_pvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(module_trait_cor)

# open a PDF device—everything that follows irá a ese archivo
# pdf(file   = paste0(date_export, "_", condition_prefix, "_module_eigengenes_pheatmap_pvals.pdf"),
#    width  = 16,
#    height = 10)

par(mar = c(5, 10, 3, 4))  # values in order: bottom, left, top, right
# Generar el heatmap
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = names(traits),
               yLabels = names(module_eigengenes),
               ySymbols = names(module_eigengenes),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# close the PDF device and write the file to disk
# dev.off()

# Save the workspace image
save.image(file = "salt_WGCNA_session1.RData")

# Save the session info
saveRDS(sessionInfo(), file = "session_info_salt_WGCNA_session1.rds")

# session_info <- readRDS("session_info_salt_WGCNA_workshop_session1.rds")
# print(session_info)
