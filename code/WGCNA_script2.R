
###########--------------------------------------------------------------------###########
#                                 WGCNA Analysis (Salt Stress-Session 2)
###########--------------------------------------------------------------------###########

# Prefix for all output files
condition_prefix <- "salt"
# Current date for file naming
date_export <- format(Sys.Date(), "%Y%m%d")

# Load necessary libraries
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

###########-----------------------0. Loading information from previus session ---------------

load("salt_WGCNA_session1.RData") # Load the WGCNA network object from the previous session
here() # Set the working directory to the project root
# i_am("code/WGCNA_script2.R") 

###########-----------------------13. Gene significance and trait relationship---------------

selected_sample <- "salt8h_1"  #<---- Define time to analyze

# Create a binary trait matrix
sample_headers <- coldata_filtered$id
traits <- as.data.frame(diag(1, nrow = length(sample_headers), ncol = length(sample_headers)))
rownames(traits) <- sample_headers
colnames(traits) <- sample_headers

# Calculate number of samples and genes
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

# Calculate gene significance and p-val
time_to_analize <- traits[[selected_sample]]
geneTraitSignificance <- cor(datExpr, time_to_analize, use = 'p')
GSPvalue <- corPvalueStudent(geneTraitSignificance, nSamples)

head(geneTraitSignificance)
# Convert to data.frame with appropriate names
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
GSPvalue <- as.data.frame(GSPvalue)

# Check result
head(geneTraitSignificance)
head(GSPvalue)


# Module colors and module of interest
moduleColors <- bwnet$colors
module <- "blue"  # Define module to be analyzed

# Get Module Names (remove "ME" prefix)
modNames <- gsub("^ME", "", names(module_eigengenes))

# Calculate Gene-Module Correlation (Module Membership) for ALL modules
geneModuleMembership <- cor(datExpr, module_eigengenes, use = 'p') # Output: genes x modules
colnames(geneModuleMembership) <- paste0("MM.", modNames)
rownames(geneModuleMembership) <- colnames(datExpr)

# Calculate p-values for Module Membership
MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
colnames(MMPvalue) <- paste0("p.MM.", modNames)


column <- match(module, modNames)
moduleGenes <- moduleColors == module

# Scatter plot to see relationship between MM and GS
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", selected_sample),
                   main = "Module membership vs. gene significance",
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


###########-----------------------14. Hub genes identification---------------

####------- Hub genes identification
# First, check chooseTopHubInEachModule function in WGCNA documentarion
?chooseTopHubInEachModule


#Function modified to get the top ten hub genes per module
chooseTopHubInEachModule_modf <- function (datExpr, colorh, omitColors = "grey", power = 2, type = "signed", top_n = 10, ...) {
  isIndex = FALSE
  modules = names(table(colorh))
  if (!(is.na(omitColors)[1])) 
    modules = modules[!is.element(modules, omitColors)]
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = 1:dim(datExpr)[2]
    isIndex = TRUE
  }
  hubs = list()
  for (m in modules) {
    adj = adjacency(datExpr[, colorh == m], power = power, 
                    type = type, ...)
    top_hubs_indices = order(rowSums(adj), decreasing = TRUE)[1:top_n]
    top_hubs = colnames(adj)[top_hubs_indices]
    hubs[[m]] = top_hubs
  }
  if (isIndex) {
    hubs = as.numeric(hubs)
    names(hubs) = modules
  }
  return(hubs)
}


as.data.frame(chooseTopHubInEachModule_modf(
  datExpr, 
  moduleColors, 
  omitColors = "grey", 
  power = soft_power_val, 
  type = "signed hybrid")) #signed


module_eigengenes = bwnet$MEs
moduleColors <- bwnet$colors

###########----------------------- 15. Network Visualization (Export to Cytoscape) ---------------

# This section handles the export of network data (TOM-based) for visualization in Cytoscape.
# Simplified version: assumes all inputs are correct and aligned.

# --- Prerequisites and Setup  ---
# datExpr: Expression data (samples x genes). Colnames are Gene IDs.
# module_colors_assigned: Vector of module assignments for each gene in datExpr (names are Gene IDs).
# bwnet$TOMFiles: Path to the saved TOM RData file.

# Load annotation data
annot_file_path <- (here("data", 'proteome_uniprot.csv'))

cre_annotations_full <- readr::read_csv(annot_file_path, show_col_types = FALSE)

annotation_map_df <- cre_annotations_full %>%
                      dplyr::select(Gene_id, Gene_names) %>%
                      dplyr::distinct(Gene_id, .keep_all = TRUE)

# Load the TOM matrix
load(bwnet$TOMFiles[1]) # Loads object, typically named 'TOM'
TOM_matrix <- if (inherits(TOM, "dist")) as.matrix(TOM) else TOM

rm(TOM)

# Align TOM_matrix with datExpr gene IDs (colnames of datExpr)
# This direct subsetting assumes datExpr colnames are a subset of or equal to TOM rownames/colnames.
# And that module_colors_assigned names perfectly match datExpr colnames.
rownames(TOM_matrix) <- colnames(datExpr)
colnames(TOM_matrix) <- colnames(datExpr)

current_gene_ids_in_datExpr <- colnames(datExpr)
TOM_aligned_matrix <- TOM_matrix[current_gene_ids_in_datExpr, current_gene_ids_in_datExpr]
module_colors_for_network <- module_colors_assigned[current_gene_ids_in_datExpr]

# Get alternative node names (gene symbols)
alt_node_names_for_network <- annotation_map_df$Gene_names[match(current_gene_ids_in_datExpr, annotation_map_df$Gene_id)]
alt_node_names_for_network[is.na(alt_node_names_for_network)] <- current_gene_ids_in_datExpr[is.na(alt_node_names_for_network)]

###########----------------------- 17. Network visualization ---------------

# --- Module visualization using igraph ---
library(igraph)

# Define genes of interest to highlight
genes_of_interest <- c("Cre14.g632176")  # Add your genes of interest here
highlight_color <- "red"  # Color for highlighting genes of interest

# Select module to visualize
module_color_to_visualize <- "blue"  # CHANGE THIS to the module you want to visualize
edge_threshold <- 0.1  # Threshold for TOM values to consider as connections

# Select genes in the module
genes_in_module <- current_gene_ids_in_datExpr[module_colors_for_network == module_color_to_visualize]

# Get TOM matrix for the module
module_tom <- TOM_aligned_matrix[genes_in_module, genes_in_module]

# Create binary adjacency matrix from TOM
adj_matrix <- module_tom > edge_threshold
diag(adj_matrix) <- FALSE  # Remove self-connections

# Create graph
g <- graph_from_adjacency_matrix(adj_matrix, 
                                mode = "undirected", 
                                weighted = TRUE)

# Add TOM values as edge weights
edge_indices <- which(upper.tri(adj_matrix) & adj_matrix, arr.ind = TRUE)
edge_weights <- module_tom[edge_indices]
E(g)$weight <- edge_weights

# Find and print top connected genes: Top 10 genes by connectivity in module
degree_dist <- sort(degree(g), decreasing = TRUE)
top_genes <- head(degree_dist, 10)
print(top_genes)

# Configure visual attributes
V(g)$size <- 5
V(g)$label.cex <- 0.5
V(g)$color <- module_color_to_visualize

genes_in_module_uniprot<- cre_annotations_full$Gene_names[match(genes_in_module, cre_annotations_full$Gene_id)]
V(g)$label <- genes_in_module_uniprot # or genes_in_module to see the Phytozome ID


# Highlight genes of interest
V(g)$color[V(g)$name %in% genes_of_interest] <- highlight_color
V(g)$size[V(g)$name %in% genes_of_interest] <- 8  # Make highlighted genes bigger

# Define layout
# Different layout options available:
    # layout_matrix <- layout_with_fr(g)      # Fruchterman-Reingold (organic, by default)
    # layout_matrix <- layout_with_kk(g)      # Kamada-Kawai (good for small networks)
    # layout_matrix <- layout_with_drl(g)     # DrL (Distributed Recursive Layout)
    # layout_matrix <- layout_in_circle(g)    # Circular simple
    # layout_matrix <- layout_as_star(g)      # Star (a central node)
    # layout_matrix <- layout_as_tree(g)      # Tree (hierarchical)
    # layout_matrix <- layout_on_grid(g)      # Grid
    # layout_matrix <- layout_on_sphere(g)    # Esfera 3D
    # layout_matrix <- layout_randomly(g)     # Random
    # layout_matrix <- layout_with_sugiyama(g) # For directed networks
    # layout_matrix <- layout_with_gem(g)      # GEM (Graph Embedder)
    # layout_matrix <- layout_with_mds(g)      # MDS (Multidimensional Scaling)

layout_matrix <- layout_with_kk(g)

# Create the plot
plot(g, 
     layout = layout_matrix,
     vertex.label.dist = 0.5,
     vertex.label.color = "black",
     edge.width = E(g)$weight * 2,
     main = paste("Module", module_color_to_visualize, "\nHighlighted genes:", 
                 paste(genes_of_interest, collapse = ", ")))

# Basic graph analysis
cat("\nSummary of the graph for module", module_color_to_visualize, ":\n")
print(summary(g))

# Degree distribution
degree_dist <- degree(g)
cat("\nDegree distribution statistics:\n")
print(summary(degree_dist))

# Print information about highlighted genes
cat("\nInformation about highlighted genes:\n")
for(gene in genes_of_interest) {
  if(gene %in% V(g)$name) {
    cat("\nGene:", gene)
    cat("\nDegree:", degree(g)[gene])
    cat("\nNeighbors:", paste(neighbors(g, gene)$name, collapse = ", "))
    # Add TOM-based connectivity information
    gene_tom <- module_tom[gene,]
    cat("\nTop 5 TOM connections:", 
        paste(names(sort(gene_tom, decreasing = TRUE)[2:6]), 
              round(sort(gene_tom, decreasing = TRUE)[2:6], 3), 
              sep = " (", collapse = "), "), ")")
  } else {
    cat("\nGene", gene, "not found in the network")
  }
}

###########----------------------- 18. Save R Session ---------------
# save.image(file = condition_prefix, "_WGCNA_workshop_session.RData"))

message("WGCNA analysis for salt stress completed.")
