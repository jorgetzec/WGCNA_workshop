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

# Were all packages loaded correctly?