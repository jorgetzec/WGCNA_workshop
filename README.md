# Gene Co-expression Network Workshop

This repository contains the materials for a workshop on Gene Co-expression Networks, with a focus on the implementation and exploration of the Weighted Gene Co-expression Network Analysis (WGCNA) methodology.

## Workshop Overview

The workshop is based on a subset of the *Chlamydomonas reinhardtii* transcriptome under salt stress conditions. The main goal is to explore this transcriptome dataset, apply WGCNA concepts, and investigate how different parameters affect the construction and interpretation of gene co-expression networks.

## Repository Structure

- [data/](data)  
  Contains raw count data, example network files, and other datasets used in the workshop.
  - [data_example1_net_edges.txt](data/data_example1_net_edges.txt): Example network edges file.
  - [data_example2_network.csv](data/data_example2_network.csv): Example network file from different tool.
  - [rawData.csv](data/rawData.csv): Raw count matrix.
  - [proteome_uniprot.csv](data/proteome_uniprot.csv): Protein annotation file.
  - [coldata.csv](data/coldata.csv): Sample metadata.

- [docs/](docs)  
  Contains course notes, slides, syllabus, a list of useful tools, and report instructions.
  - [notes_Co-expressionNetwork.md](docs/notes_Co-expressionNetwork.md): Workshop notes.
  - [Syllabus_Co-expressionNetwork.md](docs/Syllabus_Co-expressionNetwork.md): Workshop syllabus.
  - [ToolLists_Co-expressionNetwork.md](docs/ToolLists_Co-expressionNetwork.md): List of tools and resources.
  - [Co-expressionNetwork_WorkshopContent.pdf](docs/Co-expressionNetwork_WorkshopContent.pdf): Workshop slides.
  - [WGCNA_workshop_report_instructions.md](docs/WGCNA_workshop_report_instructions.md): **Instructions for preparing your research report.**
  - [Figures/](docs/Figures): Figures and images used in the workshop.

- [code/](code)  
  Contains all R scripts and code used for the workshop exercises and analyses.
  - [WGCNA_script1.R](code/WGCNA_script1.R): Main WGCNA analysis script.
  - [WGCNA_script2.R](code/WGCNA_script2.R): Additional WGCNA analysis script.
  - [data_example2_net_edges.txt](code/data_example2_net_edges.txt): Example network edges for exercises.
  - [data_example2_net_nodes.txt](code/data_example2_net_nodes.txt): Example network nodes for exercises.

## Workshop Report Instructions

Please refer to [docs/WGCNA_workshop_report_instructions.md](docs/WGCNA_workshop_report_instructions.md) for detailed guidelines on preparing your research report for this workshop.

## Getting Started

1. Clone or download (zip and unzip) this repository.
2. To use the scripts, first open the `workshop.Rproj` file located in the root folder. This will ensure the correct working directory is set.
3. Follow the scripts in the `code/` directory and refer to the notes in `docs/` for guidance. The scripts include details on the necessary R packages to install.

## About the Example Dataset

The example dataset is a subset of the *C. reinhardtii* transcriptome under salt stress. Participants will explore this dataset, apply WGCNA, and learn how to adjust parameters to construct and interpret gene co-expression networks.

---

For questions or suggestions, please open an issue or contact the workshop organizers.
