### Gene Co-expression Network (WGCNA) and Connectivity Visualization Tools

|Tool|Interface|Language/Platform|Main Features|Official Link|
|:--|:--|:--|:--|:--|
|**WGCNA (R)**|CLI|R|Built-in visualization for heatmaps and dendrograms; exports to Cytoscape/VisANT. Offers limited static plots.|[CRAN WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html)|
|**igraph**|CLI|R, Python, C/C++|Versatile for network analysis and visualization. Supports custom attributes and is high-performance for large networks.|[igraph R](https://igraph.org/r/), [igraph Python](https://igraph.org/python/)|
|**NetworkX**|CLI|Python|User-friendly Python library for building and analyzing networks. Supports WGCNA edge lists. Slower for large-scale networks.|[NetworkX PyPI](https://pypi.org/project/networkx/)|
|**Graph-tool**|CLI|Python, C++|High-performance C++ library for statistical network analysis. Visualizes WGCNA networks using edge weights. Steeper learning curve.|[Graph-tool](https://graph-tool.skewed.de/)|
|**Tabnetviz**|CLI|Python, Graphviz|Lightweight CLI tool that uses Graphviz for SVG visualizations. Maps WGCNA node attributes to visual properties.|[Tabnetviz, Python version](https://pypi.org/project/tabnetviz/)|
|**Cytoscape**|GUI|Java|Open-source platform for complex network visualization. Imports WGCNA via edge lists or R. Supports GO/KEGG terms.|[Cytoscape](https://cytoscape.org/)|
|**Graphia**|GUI|Java|Open-source tool for interactive network exploration. Visualizes WGCNA networks with customizable layouts. Handles large, complex datasets.|[Graphia](https://graphia.net/)|
|**Gephi**|GUI|Java|Open-source tool for interactive network exploration. Visualizes WGCNA networks with customizable layouts. Suitable for large-scale networks.|[Gephi](https://gephi.org/)|
|**Omics Playground**|GUI|Web-based|User-friendly platform for WGCNA analysis and visualization. Limited customization. Registration required.|[Omics Playground](https://omicsplayground.com/main/wgcna_module)|
|**PyWGCNA**|CLI|Python|Python implementation of WGCNA with visualization features (e.g., module-trait heatmaps, PPI networks). Supports functional enrichment.|[PyWGCNA PyPI](https://pypi.org/project/PyWGCNA/), [GitHub](https://github.com/SysbioBig/PyWGCNA)|

### WGCNA-Based or Similar Tools for Co-expression Network Analysis

|Tool|Interface|Language/Platform|Main Features|Official Link|
|:--|:--|:--|:--|:--|
|**WGCNA**|CLI|R|Builds co-expression networks using weighted correlations.|[CRAN WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html)|
|**BioNERO**|CLI|R|All-in-one network inference package. Uses WGCNA for gene co-expression networks (GCNs) and infers gene regulatory networks (GRNs) via CLR, GENIE3, ARACNe. Includes preprocessing, PCA, heatmaps, and topological analysis. Beginner-friendly.|[Bioconductor BioNERO](https://bioconductor.org/packages/release/bioc/html/BioNERO.html)|
|**multiWGCNA**|CLI|R|Extends WGCNA for multidimensional datasets. Builds networks at three levels and maps modules across them. Analyzes differential expression and module preservation.|[CRAN multiWGCNA](https://cran.r-project.org/web/packages/multiWGCNA/index.html)|
|**PyWGCNA**|CLI|Python|Python implementation of WGCNA. Faster than R version for large datasets. Includes module-trait correlation, hub gene detection, and functional enrichment (GO, KEGG, REACTOME).|[PyPI PyWGCNA](https://pypi.org/project/PyWGCNA/), [GitHub](https://github.com/SysbioBig/PyWGCNA)|
|**CEMiTool**|CLI|R|Infers co-expression modules using correlation-based methods similar to WGCNA. Integrates functional enrichment (GO, KEGG) and interaction data (e.g., PPI). Generates interactive HTML reports. Easier to use than WGCNA for beginners.|[Bioconductor CEMiTool](https://bioconductor.org/packages/release/bioc/html/CEMiTool.html)|
|**iterativeWGCNA**|CLI|Python, R|Extension of WGCNA wrapped in Python. Improves module robustness via iterative refinement and increases gene inclusion. Requires R WGCNA.|[PyPI iterativeWGCNA](https://pypi.org/project/iterativeWGCNA/), [GitHub](https://github.com/SysbioBig/iterativeWGCNA)|
|**GEN**|GUI|Web-based|Online platform for bulk RNA-seq analysis. Includes WGCNA for co-expression networks. Offers differential expression, enrichment, and GRN inference. No coding required. Limited customization.|[GEN Tools](https://www.gentools.cn/)|
|**Exp Omics**|GUI|Web-based|Suite of six user-friendly web applications with 28 advanced features for multi-omics expression analysis, including transcriptome analysis with WGCNA. Limited customization.|[ExpOmics](https://www.expomics.com/)|
|**Gene Expression Nebulas**|GUI|Web-based|Tools for downstream analysis of bulk and single-cell RNA-seq data, including WGCNA. Limited customization.|[Gene Expression Nebulas](https://geneexpressionnebulas.com/)|

### Gene Regulatory Network Inference Tools

|Tool|Interface|Language/Platform|Main Features|Official Link|
|:--|:--|:--|:--|:--|
|**ARACNe**|CLI|R, Python|Infers regulatory networks using mutual information. Removes indirect interactions via DPI. Suitable for large datasets.|[ARACNe GitHub](https://github.com/califano-lab/ARACNe)|
|**GENIE3**|CLI|R, Python|Uses machine learning (random forests) to infer directed regulatory interactions. High accuracy for transcription factor targets. Integrates with WGCNA outputs.|[GENIE3 GitHub](https://github.com/aertslab/GENIE3), [Bioconductor GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)|
|**SCENIC**|CLI|Python, R|Combines expression and motif analysis to infer regulons (TF + targets). Visualizes via heatmaps or t-SNE. Focused on single-cell data.|[SCENIC](https://github.com/aertslab/SCENIC)|
|**iRegulon**|GUI, CLI|Cytoscape plugin|Predicts transcription factors via motif enrichment. Integrates with Cytoscape for visualization. User-friendly GUI.|[iRegulon](https://braph.github.io/iregulon/)|
|**TRRUST**|GUI|Web-based|Curated database of TF-target interactions in humans/mice. Visualizes networks online. Limited to known interactions.|[TRRUST](https://www.grnpedia.org/trrust/)|
|**GeNeCK**|GUI|Web-based|Comprehensive online toolkit for building gene networks using expression data and optional hub gene input. No specific external link beyond site name.|[GeNeCK](http://www.genetoolkit.net/)|
