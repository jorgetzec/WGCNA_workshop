

---

# WGCNA: From Gene Expression to Functional Co-expression Networks

### Introduction to Network Biology

Biological systems are frequently represented as **networks**, which are complex sets of binary interactions or relations between different entities. Network biology applies **graph theory** tools to represent and analyze these systems. In this context, **nodes** represent entities like proteins or genes, and **edges** convey information about the links between them. Networks can have different types of edges: **undirected, directed, or weighted**.

Common types of biological networks include protein-protein interaction networks, metabolic networks, genetic interaction networks, gene/transcriptional regulatory networks, cell signaling networks, and **gene co-expression networks**. Gene co-expression networks are primarily associated with **Transcriptomics** data, such as RNA-seq or Microarrays.

### WGCNA: A Standard Tool for Co-expression Networks

**WGCNA (Weighted Gene Co-expression Network Analysis)** is a standard tool used to construct co-expression networks from gene expression data. It is "Weighted" because it uses **continuous edge weights (0 to 1)**, "Correlation-based" as it utilizes Pearson/Spearman correlations, and designed for "Scale-free topology". WGCNA is **noise-robust** due to the use of Topological Overlap Matrix (TOM) and is **scalable** for large datasets. A key limitation is that it infers **no causality** and can be sensitive to batch effects. Co-expression networks are typically **undirected** (symmetric gene relationships).

### WGCNA Workflow

The WGCNA workflow transforms expression data into networks through a series of steps:

1. **Expression Matrix (X)
		
    - **Input**: A matrix of gene expression values (e.g., RNA-seq counts).
    - **Preprocessing**: Normalization (e.g., VST, Rlog, CPM, TMM, or log2(TPM+1)) and filtering of low-variance genes.
    - **Output**: A cleaned matrix `X` where `X = xis` (expression of gene `i` in sample `s`), typically `5000 genes × 20 samples`.
    - `[insert figure: Example of an Expression Matrix]`
2. **Correlation Matrix (R)**
    
    - **Computation**: Pairwise correlations between genes are calculated to measure co-expression. **Pearson correlation** is commonly used.
    - **Formula**: `rij = cov(Xi, Xj) / (σXi * σXj)`
        - Where `Xi` and `Xj` are expression profiles of genes `i` and `j`.
    - **Output**: A symmetric `n × n` matrix `R = [rij]`, where `rij ∈ [-1, 1]`.
3. **Adjacency Matrix (A)**
    
    - **Transformation**: Correlations are transformed into **weighted connections** using a **power function (β)**. This is known as **soft-thresholding**, which highlights strong connections and reduces the impact of weak correlations.
    - **Formula**: `aij = |cor(Xi, Xj)|^β`
        - `β` is chosen to achieve **scale-free topology** in the network.
    - **Output**: A symmetric `n × n` matrix `A = [aij]`, where `aij ∈`.
    - `[insert figure: Adjacency matrix at different power values (β)]`
4. **Topological Overlap Matrix (TOM)**
    
    - **Measure of Similarity**: TOM measures the **topological similarity** between genes, considering not only direct connections but also the **number of shared neighbors**. This helps in **reducing noise** from spurious correlations.
    - **Formula**: `TOMij = (Σu aiuauj + aij) / (min(ki, kj) + 1 - aij)`
        - `Σu aiuauj` represents shared neighbors.
        - `min(ki, kj)` normalizes by connectivity.
    - **Output**: A symmetric `n × n` matrix `TOM = [TOMij]`, where `TOMij ∈`.
    - `[insert figure: Network with Topological Overlap Matrix (TOM) simulation of i-j connections]`
5. **Dissimilarity TOM (dissTOM)**
    
    - **Transformation for Clustering**: The `TOM` matrix is converted into a dissimilarity measure for easier clustering.
    - **Formula**: `dissTOM = 1 - TOMij`
    - **Output**: A symmetric `n × n` matrix `dissTOM = [dissTOMij]`, where `dissTOM ∈`.
6. **Clustering Modules**
    
    - **Hierarchical Clustering**: `dissTOM` is used as a distance measure for hierarchical clustering, constructing a dendrogram.
    - **Dynamic Tree Cutting**: This method identifies **modules** (groups of genes with low `dissTOM`) based on parameters like `minModuleSize` and `deepSplit`. Genes with low `dissTOM` (high `TOM`) are clustered together.
    - **Eigengene Calculation**: For each module, a **module eigengene (ME)** is calculated. The ME is the **first principal component (PC1)** of the expression profiles of all genes within that module. It **summarizes the module's overall expression pattern** in a single vector, capturing the dominant trend of gene expression.
    - **Module Merging**: Modules can be merged if their eigengenes are highly correlated (e.g., `>0.75`).

### Module Analysis and Trait Correlation

After module detection, WGCNA allows for the **correlation of module eigengenes with external biological traits** (e.g., disease status, treatment) to identify biologically relevant modules. Modules highly correlated with a trait are likely involved in that trait's biology.

Key concepts for assessing gene importance within modules and their relation to traits are:

- **Gene Significance (GS)**: Measures how strongly a gene's expression correlates with a specific trait.
    - **Formula**: `GSi = |cor(Xi, trait)|`
- **Module Membership (MM)**: Quantifies how central a gene is to its module by measuring the correlation between a gene's expression and its module's eigengene.
    - **Formula**: `MMi = |cor(Xi, ME_Module)|`
    - Genes with **high GS and high MM** in a trait-correlated module are considered key candidates for biological roles.
- **Hub Genes**: These are genes that are **highly connected within a module**. They are identified by selecting genes with high MM (e.g., top 10%) or high **intramodule connectivity**. Hub genes are often functionally important, such as transcription factors in biological pathways.
    - **Intramodule Connectivity Formula**: `kwithin,i = Σj∈module,j≠i aij`
    - `[insert figure: Example of hub genes in a module]`

### Biological Interpretation

Once modules and hub genes are identified, further biological insights can be gained through:

- **Gene Ontology (GO) enrichment analysis**.
- **Identification of transcription factors**.
- **KEGG pathways overrepresentation**.
- **Functional annotation** (e.g., BLAST) and **domains/motifs identification**.

### Network Visualization

Network visualization tools are used to represent the inferred networks. Common steps include:

1. **Exporting** node (genes with MM, GS) and edge (adjacency/TOM) files.
2. **Visualizing** networks, often coloring nodes by module or connectivity.
3. **Annotating** with functional terms (e.g., GO/KEGG).

Tools for WGCNA network visualization and analysis include **Cytoscape**, igraph, NetworkX, Graph-tool, Gephi, and specialized WGCNA modules in platforms like Omics Playground or PyWGCNA.

### WGCNA Applications

WGCNA is widely applied across various scientific and biological systems. Examples include:

- **Animal Science**: Identifying hub genes in bovine mammary glands during lactation.
- **Human Health (Oncology)**: Clustering differentially expressed genes (DEGs) into modules and identifying hub genes linked to prognosis in cancers like human tongue squamous cell carcinoma, oral squamous cell carcinoma, osteosarcoma, and lung adenocarcinoma. It has also been used with microbiome data in tumor contexts.
- **Comparative Genomics**: Comparing WGCNA networks across human and mouse tissues to identify conserved gene roles.
- **Plant Biology**: Identifying hub genes related to pathogen resistance in cassava and drought tolerance in maize.

### Co-expression Networks vs. Regulatory Networks

It is crucial to distinguish between **co-expression networks** and **gene regulatory networks**:

- **Co-expression Networks**: Edges are based on **correlation of expression profiles**, implying **no causality**. They are **undirected** and group genes with similar expression patterns.
- **Regulatory Networks**: Edges represent **causal interactions** (e.g., a transcription factor binding to a gene promoter). They are **directed** and model mechanistic control. Inference often relies on binding data or motif databases. Tools like ARACNe, GENIE3, and SCENIC are used for regulatory network inference.

Both network types share similarities: they represent biological relationships as graphs, identify functional modules, are data-driven, and utilize visualization tools.

---
---

# WGCNA: From Gene Expression to Functional Co-expression Networks

## Introduction to Network Biology

- Biological systems are often modeled as **networks**, representing complex interactions between entities.
- **Network biology** uses **graph theory** to analyze these systems:
    - **Nodes**: Entities like proteins or genes.
    - **Edges**: Relationships, which can be **undirected**, **directed**, or **weighted**.
- Common biological networks:
    - Protein-protein interaction networks
    - Metabolic networks
    - Genetic interaction networks
    - Gene/transcriptional regulatory networks
    - Cell signaling networks
    - **Gene co-expression networks** (linked to **transcriptomics** data, e.g., RNA-seq, microarrays)

## WGCNA: A Standard Tool for Co-expression Networks

- **WGCNA (Weighted Gene Co-expression Network Analysis)** constructs co-expression networks from gene expression data.
- Key features:
    - **Weighted**: Uses continuous edge weights (0 to 1).
    - **Correlation-based**: Employs Pearson or Spearman correlations.
    - **Scale-free topology**: Designed to mimic biological networks.
    - **Noise-robust**: Uses Topological Overlap Matrix (TOM).
    - **Scalable**: Handles large datasets.
- Limitations:
    - Infers **no causality**.
    - Sensitive to batch effects.
- Co-expression networks are **undirected** (symmetric gene relationships).

## WGCNA Workflow

The WGCNA workflow transforms expression data into networks through these steps:

### 1. Expression Matrix (X)

- **Input**: Matrix of gene expression values (e.g., RNA-seq counts).
- **Preprocessing**: Normalize (e.g., VST, Rlog, CPM, TMM, log2(TPM+1)) and filter low-variance genes.
- **Output**: Cleaned matrix `X` where `X_i,s` is the expression of gene `i` in sample `s`, e.g., 5000 genes × 20 samples (≥15 recommended).
- _[Insert figure: Example of an Expression Matrix]_

### 2. Correlation Matrix (R)

- **Computation**: Calculate pairwise Pearson correlations between genes.
- **Formula**: `rij = cov(Xi, Xj) / (σXi * σXj)`
    - `Xi`, `Xj`: Expression profiles of genes `i`, `j`.
    - `cov`: Covariance; `σXi`, `σXj`: Standard deviations.
- **Output**: Symmetric `n × n` matrix `R = [rij]`, where `rij ∈ [-1, 1]`.

### 3. Adjacency Matrix (A)

- **Transformation**: Convert correlations to weighted connections via soft-thresholding with power `β`.
- **Formula**: `aij = |cor(Xi, Xj)|^β`
    - `β`: Chosen via `pickSoftThreshold` to achieve scale-free topology (scale-free fit index r ≥ 0.8).
- **Output**: Symmetric `n × n` matrix `A = [aij]`, where `aij ∈ [0, 1]`.
- _[Insert figure: Adjacency matrix at different power values (β)]_

### 4. Topological Overlap Matrix (TOM)

- **Measure**: Quantifies topological similarity, considering direct connections and shared neighbors to reduce noise.
- **Formula**: `TOMij = (Σu aiuauj + aij) / (min(ki, kj) + 1 - aij)`
    - `Σu aiuauj`: Sum of adjacency products over all genes `u` (except `i`, `j`).
    - `ki = Σu aiu`: Connectivity of gene `i`.
    - `min(ki, kj)`: Minimum connectivity of genes `i`, `j`.
- **Output**: Symmetric `n × n` matrix `TOM = [TOMij]`, where `TOMij ∈ [0, 1]`.
- _[Insert figure: Network with TOM simulation of i-j connections]_

### 5. Dissimilarity TOM (dissTOM)

- **Transformation**: Convert TOM to dissimilarity for clustering.
- **Formula**: `dissTOMij = 1 - TOMij`
- **Output**: Symmetric `n × n` matrix `dissTOM = [dissTOMij]`, where `dissTOMij ∈ [0, 1]`.

### 6. Clustering Modules

- **Hierarchical Clustering**: Use `dissTOM` as a distance measure to build a dendrogram.
- **Dynamic Tree Cutting**: Identify modules (gene groups with low `dissTOM`, high `TOM`) using parameters like `minModuleSize` and `deepSplit`.
- **Eigengene Calculation**: Compute **module eigengene (ME)** as the first principal component (PC1) of each module’s gene expression, summarizing its expression pattern.
- **Module Merging**: Merge modules if eigengenes are highly correlated (e.g., >0.75, corresponding to `mergeCutHeight=0.25`).

## Module Analysis and Trait Correlation

- Correlate **module eigengenes** with external traits (e.g., disease, stress) to identify biologically relevant modules.
- Key metrics:
    - **Gene Significance (GS)**: Gene-trait correlation.
        - **Formula**: `GSi = |cor(Xi, trait)|`
    - **Module Membership (MM)**: Gene centrality in a module, via correlation with its eigengene.
        - **Formula**: `MMi = |cor(Xi, ME_module)|`
        - Genes with **high GS and MM** in trait-correlated modules are key candidates.
    - **Hub Genes**: Highly connected within a module, identified by high MM or intramodule connectivity.
        - **Intramodule Connectivity Formula**: `kwithin,i = Σj∈module,j≠i aij` (within-module, vs. total `ki = Σj aij`).
- _[Insert figure: Example of hub genes in a module]_

## Biological Interpretation

- Post-module analysis:
    - **Gene Ontology (GO) enrichment**
    - **Transcription factor identification**
    - **KEGG pathway overrepresentation**
    - **Functional annotation** (e.g., BLAST, domains/motifs)

## Network Visualization

- Steps:
    - **Export**: Node (genes with MM, GS) and edge (adjacency/TOM) files via `exportNetworkToCytoscape`.
    - **Visualize**: Networks colored by module or connectivity.
    - **Annotate**: With GO/KEGG terms.
- Tools: **Cytoscape**, igraph, NetworkX, Graph-tool, Gephi, Omics Playground, PyWGCNA.

## WGCNA Applications

- **Animal Science**: Hub genes in bovine mammary glands (lactation).
- **Human Health (Oncology)**: Hub genes in tongue squamous cell carcinoma, oral squamous cell carcinoma, osteosarcoma, lung adenocarcinoma; microbiome in tumors.
- **Comparative Genomics**: Conserved gene roles in human/mouse tissues.
- **Plant Biology**: Pathogen resistance in cassava, drought tolerance in maize, somatic embryogenesis in Arabidopsis.

## Co-expression Networks vs. Regulatory Networks

- **Co-expression Networks**:
    - **Undirected**, correlation-based edges (no causality).
    - Group genes with similar expression patterns.
- **Regulatory Networks**:
    - **Directed**, causal interactions (e.g., transcription factor to gene).
    - Use binding/motif data. Tools: ARACNe, GENIE3, SCENIC.
- **Similarities**: Graph-based, identify modules, data-driven, visualized similarly.


----
----


### Fórmulas Clave en WGCNA (Formato Markdown/LaTeX para Obsidian)

A continuación, se detallan las fórmulas esenciales utilizadas en el flujo de trabajo de WGCNA, completando la información de la presentación:

#### 1. Matriz de Correlación (R)

Esta matriz mide la **co-expresión** entre pares de genes. La correlación de Pearson es comúnmente utilizada.

- **Fórmula:** $\qquad r_{ij} = \frac{cov(X_i, X_j)}{\sigma_{X_i} \sigma_{X_j}}$
    
- **Explicación de los términos:**
    
    - `rij`: Coeficiente de correlación entre el gen `i` y el gen `j`.
    - `cov(Xi, Xj)`: Covarianza de los perfiles de expresión del gen `i` (`Xi`) y el gen `j` (`Xj`) a través de las muestras.
    - `σXi`: Desviación estándar del perfil de expresión del gen `i`.
    - `σXj`: Desviación estándar del perfil de expresión del gen `j`.

#### 2. Matriz de Adyacencia (A)

Esta matriz transforma las correlaciones en **conexiones ponderadas** (`aij`) utilizando una función de potencia (`β`), lo que se conoce como **soft-thresholding**. Este paso ayuda a resaltar las conexiones fuertes y reducir el impacto de las débiles, logrando una topología de red libre de escala.

- **Fórmula:** $\qquad a_{ij} = |cor(X_i, X_j)|^\beta$
    
- **Explicación de los términos:**
    
    - `aij`: La fuerza de conexión (adyacencia) entre el gen `i` y el gen `j`.
    - `cor(Xi, Xj)`: La correlación (e.g., Pearson) entre los perfiles de expresión del gen `i` y el gen `j`.
    - `β` (beta): Es el exponente de potencia (o umbral suave) elegido para la red. Se selecciona para lograr una **topología de red libre de escala**.

#### 3. Matriz de Superposición Topológica (TOM)

TOM mide la **similitud topológica** entre genes, considerando no solo las conexiones directas sino también el **número de vecinos compartidos**. Esto es crucial para **reducir el ruido** de correlaciones espurias.

- **Fórmula:** $\qquad \text{TOM}_{ij} = \frac{\sum_u a_{iu}a_{uj} + a_{ij}}{\min(k_i, k_j) + 1 - a_{ij}}$
    
- **Explicación de los términos:**
    
    - `TOMij`: El valor de superposición topológica entre el gen `i` y el gen `j`.
    - `Σu aiuauj`: Representa la suma de las adyacencias de los vecinos compartidos entre el gen `i` y el gen `j`. Esto captura el número de caminos de longitud dos que conectan `i` y `j` a través de un vecino común `u`.
    - `aij`: La adyacencia (conexión directa) entre el gen `i` y el gen `j`.
    - `ki` y `kj`: La conectividad total (o grado) de los genes `i` y `j` respectivamente. Se define como la suma de las fuerzas de conexión de un gen con todos los demás genes en la red (`ki = Σj≠i aij`).
    - `min(ki, kj)`: Normaliza la fórmula por la conectividad de los genes `i` y `j`.

#### 4. Disimilitud TOM (dissTOM)

Para facilitar la agrupación de genes en módulos, la matriz `TOM` se convierte en una medida de **disimilitud**.

- **Fórmula:** $\qquad \text{dissTOM}_{ij} = 1 - \text{TOM}_{ij}$
    
- **Explicación de los términos:**
    
    - `dissTOMij`: La medida de disimilitud entre el gen `i` y el gen `j`.
    - `TOMij`: El valor de superposición topológica entre el gen `i` y el gen `j`. Un `dissTOM` bajo (y un `TOM` alto) indica que los genes están altamente conectados y tienden a agruparse.

#### 5. Significación Génica (GS)

La Significación Génica mide cuán fuertemente la expresión de un gen (`Xi`) se **correlaciona con un rasgo biológico externo** (e.g., estado de enfermedad, tratamiento).

- **Fórmula:** $\qquad \text{GS}_i = |cor(X_i, \text{trait})|$
    
- **Explicación de los términos:**
    
    - `GSi`: La significación génica del gen `i`.
    - `cor(Xi, trait)`: La correlación entre el perfil de expresión del gen `i` (`Xi`) y el valor del rasgo de interés. El valor absoluto `|...|` indica la fuerza de la correlación, sin importar la dirección.

#### 6. Membresía de Módulo (MM)

La Membresía de Módulo cuantifica qué tan **central es un gen dentro de su módulo** al medir la correlación entre la expresión del gen y el eigengene de su módulo.

- **Fórmula:** $\qquad \text{MM}_i = |cor(X_i, \text{ME}_{\text{Module}})|$
    
- **Explicación de los términos:**
    
    - `MMi`: La membresía de módulo del gen `i`.
    - `cor(Xi, ME_Module)`: La correlación entre el perfil de expresión del gen `i` (`Xi`) y el **eigengene del módulo** al que pertenece el gen (`ME_Module`). El eigengene del módulo (`ME`) es el primer componente principal (PC1) de los perfiles de expresión de todos los genes dentro de ese módulo, resumiendo su patrón de expresión dominante.

#### 7. Conectividad Intramodular ($k_{\text{within},i}$)

Mide la **fuerza con la que un gen se conecta con otros genes dentro de su propio módulo**. Los genes con alta conectividad intramodular son considerados **genes hub** dentro del módulo.

- **Fórmula:** $\qquad k_{\text{within},i} = \sum_{j \in \text{module}, j \neq i} a_{ij}$
    
- **Explicación de los términos:**
    
    - $k_{\text{within},i}$: La conectividad intramodular del gen `i`.
    - $a_{ij}$: La adyacencia (fuerza de conexión) entre el gen `i` y el gen `j`.
    - $\sum_{j \in \text{module}, j \neq i}$: La suma se realiza sobre todos los genes `j` que pertenecen al _mismo módulo_ que el gen `i`, excluyendo el gen `i` mismo.

----
---
---

# WGCNA: From Gene Expression to Functional Co-expression Networks

### Introduction to Network Biology

Biological systems are frequently represented as **networks**, which are complex sets of binary interactions or relations between different entities. Network biology applies **graph theory** tools to represent and analyze these systems. In this context, **nodes** represent entities like proteins or genes, and **edges** convey information about the links between them. Networks can have different types of edges: **undirected, directed, or weighted**.

Common types of biological networks include protein-protein interaction networks, metabolic networks, genetic interactionENDA to genetic networks, gene/transcriptional regulatory networks, cell signaling networks, and **gene co-expression networks**. Gene co-expression networks are primarily associated with **Transcriptomics** data, such as RNA-seq or Microarrays.

### WGCNA: A Standard Tool for Co-expression Networks

**WGCNA (Weighted Gene Co-expression Network Analysis)** is a standard tool used to construct co-expression networks from gene expression data. It is "Weighted" because it uses **continuous edge weights (0 to 1)**, "Correlation-based" as it utilizes Pearson/Spearman correlations, and designed for "Scale-free topology". WGCNA is **noise-robust** due to the use of Topological Overlap Matrix (TOM) and is **scalable** for large datasets. A key limitation is that it infers **no causality** and can be sensitive to batch effects. Co-expression networks are typically **undirected** (symmetric gene relationships).

### WGCNA Workflow

The WGCNA workflow transforms expression data into networks through a series of steps:

1. **Expression Matrix (X) Processing**
    
    - **Input**: A matrix of gene expression values (e.g., RNA-seq counts).
    - **Preprocessing**: Normalization (e.g., VST, Rlog, CPM, TMM, or log2(TPM+1)) and filtering of low-variance genes.
    - **Output**: A cleaned matrix `X` where `X_i,s` (**expression** of gene `i` in sample `s`), typically `5000 genes × 20 samples`.
    - `[insert figure: Example of an Expression Matrix]`
2. **Correlation Matrix (R) calculation**
    
    - **Computation**: Pairwise correlations between genes are calculated to measure co-expression. **Pearson correlation** is commonly used.
    - **Formula**: `rij = cov(Xi, Xj) / (σXi * σXj)`
        - Where `Xi` and `Xj` are expression profiles of genes `i` and `j`.
    - **Output**: A symmetric `n × n` matrix `R = [rij]`, where `rij ∈ [-1, 1]`.
3. **Adjacency Matrix (A)**
    
    - **Transformation**: Correlations are transformed into **weighted connections** using a **power function (β)**. This is known as **soft-thresholding**, which highlights strong connections and reduces the impact of weak correlations.
    - **Formula**: `aij = |cor(Xi, Xj)|^β`
        - `β` is chosen to achieve **scale-free topology**.
    - **Output**: A symmetric `n × n` matrix `A = [aij]`, where **`aij ∈ [0, 1]`**.
    - `[insert figure: Adjacency matrix at different power values (β)]`
4. **Topological Overlap Matrix (TOM)**
    
    - **Measure of Similarity**: TOM measures the **topological similarity** between genes, considering not only direct connections but also the **number of shared neighbors**. This helps in **reducing noise** from spurious correlations.
    - **Formula**: `TOMij = (Σu aiuauj + aij) / (min(ki, kj) + 1 - aij)`
        - **`Σu aiuauj` represents shared neighbors; `ki = Σu aiu` is connectivity of gene `i`.**
        - `min(ki, kj)` normalizes by connectivity.
    - **Output**: A symmetric `n × n` matrix `TOM = [TOMij]`, where **`TOMij ∈ [0, 1]`**.
    - `[insert figure: Network with Topological Overlap Matrix (TOM) simulation of i-j connections]`
5. **Dissimilarity TOM (dissTOM)**
    
    - **Transformation for Clustering**: The `TOM` matrix is converted into a dissimilarity measure for easier clustering.
    - **Formula**: `dissTOMij = 1 - TOMij`
    - **Output**: A symmetric `n × n` matrix `dissTOM = [dissTOMij]`, where **`dissTOMij ∈ [0, 1]`**.
6. **Clustering Modules**
    
    - **Hierarchical Clustering**: `dissTOM` is used as a distance measure for hierarchical clustering, constructing a dendrogram.
    - **Dynamic Tree Cutting**: This method identifies **modules** (groups of genes with low `dissTOM`) based on parameters like `minModuleSize` and `deepSplit`. Genes with low `dissTOM` (high `TOM`) are clustered together.
    - **Eigengene Calculation**: For each module, a **module eigengene (ME)** is calculated. The ME is the **first principal component (PC1)** of the expression profiles of all genes within that module. It **summarizes the module's overall expression pattern** in a single vector, capturing the dominant trend of gene expression.
    - **Module Merging**: Modules can be merged if their eigengenes are highly correlated (**e.g., >0.75, corresponding to mergeCutHeight=0.25**).

### Module Analysis and Trait Correlation

After module detection, WGCNA allows for the **correlation of module eigengenes with external biological traits** (e.g., disease status, treatment) to identify biologically relevant modules. Modules highly correlated with a trait are likely involved in that trait's biology.

Key concepts for assessing gene importance within modules and their relation to traits are:

- **Gene Significance (GS)**: Measures how strongly a gene's expression correlates with a specific trait.
    - **Formula**: `GSi = |cor(Xi, trait)|`
- **Module Membership (MM)**: Quantifies how central a gene is to its module by measuring the correlation between a gene's expression and its module's eigengene.
    - **Formula**: `MMi = |cor(Xi, ME_module)|`
    - Genes with **high GS and high MM** in a trait-correlated module are considered key candidates for biological roles.
- **Hub Genes**: These are genes that are **highly connected within a module**. They are identified by selecting genes with high MM (e.g., top 10%) or high **intramodule connectivity**. Hub genes are often functionally important, such as transcription factors in biological pathways.
    - **Intramodule Connectivity Formula**: `kwithin,i = Σj∈module,j≠i aij` (**within-module, vs. total ki = Σj aij**)
    - `[insert figure: Example of hub genes in a module]`

### Biological Interpretation

Once modules and hub genes are identified, further biological insights can be gained through:

- **Gene Ontology (GO) enrichment analysis**.
- **Identification of transcription factors**.
- **KEGG pathways overrepresentation**.
- **Functional annotation** (e.g., BLAST) and **domains/motifs identification**.

### Network Visualization

Network visualization tools are used to represent the inferred networks. Common steps include:

1. **Exporting** node (genes with MM, GS) and edge (adjacency/TOM) files.
2. **Visualizing** networks, often coloring nodes by module or connectivity.
3. **Annotating** with functional terms (e.g., GO/KEGG).

Tools for WGCNA network visualization and analysis include **Cytoscape**, igraph, NetworkX, Graph-tool, Gephi, and specialized WGCNA modules in platforms like Omics Playground or PyWGCNA.

### WGCNA Applications

WGCNA is widely applied across various scientific and biological systems. Examples include:

- **Animal Science**: Identifying hub genes in bovine mammary glands during lactation.
- **Human Health (Oncology)**: Clustering differentially expressed genes (DEGs) into modules and identifying hub genes linked to prognosis in cancers like human tongue squamous cell carcinoma, oral squamous cell carcinoma, osteosarcoma, and lung adenocarcinoma. It has also been used with microbiome data in tumor contexts.
- **Comparative Genomics**: Comparing WGCNA networks across human and mouse tissues to identify conserved gene roles.
- **Plant Biology**: Identifying hub genes related to pathogen resistance in cassava, drought tolerance in maize, **and somatic embryogenesis in Arabidopsis**.

### Co-expression Networks vs. Regulatory Networks

It is crucial to distinguish between **co-expression networks** and **gene regulatory networks**:

- **Co-expression Networks**: Edges are based on **correlation of expression profiles**, implying **no causality**. They are **undirected** and group genes with similar expression patterns.
- **Regulatory Networks**: Edges represent **causal interactions** (e.g., a transcription factor binding to a gene promoter). They are **directed** and model mechanistic control. Inference often relies on binding data or motif databases. Tools like ARACNe, GENIE3, and SCENIC are used for regulatory network inference.

Both network types share similarities: they represent biological relationships as graphs, identify functional modules, are data-driven, and utilize visualization tools.

---
---

# WGCNA: From Gene Expression to Functional Co-expression Networks

### Introduction to Network Biology

Biological systems are frequently represented as **networks**, which are complex sets of binary interactions or relations between different entities. Network biology applies **graph theory** tools to represent and analyze these systems. In this context, **nodes** represent entities like proteins or genes, and **edges** convey information about the links between them. Networks can have different types of edges: **undirected**, **directed**, or **weighted**.

Common types of biological networks include protein-protein interaction networks, metabolic networks, genetic interaction networks, gene/transcriptional regulatory networks, cell signaling networks, and **gene co-expression networks**. Gene co-expression networks are primarily associated with **Transcriptomics** data, such as RNA-seq or Microarrays.

### WGCNA: A Standard Tool for Co-expression Networks

**WGCNA (Weighted Gene Co-expression Network Analysis)** is a standard tool used to construct co-expression networks from gene expression data. It is "Weighted" because it uses **continuous edge weights ($[0, 1]$)**, "Correlation-based" as it utilizes Pearson/Spearman correlations, and designed for "Scale-free topology". WGCNA is **noise-robust** due to the use of Topological Overlap Matrix (TOM) and is **scalable** for large datasets. A key limitation is that it infers **no causality** and can be sensitive to batch effects. Co-expression networks are typically **undirected** (symmetric gene relationships).

### WGCNA Workflow

The WGCNA workflow transforms expression data into networks through a series of steps:

1. **Expression Matrix ($X$) Preprocessing**
    
    - **Input**: A matrix of gene expression values (e.g., RNA-seq counts).
    - **Preprocessing**: Normalization (e.g., VST, Rlog, CPM, TMM, or $\log_2(\text{TPM}+1)$) and filtering of low-variance genes.
    - **Output**: A cleaned matrix $X$ where $X_{i,s}$ (expression of gene $i$ in samples $s$).
    - **Formula**:  
        $$  
        X = [X_{i,s}], \quad i=1,\dots,n, \quad s=1,\dots,m  
        $$  
        where $n$ is the number of genes, $m$ is the number of samples.
Example: 

| Gene   | Sample 1 | Sample 2 | Sample 3 | Sample 4 | Sample 5 |
| ------ | -------- | -------- | -------- | -------- | -------- |
| Gene A | 3.745401 | 9.507143 | 7.319939 | 5.986585 | 1.560186 |
| Gene B | 1.559945 | 0.580836 | 8.661761 | 6.011150 | 7.080726 |
| Gene C | 0.205845 | 9.699099 | 8.324426 | 2.123391 | 1.818250 |
| Gene D | 1.834045 | 3.042422 | 5.247564 | 4.319450 | 2.912291 |
| Gene E | 6.118529 | 1.394939 | 2.921446 | 3.663618 | 4.560700 |

2. **Correlation Matrix ($R$)**
    
    - **Computation**: Pairwise correlations between genes are calculated to measure co-expression. **Pearson correlation** is commonly used.
    - **Formula**:  
        $$  
        r_{ij} = \frac{\text{cov}(X_i, X_j)}{\sigma_{X_i} \sigma_{X_j}}  
        $$  
        where $X_i$, $X_j$ are expression profiles of genes $i$, $j$; $\text{cov}$ is covariance; $\sigma_{X_i}$, $\sigma_{X_j}$ are standard deviations.
    - **Output**: A symmetric $n \times n$ matrix $R = [r_{ij}]$, where $r_{ij} \in [-1, 1]$.
    
    Example of a correlation matrix (*R*)![[matrixCor_heatmap.png]]
3. **Adjacency Matrix ($A$)**
    
    - **Transformation**: Correlations are transformed into **weighted connections** using a **power function ($\beta$)**. This is known as **soft-thresholding**, which highlights strong connections and reduces the impact of weak correlations.
    - **Formula**:  
        $$  
        a_{ij} = \left| \text{cor}(X_i, X_j) \right|^\beta  
        $$  
        where $\beta$ is chosen to achieve **scale-free topology**.
    - **Output**: A symmetric $n \times n$ matrix $A = [a_{ij}]$, where $a_{ij} \in [0, 1]$.
    
     Exemplification of a adjiacency matrix (*A*) calculation.
     ![[matrixAdj_hetmap.png]]
4. **Topological Overlap Matrix (TOM)**
    
    - **Measure of Similarity**: TOM measures the **topological similarity** between genes, considering not only direct connections but also the **number of shared neighbors**. This helps in **reducing noise** from spurious correlations.
    - **Formula**:  
        $$  
        TOM_{ij} = \frac{\sum_{u \neq i,j} a_{iu}a_{uj} + a_{ij}}{\min(k_i, k_j) + 1 - a_{ij}}  
        $$  
        where $\sum_{u \neq i,j} a_{iu}a_{uj}$ represents shared neighbors; $k_i = \sum_u a_{iu}$ is connectivity of gene $i$.
    - **Output**: A symmetric $n \times n$ matrix $TOM = [TOM_{ij}]$, where $TOM_{ij} \in [0, 1]$.
    
    Exemplification of a Network with Topological Overlap Matrix (TOM) simulation of i-j connections
    ![[TOM_simulation 1.png]]
5. **Dissimilarity TOM (dissTOM)**
    
    - **Transformation for Clustering**: The $TOM$ matrix is converted into a dissimilarity measure for easier clustering.
    - **Formula**:  
        $$  
        \text{dissTOM}_{ij} = 1 - TOM_{ij}  
        $$
    - **Output**: A symmetric $n \times n$ matrix $\text{dissTOM} = [\text{dissTOM}_{ij}]$, where $\text{dissTOM}_{ij} \in [0, 1]$.
6. **Clustering Modules**
    
    - **Hierarchical Clustering**: $\text{dissTOM}$ is used as a distance measure for hierarchical clustering, constructing a dendrogram.
    - **Dynamic Tree Cutting**: This method identifies **modules** (groups of genes with low $\text{dissTOM}$) based on parameters like $\text{minModuleSize}$ and $\text{deepSplit}$. Genes with low $\text{dissTOM}$ (high $TOM$) are clustered together.
    - **Eigengene Calculation**: For each module, a **module eigengene (ME)** is calculated. The $ME$ is the **first principal component (PC1)** of the expression profiles of all genes within that module. It **summarizes the module's overall expression pattern** in a single vector, capturing the dominant trend of gene expression.
    - **Formula**:  
        $$  
        ME = \sum_{i \in \text{module}} w_i X_i  
        $$  
        where $w_i$ are PCA loadings, $X_i$ are standardized expressions.
    - **Module Merging**: Modules can be merged if their eigengenes are highly correlated (e.g., $>0.75$, **corresponding to $\text{mergeCutHeight}=0.25$**).
    
	Exemplification of module formation.
    ![[modules_dendro_merged.png]]

### Module Analysis and Trait Correlation

After module detection, WGCNA allows for the **correlation of module eigengenes with external biological traits** (e.g., disease status, treatment) to identify biologically relevant modules. Modules highly correlated with a trait are likely involved in that trait's biology.

![[Trait-moduleCor.png]]

Key concepts for assessing gene importance within modules and their relation to traits are:

- **Gene Significance ($GS$)**: Measures how strongly a gene's expression correlates with a specific trait.
    - **Formula**:  
        $$  
        GS_i = \left| \text{cor}(X_i, \text{trait}) \right|  
        $$
- **Module Membership ($MM$)**: Quantifies how central a gene is to its module by measuring the correlation between a gene's expression and its module's eigengene.
    - **Formula**:  
        $$  
        MM_i = \left| \text{cor}(X_i, ME_{\text{module}}) \right|  
        $$
    - Genes with **high $GS$ and high $MM$** in a trait-correlated module are considered key candidates for biological roles.
- **Hub Genes**: These are genes that are **highly connected within a module**. They are identified by selecting genes with high $MM$ (e.g., top 10%) or high **intramodule connectivity**. Hub genes are often functionally important, such as transcription factors in biological pathways.
    - **Intramodule Connectivity Formula**:  
        $$  
        k_{\text{within},i} = \sum_{j \in \text{module}, j \neq i} a_{ij}  
        $$  
        (**within-module, vs. total $k_i = \sum_j a_{ij}$**)
        
        Example of a network with hub genes identified:
		![[hub_genes.png]]

### Biological Interpretation

Once modules and hub genes are identified, further biological insights can be gained through:

- **Gene Ontology (GO) enrichment analysis**.
- **Identification of transcription factors**.
- **KEGG pathways overrepresentation**.
- **Functional annotation** (e.g., BLAST) and **domains/motifs identification**.

### Network Visualization

Network visualization tools are used to represent the inferred networks. Common steps include:

1. **Exporting** node (genes with $MM$, $GS$) and edge (adjacency/TOM) files.
2. **Visualizing** networks, often coloring nodes by module or connectivity.
3. **Annotating** with functional terms (e.g., GO/KEGG).

Tools for WGCNA network visualization and analysis include **Cytoscape**, igraph, NetworkX, Graph-tool, Gephi, and specialized WGCNA modules in platforms like Omics Playground or PyWGCNA.

### WGCNA Applications

WGCNA is widely applied across various scientific and biological systems. Examples include:

- **Animal Science**: Identifying hub genes in bovine mammary glands during lactation.
- **Human Health (Oncology)**: Clustering differentially expressed genes (DEGs) into modules and identifying hub genes linked to prognosis in cancers like human tongue squamous cell carcinoma, oral squamous cell carcinoma, osteosarcoma, and lung adenocarcinoma. It has also been used with microbiome data in tumor contexts.
- **Comparative Genomics**: Comparing WGCNA networks across human and mouse tissues to identify conserved gene roles.
- **Plant Biology**: Identifying hub genes related to pathogen resistance in cassava, drought tolerance in maize, **and somatic embryogenesis in Arabidopsis**.

### Co-expression Networks vs. Regulatory Networks

It is crucial to distinguish between **co-expression networks** and **gene regulatory networks**:

- **Co-expression Networks**: Edges are based on **correlation of expression profiles**, implying **no causality**. They are **undirected** and group genes with similar expression patterns.
- **Regulatory Networks**: Edges represent **causal interactions** (e.g., a transcription factor binding to a gene promoter). They are **directed** and model mechanistic control. Inference often relies on binding data or motif databases. Tools like ARACNe, GENIE3, and SCENIC are used for regulatory network inference.

Both network types share similarities: they represent biological relationships as graphs, identify functional modules, are data-driven, and utilize visualization tools.

---