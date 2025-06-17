# Syllabus: Gene Co-expression Network Analysis with WGCNA

### **Session 1: Introduction to Network Biology and WGCNA Fundamentals (2 hours)**

1. Introduction to Network Biology 
    1.1. Core Concepts of Networks  
    1.2. Graph Theory and Network Components  
    1.3. Network Properties and Metrics  
    1.4. Common Types of Biological Networks and Associated Omics Data  

2. Fundamentals of WGCNA (Weighted Gene Co-expression Network Analysis) 
    2.1. WGCNA as a Standard Tool  
    2.2. Features and Limitations of WGCNA  
    2.3. WGCNA Workflow: General Overview  
    2.4. Step-by-Step Network Construction: Part I  
      2.4.1. Expression Matrix (X)  
      2.4.2. Correlation Matrix (R)  
      2.4.3. Adjacency Matrix (A)  


### **Session 2: WGCNA Application, Analysis, and Tools (2 hours)**

3. Detailed Construction and Analysis of the WGCNA Network  
    3.1. Step-by-Step Network Construction: Part II  
      3.1.1. Topological Overlap Matrix (TOM)  
      3.1.2. Dissimilarity Matrix (dissTOM)  
    3.2. Module Clustering  
    3.3. Module-Trait Correlation and Key Genes  
    3.4. Hub Gene Identification  
    3.5. Functional Analysis of Modules  

4. WGCNA Applications and Tools  
    4.1. Network Visualization  
    4.2. WGCNA Applications in Various Fields  
    4.3. WGCNA and Multi-Omics Integration  
    4.4. Tools Similar to or Based on WGCNA  
    4.5. WGCNA vs. Gene Regulatory Networks  


### **Session 3: Practical WGCNA Analysis - Part I (2 hours)**

5. Data Preparation and Preprocessing
    5.1. Setting Up the R Environment
    5.2. Data Loading and Quality Control
    5.3. Data Filtering and Normalization
        5.3.1. Low Count Filtering
        5.3.2. Variance Stabilization
    5.4. Outlier Detection and Removal

6. Network Construction
    6.1. Data Transformation for WGCNA
    6.2. Soft Threshold Selection
        6.2.1. Power Selection Analysis
        6.2.2. Scale-Free Topology Assessment
    6.3. Network Construction
    6.4. Module Detection


### **Session 4: Practical WGCNA Analysis - Part II (2 hours)**

7. Module Analysis and Visualization
    7.1. Module Eigengene Analysis
        7.1.1. Eigengene Calculation
        7.1.2. Module-Trait Relationships
    7.2. Hub Gene Identification
        7.2.1. Connectivity Analysis
        7.2.2. Top Hub Genes Selection
    7.3. Module Profile Analysis
        7.3.1. Expression Pattern Visualization
        7.3.2. Gene Significance Assessment

8. Network Visualization and Export
    8.1. Module Visualization
    8.2. Cytoscape Integration
    8.3. Results Interpretation
    8.4. Results Export and Documentation
