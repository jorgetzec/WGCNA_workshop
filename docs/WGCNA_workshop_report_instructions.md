# WGCNA Workshop Report Instructions

Using the provided scripts [`WGCNA_script1.R`](../code/WGCNA_script1.R) and [`WGCNA_script2.R`](../code/WGCNA_script2.R), generate a research report that addresses the following questions and prompts. Your report should include figures, tables, and clear explanations for each point.

## Report Guidelines

1. **Filtering Criteria**
   - What was the best filtering criterion for this dataset?  
   - Is it the most recommended for this type of data? Why or why not?

2. **Gene Filtering Outcome**
   - How many genes were removed after pre-filtering?  
   - Provide before-and-after gene counts.

3. **Sample Analysis**
   - What information did the sample analysis provide?  
   - Was any sample removed from the analysis? If so, why?

4. **Network Construction Parameters**
   - What power (β) was chosen to achieve scale-free topology?  
   - What type of network was selected for the analysis (signed, unsigned, or hybrid), and why?

5. **MYB Transcription Factors**
   - In which modules are the MYB genes located?  
   - What is the expression profile of the MYB genes (upregulated, downregulated, or no major changes)?

6. **MYB Gene Metrics**
   - What are the connectivity and module membership metrics for the MYB genes within their modules?

7. **Module Exploration**
   - Choose one module to analyze in detail.
   - Identify the hub genes in this module.
   - Visualize the network of the selected module (include a network plot).

8. **Biological Interpretation**
   - What do you think is the role of the MYB transcription factor found in your selected module?  
   - Support your answer with evidence from your analysis and relevant literature if possible.

---

## Report Format

- **Introduction:** Briefly describe the dataset and the goal of the analysis.
- **Methods:** Summarize the main steps and parameters used (referencing the scripts).
- **Results:** Address each question above, including figures and tables as needed.
- **Discussion:** Interpret your findings, especially regarding the MYB transcription factors and the selected module.
- **References:** List any literature or resources you used.

