# 🧬 PTAFR Expression and Tumor Microenvironment in Carcinomas

## 📌 Abstract

Solid tumor carcinogenesis is driven by complex interactions within the tumor microenvironment (TME), which influence tumor initiation, progression, and response to treatment. The phospholipid platelet-activating factor (PAF) has been described to affect TME and growth; however, evidence so far has been limited to experimental in vitro and in vivo models. This study investigates the role of PAF receptor (PTAFR) across various carcinomas through patient data-driven analysis. We analyzed RNA sequencing and clinical data from publicly available tumor tissue samples of 8,977 carcinoma patients in the Genomic Data Commons. Our findings revealed that PTAFR overexpression correlates with reduced patient survival, increased tumor size, and enhanced lymph node metastasis. Further analysis identified 151 differentially expressed genes associated with PTAFR, many of which regulate immune responses and proto-oncogenic signaling. Histological and single-cell analyses demonstrated an immunosuppressive TME in carcinoma tissues with elevated PTAFR expression. In vivo validation indicated reduced tumor growth and increased infiltration of pro-inflammatory cells in PTAFR knockout mice. Protein modeling suggests that PTAFR forms physical complexes with immunosuppressive mediators STAT3 and PD-L2, and that pathogenic variants may disrupt the PAF-binding domain, potentially altering its signaling. Together, our findings establish PTAFR as a critical regulator of tumor progression, directly linking its expression to immunosuppressive signaling in the TME, enhanced tumor growth, and reduced patient survival, underscoring its potential as a therapeutic target.

## 👩‍🔬 Authors

**Barbara Dalmaso**¹², Ildefonso Alves da Silva-Junior², Sonia Jancar², Laura Steenpass¹, Carolina Beltrame Del Debbio², Claudia Pommerenke¹  

¹ Leibniz-Institute DSMZ – German Collection of Microorganisms and Cell Cultures GmbH, Braunschweig, Germany  
² Institute of Biomedical Sciences, University of São Paulo, São Paulo, Brazil

## 🔗 Publication

**Journal**: Experimental Cell Research (2025)  
**URL**: https://www.sciencedirect.com/science/article/abs/pii/S0014482725002435  
**DOI**: https://doi.org/10.1016/j.yexcr.2025.114643

## 💾 Data Availability

- Bulk RNA-seq and clinical data: [Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/)  
- Single-cell and histology data: Available via referenced public datasets in the manuscript  
- Processed data: Available upon request (see contact below)

## 🧰 Tools and Libraries

### 📦 Installation

You can install the necessary R packages with the following commands:

```r
install.packages(c("ggplot2", "patchwork", "survival", "survminer"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))
# Optional: for protein structure modeling
BiocManager::install("alphafoldR")
install.packages("bio3d")
```

- R version 4.3+  
- `DESeq2` – Differential gene expression analysis  
- `ggplot2`, `patchwork` – Data visualization  
- `enrichplot`, `clusterProfiler` – Functional enrichment and gene ontology  
- `survival`, `survminer` – Survival analysis  
- `BiocManager`, `org.Hs.eg.db` – Annotation and biological databases  
- `alphafoldR`, `bio3d` – Protein structure modeling and visualization

## 📁 Project Structure

```
📄 PTAFR-tumor-analysis.R # R Scripts
📄 README.md           # Project documentation  
```

## 🧠 Citation

If you use this code or data, please cite:

> Dalmaso B, da Silva-Junior IA, Jancar S, Steenpass L, Del Debbio CB, Pommerenke C.  
> "PTAFR promotes tumor growth and shapes an immunosuppressive microenvironment in carcinomas."  
> *Experimental Cell Research* (2025). https://doi.org/10.1016/j.yexcr.2025.114643

## 📧 Contact

**Email**: barbdalmaso@gmail.com

Feel free to reach out with questions, feedback, or collaboration proposals.

## 📜 License

This repository is shared for academic and non-commercial use.  
For reuse or collaboration, in-house generated data and full article are **available upon request**.
