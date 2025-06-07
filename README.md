# 🧬 PTAFR Expression and Tumor Microenvironment in Carcinomas

## 📌 Abstract

Solid tumor carcinogenesis is driven by complex interactions within the tumor microenvironment (TME), which influence tumor initiation, progression, and response to treatment. The phospholipid platelet-activating factor (PAF) has been described to affect TME and growth; however, evidence so far has been limited to experimental in vitro and in vivo models.

This study investigates the role of PAF receptor (PTAFR) across various carcinomas through patient data-driven analysis. We analyzed RNA sequencing and clinical data from publicly available tumor tissue samples of 8,977 carcinoma patients in the Genomic Data Commons.

Our findings revealed that PTAFR overexpression correlates with reduced patient survival, increased tumor size, and enhanced lymph node metastasis. Further analysis identified 151 differentially expressed genes associated with PTAFR, many of which regulate immune responses and proto-oncogenic signaling.

Histological and single-cell analyses demonstrated an immunosuppressive TME in carcinoma tissues with elevated PTAFR expression. In vivo validation indicated reduced tumor growth and increased infiltration of pro-inflammatory cells in PTAFR knockout mice. Protein modeling suggests that PTAFR forms physical complexes with immunosuppressive mediators STAT3 and PD-L2, and that pathogenic variants may disrupt the PAF-binding domain, potentially altering its signaling.

Together, our findings establish PTAFR as a critical regulator of tumor progression, directly linking its expression to immunosuppressive signaling in the TME, enhanced tumor growth, and reduced patient survival, underscoring its potential as a therapeutic target.

## 👩‍🔬 Authors

Barbara Dalmaso¹²  
Ildefonso Alves da Silva-Junior²  
Sonia Jancar²  
Laura Steenpass¹  
Carolina Beltrame Del Debbio²  
Claudia Pommerenke¹  

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

## 🔁 How to Reproduce

1. Clone this repository:
   ```bash
   git clone https://github.com/your-user/ptafr-tme-analysis.git
   cd ptafr-tme-analysis
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the analysis:
   Open and execute the notebooks in the `notebooks/` directory to replicate the core results.

## 🧰 Tools and Libraries

- Python 3.10  
- Pandas, NumPy, SciPy  
- Scanpy, Seaborn, Matplotlib  
- scikit-learn, XGBoost, LightGBM  
- Lifelines (survival analysis)  
- AlphaFold / ColabFold (protein modeling)

## 📁 Project Structure

```
📁 data/               # Processed input and metadata  
📁 notebooks/          # Analysis and visualization notebooks  
📁 figures/            # Figures and supplementary plots  
📁 structure_modeling/ # PTAFR structural predictions and interaction models  
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
