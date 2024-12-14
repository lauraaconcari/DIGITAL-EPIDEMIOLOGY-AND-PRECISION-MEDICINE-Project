# DIGITAL-EPIDEMIOLOGY-AND-PRECISION-MEDICINE-Project

## Differential Analyses of Gene Expression

### Master Degree: Data Science (Academic Year: 2024/2025)

---

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Data Sources](#data-sources)
- [Technologies Used](#technologies-used)
- [Setup Instructions](#setup-instructions)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Biological and Clinical Implications](#biological-and-clinical-implications)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

This project explores the molecular landscape of **Lung Adenocarcinoma (LUAD)** by analyzing transcriptomic and mutational data from The Cancer Genome Atlas (TCGA). The main goal is to investigate the biological mechanisms underlying LUAD progression, identify molecular subtypes, and propose potential pathways for personalized treatments. This study integrates differential gene expression analysis, network science, and pathway enrichment to provide a comprehensive understanding of LUAD.

### Objectives:
- Identify key **Differentially Expressed Genes (DEGs)** to differentiate tumor tissues from normal tissues.
- Classify patients into **distinct molecular subtypes** using co-expression and similarity networks.
- Perform pathway enrichment to link patient-specific biological mechanisms to tumor subtypes.
- Integrate mutational and transcriptomic data using **Similarity Network Fusion (SNF)**.

### Biological Insights:
To uncover LUAD subtypes, patients with high expression levels of **IL6** (a marker of the inflammatory subtype) and **SFTPC** (a marker of the terminal respiratory unit subtype) were identified. These genes represent key biological pathways:
- **IL6**: Associated with tumor-promoting inflammation and the proximal inflammatory subtype.
- **SFTPC**: Linked to terminal respiratory unit cells, representing less aggressive and more differentiated tumors.

Community and cluster analyses using Patient Similarity Networks (PSN) and the Fused Network revealed distinct patient groupings, each enriched for specific biological pathways. These findings were critical in associating pathway enrichment results with LUAD molecular subtypes.

---

## Features

- **Differential Expression Analysis**:
  - Identify DEGs between tumor and normal tissues using statistical tests.
  - Map DEGs to external gene names and annotate key markers (e.g., IL6, SFTPC, KRAS, EGFR).

- **Patient Subtyping**:
  - Use **Patient Similarity Networks (PSN)** and clustering techniques to identify patient groups.
  - Assign patients to communities based on shared molecular signatures.

- **Pathway Enrichment**:
  - Perform **Gene Ontology (GO)** enrichment analysis to link subtypes to key biological processes.
  - Visualize and interpret pathways enriched in PSN communities and Fused Network clusters.

- **Network Analysis**:
  - Construct and analyze **Co-expression Networks** to identify hubs and regulatory relationships.
  - Integrate mutational and transcriptomic data using **SNF** for a holistic view.

- **Biological and Clinical Interpretation**:
  - Identify LUAD subtypes (e.g., proliferative, inflammatory, terminal respiratory unit).
  - Propose therapeutic targets based on enriched pathways for specific subtypes.

---

## Data Sources

1. **Transcriptomic Data**:
   - Retrieved from the **TCGA-LUAD dataset** using the `TCGAbiolinks` package.
   - Includes raw gene expression profiles for tumor and normal samples.

2. **Mutational Data**:
   - Processed into binary mutation matrices for similarity analysis using `maftools`.

3. **Pathway Annotations**:
   - GO Biological Processes from the Ensembl database and `enrichR`.

---

## Technologies Used

- **Programming Language**: R
- **Libraries**:
  - `TCGAbiolinks` for data retrieval.
  - `DESeq2` for DEG analysis.
  - `igraph` for network construction.
  - `enrichR` for pathway enrichment analysis.
  - `SNFtool` for similarity network fusion.
  - `maftools` for mutational data analysis.
  - `ggplot2` for data visualization.

---

## Setup Instructions

### Prerequisites

1. Install **R (>= 4.0)** and **RStudio**.
2. Ensure the following R libraries are installed:
   ```R
   install.packages(c("BiocManager", "psych", "network", "ggplot2"))
   BiocManager::install(c("TCGAbiolinks", "DESeq2", "maftools", "enrichR", "SNFtool"))
   ```

### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/DIGITAL-EPIDEMIOLOGY-AND-PRECISION-MEDICINE-Project.git
   ```
2. Navigate to the project directory:
   ```bash
   cd DIGITAL-EPIDEMIOLOGY-AND-PRECISION-MEDICINE-Project
   ```
3. Open and execute the R script:
   ```bash
   Rscript PM_Proj_Code.R
   ```

---

## Usage

1. **Preprocess Data**:
   - Normalize and clean raw gene expression data.
2. **Perform DEG Analysis**:
   - Identify significant DEGs between tumor and normal tissues.
3. **Construct Networks**:
   - Build and analyze Patient Similarity Networks and Co-expression Networks.
4. **Pathway Enrichment**:
   - Perform GO enrichment analysis on DEGs and visualize results.
5. **Subtyping and Interpretation**:
   - Link enriched pathways to LUAD subtypes for personalized insights.

---

## Project Structure

```
DEPM-Project/
├── PM_Proj_Code.R            # Main R script
├── DEPM-1_proj_group-10.pdf  # Project report
├── README.md                 # Documentation
├── data/                     # Raw and processed data
├── results/                  # Output figures and tables
├── references/               # Literature and documentation
```

---

## Biological and Clinical Implications

This project demonstrates the power of omics data in unraveling LUAD heterogeneity. The pathway enrichment results show that communities and clusters identified from network analyses align with key biological subtypes: proliferative, inflammatory, and terminal respiratory unit. Enriched pathways such as **spindle assembly checkpoint signaling** (proliferative) and **negative regulation of cell division** (terminal respiratory unit) suggest actionable biological targets. These findings can guide the development of targeted therapies, such as inhibitors for cell cycle dysregulation in the proliferative subtype or immunomodulators for the inflammatory subtype.

Personalized prognosis is enhanced by linking patients' molecular profiles to clinical behavior. This integrative approach supports precision medicine, paving the way for better therapeutic strategies for LUAD.

---

## Contributing

Contributions are welcome! Follow these steps:

1. Fork the repository.
2. Create a feature branch:
   ```bash
   git checkout -b feature-name
   ```
3. Commit your changes:
   ```bash
   git commit -m "Add feature description"
   ```
4. Push your changes:
   ```bash
   git push origin feature-name
   ```
5. Submit a pull request for review.

---

## License

This project is licensed under the [MIT License](LICENSE). Use responsibly.

