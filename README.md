# Identifying Driver Genes from Single-Cell Gene Expression Data of Acute Myeloid Leukemia

This repository contains the code and analysis workflow for my Master’s thesis project: **“Identifying Driver Genes from Single-Cell Gene Expression Data of Acute Myeloid Leukemia.”**

The project focuses on applying **deep learning and bioinformatics methods** to single-cell RNA sequencing data from Acute Myeloid Leukemia (AML) patients. The main goal was to identify biologically relevant genes, Gene Ontology terms, and pathways that may be important in AML progression.

## Project Overview

Acute Myeloid Leukemia is a complex blood cancer with high genetic and biological heterogeneity. Single-cell RNA sequencing provides a detailed view of gene expression at the individual cell level, but the data is high-dimensional, noisy, and difficult to interpret directly.

In this project, I built a computational workflow that combines:

* Single-cell RNA-seq preprocessing
* Gene Ontology-based gene grouping
* Deep neural network-based gene expression inference
* Model evaluation using Mean Squared Error
* Bootstrapping and repeated sampling for robustness
* Gene Ontology and KEGG pathway enrichment analysis
* Biological interpretation of candidate AML-related genes and pathways

The project connects machine learning results with biological interpretation instead of treating the neural network as a black box.

## Main Research Goal

The main goal of this project was to investigate whether deep neural networks can help identify important genes and biological processes related to AML by learning from single-cell gene expression patterns.

More specifically, the project aimed to:

* Preprocess raw single-cell RNA-seq data into a machine-learning-ready format
* Train neural network models for gene expression inference
* Evaluate model performance across different training sample sizes
* Rank Gene Ontology terms using prediction error
* Identify biological processes and pathways relevant to AML
* Interpret candidate genes and pathways using GO and KEGG enrichment analysis

## Dataset

The project used publicly available single-cell RNA-seq data from Acute Myeloid Leukemia patients.

The raw sequencing files were processed using a bioinformatics pipeline before being used for machine learning analysis. The data included thousands of genes and single-cell expression profiles from AML samples.

Due to dataset size and access restrictions, raw data files may not be included directly in this repository. The repository focuses on the analysis code, workflow, and methodology.

## Methodology

The project followed a multi-step workflow.

### 1. Data Preprocessing

Raw single-cell RNA-seq data was first processed to generate a gene expression count matrix.

The preprocessing steps included:

* Processing raw FASTQ files
* Generating count matrices
* Filtering low-quality cells
* Removing cells with high mitochondrial content
* Normalizing gene expression values
* Scaling the data
* Preparing the dataset for downstream machine learning

Tools used in this stage included:

* STARsolo
* Seurat
* R
* Python

### 2. Gene Ontology Mapping

After preprocessing, genes were mapped to their corresponding Gene Ontology terms.

This step helped group genes according to biological processes. Instead of analysing genes only as isolated features, the project used GO terms to connect gene expression patterns with biological meaning.

### 3. Deep Learning Model

A deep neural network was trained to perform gene expression inference.

The model was trained using different proportions of the dataset, starting from smaller training samples and gradually increasing the amount of data. The purpose was to study how training data size affected model performance.

The model performance was evaluated using **Mean Squared Error (MSE)**.

Lower MSE values indicated that the model was better able to infer gene expression patterns for genes associated with specific GO terms.

### 4. Bootstrapping and Robustness Testing

To avoid relying on a single train-test split, I used repeated sampling and bootstrapping.

This helped test whether the results were stable across different random samples of the data. The model was trained multiple times, and the MSE values were averaged to get a more reliable performance estimate.

This was one of the important technical parts of the project because biomedical datasets can be noisy and sensitive to sampling.

### 5. GO and KEGG Enrichment Analysis

After ranking GO terms based on model performance, downstream biological analysis was performed.

This included:

* Gene Ontology enrichment analysis
* KEGG pathway enrichment analysis
* Gene-set overlap analysis
* Biological interpretation of enriched pathways

The goal was to identify biological processes and pathways that may be important in AML.

## Key Results

The analysis highlighted several biological processes and pathways relevant to AML, including:

* Translation
* Protein phosphorylation
* Immune response
* Apoptotic process
* Cell adhesion
* Protein ubiquitination
* Cell cycle-related processes
* Autophagy-related processes

KEGG pathway analysis also identified AML-related and cancer-related pathways such as:

* MAPK signaling pathway
* PI3K-Akt signaling pathway
* Ras signaling pathway
* FoxO signaling pathway
* Acute Myeloid Leukemia pathway

Several important AML-related genes appeared in the downstream analysis, including genes such as:

* FLT3
* CEBPA
* RUNX1
* KIT
* NRAS
* KRAS
* TP53-related pathways

These findings showed that the workflow could connect deep learning model outputs with biologically meaningful AML-related mechanisms.

## Technologies Used

### Programming Languages

* Python
* R

### Python Libraries

* TensorFlow
* Keras
* Pandas
* NumPy
* Scikit-learn

### R Libraries

* Seurat
* clusterProfiler
* dplyr
* ggplot2
* GOxploreR

### Bioinformatics Tools

* STARsolo
* Gene Ontology
* KEGG pathway analysis

## Repository Structure

```text
thesis/
├── data/
│   └── README.md
├── preprocessing/
│   └── single_cell_preprocessing.R
├── model/
│   └── neural_network_training.py
├── analysis/
│   ├── go_analysis.R
│   ├── kegg_analysis.R
│   └── overlap_analysis.R
├── results/
│   └── output_files/
├── figures/
│   └── thesis_figures/
├── README.md
└── requirements.txt
```

The exact file structure may differ depending on the current repository version, but the main workflow follows the same logic: preprocessing, model training, downstream analysis, and results interpretation.

## How to Run the Project

### 1. Clone the Repository

```bash
git clone https://github.com/aln1234/thesis.git
cd thesis
```

### 2. Create a Python Environment

```bash
python -m venv venv
source venv/bin/activate
```

On Windows:

```bash
venv\Scripts\activate
```

### 3. Install Python Dependencies

```bash
pip install -r requirements.txt
```

### 4. Install Required R Packages

Open R or RStudio and install the required packages:

```r
install.packages("Seurat")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
```

### 5. Run Preprocessing

Run the preprocessing scripts to prepare the single-cell RNA-seq expression matrix.

```r
source("preprocessing/single_cell_preprocessing.R")
```

### 6. Train the Deep Learning Model

```bash
python model/neural_network_training.py
```

### 7. Run Downstream Analysis

```r
source("analysis/go_analysis.R")
source("analysis/kegg_analysis.R")
source("analysis/overlap_analysis.R")
```

## Notes About Data Availability

The raw single-cell RNA-seq data files are large and may not be included directly in this repository. Users who want to reproduce the full analysis should download the original dataset from the relevant public biomedical data archive and follow the preprocessing steps described in the thesis.

Processed or example files may be included where possible to demonstrate the workflow.

## Limitations

This project was completed as a Master’s thesis, so the main focus was on developing and evaluating a computational framework rather than producing a clinically validated diagnostic model.

Some limitations include:

* The analysis was based on a limited AML dataset
* Experimental biological validation was not performed
* Model interpretation depends on GO and KEGG annotations
* Further testing on independent patient datasets would be needed
* The workflow is research-oriented and not intended for clinical use

## What I Learned

Through this project, I learned that deep learning in biomedical research is not only about model architecture. Data preprocessing, robustness testing, biological interpretation, and careful evaluation are equally important.

The project also helped me understand how machine learning outputs can be connected with domain knowledge, especially in complex biomedical problems where interpretability matters.

## Thesis Information

**Title:** Identifying Driver Genes from Single-Cell Gene Expression Data of Acute Myeloid Leukemia
**Degree:** Master’s in Computing Sciences, Data Science
**University:** Tampere University
**Author:** Albin Lamichhane
**Year:** 2024

## Author

**Albin Lamichhane**
Software Developer / AI-focused Full-Stack Developer
Email: [albinlamichhane9@gmail.com](mailto:albinlamichhane9@gmail.com)
GitHub: https://github.com/aln1234
LinkedIn: https://www.linkedin.com/in/albinlamichhane/
Portfolio: https://albinlamichhane.vercel.app/

## Disclaimer

This repository is for academic and research purposes only. The results should not be interpreted as medical advice or clinical evidence without further validation.
