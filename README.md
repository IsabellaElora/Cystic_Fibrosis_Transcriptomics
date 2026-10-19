# Cystic Fibrosis Transcriptomic Analysis

## Project Overview

This project investigates differential gene expression between cystic fibrosis (CF) and non-CF samples using RNA sequencing data and computational bioinformatics approaches.

The analysis was conducted using R/RStudio, incorporating differential expression analysis, exploratory data analysis, statistical visualisation and sensitivity analysis.

## Research Objectives

- Identify genes showing differential expression between CF and non-CF samples.
- Investigate disease-associated transcriptional changes.
- Assess the robustness of differential expression results to the exclusion of an individual sample.
- Explore selected airway-associated genes and their potential biological relevance.
- Visualise CFTR variants in relation to functional protein domains.

## Dataset

**Source:** NCBI Gene Expression Omnibus (GEO)

**Accession:** GSE154905

**Comparison:** CF versus non-CF samples under UNC medium conditions.

## Computational Methods

**R and RStudio**

- DESeq2 for differential gene expression analysis.
- Normalisation and filtering of RNA-seq count data.
- Interpretation of log2 fold changes and adjusted p-values.
- Principal component analysis (PCA) for exploratory sample assessment.
- MA plots, volcano plots and gene-level visualisations.
- Sensitivity analysis comparing primary and sample-exclusion results.
- CFTR variant retrieval and protein-domain mapping.

## Selected Findings and Visualisations

The project includes visualisations of differential gene expression, selected airway-related genes (DEFB1, IL20RB and SCGB3A1), and CFTR variants.

The sensitivity analysis examines how excluding one non-CF sample influences estimated gene expression changes.

## Research Skills Demonstrated

- R programming and bioinformatics
- RNA-seq differential expression analysis
- Statistical interpretation and data visualisation
- Exploratory data analysis and quality control
- Sensitivity analysis and reproducibility assessment
- Integration of genomic annotations with protein-domain information

## Limitations

Differential gene expression does not establish causality. Results are specific to the analysed dataset, experimental conditions and statistical design.

This repository contains an exploratory research workflow 
