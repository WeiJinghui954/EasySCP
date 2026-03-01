# EasySCP
The widespread adoption of single-cell proteomics (SCP) in biology has been limited by complex workflows and reliance on specialized instrumentation. Here we present EasySCP, a high-throughput method that integrates FACS-based single-cell sorting, an all-in-one, single-step digestion process in 384-well plates, and sensitive mass spectrometry. EasySCP identifies nearly 5,000 proteins from individual HEK293 cell. Applied to murine liver, EasySCP achieved spatially resolved proteomics profiling of hepatocytes zonation, detecting an average of 3,500 proteins per hepatocyte and uncovering zonation patterns for 3,277 out of 5,267 proteins. Building on 215 conserved zonation markers, we further developed hepatocyte spatial status score (HSS) that enables accurately reconstruct liver zonation across single-cell and multi-omics datasets. Together, our study introduces EasySCP, a broadly accessible tool for dissecting cellular heterogeneity at single-cell proteomics resolution in both healthy and disease states, effectively bridging the gap between transcriptomics and functional proteomics.

This repository contains R scripts used to reproduce all results and figures in the manuscript.

## Data Availability

Due to GitHub file size limitations, the full dataset (data/ and output/ folders) is publicly available at:

https://doi.org/10.5281/zenodo.18824044

Please download the dataset before running any scripts.

## Project Structure

After downloading the Zenodo archive, your project directory must contain the following three folders at the same level:

project/
├── Data/
├── Output/
└── R_scripts/

All three folders must coexist in the same root directory.  
The scripts in `R_scripts/` rely on relative paths and will not run if the structure is modified.

## How to Reproduce the Results

### Step 1 – Download data

Download the dataset from Zenodo and extract it into the project root directory so that:

- `Data/` contains raw and intermediate data  
- `Output/` contains processed results  

### Step 2 – Open R

We recommend using:

- R >= 4.4.1
- RStudio (optional but recommended)

### Step 3 – Set working directory

In R:

```r
setwd("path_to_project_root")
