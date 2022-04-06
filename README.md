# Project Description

This repository contains the analysis for BF528 Project 3 - Concordance of microarray and RNA-Seq differential gene expression.

This project seeks to reproduce a subset of analyses that were performed in Wang et al. 2014. The analyses done in that paper had several goals: to characterize the concordance of differential gene expression across platforms, test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s mechanism of action, and assess the mechanism of action prediction accuracy of each platform using a test set.



# Contributors

 - Data Curator

 - Programmer

Michael Peters - Analyst

 - Biologist 

# Repository Contents

## Step 1
text
```
code
```

## Step 2
text
```
code
```

## Step X
The `analyst.R` script contains all of the code needed for the analyst role. This script runs DEG analysis using _limma_, and computes concordance between samples.
To run this script, from the parent directory of repository, run:
```
module load R
Rscript analyst/analyst.R
```
