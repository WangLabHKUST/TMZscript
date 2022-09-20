# Pharmacogenomic profiling reveals molecular features of chemotherapy resistance in IDH wild type primary glioblastoma


# Directory structure
```
└─ data/:
|   └─ Supplementary File 1.xlsx: supplementary files and machine learning training datasets
|   └─ TCGA_138_zscored.xlsx: TCGA expression data with z-scored normalization
|   └─ TCGA_clinical.xlsx: TCGA clinical data 
|   └─ TCGA_testing.csv: machine learning test TCGA dataset
|   └─ pairs40_featureCounts.csv: 40 paired RNA-seq data in Fig. 3c
└─ scripts/:
|   └─ TMZ_scripts.ipynb: analysis code and scripts to reproduce figures
|   └─ TMZep.py: machine learning model
|   └─ bioNet.py: finding genes related to target
|   └─ cnape.ipynb: predicting CNA using RNA-seq data
|   └─ norm.r: normalization of RNA-seq data
|   └─ preprocess_tcga.r: merge TCGA CNA data and RNA data
```


