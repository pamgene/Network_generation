This project generates networks from UKA files based on the PCSF algorithm.

Steps to install kinograte

```
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")

BiocManager::install("org.Hs.eg.db")

devtools::install_github("IOR-Bioinformatics/PCSF")
devtools::install_github("kalganem/Kinograte")

```