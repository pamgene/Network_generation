
# Network generation from UKA with PCSF

This project generates networks from UKA files and a subset of the STRING db based on the PCSF algorithm. Then it enriches the network with pathways with EnrichR and generates kinase-to-pathway heatmaps.

See details: UKA_networks_interpretation_guide.pdf

The project was developed based on [Kinograte](https://github.com/CogDisResLab/Kinograte).

See how STRING db was downloaded and filtered [here](https://github.com/pamgene/STRING_download).



**Steps to install PCSF:**

```
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")

BiocManager::install("org.Hs.eg.db")

devtools::install_github("IOR-Bioinformatics/PCSF")

```