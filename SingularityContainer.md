---
title: "SingularityContainer"
author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "30 January, 2025"
output: 
  html_document:
    keep_md: true
    df_print: paged
    toc: true
    toc_depth: 6
  pdf_document: default
editor_options:
  chunk_output_type: console
---



#### *1. Create a Singularity container*

Work within a Singularity container to ensure reproducibility


``` bash
# Initiate virtual box Vagrant
cd $HOME/Manuscript_1_vagrant

vagrant up --provision
# vagrant reload
vagrant ssh

# Check the number of allocated cores and memory
free -h
nproc

##
#rsync -avzPhu iris-cluster:/home/users/dgerard/Singularity_containers/Manuscript_1_singularity.sif ~/Downloads/

#scp -r -P 2222 -i /Users/deborah.gerard/singularity-vm/.vagrant/machines/default/virtualbox/private_key \
#~/Downloads/Manuscript_1_singularity.sif vagrant@127.0.0.1:/home/vagrant/Manuscript_1_singularity/

##
# Check singularity version
singularity version   # 3.9.0

# Prepare the singularity definition file
mkdir Manuscript_1_singularity && cd Manuscript_1_singularity
touch Manuscript_1_singularity.def

# Build the temporary container based on Ubuntu 20.04 from Docker for development 
sudo singularity build --sandbox Manuscript_1_singularity_tmp Manuscript_1_singularity.def

# Run the container in writable mode to make changes
sudo singularity shell --writable Manuscript_1_singularity_tmp

# TEST
sudo singularity build --sandbox Manuscript_1_singularity_tmp Manuscript_1_singularity.sif

# Now that all necessay packages and library have been installed, build the container
sudo singularity build Manuscript_1_singularity.sif Manuscript_1_singularity_tmp

# Start R
singularity exec Manuscript_1_singularity.sif R

# Check R version
R.Version()
```

**Note**: R version is 4.2.3 (2023-03-15)

##### *1.1 Customize the Singularity container*

Install different R libraries


``` r
# tidyverse, ggpubr, rstatix
install.packages(c("tidyverse", "ggpubr", "rstatix"))
# plotgardener
BiocManager::install("plotgardener")

# BSgenome.Hsapiens.UCSC.hg38
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# extrafont
install.packages("extrafont")

# showtext
install.packages("showtext")

# TxDb.Hsapiens.UCSC.hg38.knownGene
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# org.Hs.eg.db
BiocManager::install("org.Hs.eg.db")

# JASPAR2020
BiocManager::install("JASPAR2020")

# TFBSTools
BiocManager::install("TFBSTools")

# ggseqlogo
install.packages("ggseqlogo")

# TxDb.Hsapiens.UCSC.hg38.refGene
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")

# LDlinkR
install.packages("LDlinkR")

# AllelicImbalance
BiocManager::install("AllelicImbalance")

# ggflowchart
install.packages("ggflowchart")

# Issue between dplyr, BiocFilrCache and biomaRt - > downgrade dbplyr
install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

# ggmanh
BiocManager::install("ggmanh")
```
