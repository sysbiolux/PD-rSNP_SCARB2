---
title: "Figure3"
author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "31 January, 2025"
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



Start R


``` bash
# Start R
singularity exec Manuscript_1_singularity.sif R

# Check R version
R.Version()
```

**Note**: R version is 4.2.3 (2023-03-15)

Load libraries and check their versions


``` r
# Load
library("tidyverse")
library("rstatix")
library("ggpubr")
library("plotgardener")
# library('BSgenome.Hsapiens.UCSC.hg38')
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.refGene")
library("org.Hs.eg.db")
library("extrafont")
library("showtext")
library("JASPAR2020")
library("TFBSTools")
library("ggseqlogo")
library("biomaRt")
library("LDlinkR")
library("AllelicImbalance")
library("png")
library("AnnotationDbi")
library("scales")
library("flowCore")
library("ggcyto")
library("ggmanh")
library("rtracklayer")
library("plyranges")
library("RColorBrewer")

# Check
packageVersion("tidyverse")  #2.0.0
packageVersion("rstatix")  #0.7.2
packageVersion("ggpubr")  #0.6.0
packageVersion("plotgardener")  #1.4.2
# packageVersion('BSgenome.Hsapiens.UCSC.hg38') #1.4.5
# packageVersion('TxDb.Hsapiens.UCSC.hg38.knownGene') #3.16.0
packageVersion("TxDb.Hsapiens.UCSC.hg38.refGene")  #3.15.0
packageVersion("org.Hs.eg.db")  #3.16.0
packageVersion("extrafont")  #0.19
packageVersion("showtext")  #0.9.7
packageVersion("JASPAR2020")  #0.99.10
packageVersion("TFBSTools")  #1.36.0
packageVersion("ggseqlogo")  #0.2
packageVersion("biomaRt")  #2.54.1
packageVersion("LDlinkR")  #1.4.0
packageVersion("AllelicImbalance")  #1.36.0
packageVersion("png")  #0.1.8
packageVersion("AnnotationDbi")  #1.60.2
# packageVersion('dbplyr') #2.5.0
packageVersion("scales")  #1.3.0
packageVersion("flowCore")  #2.10.0
packageVersion("ggcyto")  #1.26.4
packageVersion("ggmanh")  #1.2.0
packageVersion("rtracklayer")  #1.58.0
packageVersion("plyranges")  #1.18.0
packageVersion("RColorBrewer")  #1.1.3

# Import all available fonts
font_import()
```
