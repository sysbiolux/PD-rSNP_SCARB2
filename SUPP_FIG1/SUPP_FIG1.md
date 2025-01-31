---
title: "SuppFig1"
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

Plot genes of interest from the TH REP2 mCHERRY cell line RNAseq data to demonstrate that expression is similar.  
In addition, do the same for FOUNDIN-PD.

``` r
############################################ SUPPLEMENTARY FIGURE 1 ##########
############################################ Save as TIFF, 300 dpi
tiff("/home/vagrant/Manuscript_1/SUPPLEMENTARY_FIGURE1/SUPPLEMENTARY_FIGURE1.tiff",
    width = 8.27, height = 9, units = "in", res = 300, compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, height = 9, default.units = "inches", showGuides = TRUE)

#### PANEL A - text Supplementary Figure 1
plotText(label = "Supplementary Figure 1", fontsize = 14, x = 0.25, y = 0.25, just = "left",
    default.units = "inches", fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# Load data
dat.THREP2.FPKM = read_delim("/home/vagrant/epifunc/RNAseq/TH_REP2_mCHERRY/180604_HT-Rep2_fpkm.tsv",
    delim = "\t", col_names = TRUE)

# Filter for genes of interest - BAG3 and LHX1 - SCARB2, NR2C2, FAM47E
BAG3.LHX1 = c("BAG3", "LHX1")

SCARB2.NR2C2 = c("SCARB2", "NR2C2", "FAM47E")

# Reformat the matrix and take only positively sorted neurons
dat.THREP2.FPKM = dat.THREP2.FPKM %>%
    pivot_longer(starts_with("TH"), names_to = "Samples", values_to = "FPKM") %>%
    # dplyr::filter(str_detect(Samples, paste('SmNPC', 'pos', sep = '|'))) %>%
mutate(Cond = case_when(str_detect(Samples, "D15pos") ~ "Positively sorted neurons D15",
    str_detect(Samples, "D15neg") ~ "Negatively sorted neurons D15", str_detect(Samples,
        "D30pos") ~ "Positively sorted neurons D30", str_detect(Samples, "D30neg") ~
        "Negatively sorted neurons D30", TRUE ~ "smNPCs"), Cond = factor(Cond, levels = c("smNPCs",
    "Negatively sorted neurons D15", "Positively sorted neurons D15", "Negatively sorted neurons D30",
    "Positively sorted neurons D30")))

# Plot the expression of BAG3 and LHX1 in smNPC and mDAN
cdt = c("smNPCs", "Positively")

TH.REP2_BAG3.LHX1.exp = ggboxplot(dat.THREP2.FPKM %>%
    dplyr::filter(gene_name %in% BAG3.LHX1, str_detect(Cond, paste0(cdt, collapse = "|"))),
    x = "Cond", y = "FPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C"), facet.by = "gene_name", xlab = "", ylab = "Fragments Per Kilobase Million (FPKM)") +
    scale_x_discrete(labels = c("smNPC", expression("mDAN D15"), expression("mDAN D30"))) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = TH.REP2_BAG3.LHX1.exp, x = 0.5, y = 0.625, width = 3.5, height = 3.5,
    just = c("left", "top"), default.units = "inches")

# Plot the expression of SCARB2 and NR2C2
TH.REP2_SCARB2.NR2C2.exp = ggboxplot(dat.THREP2.FPKM %>%
    dplyr::filter(gene_name %in% SCARB2.NR2C2, str_detect(Cond, paste0(cdt, collapse = "|"))) %>%
    mutate(gene_name = factor(gene_name, levels = c("SCARB2", "NR2C2", "FAM47E"))),
    x = "Cond", y = "FPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C"), facet.by = "gene_name", xlab = "", ylab = "Fragments Per Kilobase Million (FPKM)") +
    scale_x_discrete(labels = c("smNPC", expression("mDAN D15"), expression("mDAN D30"))) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = TH.REP2_SCARB2.NR2C2.exp, x = 4.5, y = 0.625, width = 3.5, height = 3.5,
    just = c("left", "top"), default.units = "inches")

#### PANEL B - text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 4.5, just = "left", default.units = "inches",
    fontface = "bold")

# FOUNDIN-PD data Load
FIPD = read_delim("/home/vagrant/epifunc/PPMI_data/FOUNDIN-PD_150.9_RNAB/aggregated_expression/cpmTable.tsv",
    delim = "\t", col_names = TRUE)

# Filter for genes of interest and get the ensembl and gene names
goi_ens = AnnotationDbi::select(org.Hs.eg.db, keys = c("BAG3", "LHX1", "SCARB2",
    "NR2C2", "FAM47E"), columns = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")

# 2 LHX1 ensembl IDs came back. Check on the Ensembl website which one is
# correct
goi_ens = goi_ens %>%
    as_tibble() %>%
    dplyr::filter(!ENSEMBL %in% c("ENSG00000274577", "ENSG00000262446")) %>%
    dplyr::rename(ENSEMBL.not.dot = ENSEMBL)

# Filter the matrix of expression for my genes of interest
FIPD.fin = FIPD %>%
    dplyr::filter(str_detect(Geneid, paste(goi_ens$ENSEMBL.not.dot, collapse = "|"))) %>%
    dplyr::rename(ENSEMBL = Geneid) %>%
    mutate(ENSEMBL.not.dot = gsub("\\..*", "", ENSEMBL)) %>%
    dplyr::select(ENSEMBL, ENSEMBL.not.dot, everything()) %>%
    pivot_longer(cols = starts_with("RNAB"), names_to = "Samples", values_to = "CPM") %>%
    mutate(Day = str_extract(Samples, "da+\\d+")) %>%
    left_join(goi_ens, by = "ENSEMBL.not.dot") %>%
    dplyr::select(SYMBOL, Day, CPM)

# Plot expression of BAG3 and LHX1
FIPD.exp_BAG3.LHX1 = ggboxplot(FIPD.fin %>%
    dplyr::filter(SYMBOL %in% BAG3.LHX1), x = "Day", y = "CPM", add = "jitter", fill = "Day",
    palette = c("#482677FF", "#2D708EFF", "#29AF7FFF"), facet.by = "SYMBOL", xlab = "",
    ylab = "Counts per Million (CPM)", add.params = list(size = 0.1)) + scale_x_discrete(labels = c("Day 0",
    "Day 25", "Day 65")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = FIPD.exp_BAG3.LHX1, x = 0.5, y = 4.75, width = 3.5, height = 3.5, just = c("left",
    "top"), default.units = "inches")

# Plot expression of SCARB2 alone
FIPD.exp_SCARB2 = ggboxplot(FIPD.fin %>%
    dplyr::filter(SYMBOL %in% "SCARB2"), x = "Day", y = "CPM", add = "jitter", fill = "Day",
    palette = c("#482677FF", "#2D708EFF", "#29AF7FFF"), facet.by = "SYMBOL", xlab = "",
    ylab = "Counts per Million (CPM)", add.params = list(size = 0.1)) + scale_x_discrete(labels = c("Day 0",
    "Day 25", "Day 65")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = FIPD.exp_SCARB2, x = 4.25, y = 4.75, width = 1.5, height = 3.5, just = c("left",
    "top"), default.units = "inches")

# Plot expression of NR2C2 and FAM47E
FIPD.exp_NR2C2.FAM47E = ggboxplot(FIPD.fin %>%
    dplyr::filter(SYMBOL %in% c("NR2C2", "FAM47E")) %>%
    mutate(SYMBOL = factor(SYMBOL, levels = c("NR2C2", "FAM47E"))), x = "Day", y = "CPM",
    add = "jitter", fill = "Day", palette = c("#482677FF", "#2D708EFF", "#29AF7FFF"),
    facet.by = "SYMBOL", xlab = "", ylab = "Counts per Million (CPM)", add.params = list(size = 0.1)) +
    scale_x_discrete(labels = c("Day 0", "Day 25", "Day 65")) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = FIPD.exp_NR2C2.FAM47E, x = 6, y = 4.75, width = 2, height = 3.5, just = c("left",
    "top"), default.units = "inches")

pageGuideHide()
dev.off()
```
