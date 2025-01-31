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

### FIGURE 3

``` r
########## FIGURE 3 ##########
tiff("/home/vagrant/Manuscript_1/FIGURE3/FIGURE3.tiff", width = 8.27, height = 11.67,
    units = "in", res = 300, compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, height = 11.67, default.units = "inches", showGuides = TRUE)

# text Figure 3
plotText(label = "Figure 3", fontsize = 14, fontface = "bold", fontfamily = "Helvetica",
    x = 0.25, y = 0.25, just = "left", default.units = "inches")

#### PANEL A - Scheme of transduction for knockdown and luciferase text A
plotText(label = "A", fontsize = 12, fontface = "bold", fontfamily = "Helvetica",
    x = 0.25, y = 0.5, just = "left", default.units = "inches")

# Fig3A
fig3A = readPNG("/home/vagrant/Manuscript_1/FIGURE3/Fig3A_BIORENDER.png")

plotRaster(image = fig3A, x = 0.25, y = 1, width = 8, height = 1, just = "left",
    interpolate = TRUE)

#### PANEL B - luciferase barplots with statistics text B
plotText(label = "B", fontsize = 12, fontface = "bold", fontfamily = "Helvetica",
    x = 0.25, y = 2, just = "left", default.units = "inches")

Gluc_signal = tibble(Sample = c("BAG3-WT", "BAG3-MUT", "IDUA-WT", "IDUA-MUT", "SCARB2-WT",
    "SCARB2-MUT", "miniCMV", "Neg_CTRL"), GLUC_n1_rep1 = c(319, 240, 154, 196, 379,
    211, 431, 33), GLUC_n1_rep2 = c(323, 237, 161, 205, 371, 224, 459, 33), GLUC_n2_rep1 = c(201,
    143, 121, 128, 278, 172, 253, 31), GLUC_n2_rep2 = c(191, 149, 119, 126, 278,
    162, 249, 31), GLUC_n3_rep1 = c(249, 146, 136, 142, 223, 175, 212, 32), GLUC_n3_rep2 = c(203,
    152, 121, 125, 209, 168, 190, 33), GLUC_n4_rep1 = c(348, 333, 220, 185, 396,
    244, 319, 31), GLUC_n4_rep2 = c(322, 306, 229, 186, 369, 227, 294, 32), SEAP_n1_rep1 = c(72,
    78, 54, 67, 99, 75, 128, 32), SEAP_n1_rep2 = c(69, 70, 49, 66, 97, 79, 132, 29),
    SEAP_n2_rep1 = c(50, 52, 45, 45, 76, 64, 80, 32), SEAP_n2_rep2 = c(50, 52, 45,
        51, 85, 63, 80, 32), SEAP_n3_rep1 = c(52, 50, 43, 50, 66, 55, 68, 32), SEAP_n3_rep2 = c(49,
        53, 48, 46, 66, 62, 69, 31), SEAP_n4_rep1 = c(64, 75, 58, 54, 103, 75, 111,
        33), SEAP_n4_rep2 = c(61, 74, 62, 56, 108, 77, 111, 31)) %>%
    pivot_longer(!Sample, names_to = "Meas", values_to = "Abs_val") %>%
    mutate(Meas2 = gsub("^([^_]*_[^_]*)_.*", "\\1", Meas)) %>%
    group_by(Sample, Meas2) %>%
    mutate(AVG = mean(Abs_val, na.rm = TRUE))

# Calculate ratio Gluc over SEAP
Gluc_signal_ratio = Gluc_signal %>%
    dplyr::select(-Meas, -Abs_val) %>%
    distinct(Meas2, .keep_all = TRUE) %>%
    pivot_wider(names_from = "Meas2", values_from = "AVG") %>%
    mutate(Ratio_within_bio_rep1 = GLUC_n1/SEAP_n1, Ratio_within_bio_rep2 = GLUC_n2/SEAP_n2,
        Ratio_within_bio_rep3 = GLUC_n3/SEAP_n3, Ratio_within_bio_rep4 = GLUC_n4/SEAP_n4) %>%
    dplyr::select(Sample, Ratio_within_bio_rep1:Ratio_within_bio_rep4) %>%
    pivot_longer(cols = c("Ratio_within_bio_rep1", "Ratio_within_bio_rep2", "Ratio_within_bio_rep3",
        "Ratio_within_bio_rep4"), names_to = "Ratio", values_to = "value") %>%
    mutate(Genes = gsub("(.*)-.*", "\\1", Sample))


# Calculate statistics - T-test
GLUC_stat = compare_means(value ~ Sample, data = Gluc_signal_ratio %>%
    ungroup(), method = "t.test")

# Plot for BAG3
GLUC_BAG3 = ggbarplot(Gluc_signal_ratio %>%
    dplyr::filter(Sample %in% c("BAG3-WT", "BAG3-MUT", "miniCMV", "Neg_CTRL")), x = "Sample",
    y = "value", width = 0.25, add = c("mean_se", "jitter"), fill = "grey", xlab = "",
    ylab = "Ratio Gluc/SEAP", position = position_dodge(0.9)) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9, face = "bold", family = "Helvetica"),
    axis.text.y = element_text(size = 9, face = "bold"), axis.title.y = element_text(size = 10,
        face = "bold")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
    scale_x_discrete(labels = c("BAG3 other allele", "BAG3 effect allele", "Positive control",
        "Negative control")) + stat_pvalue_manual(GLUC_stat %>%
    dplyr::filter(group1 == "BAG3-WT", group2 == "BAG3-MUT"), label = "p = {p.format}",
    y.position = 5.6, size = 3)

# Place the plot
plotGG(plot = GLUC_BAG3, x = 0.5, y = 2.25, width = 3.5, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Plot for SCARB2
GLUC_SCARB2 = ggbarplot(Gluc_signal_ratio %>%
    dplyr::filter(Sample %in% c("SCARB2-WT", "SCARB2-MUT", "miniCMV", "Neg_CTRL")),
    x = "Sample", y = "value", width = 0.25, add = c("mean_se", "jitter"), fill = "grey",
    xlab = "", ylab = "Ratio Gluc/SEAP", position = position_dodge(0.9)) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9, face = "bold"), axis.text.y = element_text(size = 9,
        face = "bold"), axis.title.y = element_text(size = 10, face = "bold")) +
    scale_x_discrete(labels = c("SCARB2 other allele", "SCARB2 effect allele", "Positive control",
        "Negative control")) + scale_y_continuous(expand = c(0, 0), limits = c(0,
    6)) + stat_pvalue_manual(GLUC_stat %>%
    dplyr::filter(group1 == "SCARB2-WT", group2 == "SCARB2-MUT"), label = "p = {p.format}",
    y.position = 4, size = 3)

# Place the plot
plotGG(plot = GLUC_SCARB2, x = 4.5, y = 2.25, width = 3.5, height = 4.25, just = c("left",
    "top"), default.units = "inches")

#### PANEL C - knockdown barplots with statistics text C
plotText(label = "C", fontsize = 12, fontfamily = "Helvetica", x = 0.25, y = 6.75,
    just = "left", default.units = "inches", fontface = "bold")

# Knockdown of LHX1 and quantification of LHX1 and BAG3 expression 3 days post
# transduction
KD.BAG3.LHX1.3days = tibble(Sple = rep("shLHX1", 6), Gene = c(rep("LHX1", 3), rep("BAG3",
    3)), Exp = c(0.411253539, 0.359459217, 0.367190337, 0.972688658, 0.965902853,
    0.88813445))

# Stats
KD.BAG3.LHX1.3days.stats = tibble(Sple = c(rep("shLHX1", 6), rep("shCTRL", 6)), Gene = rep(c(rep("LHX1",
    3), rep("BAG3", 3)), 2), Exp = c(9.5167, 8.39675, 8.1855, 9.13605, 8.33795, 8.76415,
    8.2348, 6.92065, 6.7401, 9.0961, 8.2879, 8.593))

stat_LHX1.KD = compare_means(Exp ~ Sple, data = KD.BAG3.LHX1.3days.stats, group.by = "Gene",
    method = "t.test", paired = TRUE)
# Plot
pl.KD_BAG3_LHX1 = KD.BAG3.LHX1.3days %>%
    mutate(Gene = factor(Gene, levels = c("LHX1", "BAG3"))) %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("Gene"), position = position_dodge(0.9),
        palette = c("#440154FF", "#287C8EFF"), add = c("mean_se", "jitter"), xlab = "",
        ylab = "Relative expression to shSCRAMBLE", title = "shLHX1 - 3 days") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold.italic", size = 10, angle = 90,
    family = "Helvetica"), axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none", plot.title = element_text(face = "bold", size = 10,
        hjust = 0.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
    stat_pvalue_manual(stat_LHX1.KD, label = "p = {p.format}", x = "Gene", y.position = c(0.5,
        1.2), size = 3, hide.ns = FALSE)

plotGG(plot = pl.KD_BAG3_LHX1, x = 0.5, y = 7, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Knockdown of NR2C2 and quantification of NR2C2 and SCARB2 expression 3 days
# post transduction
KD.SCARB2.NR2C2.3days = tibble(Sple = rep("shNR2C2", 6), Gene = c(rep("NR2C2", 3),
    rep("SCARB2", 3)), Exp = c(0.345210918, 0.512864636, 0.548665969, 1.239922544,
    1.130295504, 1.170885498))

# Stats
KD.SCARB2.NR2C2.3days.stats = tibble(Sple = c(rep("shNR2C2", 6), rep("shCTRL", 6)),
    Gene = rep(c(rep("NR2C2", 3), rep("SCARB2", 3)), 2), Exp = c(9.43575, 9.3043,
        9.263, 5.6841, 5.92295, 5.9289, 7.9013, 8.34095, 8.397, 5.99435, 6.09965,
        6.1565))

stat_NR2C2.KD = compare_means(Exp ~ Sple, data = KD.SCARB2.NR2C2.3days.stats, group.by = "Gene",
    method = "t.test", paired = TRUE)

# Plot
pl.KD_SCARB2_NR2C2 = KD.SCARB2.NR2C2.3days %>%
    mutate(Gene = factor(Gene, levels = c("NR2C2", "SCARB2"))) %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("Gene"), position = position_dodge(0.9),
        palette = c("#440154FF", "#287C8EFF"), add = c("mean_se", "jitter"), xlab = "",
        ylab = "Relative expression to shSCRAMBLE", title = "shNR2C2 - 3 days") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold.italic", size = 10, angle = 90,
    family = "Helvetica"), axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none", plot.title = element_text(face = "bold", size = 10,
        hjust = 0.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
    stat_pvalue_manual(stat_NR2C2.KD, label = "p = {p.format}", x = "Gene", y.position = c(0.65,
        1.5), size = 3, hide.ns = FALSE)

plotGG(plot = pl.KD_SCARB2_NR2C2, x = 3, y = 7, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Knockdown of NR2C2 and quantification of NR2C2 and SCARB2 expression 6 days
# post transduction
KD.SCARB2.NR2C2.6days = tibble(Sple = rep("shNR2C2", 6), Gene = c(rep("NR2C2", 3),
    rep("SCARB2", 3)), Exp = c(0.662710186, 0.456029798, 0.787662775, 1.410640142,
    1.425976296, 1.332143962))

# Stats
KD.SCARB2.NR2C2.6days.stats = tibble(Sple = c(rep("shNR2C2", 6), rep("shCTRL", 6)),
    Gene = rep(c(rep("NR2C2", 3), rep("SCARB2", 3)), 2), Exp = c(8.34355, 9.0917,
        8.32445, 5.3595, 5.19535, 6.5266, 7.75, 7.9589, 7.9801, 5.85585, 5.7073,
        6.94035))

stat_NR2C2.KD.6days = compare_means(Exp ~ Sple, data = KD.SCARB2.NR2C2.6days.stats,
    group.by = "Gene", method = "t.test", paired = TRUE)

# Plot
pl.KD_SCARB2_NR2C2_6days = KD.SCARB2.NR2C2.6days %>%
    mutate(Gene = factor(Gene, levels = c("NR2C2", "SCARB2"))) %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("Gene"), position = position_dodge(0.9),
        palette = c("#440154FF", "#287C8EFF"), add = c("mean_se", "jitter"), xlab = "",
        ylab = "Relative expression to shSCRAMBLE", title = "shNR2C2 - 6 days") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold.italic", size = 10, angle = 90,
    family = "Helvetica"), axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none", plot.title = element_text(face = "bold", size = 10,
        hjust = 0.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
    stat_pvalue_manual(stat_NR2C2.KD.6days, label = "p = {p.format}", x = "Gene",
        y.position = c(0.85, 1.5), size = 3, hide.ns = FALSE)

plotGG(plot = pl.KD_SCARB2_NR2C2_6days, x = 5.5, y = 7, width = 2, height = 4.25,
    just = c("left", "top"), default.units = "inches")

pageGuideHide()
dev.off()
```
