---
title: "SUPP_FIG2"
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

IHEC consortium has one sample of substantia nigra for H3K27me3, H3K9me3, H3K4me1, H3K4me3 and H3K36me3 from a healthy brain for BAG3, IDUA and SCARB2 loci

``` r
############################################ SUPPLEMENTARY FIGURE 2 ##########
############################################ Save as TIFF, 300 dpi
tiff("/home/vagrant/Manuscript_1/SUPPLEMENTARY_FIGURE2/SUPPLEMENTARY_FIGURE2.tiff",
    width = 8.27, height = 6.5, units = "in", res = 300, compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, height = 6.5, default.units = "inches", showGuides = TRUE)

#### PANEL A - text Supplementary Figure 2
plotText(label = "Supplementary Figure 2", fontsize = 14, x = 0.25, y = 0.25, just = "left",
    default.units = "inches", fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# ChIP-seq signal BAG3 path where bigwig files are
bw.IHEC.path = "/home/vagrant/Documents/IHEC/ChIP_Seq/BIGWIG/"

# Reads BIGWIG file at BAG3 location - H3K4me1
bw.BAG3.H3K4me1 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6a1889bf-7312-4e41-b1cc-e2f3be7a8076.fc.signal.bigwig"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Reads BIGWIG file at BAG3 location - H3K4me3
bw.BAG3.H3K4me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.4cb9300f-9097-4374-8f2a-4afd8ff855fb.fc.signal.bigwig"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Reads BIGWIG file at BAG3 location - H3K36me3
bw.BAG3.H3K36me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6586e624-852f-4cac-83df-2325cee19699.fc.signal.bigwig"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Reads BIGWIG file at BAG3 location - H3K9me3
bw.BAG3.H3K9me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.e8fe7664-906e-4d87-8f3a-6ce88dbeadc6.fc.signal.bigwig"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Reads BIGWIG file at BAG3 location - H3K27me3
bw.BAG3.H3K27me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.5f5024a9-2e3f-4527-82ec-0ed04bf4a848.fc.signal.bigwig"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Add a scale next to the bigwig files - check the maximum to choose
scale.IHEC.BAG3.max = max(c(bw.BAG3.H3K4me1$score, bw.BAG3.H3K4me3$score, bw.BAG3.H3K36me3$score,
    bw.BAG3.H3K9me3$score, bw.BAG3.H3K27me3$score)) %>%
    round(., digits = 1)

print(scale.IHEC.BAG3.max)

# Define parameters for the regions
region.p.IHEC.BAG3 = pgParams(chrom = "chr10", chromstart = 119611405, chromend = 119691405,
    assembly = "hg38", range = c(0, scale.IHEC.BAG3.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr10", chromstart = 119611405, chromend = 119691405, assembly = "hg38",
    x = 0.75, y = 0.625, length = 4, default.units = "inches", scale = "Mb")

# ChIPseq signal H3K4me1
H3K4me1_BAG3 = plotSignal(data = bw.BAG3.H3K4me1, params = region.p.IHEC.BAG3, fill = "#b5de2b",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 0.875, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K4me1_BAG3, at = c(0, scale.IHEC.BAG3.max), fontsize = 6)

# ChIPseq signal H3K4me3
H3K4me3_BAG3 = plotSignal(data = bw.BAG3.H3K4me3, params = region.p.IHEC.BAG3, fill = "#35b779",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 1.25, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K4me3_BAG3, at = c(0, scale.IHEC.BAG3.max), fontsize = 6)

# ChIPseq signal H3K36me3
H3K36me3_BAG3 = plotSignal(data = bw.BAG3.H3K36me3, params = region.p.IHEC.BAG3,
    fill = "#26828e", alpha = 0.7, linecolor = NA, x = 0.75, y = 1.625, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K36me3_BAG3, at = c(0, scale.IHEC.BAG3.max), fontsize = 6)

# ChIPseq signal H3K9me3
H3K9me3_BAG3 = plotSignal(data = bw.BAG3.H3K9me3, params = region.p.IHEC.BAG3, fill = "#3e4989",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 2, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K9me3_BAG3, at = c(0, scale.IHEC.BAG3.max), fontsize = 6)

# ChIPseq signal H3K27me3
H3K27me3_BAG3 = plotSignal(data = bw.BAG3.H3K27me3, params = region.p.IHEC.BAG3,
    fill = "#440154", alpha = 0.7, linecolor = NA, x = 0.75, y = 2.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K27me3_BAG3, at = c(0, scale.IHEC.BAG3.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr10", chromstart = 119611405, chromend = 119691405, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 2.75, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me1", fontsize = 6, fontcolor = "#b5de2b", x = 0.75, y = 0.875,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me3", fontsize = 6, fontcolor = "#35b779", x = 0.75, y = 1.25,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K36me3", fontsize = 6, fontcolor = "#26828e", x = 0.75, y = 1.625,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K9me3", fontsize = 6, fontcolor = "#3e4989", x = 0.75, y = 2,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K27me3", fontsize = 6, fontcolor = "#440154", x = 0.75, y = 2.375,
    just = c("left", "top"), default.units = "inches")

#### PANEL B - text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 3.25, just = "left", default.units = "inches",
    fontface = "bold")

# Reads BIGWIG file at SCARB2 location - H3K4me1
bw.SCARB2.H3K4me1 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6a1889bf-7312-4e41-b1cc-e2f3be7a8076.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Reads BIGWIG file at SCARB2 location - H3K4me3
bw.SCARB2.H3K4me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.4cb9300f-9097-4374-8f2a-4afd8ff855fb.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Reads BIGWIG file at SCARB2 location - H3K36me3
bw.SCARB2.H3K36me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6586e624-852f-4cac-83df-2325cee19699.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Reads BIGWIG file at SCARB2 location - H3K9me3
bw.SCARB2.H3K9me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.e8fe7664-906e-4d87-8f3a-6ce88dbeadc6.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Reads BIGWIG file at SCARB2 location - H3K27me3
bw.SCARB2.H3K27me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.5f5024a9-2e3f-4527-82ec-0ed04bf4a848.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Add a scale next to the bigwig files - check the maximum to choose
scale.IHEC.SCARB2.max = max(c(bw.SCARB2.H3K4me1$score, bw.SCARB2.H3K4me3$score, bw.SCARB2.H3K36me3$score,
    bw.SCARB2.H3K9me3$score, bw.SCARB2.H3K27me3$score)) %>%
    round(., digits = 1)

print(scale.IHEC.SCARB2.max)

# Define parameters for the regions
region.p.IHEC.SCARB2 = pgParams(chrom = "chr4", chromstart = 76173717, chromend = 76253717,
    assembly = "hg38", range = c(0, scale.IHEC.SCARB2.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr4", chromstart = 76173717, chromend = 76253717, assembly = "hg38",
    x = 0.75, y = 3.5, length = 4, default.units = "inches", scale = "Mb")

# ChIPseq signal H3K4me1
H3K4me1_SCARB2 = plotSignal(data = bw.SCARB2.H3K4me1, params = region.p.IHEC.SCARB2,
    fill = "#b5de2b", alpha = 0.7, linecolor = NA, x = 0.75, y = 3.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K4me1_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K4me3
H3K4me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K4me3, params = region.p.IHEC.SCARB2,
    fill = "#35b779", alpha = 0.7, linecolor = NA, x = 0.75, y = 4.125, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K4me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K36me3
H3K36me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K36me3, params = region.p.IHEC.SCARB2,
    fill = "#26828e", alpha = 0.7, linecolor = NA, x = 0.75, y = 4.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K36me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K9me3
H3K9me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K9me3, params = region.p.IHEC.SCARB2,
    fill = "#3e4989", alpha = 0.7, linecolor = NA, x = 0.75, y = 4.875, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K9me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K27me3
H3K27me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K27me3, params = region.p.IHEC.SCARB2,
    fill = "#440154", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.25, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K27me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr4", chromstart = 76173717, chromend = 76253717, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 5.625, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me1", fontsize = 6, fontcolor = "#b5de2b", x = 0.75, y = 3.75,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me3", fontsize = 6, fontcolor = "#35b779", x = 0.75, y = 4.125,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K36me3", fontsize = 6, fontcolor = "#26828e", x = 0.75, y = 4.5,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K9me3", fontsize = 6, fontcolor = "#3e4989", x = 0.75, y = 4.875,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K27me3", fontsize = 6, fontcolor = "#440154", x = 0.75, y = 5.25,
    just = c("left", "top"), default.units = "inches")

pageGuideHide()
dev.off()
```
