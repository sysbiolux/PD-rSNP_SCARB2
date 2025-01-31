---
title: "Figure4"
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


``` r
########## FIGURE 4 ########## Save as TIFF, 300 dpi
tiff("/home/vagrant/Manuscript_1/FIGURE4/FIGURE4.tiff", width = 8.27, height = 11.67,
    units = "in", res = 300, compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, height = 11.67, default.units = "inches", showGuides = TRUE)

plotText(label = "Figure 4", fontsize = 14, fontface = "bold", fontfamily = "Helvetica",
    x = 0.25, y = 0.25, just = "left", default.units = "inches")

#### PANEL A - Chromatin imbalance in TH REP1 mCHERRY cell line at rs1465922
#### text A
plotText(label = "A", fontsize = 12, fontface = "bold", fontfamily = "Helvetica",
    x = 0.25, y = 0.5, just = "left", default.units = "inches")


# Chromatin imbalance values as calculated by the AllelicImbalance R package
countList = readRDS("/home/vagrant/Manuscript_1/FIGURE4/SCARB2_het_counts_THREP1_ATACseq.rds")

# Extract and reformat for positively sorted TH population
AI_ATAC = countList$chr4_76213717 %>%
    as_tibble() %>%
    mutate(Sample = rownames(countList$chr4_76213717)) %>%
    pivot_longer(cols = c(A, C, G, T), names_to = "Base", values_to = "Reads") %>%
    mutate(Sample.n = case_when(str_detect(Sample, "astrocytes") ~ "Astrocytes",
        str_detect(Sample, "D15_neg") ~ "Negative DAN D15", str_detect(Sample, "D15_pos") ~
            "mDAN D15", str_detect(Sample, "D30_pos") ~ "mDAN D30", str_detect(Sample,
            "D50_neg") ~ "Negative DAN D50", str_detect(Sample, "D50_pos") ~ "mDAN D50",
        TRUE ~ "smNPCs")) %>%
    dplyr::filter(!str_detect(Sample, "neg"), Base %in% c("A", "G")) %>%
    mutate(Sample.n = factor(Sample.n, levels = c("smNPCs", "mDAN D15", "mDAN D30",
        "mDAN D50", "Astrocytes")), Allele = case_when(Base == "A" ~ "PD allele",
        TRUE ~ "Other allele"))
# Plot
chrom_AI = ggbarplot(AI_ATAC, x = "Sample.n", y = "Reads", fill = "Allele", position = position_dodge(0.7),
    xlab = "", ylab = "Reads", label = format(c(paste0("p = ", 0.3), "", paste0("p = ",
        0.3), "", paste0("p = ", 0.003), "", paste0("p = ", 3e-07), "", paste0("p = ",
        0.2), ""), scientific = TRUE), label.pos = "out", lab.size = 3, lab.vjust = 0.5,
    lab.hjust = -0.5, orientation = "horiz") + theme(legend.position = "right", axis.text.x = element_text(angle = 45,
    hjust = 1, face = "bold", size = 10), axis.text.y = element_text(face = "bold",
    size = 10), axis.title.y = element_text(face = "bold", size = 10), axis.title.x = element_text(face = "bold",
    size = 10), legend.text = element_text(face = "bold", size = 10), legend.title = element_text(face = "bold",
    size = 10), plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = c("#80A44F",
    "#427F96")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 330))

# Place the plot
plotGG(plot = chrom_AI, x = 0.5, y = 0.75, width = 7, height = 4, just = c("left",
    "top"), default.units = "inches")

# text B
plotText(label = "B", fontsize = 12, fontfamily = "Helvetica", x = 0.25, y = 5, just = "left",
    default.units = "inches", fontface = "bold")

# Load GTEX data (v8)
GTEX_SNPs_interest = read_delim("/home/vagrant/epifunc/eQTpLot/20230117_GTEX_data.csv",
    delim = ";", col_names = TRUE)

# Take the SNPs of interest
GTEX_SNPs_SCARB2_brain = GTEX_SNPs_interest %>%
    mutate(GTEX_sig = if_else(P.Value < 0.05, "*", "")) %>%
    dplyr::filter(str_detect(Tissue, "Brain"), Gene.Symbol == "SCARB2") %>%
    mutate(SNP.Id = factor(SNP.Id, levels = c("rs6825004", "rs7697073", "rs11547135",
        "rs6812193", "rs1465922")))

# Plot
GTEX_eqtl = ggplot(data = GTEX_SNPs_SCARB2_brain, mapping = aes(x = SNP.Id, y = Tissue,
    fill = NES)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 60,
    hjust = 1, face = "bold", size = 10), axis.title.x = element_blank(), axis.text.y = element_text(face = "bold",
    size = 10), axis.title.y = element_blank(), legend.text = element_text(face = "bold",
    size = 10), legend.title = element_text(face = "bold", size = 10)) + scale_fill_gradient2(low = muted("blue"),
    mid = "white", high = muted("red"), midpoint = 0, name = "eQTL \neffect size") +
    geom_text(aes(label = GTEX_sig))

# Palce the plot
plotGG(plot = GTEX_eqtl, x = 0.5, y = 5.25, width = 7, height = 6, just = c("left",
    "top"), default.units = "inches")


pageGuideHide()
dev.off()
```
