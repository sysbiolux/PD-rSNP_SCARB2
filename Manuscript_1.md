---
title: "Manuscript_1"
author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "06 February, 2024"
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


```bash
# Initiate virtual box Vagrant
cd $HOME/Manuscript_1_vagrant

vagrant up --provision
# vagrant reload
vagrant ssh

# Check the number of allocated cores and memory
free -h
nproc

# Check singularity version
singularity version   # 3.9.0

# Prepare the singularity definition file
mkdir Manuscript_1_singularity && cd Manuscript_1_singularity
touch Manuscript_1_singularity.def

# Build the temporary container based on Ubuntu 20.04 from Docker for development 
sudo singularity build --sandbox Manuscript_1_singularity_tmp Manuscript_1_singularity.def

# Run the container in writable mode to make changes
sudo singularity shell --writable Manuscript_1_singularity_tmp

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


```r
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
```

Load libraries and check their versions


```r
# Load
library("tidyverse")
library("rstatix")
library("ggpubr")
library("plotgardener")
# library('BSgenome.Hsapiens.UCSC.hg38')
# library('TxDb.Hsapiens.UCSC.hg38.knownGene')
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
packageVersion("showtext")  #0.9.6
packageVersion("JASPAR2020")  #0.99.10
packageVersion("TFBSTools")  #1.36.0
packageVersion("ggseqlogo")  #0.1
packageVersion("biomaRt")  #2.54.1
packageVersion("LDlinkR")  #1.3.0
packageVersion("AllelicImbalance")  #1.36.0
packageVersion("png")  #0.1.8
packageVersion("AnnotationDbi")  #1.60.2
packageVersion("dbplyr")  #2.3.4
```

# FIGURES
### FIGURE 1

```r
########## FIGURE 1 ########## Save as a PDF
pdf("/home/vagrant/Manuscript_1/FIGURE1/FIGURE1.pdf", width = 8.3, height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = FALSE)

#### PANEL A - text Figure 1
plotText(label = "Figure 1", fontsize = 14, x = 0.25, y = 0.25, just = "left", default.units = "inches",
    fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 1.5, just = "left", default.units = "inches",
    fontface = "bold")

# text C
plotText(label = "C", fontsize = 12, x = 0.25, y = 4.5, just = "left", default.units = "inches",
    fontface = "bold")

# text D
plotText(label = "D", fontsize = 12, x = 0.25, y = 8.25, just = "left", default.units = "inches",
    fontface = "bold")

################### Figure 1A ####
Fig1A = readPNG("/home/vagrant/Manuscript_1/FIGURE1/Fig1A.png")

# Plot Figure 1A
plotRaster(image = Fig1A, x = 0.375, y = 0.625, width = 7.55, height = 1, just = c("left",
    "top"), interpolate = FALSE)

################### Figure 1B #### Load data obtained from Jochen
RNAseq_odd = read_delim("/home/vagrant/Manuscript_1/FIGURE1/REFORMAT_tpm_gene_sets.gsa.txt",
    delim = "\t", col_names = TRUE)

# RNAseq_odd =
# read_delim('/Volumes/deborah.gerard/Documents/Manuscript_1/FIGURE1/REFORMAT_tpm_gene_sets.gsa.txt',
# delim = '\t', col_names = TRUE)

# Calculate odd ratio
RNAseq_odd = RNAseq_odd %>%
    dplyr::filter(!VARIABLE %in% c("microglia")) %>%
    mutate(odd_R = exp(BETA)) %>%
    mutate(VARIABLE = factor(VARIABLE, levels = c("smNPC", "D15negsort", "D15possort",
        "D30postsort", "D50negsort", "D50possort", "astrocytes")))

# Plot
fig1B = ggplot(RNAseq_odd, aes(x = odd_R, y = VARIABLE)) + geom_point(shape = 16,
    size = 3) + geom_errorbarh(aes(xmin = odd_R - SE, xmax = odd_R + SE), height = 0.25) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") + theme_classic() +
    theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold")) +
    scale_y_discrete(labels = c(smNPC = "smNPC", D15negsort = "non-mDAN D15", D15possort = "mDAN D15",
        D30postsort = "mDAN D30", D50negsort = "non-mDAN D50", D50possort = "mDAN D50",
        astrocytes = "Astrocytes")) + geom_text(label = round(RNAseq_odd$P, digits = 3),
    nudge_y = 0.4) + xlab("Odds ratio (95% Confidence Interval)") + ylab("")

# Place
plotGG(plot = fig1B, x = 0.375, y = 1.75, width = 3, height = 2.5, just = c("left",
    "top"), default.units = "inches")

################### Figure 1C ####
Fig1C_pip = readPNG("/home/vagrant/Manuscript_1/FIGURE1/Fig1C_pipeline.png")

Fig1C_tab1.1 = readPNG("/home/vagrant/Manuscript_1/FIGURE1/Fig1C_SNPs_tab1.png/Fig1C_SNPs_tab1.png-1.png")
Fig1C_tab1.2 = readPNG("/home/vagrant/Manuscript_1/FIGURE1/Fig1C_SNPs_tab1.png/Fig1C_SNPs_tab1.png-2.png")

# Plot Figure 1C
plotRaster(image = Fig1C_pip, x = 0.375, y = 5.25, width = 3, height = 2.5, just = c("left",
    "top"), interpolate = FALSE)

plotRaster(image = Fig1C_tab1.1, x = 3.375, y = 5.25, width = 2.3, height = 2.5,
    just = c("left", "top"), interpolate = FALSE)

plotRaster(image = Fig1C_tab1.2, x = 5.8, y = 5.25, width = 2.3, height = 2.5, just = c("left",
    "top"), interpolate = FALSE)

# diag_edge_data = tibble(from = c('~17x10^6\nPD GWAS SNPs', '~3400
# significant\nPD GWAS SNPs', '28 significant\nPD rSNPs'), to = c('~3400
# significant\nPD GWAS SNPs', '28 significant\nPD rSNPs', '3 significant\nPD
# rSNPs')) diag_node_data = tibble(name = c('~17x10^6\nPD GWAS SNPs', '~3400
# significant\nPD GWAS SNPs', '28 significant\nPD rSNPs', '3 significant\nPD
# rSNPs'), label = name)

# Diagramme fig1C = ggflowchart(diag_edge_data, diag_node_data, text_size =
# 3.0, arrow_size = 0.2)

# Place the diagramme plotGG(plot = fig1C, x = 0.5, y = 4.5, width = 1.5,
# height = 3.0, just = c('left', 'top'), default.units = 'inches')

################### Figure 1D ####

# Locus plot Load the Nalls et al. GWAS data
nalls_allSNPs.HG38.GR = readRDS("/home/vagrant/Manuscript_1/FIGURE1/nalls_allSNPs.HG38.GR.rds")

# nalls_allSNPs.HG38.GR = readRDS(
# '/Volumes/deborah.gerard/Documents/Manuscript_1/FIGURE1/nalls_allSNPs.HG38.GR.rds')

# Here our SNPs of interest and their location
snp = GRanges(seqnames = c("chr10", "chr4", "chr4"), ranges = IRanges(start = c(119651404,
    987107, 76213716), end = c(119651405, 987108, 76213717)))

# Retrieve them from Nalls et al data
snp_Nalls = subsetByOverlaps(nalls_allSNPs.HG38.GR, snp) %>%
    as_tibble() %>%
    dplyr::select(-width, -strand) %>%
    mutate(rsID = c("rs144814361", "rs1465922", "rs11248061"))

# Visualise the genomic region per SNP (+- 400100 bp) rs144814361
r = 40100

# Filter for the SNP in BAG3 regulatory region
Nalls_sel_reg.BAG3 = nalls_allSNPs.HG38.GR %>%
    as_tibble() %>%
    dplyr::filter(seqnames == "chr10", between(base_pair_position, 119651405 - r,
        119651405 + r)) %>%
    mutate(minus_log_10_p = -log10(p), GWAS_p_sig = if_else(minus_log_10_p > -log10(5e-08),
        "yes", "no"), rsID = if_else(seqnames == "chr10" & base_pair_position ==
        119651405, "rs144814361", ""))

# Get SNPs info from Ensembl
snp.ensembl = useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

# SNP info
query.snp = getBM(c("ensembl_gene_stable_id", "refsnp_id", "chr_name", "chrom_start",
    "chrom_end", "minor_allele", "minor_allele_freq"), filters = "snp_filter", values = c("rs144814361",
    "rs1465922", "rs11248061"), mart = snp.ensembl)

# Get rsID of SNPs that are significant and append
position_chr10 = Nalls_sel_reg.BAG3 %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    filter(GWAS_p_sig != "no") %>%
    dplyr::select(position) %>%
    pull()


rsID_sig_SNPs.chr10 = getBM(c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
    "allele", "minor_allele"), filters = "chromosomal_region", values = position_chr10,
    mart = snp.ensembl)

rsID_sig_SNPs.chr10 = rsID_sig_SNPs.chr10 %>%
    unite(col = "position", chr_name:chrom_end, sep = ":", remove = FALSE)

Nalls_sel_reg.BAG3 = Nalls_sel_reg.BAG3 %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    left_join(rsID_sig_SNPs.chr10, by = "position")

# Calculate R2 to check if SNPs are in LD
BAG3_SNPs_LD = Nalls_sel_reg.BAG3 %>%
    dplyr::select(refsnp_id) %>%
    drop_na()

BAG3_R2 = LDproxy(snp = "rs144814361", pop = "EUR", r2d = "r2", token = "5bba40af90d7",
    genome_build = "grch38_high_coverage")

BAG3_R2 = BAG3_R2 %>%
    as_tibble() %>%
    dplyr::select(RS_Number, Coord, R2) %>%
    dplyr::rename(refsnp_id = RS_Number)

# Plot
BAG3.pl = Nalls_sel_reg.BAG3 %>%
    unite("Coord", seqnames, base_pair_position, sep = ":", remove = FALSE) %>%
    left_join(BAG3_R2, by = "Coord") %>%
    dplyr::select(seqnames, base_pair_position, p, refsnp_id.x, R2) %>%
    dplyr::rename(chrom = seqnames, pos = base_pair_position, snp = refsnp_id.x,
        LD = R2) %>%
    mutate(chrom = as.character(chrom)) %>%
    group_by(LD) %>%
    mutate(LD_grp = cut(LD, seq(0, 1, 0.2)))

BAG3.pl$LD_grp = addNA(BAG3.pl$LD_grp)

BAG3_leadSNP = "rs144814361"

BAG3_MP = plotManhattan(data = BAG3.pl, chrom = "chr10", chromstart = 119611405,
    chromend = 119691405, assembly = "hg38", fill = colorby("LD_grp", palette = colorRampPalette(c("#35b779",
        "#21918c", "#31688e", "#443983", "#440154", "#CCCCCC"))), trans = "-log10",
    sigLine = TRUE, col = "grey", lty = 2, range = c(0, 16), leadSNP = list(snp = BAG3_leadSNP,
        pch = 18, cex = 0.75, fill = "#990066", fontsize = 8), x = 0.5, y = 8.5,
    width = 1.5, height = 2.5, just = c("left", "top"), default.units = "inches")

## Annotate genome label
annoGenomeLabel(plot = BAG3_MP, x = 0.5, y = 11, fontsize = 8, scale = "Mb", just = c("left",
    "top"), default.units = "inches")

## Annotate y-axis
annoYaxis(plot = BAG3_MP, at = c(0, 2, 4, 6, 8, 10, 12), axisLine = TRUE, fontsize = 8)

## Plot y-axis label
plotText(label = "-log10(p-value)", x = 0.125, y = 10.125, rot = 90, fontsize = 8,
    fontface = "bold", just = "center", default.units = "inches")

# Plot gene BAG3
plotGenes(chrom = "chr10", chromstart = 119611405, chromend = 119691405, x = 0.5,
    y = 8.5, width = 1.5, height = 0.5, just = c("left", "top"), default.units = "inches",
    fontcolor = "black", fill = "black", assembly = assembly(Genome = "hg38refGene",
        TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"))

# Filter for the SNP in IDUA regulatory region
Nalls_sel_reg.IDUA = nalls_allSNPs.HG38.GR %>%
    as_tibble() %>%
    filter(seqnames == "chr4", between(base_pair_position, 987108 - r, 987108 + r)) %>%
    mutate(minus_log_10_p = -log10(p), GWAS_p_sig = if_else(minus_log_10_p > -log10(5e-08),
        "yes", "no"), rsID = if_else(seqnames == "chr4" & base_pair_position == 987108,
        "rs11248061", ""))

# Get rsID of SNPs that significant and append
position_chr4.IDUA = Nalls_sel_reg.IDUA %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    filter(GWAS_p_sig != "no") %>%
    dplyr::select(position) %>%
    pull()

rsID_sig_SNPs.chr4.IDUA = getBM(c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
    "allele", "minor_allele", "minor_allele_freq"), filters = "chromosomal_region",
    values = position_chr4.IDUA, mart = snp.ensembl)

rsID_sig_SNPs.chr4.IDUA = rsID_sig_SNPs.chr4.IDUA %>%
    unite(col = "position", chr_name:chrom_end, sep = ":", remove = FALSE)

Nalls_sel_reg.IDUA = Nalls_sel_reg.IDUA %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    left_join(rsID_sig_SNPs.chr4.IDUA, by = "position")

# Calculate R2 to check if SNPs are in LD
IDUA_SNPs_LD = Nalls_sel_reg.IDUA %>%
    dplyr::select(refsnp_id) %>%
    drop_na()

IDUA_SNPs_LD = IDUA_SNPs_LD %>%
    filter(!str_detect(refsnp_id, "^CR"))

IDUA_R2 = LDproxy(snp = "rs11248061", pop = "EUR", r2d = "r2", token = "5bba40af90d7",
    genome_build = "grch38_high_coverage")

IDUA_R2 = IDUA_R2 %>%
    as_tibble() %>%
    dplyr::select(RS_Number, Coord, R2) %>%
    dplyr::rename(refsnp_id = RS_Number)


# Plot
IDUA.pl = Nalls_sel_reg.IDUA %>%
    unite("Coord", seqnames, base_pair_position, sep = ":", remove = FALSE) %>%
    distinct(Coord, .keep_all = TRUE) %>%
    left_join(IDUA_R2, by = "Coord") %>%
    dplyr::select(seqnames, base_pair_position, p, refsnp_id.x, R2) %>%
    dplyr::rename(chrom = seqnames, pos = base_pair_position, snp = refsnp_id.x,
        LD = R2) %>%
    mutate(chrom = as.character(chrom)) %>%
    group_by(LD) %>%
    mutate(LD_grp = cut(LD, seq(0, 1, 0.2)))

IDUA.pl$LD_grp = addNA(IDUA.pl$LD_grp)

IDUA_leadSNP = "rs11248061"

IDUA_MP = plotManhattan(data = IDUA.pl, chrom = "chr4", chromstart = 947108, chromend = 1027108,
    assembly = "hg38", fill = colorby("LD_grp", palette = colorRampPalette(c("#35b779",
        "#21918c", "#31688e", "#443983", "#440154", "#CCCCCC"))), trans = "-log10",
    sigLine = TRUE, col = "grey", lty = 2, range = c(0, 16), leadSNP = list(snp = IDUA_leadSNP,
        pch = 18, cex = 0.75, fill = "#990066", fontsize = 8), x = 2.5, y = 8.5,
    width = 1.5, height = 2.5, just = c("left", "top"), default.units = "inches")

## Annotate genome label
annoGenomeLabel(plot = IDUA_MP, x = 2.5, y = 11, fontsize = 8, scale = "Mb", just = c("left",
    "top"), default.units = "inches")

## Annotate y-axis
annoYaxis(plot = IDUA_MP, at = c(0, 2, 4, 6, 8, 10, 12), axisLine = TRUE, fontsize = 8)

## Plot y-axis label
plotText(label = "-log10(p-value)", x = 2.125, y = 10.125, rot = 90, fontsize = 8,
    fontface = "bold", just = "center", default.units = "inches")

# Plot gene IDUA
plotGenes(chrom = "chr4", chromstart = 947108, chromend = 1027108, x = 2.5, y = 8.5,
    width = 1.5, height = 0.5, just = c("left", "top"), default.units = "inches",
    fontcolor = "black", fill = "black", assembly = assembly(Genome = "hg38refGene",
        TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"))

# Filter for the SNP in SCARB2 regulatory region
Nalls_sel_reg.SCARB2 = nalls_allSNPs.HG38.GR %>%
    as_tibble() %>%
    filter(seqnames == "chr4", between(base_pair_position, 76213717 - r, 76213717 +
        r)) %>%
    mutate(minus_log_10_p = -log10(p), GWAS_p_sig = if_else(minus_log_10_p > -log10(5e-08),
        "yes", "no"), rsID = if_else(seqnames == "chr4" & base_pair_position == 76213717,
        "rs1465922", ""))

# Get rsID of SNPs that significant and append
position_chr4.SCARB2 = Nalls_sel_reg.SCARB2 %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    filter(GWAS_p_sig != "no") %>%
    dplyr::select(position) %>%
    pull()

rsID_sig_SNPs.chr4.SCARB2 = getBM(c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
    "allele", "minor_allele", "minor_allele_freq"), filters = "chromosomal_region",
    values = position_chr4.SCARB2, mart = snp.ensembl)

rsID_sig_SNPs.chr4.SCARB2 = rsID_sig_SNPs.chr4.SCARB2 %>%
    unite(col = "position", chr_name:chrom_end, sep = ":", remove = FALSE)

Nalls_sel_reg.SCARB2 = Nalls_sel_reg.SCARB2 %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    unite(col = "position", c(chr, end, base_pair_position), sep = ":", remove = FALSE) %>%
    left_join(rsID_sig_SNPs.chr4.SCARB2, by = "position")

# Calculate R2 to check if SNPs are in LD
SCARB2_SNPs_LD = Nalls_sel_reg.SCARB2 %>%
    dplyr::select(refsnp_id) %>%
    drop_na()

SCARB2_R2 = LDproxy(snp = "rs1465922", pop = "EUR", r2d = "r2", token = "5bba40af90d7",
    genome_build = "grch38_high_coverage")

SCARB2_R2 = SCARB2_R2 %>%
    as_tibble() %>%
    dplyr::select(RS_Number, Coord, R2) %>%
    dplyr::rename(refsnp_id = RS_Number)

# Plot
SCARB2.pl = Nalls_sel_reg.SCARB2 %>%
    unite("Coord", seqnames, base_pair_position, sep = ":", remove = FALSE) %>%
    left_join(SCARB2_R2, by = "Coord") %>%
    dplyr::select(seqnames, base_pair_position, p, refsnp_id.x, R2) %>%
    dplyr::rename(chrom = seqnames, pos = base_pair_position, snp = refsnp_id.x,
        LD = R2) %>%
    mutate(chrom = as.character(chrom)) %>%
    group_by(LD) %>%
    mutate(LD_grp = cut(LD, seq(0, 1, 0.2)))

SCARB2.pl$LD_grp = addNA(SCARB2.pl$LD_grp)

SCARB2_leadSNP = "rs1465922"

SCARB2_MP = plotManhattan(data = SCARB2.pl, chrom = "chr4", chromstart = 76173717,
    chromend = 76253717, assembly = "hg38", fill = colorby("LD_grp", palette = colorRampPalette(c("#35b779",
        "#21918c", "#31688e", "#443983", "#440154", "#CCCCCC"))), trans = "-log10",
    sigLine = TRUE, col = "grey", lty = 2, range = c(0, 16), leadSNP = list(snp = SCARB2_leadSNP,
        pch = 18, cex = 0.75, fill = "#990066", fontsize = 8), x = 4.5, y = 8.5,
    width = 1.5, height = 2.5, just = c("left", "top"), default.units = "inches")

## Annotate genome label
annoGenomeLabel(plot = SCARB2_MP, x = 4.5, y = 11, fontsize = 8, scale = "Mb", just = c("left",
    "top"), default.units = "inches")

## Annotate y-axis
annoYaxis(plot = SCARB2_MP, at = c(0, 2, 4, 6, 8, 10, 12), axisLine = TRUE, fontsize = 8)

## Plot y-axis label
plotText(label = "-log10(p-value)", x = 4.125, y = 10.125, rot = 90, fontsize = 8,
    fontface = "bold", just = "center", default.units = "inches")

# Plot gene SCARB2
plotGenes(chrom = "chr4", chromstart = 76173717, chromend = 76253717, x = 4.5, y = 8.5,
    width = 1.5, height = 0.5, just = c("left", "top"), default.units = "inches",
    fontcolor = "black", fill = "black", assembly = assembly(Genome = "hg38refGene",
        TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"))

# Legend LD score
plotLegend(legend = c("LD reference SNP", paste("0", "<", "r^2", "<= 0.2"), paste("0.2",
    "<", "r^2", "<= 0.4"), paste("0.4", "<", "r^2", "<= 0.6"), paste("0.6", "<",
    "r^2", "<= 0.8"), paste("0.8", "<", "r^2", "<= 1.0"), "no LD data"), fill = c("#990066",
    "#35b779", "#21918c", "#31688e", "#443983", "#440154", "#CCCCCC"), cex = 0.75,
    pch = c(18, 19, 19, 19, 19, 19, 19), border = FALSE, x = 6.5, y = 9.25, width = 1,
    height = 1, just = c("left", "top"), default.units = "inches")

dev.off()
```

### FIGURE 2
Create the page layout that will contain all the necessary plots


```r
########## FIGURE 2 ##########

# Save as a PDF
pdf("/home/vagrant/Manuscript_1/FIGURE2/FIGURE2.pdf", width = 8.3, height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = FALSE)


#### PANEL A - BAG3 locus text Figure 2
plotText(label = "Figure 2", fontsize = 14, x = 0.25, y = 0.25, just = "left", default.units = "inches",
    fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# ATACseq signal BAG3 path where bigwig files are
bw.path = "/home/vagrant/epifunc/"

# 1st bio replicate smNPC
bw.smNPC_1 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_I.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 2nd bio replicate smNPC
bw.smNPC_2 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_II.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 3rd bio replicate smNPC
bw.smNPC_3 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_III.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 1st bio replicate TH positive neurons D15
bw.TH.pos_D15_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_I.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 2nd bio replicate TH positive neurons D15
bw.TH.pos_D15_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_II.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 3rd bio replicate TH positive neurons D15
bw.TH.pos_D15_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_III.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 1st bio replicate TH positive neurons D30
bw.TH.pos_D30_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_I.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 2nd bio replicate TH positive neurons D30
bw.TH.pos_D30_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_II.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 3rd bio replicate TH positive neurons D30
bw.TH.pos_D30_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_III.bw"), chrom = "chr10",
    chromstart = 119611305, chromend = 119691505)

# 1st bio replicate TH positive neurons D50
bw.TH.pos_D50_1 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# 2nd bio replicate TH positive neurons D50
bw.TH.pos_D50_2 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# 3rd bio replicate TH positive neurons D50
bw.TH.pos_D50_3 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
    chrom = "chr10", chromstart = 119611305, chromend = 119691505)

# Add a scale next to the bigwig files - check the maximum to choose
scale.BAG3.max = max(c(bw.smNPC_1$score, bw.smNPC_2$score, bw.smNPC_3$score, bw.TH.pos_D15_1$score,
    bw.TH.pos_D15_2$score, bw.TH.pos_D15_3$score, bw.TH.pos_D30_1$score, bw.TH.pos_D30_2$score,
    bw.TH.pos_D30_3$score, bw.TH.pos_D50_1$score, bw.TH.pos_D50_2$score, bw.TH.pos_D50_3$score))

print(scale.BAG3.max)

# Define parameters for the regions
region.p.BAG3 = pgParams(chrom = "chr10", chromstart = 119611405, chromend = 119691405,
    assembly = "hg38", range = c(0, scale.BAG3.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr10", chromstart = 119611405, chromend = 119691405, assembly = "hg38",
    x = 0.75, y = 0.625, length = 4, default.units = "inches", scale = "Mb")

# ATACseq signal smNPC_I and scale
ATAC_smNPC_I = plotSignal(data = bw.smNPC_1, params = region.p.BAG3, fill = "#D3D3D3",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 0.875, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_I, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal smNPC_II
ATAC_smNPC_II = plotSignal(data = bw.smNPC_2, params = region.p.BAG3, fill = "#C0C0C0",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 0.875, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_II, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal smNPC_III
ATAC_smNPC_III = plotSignal(data = bw.smNPC_3, params = region.p.BAG3, fill = "#A9A9A9",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 0.875, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 I
ATAC_TH.pos_D15_I = plotSignal(data = bw.TH.pos_D15_1, params = region.p.BAG3, fill = "#B22222",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 1.25, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_I, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 II
ATAC_TH.pos_D15_II = plotSignal(data = bw.TH.pos_D15_2, params = region.p.BAG3, fill = "#A52A2A",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 1.25, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_II, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 III
ATAC_TH.pos_D15_III = plotSignal(data = bw.TH.pos_D15_3, params = region.p.BAG3,
    fill = "#8B0000", alpha = 0.7, linecolor = NA, x = 0.75, y = 1.25, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_III, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 I
ATAC_TH.pos_D30_I = plotSignal(data = bw.TH.pos_D30_1, params = region.p.BAG3, fill = "#FF6347",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 1.625, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_I, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 II
ATAC_TH.pos_D30_II = plotSignal(data = bw.TH.pos_D30_2, params = region.p.BAG3, fill = "#FF0000",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 1.625, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_II, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 III
ATAC_TH.pos_D30_III = plotSignal(data = bw.TH.pos_D30_3, params = region.p.BAG3,
    fill = "#DC143C", alpha = 0.7, linecolor = NA, x = 0.75, y = 1.625, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_III, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 I
ATAC_TH.pos_D50_I = plotSignal(data = bw.TH.pos_D50_1, params = region.p.BAG3, fill = "#FFA07A",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 2, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_I, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 II
ATAC_TH.pos_D50_II = plotSignal(data = bw.TH.pos_D50_2, params = region.p.BAG3, fill = "#FA8072",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 2, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_II, at = c(0, scale.BAG3.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 III
ATAC_TH.pos_D50_III = plotSignal(data = bw.TH.pos_D50_3, params = region.p.BAG3,
    fill = "#E9967A", alpha = 0.7, linecolor = NA, x = 0.75, y = 2, width = 4, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_III, at = c(0, scale.BAG3.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr10", chromstart = 119611405, chromend = 119691405, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 2.375, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "smNPC", fontsize = 6, fontcolor = "#A9A9A9", x = 0.75, y = 0.875,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D15", fontsize = 6, fontcolor = "#B22222", x = 0.75, y = 1.25,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D30", fontsize = 6, fontcolor = "#DC143C", x = 0.75, y = 1.625,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D50", fontsize = 6, fontcolor = "#E9967A", x = 0.75, y = 2,
    just = c("left", "top"), default.units = "inches")

#### PANEL B - BAG3 and LHX1 expression (RNAseq) text B
plotText(label = "B", fontsize = 12, x = 5, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# Load RPKM matrix obtained from Borja and rename columns dat.RPKM =
# read_delim('/home/vagrant/epifunc/RNAseq/RPKM_genename', delim = '\t',
# col_names = TRUE) dat.RPKM = dat.RPKM %>% dplyr::select(gene_id, gene_name,
# everything())

# RPKM data - Remove the negative sorted neurons and astrocytes dat.RPKM.filt =
# dat.RPKM %>% dplyr::select(!contains('negsort'), -starts_with('ASTRO')) %>%
# pivot_longer(cols = ends_with('bam'), names_to = 'Sample', values_to =
# 'RPKM') %>% mutate(Cond = case_when(str_detect(Sample, 'D15_possort') ~
# 'Positively sorted neurons D15', str_detect(Sample, 'D30_possort') ~
# 'Positively sorted neurons D30', str_detect(Sample, 'D50_possort') ~
# 'Positively sorted neurons D50', TRUE ~ 'smNPCs'), Cond = factor(Cond, levels
# = c('smNPCs', 'Positively sorted neurons D15', 'Positively sorted neurons
# D30', 'Positively sorted neurons D50')), gene_id = gsub('\\..*', '',
# gene_id))

# Save the matrix of expression write_rds(dat.RPKM.filt,
# '/home/vagrant/Manuscript_1/FIGURE2/mat_RPKM.rds')

# Load it for making the plots
dat.RPKM.filt = read_rds("/home/vagrant/Manuscript_1/FIGURE2/mat_RPKM.rds")

# Expression for BAG3 and LHX1
BAG3.LHX1.exp = ggboxplot(dat.RPKM.filt %>%
    dplyr::filter(gene_name %in% c("BAG3", "LHX1")), x = "Cond", y = "RPKM", add = "jitter",
    fill = "Cond", palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"), facet.by = "gene_name",
    xlab = "", ylab = "Reads Per Kilobase Million (RPKM)") + scale_x_discrete(labels = c("smNPC",
    expression("mDAN D15"), expression("mDAN D30"), expression("mDAN D50"))) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))


# Place the plot
plotGG(plot = BAG3.LHX1.exp, x = 5.25, y = 0.625, width = 2, height = 3.5, just = c("left",
    "top"), default.units = "inches")

# Add the transcription factor binding motif (TFBS) for LHX1
pfm.LHX1 = getMatrixByID(JASPAR2020, ID = "MA1518.1")
LHX1 = new.env()

LHX1$LHX1 = pfm.LHX1@profileMatrix

LHX1 = as.list(LHX1)
LHX1_TFBS = ggseqlogo(LHX1)

# Place the plot
plotGG(plot = LHX1_TFBS, x = 1.7125, y = 2.625, width = 1.5, height = 1.5, just = c("left",
    "top"), default.units = "inches")

# Add the SNP position that zoom in the TFBS
SNP_param = pgParams(chrom = "chr10", chromstart = 119651405, chromend = 119651405,
    assembly = "hg38")

annoHighlight(plot = ATAC_smNPC_I, params = SNP_param, fill = "#404788FF", y = 0.875,
    height = 1.75, just = c("left", "top"), default.units = "inches")

annoZoomLines(plot = ATAC_smNPC_I, params = SNP_param, y0 = 2.625, x1 = c(2.25, 3.25),
    y1 = 2.75, default.units = "inches")

# text C
plotText(label = "C", fontsize = 12, x = 0.25, y = 4.25, just = "left", default.units = "inches",
    fontface = "bold")

# ATACseq signal IDUA 1st bio replicate smNPC
bw.IDUA.smNPC_1 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_I.bw"), chrom = "chr4",
    chromstart = 947008, chromend = 1027208)

# 2nd bio replicate smNPC
bw.IDUA.smNPC_2 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_II.bw"), chrom = "chr4",
    chromstart = 947008, chromend = 1027208)

# 3rd bio replicate smNPC
bw.IDUA.smNPC_3 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_III.bw"), chrom = "chr4",
    chromstart = 947008, chromend = 1027208)

# 1st bio replicate TH positive neurons D15
bw.IDUA.TH.pos_D15_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_I.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 2nd bio replicate TH positive neurons D15
bw.IDUA.TH.pos_D15_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_II.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 3rd bio replicate TH positive neurons D15
bw.IDUA.TH.pos_D15_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_III.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 1st bio replicate TH positive neurons D30
bw.IDUA.TH.pos_D30_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_I.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 2nd bio replicate TH positive neurons D30
bw.IDUA.TH.pos_D30_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_II.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 3rd bio replicate TH positive neurons D30
bw.IDUA.TH.pos_D30_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_III.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 1st bio replicate TH positive neurons D50
bw.IDUA.TH.pos_D50_1 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 2nd bio replicate TH positive neurons D50
bw.IDUA.TH.pos_D50_2 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# 3rd bio replicate TH positive neurons D50
bw.IDUA.TH.pos_D50_3 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Add a scale next to the bigwig files - check the maximum to choose
scale.IDUA.max = max(c(bw.IDUA.smNPC_1$score, bw.IDUA.smNPC_2$score, bw.IDUA.smNPC_3$score,
    bw.IDUA.TH.pos_D15_1$score, bw.IDUA.TH.pos_D15_2$score, bw.IDUA.TH.pos_D15_3$score,
    bw.IDUA.TH.pos_D30_1$score, bw.IDUA.TH.pos_D30_2$score, bw.IDUA.TH.pos_D30_3$score,
    bw.IDUA.TH.pos_D50_1$score, bw.IDUA.TH.pos_D50_2$score, bw.IDUA.TH.pos_D50_3$score))

print(scale.IDUA.max)

# Define parameters for the regions
region.p.IDUA = pgParams(chrom = "chr4", chromstart = 947108, chromend = 1027108,
    assembly = "hg38", range = c(0, scale.IDUA.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr4", chromstart = 947108, chromend = 1027108, assembly = "hg38",
    x = 0.75, y = 4.375, length = 4, default.units = "inches", scale = "Mb")

# ATACseq signal smNPC_I
ATAC_smNPC_I_IDUA = plotSignal(data = bw.IDUA.smNPC_1, params = region.p.IDUA, fill = "#D3D3D3",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 4.625, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_I_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal smNPC_II
ATAC_smNPC_II_IDUA = plotSignal(data = bw.IDUA.smNPC_2, params = region.p.IDUA, fill = "#C0C0C0",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 4.625, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_II_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# Add the SNP position that zoom in the TFBS
SNP_param_IDUA = pgParams(chrom = "chr4", chromstart = 987108, chromend = 987108,
    assembly = "hg38")

annoHighlight(plot = ATAC_smNPC_I_IDUA, params = SNP_param_IDUA, fill = "#404788FF",
    y = 4.75, height = 1.875, just = c("left", "top"), default.units = "inches")

annoZoomLines(plot = ATAC_smNPC_I_IDUA, params = SNP_param_IDUA, y0 = 6.625, x1 = c(2.25,
    3.25), y1 = 6.75, default.units = "inches")

# ATACseq signal smNPC_III
ATAC_smNPC_III_IDUA = plotSignal(data = bw.IDUA.smNPC_3, params = region.p.IDUA,
    fill = "#A9A9A9", alpha = 0.7, linecolor = NA, x = 0.75, y = 4.625, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 I
ATAC_TH.pos_D15_I_IDUA = plotSignal(data = bw.IDUA.TH.pos_D15_1, params = region.p.IDUA,
    fill = "#B22222", alpha = 0.7, linecolor = NA, x = 0.75, y = 5, width = 4, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_I_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 II
ATAC_TH.pos_D15_II_IDUA = plotSignal(data = bw.IDUA.TH.pos_D15_2, params = region.p.IDUA,
    fill = "#A52A2A", alpha = 0.7, linecolor = NA, x = 0.75, y = 5, width = 4, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_II_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 III
ATAC_TH.pos_D15_III_IDUA = plotSignal(data = bw.IDUA.TH.pos_D15_3, params = region.p.IDUA,
    fill = "#8B0000", alpha = 0.7, linecolor = NA, x = 0.75, y = 5, width = 4, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_III_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 I
ATAC_TH.pos_D30_I_IDUA = plotSignal(data = bw.IDUA.TH.pos_D30_1, params = region.p.IDUA,
    fill = "#FF6347", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_I_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 II
ATAC_TH.pos_D30_II_IDUA = plotSignal(data = bw.IDUA.TH.pos_D30_2, params = region.p.IDUA,
    fill = "#FF0000", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_II_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 III
ATAC_TH.pos_D30_III_IDUA = plotSignal(data = bw.IDUA.TH.pos_D30_3, params = region.p.IDUA,
    fill = "#DC143C", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_III_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 I
ATAC_TH.pos_D50_I_IDUA = plotSignal(data = bw.IDUA.TH.pos_D50_1, params = region.p.IDUA,
    fill = "#FFA07A", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_I_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 II
ATAC_TH.pos_D50_II_IDUA = plotSignal(data = bw.IDUA.TH.pos_D50_2, params = region.p.IDUA,
    fill = "#FA8072", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_II_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 III
ATAC_TH.pos_D50_III_IDUA = plotSignal(data = bw.IDUA.TH.pos_D50_3, params = region.p.IDUA,
    fill = "#E9967A", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_III_IDUA, at = c(0, scale.IDUA.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr4", chromstart = 947108, chromend = 1027108, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 6.125, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "smNPC", fontsize = 6, fontcolor = "#A9A9A9", x = 0.75, y = 4.625,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D15", fontsize = 6, fontcolor = "#B22222", x = 0.75, y = 5,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D30", fontsize = 6, fontcolor = "#DC143C", x = 0.75, y = 5.375,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D50", fontsize = 6, fontcolor = "#E9967A", x = 0.75, y = 5.75,
    just = c("left", "top"), default.units = "inches")

#### PANEL D - IDUA and ZBTB14 expression (RNAseq) text D
plotText(label = "D", fontsize = 12, x = 5, y = 4.25, just = "left", default.units = "inches",
    fontface = "bold")

IDUA.ZBTB14.exp = ggboxplot(dat.RPKM.filt %>%
    dplyr::filter(gene_name %in% c("IDUA", "ZBTB14", "SLC26A1")) %>%
    mutate(gene_name = factor(gene_name, levels = c("IDUA", "ZBTB14", "SLC26A1"))),
    x = "Cond", y = "RPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C", "#E9967A"), facet.by = "gene_name", xlab = "", ylab = "Reads Per Kilobase Million (RPKM)") +
    scale_x_discrete(labels = c("smNPC", "mDAN D15", "mDAN D30", "mDAN D50")) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))


# Place the plot
plotGG(plot = IDUA.ZBTB14.exp, x = 5.25, y = 4.25, width = 3, height = 3.5, just = c("left",
    "top"), default.units = "inches")

# Add the transcription factor binding motif (TFBS) for ZBTB14
pfm.ZBTB14 = getMatrixByID(JASPAR2020, ID = "MA1650.1")
ZBTB14 = new.env()

ZBTB14$ZBTB14 = pfm.ZBTB14@profileMatrix

ZBTB14 = as.list(ZBTB14)
ZBTB14_TFBS = ggseqlogo(ZBTB14)

# Place the plot
plotGG(plot = ZBTB14_TFBS, x = 1.7125, y = 6.625, width = 1.5, height = 1.5, just = c("left",
    "top"), default.units = "inches")

# text E
plotText(label = "E", fontsize = 12, x = 0.25, y = 8, just = "left", default.units = "inches",
    fontface = "bold")

# ATACseq signal SCARB2 1st bio replicate smNPC
bw.SCARB2.smNPC_1 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_I.bw"), chrom = "chr4",
    chromstart = 76173617, chromend = 76253817)

# 2nd bio replicate smNPC
bw.SCARB2.smNPC_2 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_II.bw"), chrom = "chr4",
    chromstart = 76173617, chromend = 76253817)

# 3rd bio replicate smNPC
bw.SCARB2.smNPC_3 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_III.bw"), chrom = "chr4",
    chromstart = 76173617, chromend = 76253817)

# 1st bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_I.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 2nd bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_II.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 3rd bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_III.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 1st bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_1 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_I.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 2nd bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_2 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_II.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 3rd bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_3 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_III.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 1st bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_1 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 2nd bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_2 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# 3rd bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_3 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
    chrom = "chr4", chromstart = 76173617, chromend = 76253817)

# Add a scale next to the bigwig files - check the maximum to choose
scale.SCARB2.max = max(c(bw.SCARB2.smNPC_1$score, bw.SCARB2.smNPC_2$score, bw.SCARB2.smNPC_3$score,
    bw.SCARB2.TH.pos_D15_1$score, bw.SCARB2.TH.pos_D15_2$score, bw.SCARB2.TH.pos_D15_3$score,
    bw.SCARB2.TH.pos_D30_1$score, bw.SCARB2.TH.pos_D30_2$score, bw.SCARB2.TH.pos_D30_3$score,
    bw.SCARB2.TH.pos_D50_1$score, bw.SCARB2.TH.pos_D50_2$score, bw.SCARB2.TH.pos_D50_3$score))

print(scale.SCARB2.max)

# Define parameters for the regions
region.p.SCARB2 = pgParams(chrom = "chr4", chromstart = 76173717, chromend = 76253717,
    assembly = "hg38", range = c(0, scale.SCARB2.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr4", chromstart = 76173717, chromend = 76253717, assembly = "hg38",
    x = 0.75, y = 8.125, length = 4, default.units = "inches", scale = "Mb")

# ATACseq signal smNPC_I
ATAC_smNPC_I_SCARB2 = plotSignal(data = bw.SCARB2.smNPC_1, params = region.p.SCARB2,
    fill = "#D3D3D3", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_I_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal smNPC_II
ATAC_smNPC_II_SCARB2 = plotSignal(data = bw.SCARB2.smNPC_2, params = region.p.SCARB2,
    fill = "#C0C0C0", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_II_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal smNPC_III
ATAC_smNPC_III_SCARB2 = plotSignal(data = bw.SCARB2.smNPC_3, params = region.p.SCARB2,
    fill = "#A9A9A9", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.375, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# Add the SNP position that zoom in the TFBS
SNP_param_SCARB2 = pgParams(chrom = "chr4", chromstart = 76213717, chromend = 76213717,
    assembly = "hg38")

annoHighlight(plot = ATAC_smNPC_I_SCARB2, params = SNP_param_SCARB2, fill = "#404788FF",
    y = 8.375, height = 1.875, just = c("left", "top"), default.units = "inches")

annoZoomLines(plot = ATAC_smNPC_I_SCARB2, params = SNP_param_SCARB2, y0 = 10.25,
    x1 = c(2.25, 3.25), y1 = 10.5, default.units = "inches")

# ATACseq signal TH positive neurons D15 I
ATAC_TH.pos_D15_I_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D15_1, params = region.p.SCARB2,
    fill = "#B22222", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_I_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 II
ATAC_TH.pos_D15_II_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D15_2, params = region.p.SCARB2,
    fill = "#A52A2A", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_II_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D15 III
ATAC_TH.pos_D15_III_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D15_3, params = region.p.SCARB2,
    fill = "#8B0000", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_III_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 I
ATAC_TH.pos_D30_I_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D30_1, params = region.p.SCARB2,
    fill = "#FF6347", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.125, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_I_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 II
ATAC_TH.pos_D30_II_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D30_2, params = region.p.SCARB2,
    fill = "#FF0000", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.125, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_II_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D30 III
ATAC_TH.pos_D30_III_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D30_3, params = region.p.SCARB2,
    fill = "#DC143C", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.125, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_III_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 I
ATAC_TH.pos_D50_I_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D50_1, params = region.p.SCARB2,
    fill = "#FFA07A", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_I_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 II
ATAC_TH.pos_D50_II_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D50_2, params = region.p.SCARB2,
    fill = "#FA8072", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_II_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# ATACseq signal TH positive neurons D50 III
ATAC_TH.pos_D50_III_SCARB2 = plotSignal(data = bw.SCARB2.TH.pos_D50_3, params = region.p.SCARB2,
    fill = "#E9967A", alpha = 0.7, linecolor = NA, x = 0.75, y = 9.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_III_SCARB2, at = c(0, scale.SCARB2.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr4", chromstart = 76173717, chromend = 76253717, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 9.875, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "smNPC", fontsize = 6, fontcolor = "#A9A9A9", x = 0.75, y = 8.375,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D15", fontsize = 6, fontcolor = "#B22222", x = 0.75, y = 8.75,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D30", fontsize = 6, fontcolor = "#DC143C", x = 0.75, y = 9.125,
    just = c("left", "top"), default.units = "inches")

plotText(label = "mDAN D50", fontsize = 6, fontcolor = "#E9967A", x = 0.75, y = 9.5,
    just = c("left", "top"), default.units = "inches")

# text F
plotText(label = "F", fontsize = 12, x = 5, y = 8, just = "left", default.units = "inches",
    fontface = "bold")

SCARB2.NR2C2.exp = ggboxplot(dat.RPKM.filt %>%
    dplyr::filter(gene_name %in% c("SCARB2", "NR2C2", "FAM47E")) %>%
    mutate(gene_name = factor(gene_name, levels = c("SCARB2", "NR2C2", "FAM47E"))),
    x = "Cond", y = "RPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C", "#E9967A"), facet.by = "gene_name", xlab = "", ylab = "Reads Per Kilobase Million (RPKM)") +
    scale_x_discrete(labels = c("smNPC", "mDAN D15", "mDAN D30", "mDAN D50")) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))


# Place the plot
plotGG(plot = SCARB2.NR2C2.exp, x = 5.25, y = 7.875, width = 3, height = 3.5, just = c("left",
    "top"), default.units = "inches")

# Add the transcription factor binding motif (TFBS) for NR2C2
pfm.NR2C2 = getMatrixByID(JASPAR2020, ID = "MA1536.1")

# Reverse to position weight matrix as the SNP is on the reverse strand
pfm.NR2C2 = reverseComplement(pfm.NR2C2)

NR2C2 = new.env()

NR2C2$NR2C2 = pfm.NR2C2@profileMatrix

NR2C2 = as.list(NR2C2)
NR2C2_TFBS = ggseqlogo(NR2C2)

# Place the plot
plotGG(plot = NR2C2_TFBS, x = 1.6875, y = 10.375, width = 1.5, height = 1.5, just = c("left",
    "top"), default.units = "inches")

dev.off()
```

### FIGURE 3

```r
########## FIGURE 3 ########## Save as a PDF
pdf("/home/vagrant/Manuscript_1/FIGURE3/FIGURE3.pdf", width = 8.3, height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = FALSE)

# text Figure 3
plotText(label = "Figure 3", fontsize = 14, x = 0.25, y = 0.25, just = "left", default.units = "inches",
    fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")


# Panel B text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 2, just = "left", default.units = "inches",
    fontface = "bold")

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


# Plot for BAG3
GLUC_BAG3 = ggbarplot(Gluc_signal_ratio %>%
    filter(Sample %in% c("BAG3-WT", "BAG3-MUT", "miniCMV", "Neg_CTRL")), x = "Sample",
    y = "value", width = 0.25, add = c("mean_se", "jitter"), fill = "grey", xlab = "",
    ylab = "Ratio Gluc/SEAP", position = position_dodge(0.9)) + ylim(0, 6) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10,
        face = "bold")) + scale_x_discrete(labels = c("BAG3-reference allele", "BAG3-PD associated allele",
    "miniCMV", "Negative CTRL")) + scale_y_continuous(expand = c(0, 0), limits = c(0,
    6)) + stat_compare_means(comparisons = list(c("BAG3-WT", "BAG3-MUT")), method = "t.test",
    size = 4, label = "p.format") + stat_compare_means(comparisons = list(c("BAG3-WT",
    "miniCMV")), method = "t.test", size = 4, label = "p.format")

plotGG(plot = GLUC_BAG3, x = 0.5, y = 2.25, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Plot for IDUA
GLUC_IDUA = ggbarplot(Gluc_signal_ratio %>%
    filter(Sample %in% c("IDUA-WT", "IDUA-MUT", "miniCMV", "Neg_CTRL")), x = "Sample",
    y = "value", width = 0.25, add = c("mean_se", "jitter"), fill = "grey", xlab = "",
    ylab = "Ratio Gluc/SEAP", position = position_dodge(0.9)) + ylim(0, 6) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10,
        face = "bold")) + scale_x_discrete(labels = c("IDUA-reference allele", "IDUA-PD associated allele",
    "miniCMV", "Negative CTRL")) + scale_y_continuous(expand = c(0, 0), limits = c(0,
    6)) + stat_compare_means(comparisons = list(c("IDUA-WT", "IDUA-MUT")), method = "t.test",
    size = 4, label = "p.format")

plotGG(plot = GLUC_IDUA, x = 3, y = 2.25, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Plot for SCARB2
GLUC_SCARB2 = ggbarplot(Gluc_signal_ratio %>%
    filter(Sample %in% c("SCARB2-WT", "SCARB2-MUT", "miniCMV", "Neg_CTRL")), x = "Sample",
    y = "value", width = 0.25, add = c("mean_se", "jitter"), fill = "grey", xlab = "",
    ylab = "Ratio Gluc/SEAP", position = position_dodge(0.9)) + ylim(0, 6) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10,
        face = "bold")) + scale_x_discrete(labels = c("SCARB2-reference allele",
    "SCARB2-PD associated allele", "miniCMV", "Negative CTRL")) + scale_y_continuous(expand = c(0,
    0), limits = c(0, 6)) + stat_compare_means(comparisons = list(c("SCARB2-WT",
    "SCARB2-MUT")), method = "t.test", size = 4, label = "p.format") + stat_compare_means(comparisons = list(c("SCARB2-WT",
    "miniCMV")), method = "t.test", size = 4, label = "p.format")

plotGG(plot = GLUC_SCARB2, x = 5.5, y = 2.25, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

# Panel C
plotText(label = "C", fontsize = 12, x = 0.25, y = 6.75, just = "left", default.units = "inches",
    fontface = "bold")

#### Knocdown ###
KD.BAG3.LHX1.3days = tibble(Sple = c(rep("shBAG3", 6), rep("shLHX1", 6)), Rep = rep(1:3,
    4), Gene = rep(c(rep("BAG3", 3), rep("LHX1", 3)), 2), Exp = c(0.446103385, 0.368235341,
    0.833613184, 0.921655819, 0.698775049, 0.766628746, 0.972688658, 0.965902853,
    0.88813445, 0.411253539, 0.359459217, 0.367190337))

###
test = tibble(val = c(9.13605, 8.33795, 8.76415, 9.0961, 8.2879, 8.593), Sample = c(rep("shLHX1",
    3), rep("shCTRL", 3)))

t_test(val ~ Sample, data = test, paired = TRUE)

###

pl.KD_BAG3_LHX1 = KD.BAG3.LHX1.3days %>%
    filter(Sple != "shBAG3") %>%
    mutate(Gene = factor(Gene, levels = c("LHX1", "BAG3"))) %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("Gene"), position = position_dodge(0.9),
        palette = c("#440154FF", "#287C8EFF"), add = c("mean_se", "jitter"), xlab = "",
        ylab = "Relative expression to shSCRAMBLE") + ylim(0, 1.5) + geom_hline(yintercept = 1,
    lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold", size = 10, angle = 90),
    axis.text.y = element_text(face = "bold", size = 10), legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) + stat_compare_means(comparisons = list(c("SCARB2-WT",
    "SCARB2-MUT")), method = "t.test", size = 4, label = "p.format")

plotGG(plot = pl.KD_BAG3_LHX1, x = 0.5, y = 7, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

KD.IDUA.ZBTB14.3days = tibble(Sple = c(rep("shIDUA", 6), rep("shZBTB14", 6)), Rep = rep(1:3,
    4), Gene = rep(c(rep("IDUA", 3), rep("ZBTB14", 3)), 2), Exp = c(0.580734239,
    0.722114742, 0.733236241, 0.743214215, 0.608719716, 0.676822765, 0.663813558,
    1.259979272, 1.774054093, 0.334134302, 0.316834205, 0.294206291))

pl.KD_IDUA_ZBTB14 = KD.IDUA.ZBTB14.3days %>%
    filter(Sple != "shIDUA") %>%
    mutate(Gene = factor(Gene, levels = c("ZBTB14", "IDUA"))) %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("Gene"), position = position_dodge(0.9),
        palette = c("#440154FF", "#287C8EFF"), add = c("mean_se"), xlab = "", ylab = "Relative expression to shSCRAMBLE") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold", size = 10, angle = 90),
    axis.text.y = element_text(face = "bold", size = 10), legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))

plotGG(plot = pl.KD_IDUA_ZBTB14, x = 3, y = 7, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

KD.SCARB2.NR2C2 = tibble(Sple = c(rep("shSCARB2", 12), rep("shNR2C2", 12)), Rep = rep(1:3,
    8), days = c(rep("3 days", 6), rep("6 days", 6), rep("3 days", 6), rep("6 days",
    6)), Gene = rep(c(rep("SCARB2", 3), rep("NR2C2", 3)), 4), Exp = c(0.751138057,
    0.592217899, 0.571232513, 1.135556902, 1.207647995, 0.955349154, 0.64510236,
    0.470902406, 0.792207355, 1.084251922, 0.999307093, 0.915638481, 1.239922544,
    1.130295504, 1.170885498, 0.345210918, 0.512864636, 0.548665969, 1.410640142,
    1.425976296, 1.332143962, 0.662710186, 0.456029798, 0.787662775))



pl.KD_SCARB2_NR2C2 = KD.SCARB2.NR2C2 %>%
    filter(Sple != "shSCARB2") %>%
    mutate(Gene = factor(Gene, levels = c("NR2C2", "SCARB2"))) %>%
    # filter(days == '3 days') %>% group_by(Sple, days, Gene) %>% mutate(AVG =
    # mean(Exp), SEM = sd(Exp)/sqrt(3)) %>% distinct(Sple, .keep_all = TRUE)
    # %>%
ggbarplot(., x = "Gene", y = "Exp", fill = c("Sple"), facet.by = "days", position = position_dodge(0.9),
    palette = c("#440154FF", "#287C8EFF"), add = c("mean_se"), xlab = "", ylab = "Relative expression to shSCRAMBLE") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold", size = 10, angle = 90),
    axis.text.y = element_text(face = "bold", size = 10), legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +  scale_y_continuous(expand
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +  = c(0, 0),
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +  limits = c(0,
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +  1.5)) +
plotGG(plot = pl.KD_SCARB2_NR2C2, x = 5.5, y = 7, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

dev.off()
```

### FIGURE 4

```r
# Save as a PDF
pdf("/home/vagrant/Manuscript_1/FIGURE4/FIGURE4.pdf",
    width = 8.3, 
    height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, 
           height = 11.7, 
           default.units = "inches",
           showGuides = FALSE)

plotText(label = "Figure 4", 
         fontsize = 14,
         x = 0.25, 
         y = 0.25, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")
# text A
plotText(label = "A", 
         fontsize = 12,
         x = 0.25, 
         y = 0.5, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

#FIG 4A
AI_RNA = tibble(Samples = 
                  c("smNPCs", "Neurons positively sorted D15", "Neurons positively sorted D30",
                    "Neurons positively sorted D50", "Astrocytes"),
                A = c(5, 23, 77, 56, 30),
                G = c(8, 4, 62, 34, 47)) %>% 
  gather(-Samples, key = "Base", value = "Reads")

AI_ATAC = tibble(Samples = 
                  c("smNPCs", "mDAN D15", "mDAN D30",
                    "mDAN D50", "Astrocytes"),
                A = c(31, 195, 202, 118, 155),
                G = c(39, 263, 264, 196, 179)) %>% 
  gather(-Samples, key = "Base", value = "Reads")

AI = AI_ATAC %>% 
  bind_rows(AI_RNA) %>% 
  mutate(Cond = rep(c("ATAC", "RNA"), each = 10))

# Good looking plot
chrom_AI = ggbarplot(AI_ATAC, x = "Samples", y = "Reads",
          fill = "Base",
          position = position_dodge(0.7),
          xlab = "", ylab = "Reads",
          #title = "Allelic imbalance at chr4:76213717 - SNP in SCARB2 promoter",
          label = c(0.3, 0.001, 0.004, 1e-5, 0.2,
                    rep("", 5)),
          label.pos = "out",
          lab.size = 3) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1.0, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10)) +
        # title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#80A44F", "#427F96"))

plotGG(plot = chrom_AI, 
       x = 0.5, 
       y = 0.75, 
       width = 6, 
       height = 4,
       just = c("left", "top"), 
       default.units = "inches")

# text B
plotText(label = "B", 
         fontsize = 12,
         x = 0.25, 
         y = 5.0, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

GTEX_SNPs_interest = read_delim("/home/vagrant/epifunc/eQTpLot/20230117_GTEX_data.csv",
               delim = ";",
               col_names = TRUE)

GTEX_SNPs_SCARB2_brain = GTEX_SNPs_interest %>% 
  mutate(GTEX_sig = if_else(P.Value < 0.05, "*", "")) %>% 
  filter(str_detect(Tissue, "Brain"),
         Gene.Symbol == "SCARB2")

GTEX_eqtl = ggplot(data = GTEX_SNPs_SCARB2_brain, mapping = aes(x = SNP.Id,
                                                y = Tissue,
                                                fill = NES)) + 
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold", size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10)) +
  scale_fill_gradient2(
    low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0,
    name = "eQTL \neffect size") +
  geom_text(aes(label = GTEX_sig)) 

plotGG(plot = GTEX_eqtl, 
       x = 0.5, 
       y = 5.25, 
       width = 7, 
       height = 6,
       just = c("left", "top"), 
       default.units = "inches")

dev.off()
```

### FIGURE 5

```r
########## FIGURE 5 ########## Save as a PDF
pdf("/home/vagrant/Manuscript_1/FIGURE5/FIGURE5.pdf", width = 8.3, height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = FALSE)

#### PANEL A - text Figure 5
plotText(label = "Figure 5", fontsize = 14, x = 0.25, y = 0.25, just = "left", default.units = "inches",
    fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

KD_TH_lev = tibble(Sple = rep("shNR2C2", 3), Rep = 1:3, days = "6 days", Gene = "TH",
    Exp = c(1.705920009, 2.196566231, 1.107610012))

# TH levels in knockdown of NR2C2 6 days#
TH_lev_KD = KD_TH_lev %>%
    ggbarplot(., x = "Gene", y = "Exp", fill = c("days"), facet.by = "Sple", position = position_dodge(0.9),
        palette = c("#287C8EFF"), add = c("mean_se"), xlab = "", ylab = "Relative expression to shSCRAMBLE") +
    geom_hline(yintercept = 1, lty = "dashed", color = "black") + theme(axis.title.y = element_text(face = "bold",
    size = 10), axis.text.x = element_text(face = "bold", size = 10, angle = 90),
    axis.text.y = element_text(face = "bold", size = 10), legend.position = "none")

plotGG(plot = TH_lev_KD, x = 0.5, y = 0.5, width = 2, height = 4.25, just = c("left",
    "top"), default.units = "inches")

### DEXOMAG - 3 biological repliactes ###
GCase = read_excel("/Volumes/deborah.gerard/Documents/Wetlab/epifunc/DEXOMAG_shSCARB2_shNR2C2_shCTRL_N1_N2_N3/20231214_DEXOMAG_KD_SCARB2_NR2C2_SCRAMBLE_N1_N2_N3_ANALYSIS.xlsx")

GCase %>%
    dplyr::select(GCase, `GCase MEAN`, ...10, ...11) %>%
    drop_na() %>%
    dplyr::rename(Sample = GCase, N1 = `GCase MEAN`, N2 = ...10, N3 = ...11) %>%
    pivot_longer(N1:N3, names_to = "Rep.bio", values_to = "norm.fluo.value") %>%
    mutate(norm.fluo.value = as.numeric(norm.fluo.value)) %>%
    # dplyr::filter(Sample != 'Whole Gcase NT') %>%
ggbarplot(., x = "Sample", y = "norm.fluo.value", add = c("mean_se", "jitter"), xlab = "",
    ylab = "GCase activity (normalised fluorescence values)") + scale_x_discrete(labels = c("SCARB2 KD",
    "NR2C2 KD", "SCRAMBLE", "Not transduced + dexomag", "Not transduced - dexomag",
    "Whole cell lysate")) + theme(axis.text.x = element_text(angle = 45, face = "bold",
    hjust = 1), axis.text.y = element_text(face = "bold"), axis.title.y = element_text(face = "bold")) +
    stat_compare_means(comparisons = list(c("shSCARB2", "shSCRAMBLE")), method = "t.test",
        size = 4, label = "p.format") + stat_compare_means(comparisons = list(c("shNR2C2",
    "shSCRAMBLE")), method = "t.test", size = 4, label = "p.format")


dev.off()
```

### SUPPLEMENTARY FIGURE 1
IHEC consortium has one sample of substantia nigra for H3K27me3, H3K9me3, H3K4me1, H3K4me3 and H3K36me3 from a healthy brain for BAG3, IDUA and SCARB2 loci

```r
############################################ SUPPLEMENTARY FIGURE 1 ##########

# Save as a PDF
pdf("/home/vagrant/Manuscript_1/SUPPLEMENTARY_FIGURE1/SUPPLEMENTARY_FIGURE1.pdf",
    width = 8.3, height = 11.7)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = TRUE)

#### PANEL A - text Supplementary Figure 1
plotText(label = "Supplementary Figure 1", fontsize = 14, x = 0.25, y = 0.25, just = "left",
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

# Reads BIGWIG file at IDUA location - H3K4me1
bw.IDUA.H3K4me1 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6a1889bf-7312-4e41-b1cc-e2f3be7a8076.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Reads BIGWIG file at IDUA location - H3K4me3
bw.IDUA.H3K4me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.4cb9300f-9097-4374-8f2a-4afd8ff855fb.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Reads BIGWIG file at IDUA location - H3K36me3
bw.IDUA.H3K36me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.6586e624-852f-4cac-83df-2325cee19699.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Reads BIGWIG file at IDUA location - H3K9me3
bw.IDUA.H3K9me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.e8fe7664-906e-4d87-8f3a-6ce88dbeadc6.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Reads BIGWIG file at IDUA location - H3K27me3
bw.IDUA.H3K27me3 = readBigwig(file = paste0(bw.IHEC.path, "ihec.chipseq.ihec-chipseq-containerv1.1.4.IHECRE00000986.6.5f5024a9-2e3f-4527-82ec-0ed04bf4a848.fc.signal.bigwig"),
    chrom = "chr4", chromstart = 947008, chromend = 1027208)

# Add a scale next to the bigwig files - check the maximum to choose
scale.IHEC.IDUA.max = max(c(bw.IDUA.H3K4me1$score, bw.IDUA.H3K4me3$score, bw.IDUA.H3K36me3$score,
    bw.IDUA.H3K9me3$score, bw.IDUA.H3K27me3$score)) %>%
    round(., digits = 1)

print(scale.IHEC.IDUA.max)

# Define parameters for the regions
region.p.IHEC.IDUA = pgParams(chrom = "chr4", chromstart = 947108, chromend = 1027108,
    assembly = "hg38", range = c(0, scale.IHEC.IDUA.max))

# Add the genomic label
plotGenomeLabel(chrom = "chr4", chromstart = 947108, chromend = 1027108, assembly = "hg38",
    x = 0.75, y = 3.5, length = 4, default.units = "inches", scale = "Mb")

# ChIPseq signal H3K4me1
H3K4me1_IDUA = plotSignal(data = bw.IDUA.H3K4me1, params = region.p.IHEC.IDUA, fill = "#b5de2b",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 3.75, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K4me1_IDUA, at = c(0, scale.IHEC.IDUA.max), fontsize = 6)

# ChIPseq signal H3K4me3
H3K4me3_IDUA = plotSignal(data = bw.IDUA.H3K4me3, params = region.p.IHEC.IDUA, fill = "#35b779",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 4.125, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K4me3_IDUA, at = c(0, scale.IHEC.IDUA.max), fontsize = 6)

# ChIPseq signal H3K36me3
H3K36me3_IDUA = plotSignal(data = bw.IDUA.H3K36me3, params = region.p.IHEC.IDUA,
    fill = "#26828e", alpha = 0.7, linecolor = NA, x = 0.75, y = 4.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K36me3_IDUA, at = c(0, scale.IHEC.IDUA.max), fontsize = 6)

# ChIPseq signal H3K9me3
H3K9me3_IDUA = plotSignal(data = bw.IDUA.H3K9me3, params = region.p.IHEC.IDUA, fill = "#3e4989",
    alpha = 0.7, linecolor = NA, x = 0.75, y = 4.875, width = 4, height = 0.25, just = c("left",
        "top"), default.units = "inches")

annoYaxis(plot = H3K9me3_IDUA, at = c(0, scale.IHEC.IDUA.max), fontsize = 6)

# ChIPseq signal H3K27me3
H3K27me3_IDUA = plotSignal(data = bw.IDUA.H3K27me3, params = region.p.IHEC.IDUA,
    fill = "#440154", alpha = 0.7, linecolor = NA, x = 0.75, y = 5.25, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K27me3_IDUA, at = c(0, scale.IHEC.IDUA.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr4", chromstart = 947108, chromend = 1027108, assembly = assembly(Genome = "hg38refGene",
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

# ChIP-seq signal SCARB2 Reads BIGWIG file at SCARB2 location - H3K4me1
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
    x = 0.75, y = 6.5, length = 4, default.units = "inches", scale = "Mb")

# ChIPseq signal H3K4me1
H3K4me1_SCARB2 = plotSignal(data = bw.SCARB2.H3K4me1, params = region.p.IHEC.SCARB2,
    fill = "#b5de2b", alpha = 0.7, linecolor = NA, x = 0.75, y = 6.75, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K4me1_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K4me3
H3K4me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K4me3, params = region.p.IHEC.SCARB2,
    fill = "#35b779", alpha = 0.7, linecolor = NA, x = 0.75, y = 7.125, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K4me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K36me3
H3K36me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K36me3, params = region.p.IHEC.SCARB2,
    fill = "#26828e", alpha = 0.7, linecolor = NA, x = 0.75, y = 7.5, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K36me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K9me3
H3K9me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K9me3, params = region.p.IHEC.SCARB2,
    fill = "#3e4989", alpha = 0.7, linecolor = NA, x = 0.75, y = 7.875, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K9me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# ChIPseq signal H3K27me3
H3K27me3_SCARB2 = plotSignal(data = bw.SCARB2.H3K27me3, params = region.p.IHEC.SCARB2,
    fill = "#440154", alpha = 0.7, linecolor = NA, x = 0.75, y = 8.25, width = 4,
    height = 0.25, just = c("left", "top"), default.units = "inches")

annoYaxis(plot = H3K27me3_SCARB2, at = c(0, scale.IHEC.SCARB2.max), fontsize = 6)

# Gene track
plotGenes(chrom = "chr4", chromstart = 76173717, chromend = 76253717, assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), fontcolor = "black",
    fill = "black", x = 0.75, y = 8.625, width = 4, height = 0.375, just = c("left",
        "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me1", fontsize = 6, fontcolor = "#b5de2b", x = 0.75, y = 6.75,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K4me3", fontsize = 6, fontcolor = "#35b779", x = 0.75, y = 7.125,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K36me3", fontsize = 6, fontcolor = "#26828e", x = 0.75, y = 7.5,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K9me3", fontsize = 6, fontcolor = "#3e4989", x = 0.75, y = 7.875,
    just = c("left", "top"), default.units = "inches")

# Add the samples name
plotText(label = "H3K27me3", fontsize = 6, fontcolor = "#440154", x = 0.75, y = 8.25,
    just = c("left", "top"), default.units = "inches")

#### PANEL B - text Supplementary Figure 1 text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 8.5, just = "left", default.units = "inches",
    fontface = "bold")

# Load IHEC RNA-seq data (TPM)
IHEC.exp = read_delim("/home/vagrant/Documents/IHEC/RNA_seq/genes_tpm.csv", delim = ",",
    col_names = TRUE)

IHEC.exp %>%
    mutate(id_col = gsub("\\..*", "", id_col))

dev.off()
```

### SUPPLEMENTARY FIGURE 2
Plot genes of interest from the TH REP2 mCHERRY cell line RNAseq data to demonstrate that expression is similar

```r
############################################ SUPPLEMENTARY FIGURE 2 ##########

# Save as a PDF
pdf("/home/vagrant/Manuscript_1/SUPPLEMENTARY_FIGURE2/SUPPLEMENTARY_FIGURE2.pdf",
    width = 13, height = 15)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = FALSE)

#### PANEL A - text Supplementary Figure 2
plotText(label = "Supplementary Figure 2", fontsize = 14, x = 0.25, y = 0.25, just = "left",
    default.units = "inches", fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# Load data
dat.THREP2.FPKM = read_delim("/home/vagrant/epifunc/RNAseq/TH_REP2_mCHERRY/180604_HT-Rep2_fpkm.tsv",
    delim = "\t", col_names = TRUE)

dat.THREP2.FPKM %>%
    dplyr::filter(gene_name == "SCARB2") %>%
    summary()

# Filter for genes of interest - BAG3 and LHX1
BAG3.LHX1 = c("BAG3", "LHX1")

IDUA.ZBTB14 = c("IDUA", "ZBTB14", "SLC26A1")

SCARB2.NR2C2 = c("SCARB2", "NR2C2", "FAM47E")

# Reformat the matrix and take only positively sorted neurons
dat.THREP2.FPKM = dat.THREP2.FPKM %>%
    pivot_longer(starts_with("TH"), names_to = "Samples", values_to = "FPKM") %>%
    dplyr::filter(str_detect(Samples, paste("SmNPC", "pos", sep = "|"))) %>%
    mutate(Cond = case_when(str_detect(Samples, "D15pos") ~ "Positively sorted neurons D15",
        str_detect(Samples, "D30pos") ~ "Positively sorted neurons D30", TRUE ~ "smNPCs"),
        Cond = factor(Cond, levels = c("smNPCs", "Positively sorted neurons D15",
            "Positively sorted neurons D30")))


# Plot the expression of BAG3 and LHX1
TH.REP2_BAG3.LHX1.exp = ggboxplot(dat.THREP2.FPKM %>%
    dplyr::filter(gene_name %in% BAG3.LHX1), x = "Cond", y = "FPKM", add = "jitter",
    fill = "Cond", palette = c("#A9A9A9", "#B22222", "#DC143C"), facet.by = "gene_name",
    xlab = "", ylab = "Fragments Per Kilobase Million (FPKM)") + scale_x_discrete(labels = c("smNPC",
    expression("mDAN D15"), expression("mDAN D30"))) + theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = TH.REP2_BAG3.LHX1.exp, x = 0.5, y = 0.625, width = 2, height = 3.5,
    just = c("left", "top"), default.units = "inches")

# Plot the expression of IDUA and ZBTB14
TH.REP2_IDUA.ZBTB14.exp = ggboxplot(dat.THREP2.FPKM %>%
    dplyr::filter(gene_name %in% IDUA.ZBTB14) %>%
    mutate(gene_name = factor(gene_name, levels = c("IDUA", "ZBTB14", "SLC26A1"))),
    x = "Cond", y = "FPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C"), facet.by = "gene_name", xlab = "", ylab = "Fragments Per Kilobase Million (FPKM)") +
    scale_x_discrete(labels = c("smNPC", expression("mDAN D15"), expression("mDAN D30"))) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = TH.REP2_IDUA.ZBTB14.exp, x = 2.75, y = 0.625, width = 2.5, height = 3.5,
    just = c("left", "top"), default.units = "inches")

# Plot the expression of SCARB2 and NR2C2
TH.REP2_SCARB2.NR2C2.exp = ggboxplot(dat.THREP2.FPKM %>%
    dplyr::filter(gene_name %in% SCARB2.NR2C2) %>%
    mutate(gene_name = factor(gene_name, levels = c("SCARB2", "NR2C2", "FAM47E"))),
    x = "Cond", y = "FPKM", add = "jitter", fill = "Cond", palette = c("#A9A9A9",
        "#B22222", "#DC143C"), facet.by = "gene_name", xlab = "", ylab = "Fragments Per Kilobase Million (FPKM)") +
    scale_x_discrete(labels = c("smNPC", expression("mDAN D15"), expression("mDAN D30"))) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))

# Place the plot
plotGG(plot = TH.REP2_SCARB2.NR2C2.exp, x = 5.5, y = 0.625, width = 2.5, height = 3.5,
    just = c("left", "top"), default.units = "inches")

#### PANEL B - text B
plotText(label = "B", fontsize = 12, x = 0.25, y = 4.5, just = "left", default.units = "inches",
    fontface = "bold")

# FOUNDIN-PD data Load
FIPD = read_delim("/home/vagrant/epifunc/PPMI_data/FOUNDIN-PD_150.9_RNAB/aggregated_expression/cpmTable.tsv",
    delim = "\t", col_names = TRUE)

# FIPD =
# read_delim('/home/vagrant/epifunc/PPMI_data/FOUNDIN-PD_150.9_RNAB/aggregated_expression/countTable.tsv',
# delim = '\t', col_names = TRUE)

# Filter for genes of interest and get the ensembl and gene names
goi_ens = select(org.Hs.eg.db, keys = c("BAG3", "LHX1", "IDUA", "ZBTB14", "SLC26A1",
    "SCARB2", "NR2C2", "FAM47E"), columns = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")

# 2 LHX1 ensembl IDs came back. Check on the Ensembl website which one is
# correct
goi_ens = goi_ens %>%
    as_tibble() %>%
    dplyr::filter(ENSEMBL != "ENSG00000274577") %>%
    dplyr::rename(ENSEMBL.not.dot = ENSEMBL)

# Load the information abour iPSCs cell lines
info_iPSCs = read_delim("/home/vagrant/epifunc/PPMI_data/Participant_Status_02Jun2023.csv",
    delim = ",", col_names = TRUE)

info_iPSCs = info_iPSCs %>%
    mutate(PATNO = as.character(PATNO))

# Filter the matrix of expression for my genes of interest
FIPD.fin = FIPD %>%
    # mutate(Geneid = str_replace(Geneid, '\\..', '')) %>%
dplyr::filter(str_detect(Geneid, paste(goi_ens$ENSEMBL.not.dot, collapse = "|"))) %>%
    dplyr::rename(ENSEMBL = Geneid) %>%
    mutate(ENSEMBL.not.dot = gsub("\\..*", "", ENSEMBL)) %>%
    dplyr::select(ENSEMBL, ENSEMBL.not.dot, everything()) %>%
    # left_join(goi_ens, by = 'ENSEMBL') %>%
pivot_longer(cols = starts_with("RNAB"), names_to = "Samples", values_to = "CPM") %>%
    mutate(PATNO = str_extract(Samples, "[^RNAB_PPMI_]+"), Day = str_extract(Samples,
        "da+\\d+")) %>%
    left_join(info_iPSCs, by = "PATNO") %>%
    left_join(goi_ens, by = "ENSEMBL.not.dot") %>%
    dplyr::select(ENSEMBL:Day, COHORT_DEFINITION, SYMBOL) %>%
    mutate(Condition = case_when(COHORT_DEFINITION == "Healthy Control" ~ "Healthy\ncontrol",
        COHORT_DEFINITION == "Parkinson's Disease" ~ "Parkinson's\ndisease", TRUE ~
            COHORT_DEFINITION)) %>%
    # group_by(SYMBOL, Day, COHORT_DEFINITION) %>%
group_by(ENSEMBL, Day) %>%
    mutate(AVG_CPM = mean(CPM, na.rm = TRUE))

# Plot expression of all putative target genes and TFs except SCARB2 (will be
# plot alone)
FIPD.exp.pl = ggboxplot(FIPD.fin %>%
    dplyr::filter(SYMBOL != "SCARB2"), x = "Day", y = "AVG_CPM", add = "jitter",
    fill = "Day", palette = c("#482677FF", "#2D708EFF", "#29AF7FFF"), facet.by = "SYMBOL",
    xlab = "", ylab = "CPM") + scale_x_discrete(labels = c("Day 0", "Day 25", "Day 65")) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 6))

# Place the plot
plotGG(plot = FIPD.exp.pl, x = 0.5, y = 4.75, width = 6, height = 4.5, just = c("left",
    "top"), default.units = "inches")

# Plot expression of SCARB2
FIPD.exp.SCARB2.pl = ggboxplot(FIPD.fin %>%
    dplyr::filter(SYMBOL == "SCARB2"), x = "Day", y = "AVG_CPM", add = "jitter",
    fill = "Day", palette = c("#482677FF", "#2D708EFF", "#29AF7FFF"), facet.by = "SYMBOL",
    xlab = "", ylab = "CPM") + scale_x_discrete(labels = c("Day 0", "Day 25", "Day 65")) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 6))

# Place the plot
plotGG(plot = FIPD.exp.SCARB2.pl, x = 0.5, y = 9.75, width = 6, height = 3.5, just = c("left",
    "top"), default.units = "inches")

dev.off()
```

### SUPPLEMENTARY FIGURE 3
Plot FACS results for experiments where knockdown were performed at day 9 and samples collected 3 and 6 days post-transduction

```r
############################################ SUPPLEMENTARY FIGURE 3 ##########

# Save as a PDF
pdf("/home/vagrant/Manuscript_1/SUPPLEMENTARY_FIGURE3/SUPPLEMENTARY_FIGURE3.pdf",
    width = 13, height = 15)

# Create a A4 blank page
pageCreate(width = 8.3, height = 11.7, default.units = "inches", showGuides = TRUE)

#### PANEL A - text Supplementary Figure 3
plotText(label = "Supplementary Figure 3", fontsize = 14, x = 0.25, y = 0.25, just = "left",
    default.units = "inches", fontface = "bold")
# text A
plotText(label = "A", fontsize = 12, x = 0.25, y = 0.5, just = "left", default.units = "inches",
    fontface = "bold")

# Read the facs files (.fcs) facs.path =
# '/home/vagrant/Documents/Wetlab/epifunc/'
facs.path = "/Volumes/deborah.gerard/Documents/Wetlab/epifunc/"

facs_6days_PT_N1.fs = read.flowSet(path = paste0(facs.path, "2022.11_DG_KD_SCARB2_NR2C2_1/FACS/RAW_DATA/NEW_EXPORT_06122022_1624/2022.11_DG_KD_SCARB2_NR2C2_1/"),
    pattern = ".fcs", alter.names = TRUE, truncate_max_range = FALSE)

# Check what channels are present
colnames(facs_6days_PT_N1.fs)

# Check sample
pData(facs_6days_PT_N1.fs)

# Rename some for clarity
colnames(facs_6days_PT_N1.fs)[colnames(facs_6days_PT_N1.fs) == "FITC.A"] = "GFP"
colnames(facs_6days_PT_N1.fs)[colnames(facs_6days_PT_N1.fs) == "mCherry..A"] = "mCHERRY"
colnames(facs_6days_PT_N1.fs)[colnames(facs_6days_PT_N1.fs) == "SYTOX_BLUE..A"] = "SYTOX BLUE"

# To be able to add gates, the object must be transformed as a GatingSet object
facs_6days_PT_N1.gs = GatingSet(facs_6days_PT_N1.fs)

# Define neurons gate
dead.neur.gate = polygonGate(filterId = "Neurons", FSC.A = c(50000, 2e+05, 2e+05,
    50000), SSC.A = c(0, 125000, 2e+05, 1e+05))

# Plot SYTOX BLUE positive CTRL
ggcyto(facs_6days_PT_N1.gs[[1]], aes(x = FSC.A, y = SSC.A), subset = "root") + geom_hex(bins = 200) +
    geom_gate(dead.neur.gate) + ggcyto_par_set(limits = "instrument")

# The gate for dead neurons has been defined -> add the gate to the gateset and
# recompute
add(facs_6days_PT_N1.gs, dead.neur.gate)

recompute(facs_6days_PT_N1.gs)

dev.off()
```

SCARB2 snp (rs1465922) for samples with whole genome sequencing data (75 cases)

```r
# Load data
SCARB2.SNP.WGS.75 = read_excel("/Volumes/deborah.gerard/Documents/epifunc/Final.WGS.SCARB2.OK_from_Sinthu.xlsx")

# Filter for rs1465922
SCARB2.SNP.WGS.75 %>%
    filter(avsnp150 == "rs1465922") %>%
    select(samples) %>%
    separate(samples, sep = ", ", into = )
```



```r
sessionInfo()
```

```
## R version 4.3.1 (2023-06-16)
## Platform: x86_64-apple-darwin20 (64-bit)
## Running under: macOS Sonoma 14.2.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Luxembourg
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.33     R6_2.5.1          fastmap_1.1.1     xfun_0.41        
##  [5] cachem_1.0.8      knitr_1.45        htmltools_0.5.7   rmarkdown_2.25   
##  [9] cli_3.6.1         sass_0.4.7        jquerylib_0.1.4   compiler_4.3.1   
## [13] rstudioapi_0.15.0 tools_4.3.1       evaluate_0.23     bslib_0.5.1      
## [17] yaml_2.3.7        formatR_1.14      rlang_1.1.2       jsonlite_1.8.7
```

