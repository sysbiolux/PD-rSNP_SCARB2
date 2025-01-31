---
title: "Figure2"
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

The SNEEP analysis will be run again by Nina with the PD GWAS SNPs that been filtered only for the GWAS p-value of 5e-08

``` r
# This is the file obtained by Jochen
nalls_allSNPs = vroom(
  "/home/vagrant/epifunc/SNPs_overlap/nallsEtAl2019_hg38_full.tsv")

# Select the coordinates, the alleles and the p-value columns. Filter for the GWAS p-value
nalls_allSNPs %>% 
  dplyr::select(CHR_hg38:END_hg38, A1:A2, p) %>% 
  dplyr::filter(p < 5e-08) %>% # GWAS p-value threshold
  write_delim("/home/vagrant/epifunc/SNEEP2024/input/Nalls_PD_SNPs_GWASpval_filtered.txt",
              delim = "\t",
              col_names = TRUE)

# Provide the RPKM data from the RNAseq to Nina
dat.RPKM = read_delim("/home/vagrant/epifunc/RNAseq/RPKM_genename",
                      delim = "\t", 
                      col_names = TRUE)

# Provide the file as a text file
dat.RPKM %>% 
  write_delim(.,
              "/home/vagrant/epifunc/SNEEP2024/input/RNAseq/All_samples_RPKM.txt",
              delim = "\t",
              col_names = TRUE)

# Check the hg38 coordinates 
nalls_allSNPs %>% 
  dplyr::filter(CHR_hg38 == "chr4",
                END_hg38 == "76213717")
```

For the SNEEP analysis, retrieve also the `narrowPeak` files and send them to Nina

``` bash
# narrowPeak files for smNPC, mDAN neurons day 15 - day 30 - day 50
rsync -avzPhu --no-p iris-cluster:/work/projects/epifunc/manuscript_sample_submission/ATAC/peaks/*narrowPeak /Volumes/deborah.gerard/Documents/epifunc/SNEEP2024/input/ATAC/narrowPeaks/

# narrowPeak files for astrocytes and non-mDAN neurons at day 15 and day 50 (day 30 data does not exist)
# According to the epifunc gitlab (https://git-r3lab.uni.lu/jochen.ohnmacht/epifunc) they are in /work/projects/epifunc/genrich_replicates_opt. According the epifunc gitlab, the filtered.narrowPeak files were used for HINT-ATAC and EPIC-DREM
rsync -avzPhu --no-p iris-cluster:/work/projects/epifunc/genrich_replicates_opt/*neg.filtered.narrowPeak /Volumes/deborah.gerard/Documents/epifunc/SNEEP2024/input/ATAC/narrowPeaks/

rsync -avzPhu --no-p iris-cluster:/work/projects/epifunc/genrich_replicates_opt/astrocytes.filtered.narrowPeak /Volumes/deborah.gerard/Documents/epifunc/SNEEP2024/input/ATAC/narrowPeaks/
```

Process the results of the SNEEP pipeline

``` r
## Take the smNPCs, TH+ neurons at day 15 - day 30 - day 50 Load the data smNPC
SNEEP.smNPC = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/parkinsonDisease_project_Deborah_Lasse/sneep_witATAC_mergedlowC_2409/result_PD_smNPV_lowC_ATAC_merged_2409.txt",
    delim = "\t", col_names = TRUE)

# Add a column specifying the sample name
SNEEP.smNPC = SNEEP.smNPC %>%
    mutate(Sample = "smNPC")

# Unique rSNPs
SNEEP.smNPC %>%
    distinct(SNP_position) %>%
    write_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/uniq_rSNPs_smNPC.txt",
        delim = "\t", col_names = FALSE)

# How many unique target genes associated to those 91 PD-rSNPs found in smNPC
TG_smNPC = SNEEP.smNPC %>%
    dplyr::select(SNP_position, geneNames) %>%
    dplyr::filter(geneNames != ".") %>%
    group_by(SNP_position, geneNames) %>%
    distinct(SNP_position, .keep_all = TRUE) %>%
    separate_rows(geneNames, sep = ",") %>%
    ungroup() %>%
    distinct(geneNames)

write_delim(TG_smNPC, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/enrichR/distinct_target_genes_smNPC.txt",
    delim = "\t", col_names = FALSE)

# 293 target genes are associated to 91 PD-rSNPs in smNPC

# TH+ neurons day 15
SNEEP.THposD15 = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/parkinsonDisease_project_Deborah_Lasse/sneep_witATAC_mergedlowC_2409/result_PD_mDAN_day15_lowC_ATAC_N1_2409.txt",
    delim = "\t", col_names = TRUE)

SNEEP.THposD15 = SNEEP.THposD15 %>%
    mutate(Sample = "THpos_neurons_D15")

# Unique rSNPs
SNEEP.THposD15 %>%
    distinct(SNP_position) %>%
    write_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/uniq_rSNPs_mDAND15.txt",
        delim = "\t", col_names = FALSE)

# How many unique target genes associated to those 180 PD-rSNPs found in mDAN
# D15
TG_mDAN.D15 = SNEEP.THposD15 %>%
    dplyr::select(SNP_position, geneNames) %>%
    dplyr::filter(geneNames != ".") %>%
    group_by(SNP_position, geneNames) %>%
    distinct(SNP_position, .keep_all = TRUE) %>%
    separate_rows(geneNames, sep = ",") %>%
    ungroup() %>%
    distinct(geneNames)

write_delim(TG_mDAN.D15, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/enrichR/distinct_target_genes_mDAN_D15.txt",
    delim = "\t", col_names = FALSE)

# 308 target genes are associated to 180 PD-rSNPs in mDAN D15

# TH+ neurons day 30
SNEEP.THposD30 = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/parkinsonDisease_project_Deborah_Lasse/sneep_witATAC_mergedlowC_2409/result_PD_mDAN_day30_lowC_ATAC_N1_2409.txt",
    delim = "\t", col_names = TRUE)

SNEEP.THposD30 = SNEEP.THposD30 %>%
    mutate(Sample = "THpos_neurons_D30")

# Unique rSNPs
SNEEP.THposD30 %>%
    distinct(SNP_position) %>%
    write_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/uniq_rSNPs_mDAND30.txt",
        delim = "\t", col_names = FALSE)

# How many unique target genes associated to those 181 PD-rSNPs found in mDAN
# D30
TG_mDAN.D30 = SNEEP.THposD30 %>%
    dplyr::select(SNP_position, geneNames) %>%
    dplyr::filter(geneNames != ".") %>%
    group_by(SNP_position, geneNames) %>%
    distinct(SNP_position, .keep_all = TRUE) %>%
    separate_rows(geneNames, sep = ",") %>%
    ungroup() %>%
    distinct(geneNames)

write_delim(TG_mDAN.D30, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/enrichR/distinct_target_genes_mDAN_D30.txt",
    delim = "\t", col_names = FALSE)

# 301 target genes are associated to 181 PD-rSNPs in mDAN D30

# TH+ neurons day 50
SNEEP.THposD50 = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/parkinsonDisease_project_Deborah_Lasse/sneep_witATAC_mergedlowC_2409/result_PD_mDAN_day50_lowC_ATAC_N1_2409.txt",
    delim = "\t", col_names = TRUE)


SNEEP.THposD50 = SNEEP.THposD50 %>%
    mutate(Sample = "THpos_neurons_D50")

# Unique rSNPs
SNEEP.THposD50 %>%
    distinct(SNP_position) %>%
    write_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/uniq_rSNPs_mDAND50.txt",
        delim = "\t", col_names = FALSE)

# How many unique target genes associated to those 160 PD-rSNPs found in mDAN
# D50
TG_mDAN.D50 = SNEEP.THposD50 %>%
    dplyr::select(SNP_position, geneNames) %>%
    dplyr::filter(geneNames != ".") %>%
    group_by(SNP_position, geneNames) %>%
    distinct(SNP_position, .keep_all = TRUE) %>%
    separate_rows(geneNames, sep = ",") %>%
    ungroup() %>%
    distinct(geneNames)

write_delim(TG_mDAN.D50, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/enrichR/distinct_target_genes_mDAN_D50.txt",
    delim = "\t", col_names = FALSE)

# 300 target genes are associated to 160 PD-rSNPs in mDAN D50

# Combine target genes from mDAN D15, D30 and D50 together
TG_distinct_mDAN = TG_mDAN.D15 %>%
    bind_rows(TG_mDAN.D30) %>%
    bind_rows(TG_mDAN.D50) %>%
    distinct(geneNames)

write_delim(TG_distinct_mDAN, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/enrichR/distinct_target_genes_mDAN_all_together.txt",
    delim = "\t", col_names = FALSE)

# Combine them all
SNEEP_all = SNEEP.smNPC %>%
    bind_rows(SNEEP.THposD15) %>%
    bind_rows(SNEEP.THposD30) %>%
    bind_rows(SNEEP.THposD50)

# And save them
write_delim(SNEEP_all, "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/SNEEP_all_res_comb.txt",
    delim = "\t", col_names = TRUE)

# Load when needed
SNEEP_all = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/SNEEP_all_res_comb.txt",
    delim = "\t", col_names = TRUE)

# How many unique rSNPs per samples
SNEEP_all %>%
    group_by(Sample) %>%
    summarise(dplyr::n())

SNEEP_all %>%
    group_by(Sample) %>%
    distinct(SNP_position) %>%
    summarise(dplyr::n())

# How many unique rSNPs?
SNEEP_all %>%
    distinct(SNP_position) %>%
    write_delim(., "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/SNEEP_uniq_rSNPs.txt",
        delim = "\t", col_names = TRUE)

# How many SNPs among the 254 PD-rSNPs are creating or disturbing TFBS
SNEEP_all %>%
    distinct(SNP_position, .keep_all = TRUE) %>%
    dplyr::select(SNP_position, log_pvalueBindAffVar1_pvalueBindAffVar2) %>%
    mutate(pos = sum(log_pvalueBindAffVar1_pvalueBindAffVar2 > 0), neg = sum(log_pvalueBindAffVar1_pvalueBindAffVar2 <
        0))

# How many SNPs among the 254 PD-rSNPs are associated to MAPT or SNCA MAPT
SNEEP_all %>%
    # group_by(SNP_position, geneNames) %>%
dplyr::select(SNP_position, TF, geneNames) %>%
    # distinct(SNP_position, .keep_all = TRUE) %>%
dplyr::filter(str_detect(geneNames, "MAPT")) %>%
    distinct(SNP_position, .keep_all = TRUE)

# 34 out 254 PD-rSNPs are associated to MAPT SNCA
SNEEP_all %>%
    # group_by(SNP_position, geneNames) %>%
dplyr::select(SNP_position, TF, geneNames) %>%
    # distinct(SNP_position, .keep_all = TRUE) %>%
dplyr::filter(str_detect(geneNames, "SNCA")) %>%
    distinct(SNP_position, .keep_all = TRUE)

# 7 out 254 PD-rSNPs are associated to SNCA

# How many chromosomes have at least one of the 254 PD-rSNPs
SNEEP_all %>%
    distinct(SNP_position) %>%
    separate(SNP_position, into = c("chr", "tmp"), sep = ":") %>%
    group_by(chr) %>%
    summarise(chr_nb = n())
```

Now that unique set of rSNPs is available, filter this number by computing correlation between the expression of the TFs and their predicted target genes

``` r
# Load RNAse data
dat.RPKM.filt = read_rds("/home/vagrant/Manuscript_1/FIGURE2/mat_RPKM.rds")

# Extract the TFBS predicted to be altered by the 254 rSNPs and plot their expression as well a and s the expression of their putative target genes
## As some TFs act as dimers, extract the expression of each and the lowest one (of each pair) will be the rate limiting TF
TF_dimer_exp = SNEEP_all %>% 
  dplyr::select(SNP_position:var2, TF, geneNames, Sample) %>% 
  group_by(SNP_position, Sample) %>% 
  dplyr::filter(geneNames != ".") %>% # Filter out rSNPs where there are no putative target genes
  mutate(TF = gsub("\\(.*", "\\1", TF)) %>% # remove the JASPAR ID next to some TFs
  dplyr::filter(str_detect(TF, "::")) %>%  # Filter out the TF dimers and check if one of the TFs is low expressed
  dplyr::select(SNP_position, TF) %>% 
  mutate(TF_dimer = TF) %>% 
  separate_rows(TF, sep = "::") %>% 
  ungroup() %>% 
  dplyr::select(-Sample)

## Plot
dat_2_plot = dat.RPKM.filt %>% 
  dplyr::filter(gene_name %in% TF_dimer_exp$TF) %>% 
  dplyr::rename(TF = gene_name) %>% 
  left_join(TF_dimer_exp,
            by = "TF") %>% 
  mutate(TF = as.factor(TF))

pdf("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/TF_dimer_Ind_exp.pdf",
    width = 15,
    height = 20)
ggboxplot(dat_2_plot,
          x = "Cond", 
          y = "RPKM", 
          add = "jitter", 
          color = "TF", 
          #palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"),
          facet.by = "TF_dimer",
          xlab = "", 
          ylab = "Reads Per Kilobase Million (RPKM)") +
  geom_hline(yintercept = 5,
             linetype = "dashed") +  # 5 RPKM expression threshold
  scale_x_discrete(labels = c("smNPC",
                              expression("mDAN D15"),
                              expression("mDAN D30"),
                              expression("mDAN D50"))) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90))

dev.off()

## Compute Pearson correlation to filter rSNPs
# Tibble for TF expression
TF_exp = SNEEP_all %>%
  dplyr::select(SNP_position, TF, geneNames) %>%
  dplyr::filter(geneNames != ".") %>% # no target genes idenditifed
  mutate(TF = gsub("\\(.*", "\\1", TF)) %>%
  separate_rows(TF, sep = "::") %>% # Separate TF dimers into single TF
  separate_rows(geneNames, sep = ",") %>% # Separate target genes into single target genes
  unite(col = "Item", c("SNP_position", "TF", "geneNames"), remove = FALSE) %>%
  dplyr::select(Item, TF) %>%
  dplyr::rename(gene_name = TF) %>%
  dplyr::filter(gene_name %in% dat.RPKM.filt$gene_name) %>%
  left_join(dat.RPKM.filt,
            by = "gene_name") %>%
  group_by(gene_name) %>% 
  dplyr::filter(any(RPKM > 1)) %>% # Keep only TF that have RPKM > 1 in all samples
  dplyr::select(-gene_id, -Cond) %>%
  dplyr::rename(TF = gene_name,
                RPKM.TF = RPKM)

# Tibble for TG expression  
TG_exp = SNEEP_all %>%
  dplyr::select(SNP_position, TF, geneNames) %>%
  dplyr::filter(geneNames != ".") %>% # no target genes idenditifed
  mutate(TF = gsub("\\(.*", "\\1", TF)) %>%
  separate_rows(TF, sep = "::") %>% # Separate TF dimers into single TF
  separate_rows(geneNames, sep = ",") %>% # Separate target genes into single target genes
  unite(col = "Item", c("SNP_position", "TF", "geneNames"), remove = FALSE) %>%
  dplyr::select(Item, geneNames) %>%
  separate_rows(geneNames, sep = ",") %>%
  dplyr::rename(gene_name = geneNames) %>%
  dplyr::filter(gene_name %in% dat.RPKM.filt$gene_name) %>%
  left_join(dat.RPKM.filt,
            by = "gene_name") %>%
  group_by(gene_name) %>% 
  dplyr::filter(any(RPKM > 1)) %>% # Keep only target genes that have RPKM > 1 in all samples
  dplyr::select(-gene_id, -Cond) %>%
  dplyr::rename(target_gene = gene_name,
                RPKM.TG = RPKM)

# Join both target gene and TF tables. Filter out TF or target genes with expression below 1 RPKM
TF_TG_mat = merge(TF_exp, TG_exp) %>% 
  as_tibble() %>% 
  group_by(TF, target_gene) %>% 
  distinct(Sample, .keep_all = TRUE)

# Pearson correlation
TF_TG_mat_Pearson_cor = cor_test(TF_TG_mat %>% 
                                   group_by(TF, target_gene),
                                 vars = c("RPKM.TF","RPKM.TG"),
                                 use = "everything")

write_delim(TF_TG_mat_Pearson_cor,
            "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/TF_TG_exp_pearsonCor.txt",
            delim = "\t",
            col_names = TRUE)

# Filter the correlation matrix to only keep TF-target gene pairs with a strong correlation/anti-correlation
TF_TG_mat_Pearson_cor.filt = TF_TG_mat_Pearson_cor %>% 
  dplyr::filter(cor > 0.5 | cor < -0.5)

unique(TF_TG_mat_Pearson_cor.filt$target_gene)
unique(TF_TG_mat_Pearson_cor.filt$TF)

write_delim(TF_TG_mat_Pearson_cor.filt,
            "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/TF_TG_exp_pearsonCor_abs0.5.txt",
            delim = "\t",
            col_names = TRUE)

# Reassociate the rSNPs to which the TF-target genes pairs belongs and make a table
tbl.TF.TG = TF_TG_mat_Pearson_cor.filt %>% 
  unite(col = "TF_TG_pair", TF, target_gene, sep = "_", remove = FALSE)

tbl.TF.TG.FINAL = TF_TG_mat %>% 
  ungroup() %>% 
  dplyr::filter(str_detect(Item, paste0(tbl.TF.TG$TF_TG_pair, collapse = "|"))) %>% 
  distinct(Item, .keep_all = TRUE) %>% 
  dplyr::select(Item, TF, target_gene) %>% 
  mutate(rSNPs = gsub("\\_.*", "\\1", Item)) %>% 
  dplyr::select(rSNPs, TF, target_gene)

# Select the coordinates of the rSNPs ans associate back rsID from UCSC using all snps from db155
tbl.TF.TG.FINAL %>% 
  dplyr::select(rSNPs) %>% 
  distinct(rSNPs) %>% 
  write_delim(.,
              "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/UCSC_coordinates_rSNPs.txt",
              delim = "\t",
              col_names = FALSE)

# Load the newly created file
rsID_PD_rSNP = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/UCSC_coordinates_rSNPs_rsID.txt",
              delim = ",",
              col_names = TRUE)

# Load the newly created file
rsID_PD_rSNP.filt = rsID_PD_rSNP %>% 
  dplyr::filter(class == "snv") %>% 
  unite(col = "tmp", c(`#"chrom"`, chromStart), sep = ":", remove = TRUE) %>% 
  unite(col = "rSNPs", c("tmp", "chromEnd"), sep = "-", remove = FALSE) %>% 
  dplyr::select(rSNPs, name, ref, alts) %>% 
  distinct(name, .keep_all = TRUE) # Remove doublets

# Filter the final table
tbl.TF.TG.FINAL.rsID = tbl.TF.TG.FINAL %>% 
  #dplyr::select(rSNPs) %>% 
  #distinct(rSNPs) %>% 
  dplyr::filter(rSNPs %in% rsID_PD_rSNP.filt$rSNPs) %>% 
  left_join(rsID_PD_rSNP.filt,
            by = "rSNPs") %>% 
  dplyr::select(rSNPs, name:alts, TF, target_gene) %>% 
  separate(rSNPs, into = c("chr", "tmp"), sep = ":", remove = TRUE) %>% 
  separate(tmp, into = c("start", "end"), sep = "-", remove = TRUE) %>% 
  dplyr::rename(rsID = name)

# Extract the unique rsIDs and get the MAF (minor allele frequencies)
snp2MAF = tbl.TF.TG.FINAL.rsID %>% 
  dplyr::select(rsID) %>% 
  distinct(rsID) %>% 
  pull()

# Get allele frequencies
LDproxy_batch(snp = snp2MAF,
              pop = "EUR",
              r2d = "r2",
              token = "5bba40af90d7",
              append = TRUE,
              genome_build = "grch38_high_coverage")

# The data has been automaticcaly saved -> load it
MAF.snp = read_delim("/home/vagrant/Manuscript_1/TABLES/combined_query_snp_list_grch38_high_coverage.txt",
                     delim = "\t",
                     col_names = TRUE)

# The linkage desiquilibrium linkage value has been also calculated. Filter only the snp that are in linkage desiquilibrium with themselves
MAF.snp %>% 
  dplyr::select(RS_Number:Distance) %>% 
  dplyr::rename(query_snp = RS_Number,
                rsID = Coord,
                Coord = Alleles,
                Allele = MAF,
                MAF = Distance) %>% 
  mutate(same.snp = if_else(query_snp == rsID, 1, 0)) %>% 
  # dplyr::filter(same.snp == 1,
  #               rsID %in% snp2MAF) %>% 
  dplyr::filter(same.snp == 1) %>% 
  #left_join(tbl.TF.TG.FINAL.rsID, by = "rsID") %>% 
  full_join(tbl.TF.TG.FINAL.rsID, by = "rsID") %>% 
  #distinct(rsID, .keep_all = TRUE) %>% 
  write_delim(.,
            "/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/final_tbl_PD_rSNPs_with_MAF.txt",
            delim = "\t",
            col_names = TRUE)
```

Create the page layout that will contain all the necessary plots


``` r
####################
#### Figure 2  ####
###################
# Save as TIFF, 300 ppi
tiff("/home/vagrant/Manuscript_1/FIGURE2/FIGURE2.tiff",
     width = 8.27,
     height = 11.67,
     units = "in",
     res = 300,
     compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, 
           height = 11.67, 
           default.units = "inches",
           showGuides = TRUE)

# text Figure 2
plotText(label = "Figure 2", 
         fontsize = 14,
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 0.25, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

#### PANEL A - SNEEP PIPELINE ####
# text A
plotText(label = "A", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 0.5, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

# Load the picture prepared in Biorender
Fig2A_pip = readPNG("/home/vagrant/Manuscript_1/FIGURE2/Fig2A_FlowChart.png")

# Plot Figure 2A - Pipeline
plotRaster(image = Fig2A_pip,
           x = 0.5,
           y = 0.625,
           width = 3.5,
           height = 3.0,
           just = c("left", 
                    "top"),
           interpolate = TRUE)

#### PANEL B - Pizza plot of 50 PD-rSNPs (among 54) that shows allelic imbalance #### 
# text B
plotText(label = "B", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 4.125, 
         y = 0.5, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

# Since there are 254 unique rSNPs but we do not have about them being heterozygous or not, check if there is also chromatin allelic imbalance
# hg38
# Load the 254 PD rSNPs coordinates
PD_rSNPs_254 = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/SNEEP_uniq_rSNPs.txt",
                          delim = "\t",
                          col_names = TRUE)

PD_rSNPs_254_fmt = PD_rSNPs_254 %>% 
  separate(SNP_position, into = c("chr", "tmp"), sep = ":", remove = FALSE) %>% 
  separate(tmp, into = c("start", "end"), sep = "-", remove = TRUE) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))
  
# Make a Granges object of the 254 PD rSNPs (width has to be one)
searchArea = GRanges(seqnames = PD_rSNPs_254_fmt$chr, 
                     ranges = IRanges(start = PD_rSNPs_254_fmt$end, 
                                      end = PD_rSNPs_254_fmt$end),
                     name = PD_rSNPs_254_fmt$SNP_position)

# Path where the ATACseq BAM files are
pathToFiles = "/home/vagrant/Manuscript_1/FIGURE4/"

# Import a specific region from the BAM files (based on the search area above)
reads = impBamGAL(pathToFiles, 
                  searchArea, 
                  verbose = TRUE)

# Count the number of reads for the 2 different alleles at that position
countList = getAlleleCounts(reads, 
                           searchArea, 
                            verbose = TRUE)

# Save as a rds file for loading later if necessary
saveRDS(countList, 
        file = "/home/vagrant/Manuscript_1/FIGURE4/PD_rSNPs_254_het_counts_THREP1_ATACseq.rds")

# Load the new results of the 254 rSNPs
countList = read_rds("/home/vagrant/Manuscript_1/FIGURE4/PD_rSNPs_254_het_counts_THREP1_ATACseq.rds")

# Make an ASEset object to be able to plot
countList.ASE = ASEsetFromCountList(searchArea, countList)

# And plot one graph per pdf file
for (i in 1:length(countList.ASE)){
  pdf(paste0("/home/vagrant/Manuscript_1/FIGURE4/20241111_AlellicImbalance_254rSNPs_", 
             rownames(countList.ASE[i]),
             ".pdf"))
  barplot(countList.ASE[i])
  dev.off()
}

# Now check how many are truly heterozygous thanks to the whole genome sequencing data (WGS)
# Use table browser to retrieve according to this link http://genome.ucsc.edu/FAQ/FAQreleases.html#snpConversion
searchArea %>% 
  as_tibble() %>% 
  dplyr::select(-start) %>% # Not the real start when converting a GRanges to a tibble 
  mutate(start = end - 1) %>% 
  dplyr::select(seqnames, start, end, name) %>% 
  dplyr::rename(rSNPs = name) %>% 
  inner_join(PD_rSNPs_254_rsID.filt,
             by = "rSNPs") %>% 
  dplyr::select(name) %>% 
  write_delim(.,
              "/home/vagrant/Manuscript_1/FIGURE4/rsIDs_254PDrSNPs.txt",
              delim = "\t",
              col_names = FALSE)

# Load now the coordinates that are in hg19
PD_rSNPs_254_hg19 = read_delim("/home/vagrant/Manuscript_1/FIGURE4/rsIDs_254PDrSNPs_hg19_coord.txt",
                               delim = "\t",
                               col_names = TRUE)

# Only keep SNPs that are on canonical chromosomes
PD_rSNPs_254_hg19.GR = PD_rSNPs_254_hg19 %>% 
  dplyr::filter(!str_detect(`#chrom`, "_")) %>% 
  dplyr::select(`#chrom`:name) %>% 
  dplyr::rename(chr = `#chrom`) %>% 
  makeGRangesFromDataFrame(.,
                           keep.extra.columns = TRUE,
                           seqnames.field = "chr",
                           start.field = "chromStart",
                           end.field = "chromEnd",
                           starts.in.df.are.0based = TRUE)

# Load WGS data of TH REP1 mCHERRY cell line (hg19 genome version)
dat.SNPs.only.GR = read_rds("/home/vagrant/epifunc/SNPs_overlap/E19D017a78.merge.qc.SNP_INDEL.hg19.annotated.rds")

dat.SNPs.only.GR = dat.SNPs.only.GR %>% 
  dplyr::select(CHR:POS, REF:FILTER, E19D017a78) %>% 
  mutate(START = POS - 1) %>% 
  makeGRangesFromDataFrame(.,
                           keep.extra.columns = TRUE,
                           seqnames.field = "CHR",
                           start.field = "START",
                           end.field = "POS",
                           starts.in.df.are.0based = TRUE)

seqlevelsStyle(dat.SNPs.only.GR) = "UCSC"
seqlevelsStyle(dat.SNPs.only.GR)

# Overlap
hits = findOverlaps(PD_rSNPs_254_hg19.GR, dat.SNPs.only.GR)
rSNP.254_in_E19D017a78 = PD_rSNPs_254_hg19.GR[queryHits(hits)]
E19D017a78_w_rSNP.254 = dat.SNPs.only.GR[subjectHits(hits)]

mcols(rSNP.254_in_E19D017a78) = cbind(mcols(rSNP.254_in_E19D017a78),
                                 mcols(E19D017a78_w_rSNP.254))

rSNP.254_in_E19D017a78 = as.data.frame(rSNP.254_in_E19D017a78) %>% 
  as_tibble()

# Out of 254 PD-rSNPs, we have information about 242

# Save
write_delim(rSNP.254_in_E19D017a78,
            "/home/vagrant/Manuscript_1/FIGURE4/Genotype_of_254PDrSNPs_in_THREP1mCHERRY.txt",
            delim = "\t",
            col_names = TRUE)

# Load later
rSNP.254_in_E19D017a78 = read_delim("/home/vagrant/Manuscript_1/FIGURE4/Genotype_of_254PDrSNPs_in_THREP1mCHERRY.txt",
                                    delim = "\t",
                                    col_names = TRUE)

# Check how many are "unknown" (./.:.)
rSNP.254_in_E19D017a78 %>% 
  mutate(Genotype = gsub("\\:.*", "\\1", E19D017a78)) %>% 
  group_by(Genotype) %>% 
  summarise(GT.num = n()) # 54 PD-rSNPs are heterozygous in TH REP1 mCHERRY cell line, 24 are homozygous for the alternative allele and 164 are                                  unknown

# Check if the 54 PD-rSNPs that are heterozygous in TH REP1 mCHERRY cell line also lead to chromatin allelic imbalance
# rsIDs of the 54 heterozygous PD-rSNPs 
rSNP.254_in_E19D017a78_0.1 = rSNP.254_in_E19D017a78 %>% 
  mutate(Genotype = gsub("\\:.*", "\\1", E19D017a78),
         Genotype2 = gsub("\\|", "\\/", Genotype)) %>% 
  dplyr::filter(str_detect(Genotype2, "0/1")) %>% 
  #unite(col = "toInves", seqnames, end, sep = "_") %>% 
  dplyr::select(name)

# Get their coordinates
hetero_54_rSNP_in_E19D017a78.hg38 = searchArea %>% 
  as_tibble() %>% 
  dplyr::select(-start) %>% # Not the real start when converting a GRanges to a tibble 
  mutate(start = end - 1) %>% 
  dplyr::select(seqnames, start, end, name) %>% 
  dplyr::rename(rSNPs = name) %>% 
  inner_join(PD_rSNPs_254_rsID.filt,
             by = "rSNPs") %>% 
  dplyr::filter(name %in% rSNP.254_in_E19D017a78_0.1$name) %>% 
  unite(col = "coordToInves", seqnames, end, sep = "_")

# Extract the 54 heterozygous PD-rSNPs from thre ASet object
PDrSNPs_54_het = countList[intersect(names(countList), hetero_54_rSNP_in_E19D017a78.hg38$coordToInves)]

searchArea.0.1 = searchArea %>% 
  as_tibble() %>% 
  dplyr::filter(name %in% hetero_54_rSNP_in_E19D017a78.hg38$rSNPs)

searchArea.0.1 = searchArea.0.1 %>% 
  makeGRangesFromDataFrame(.,
                           keep.extra.columns = TRUE,
                           seqnames.field = "seqnames",
                           start.field = "start",
                           end.field = "end",
                           starts.in.df.are.0based = FALSE)

# And make it a ASet object back and proceed with a chi-squared test
PDrSNPs_54_het.ASE = ASEsetFromCountList(searchArea.0.1, PDrSNPs_54_het)
ref(PDrSNPs_54_het.ASE) = c("G", "T", "C")
PDrSNPs_54_het.ASE.chisq = as_tibble(chisq.test(PDrSNPs_54_het.ASE[,1:7], "+")) %>% 
  mutate(Sample = colnames(PDrSNPs_54_het.ASE)) %>% 
  dplyr::select(Sample, everything())

# Consider that a PD-rSNP lead to allelic imbalance if at least one sample (smNPC, mDAN-D15, mDAN-D30, mDAN-D50) shows allelic imbalance
PDrSNPs_54_het.ASE.chisq %>% 
  dplyr::filter(!str_detect(Sample, paste0(c("astrocyte", "neg"), collapse = "|"))) %>% 
  #pivot_longer(!Sample, names_to = "PDrSNPs", values_to = "pvalue") %>% 
  #group_by(Sample) %>% 
  #dplyr::filter(if_any(where(is.numeric), ~ . < 0.05))
  select_if(~any(. < 0.05)) %>% 
  mutate(AI_yes = ncol(.)) %>% 
  View()
# 46 among 54 PD-rSNPs show allelic imbalance

test = tibble(AI = c("yes", "no"),
       value = c(46, 8),
       percent = round((value/sum(value))*100, digits = 1))

pizza.AI = ggplot(test,
       aes(x = "", y = value, fill = AI)) +
  geom_bar(width = 1, stat = "identity", alpha = 0.5) +
    #theme_void() +
    scale_fill_manual(values = c("#4682BD", "#64A7D1")) +
    coord_polar(theta = "y", start = 2.125, clip = "off") +
    geom_text(aes(
        x = c(1.0, 1.2),
        y = c(21.5, 50),
        label = paste0(percent, "%")),
        size = 5.0, 
        color = "black",
        fontface = "bold") +
  theme_void() +
  theme(legend.position = "none")

plotGG(plot = pizza.AI, 
       x = 3.0, 
       y = 0.5, 
       width = 6.0, 
       height = 1.75,
       just = c("left", "top"), 
       default.units = "inches")
  
# And plot one graph per pdf file
for (i in 1:length(PDrSNPs_54_het.ASE)){
  pdf(paste0("/home/vagrant/Manuscript_1/FIGURE4/20241115_AlellicImbalance_54HeterorSNPs_",
             rownames(PDrSNPs_54_het.ASE[i]),
             ".pdf"))
  barplot(PDrSNPs_54_het.ASE[i])
  dev.off()
}

#### PANEL C - Manhattan plot #### 
# text C
plotText(label = "C", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 4.125, 
         y = 2.125, 
         just = "left", 
         default.units = "inches",
         fontface = "bold") 

# Load the Nalls et al. GWAS data
nalls_allSNPs.HG38.GR = readRDS("/home/vagrant/Manuscript_1/FIGURE1/nalls_allSNPs.HG38.GR.rds")

# Load the rsIDS FOR THE 254 significant SNPs
rsID.254 = read_delim("/home/vagrant/epifunc/SNEEP2024/RESULTS_SEPTEMBER_2024/output/SNEEP_uniq_rSNPs.txt",
                     delim = "\t",
                     col_names = TRUE)

# Extract the chr and the position
rsID.254_2 = rsID.254 %>% 
  separate(col = "SNP_position", 
           into = c("chr", "pos"), 
           sep = ":",
           remove = FALSE) %>% 
  separate(col = "pos",
           into = c("start", "end"),
           sep = "-",
           remove = TRUE) %>% 
  unite(col = pos, 
        c("chr", "end"), 
        sep = "_", 
        remove = TRUE) %>% 
  dplyr::select(-start)

# Select the SNPs position from nalls et al and associate rsIDs
nalls.manH = nalls_allSNPs.HG38.GR %>% 
  as_tibble() %>% 
  dplyr::select(seqnames, base_pair_position, p) %>% 
  unite(col = "pos", seqnames, base_pair_position, sep = "_", remove = FALSE) %>%
  left_join(rsID.254_2, by = "pos") %>% 
  separate(pos, into = c("chrom", "pos"), sep = "_", remove = TRUE) %>% 
  dplyr::select(chrom, pos, p, SNP_position)

# Save
write_rds(nalls.manH,
          "/home/vagrant/Manuscript_1/FIGURE1/20241107_254_manH.rds")

# Load for later 
nalls.manH = read_rds("/home/vagrant/Manuscript_1/FIGURE1/20241107_254_manH.rds")

# Rename columns
nalls.manH = nalls.manH %>% 
  dplyr::rename(snp = SNP_position) %>% 
  mutate(pos = as.numeric(pos))

# Manhattan plot - highlight the 254 SNPs with rsID
# Get the rsID of the 254 rSNPs first
nalls.manH %>% 
  dplyr::filter(!is.na(snp),
                p <= 5e-08) %>% 
  dplyr::select(snp) %>% 
  write_delim(.,
              "/home/vagrant/Manuscript_1/FIGURE1/20241107_254PDrSNPs_coord.txt",
              col_names = FALSE,
              delim = "\t")

# Load the files with the associated rsIDs
PD_rSNPs_254_rsID = read_delim("/home/vagrant/Manuscript_1/FIGURE1/20241107_254PDrSNPs_rsID.txt",
                               delim = "\t",
                               col_names = TRUE)

# Only keep snv
PD_rSNPs_254_rsID.filt = PD_rSNPs_254_rsID %>% 
  filter(class == "snv") %>% 
  dplyr::select(`#chrom`:name) %>% 
  distinct(name, .keep_all = TRUE) %>% 
  unite(col = "tmp", `#chrom`, chromStart, sep = ":", remove = FALSE) %>% 
  unite(col = "rSNPs", tmp, chromEnd, sep = "-", remove = TRUE)

# Associated then the significant GWAS pvalue
snp_2_high = nalls.manH %>% 
  dplyr::filter(!is.na(snp),
                p <= 5e-08) %>% 
  dplyr::rename(rSNPs = snp) %>% 
  left_join(PD_rSNPs_254_rsID.filt,
            by = "rSNPs") %>% 
  dplyr::select(-rSNPs, -`#chrom`, -chromStart) %>% 
  dplyr::rename(snp = name)
  
# Get the full set of SNPs and associate rsID to the one I want to highlight
# snp_to_pl = nalls.manH %>% 
#   mutate(snp = NA_character_) %>% 
#   bind_rows(snp_2_high)

snp_to_pl = nalls.manH %>%
  mutate(snp = NA_character_) %>%
  bind_rows(snp_2_high) %>%
  mutate(pbis = -log10(p),
         chrom = factor(chrom,
                        levels = paste0("chr",
                                        1:22)),
         color = if_else(is.na(snp), "no_high", "high"),
         
         label = if_else(snp %in% c("rs1465922", "rs144814361"), snp, NA_character_))

highlight_colormap = c("no_high" = adjustcolor( "grey", 
                                                alpha.f = 0.2), 
                       "high" = "#6600FF")

no_high = manhattan_data_preprocess(snp_to_pl,
                                    pval.colname = "p",
                                    chr.colname = "chrom",
                                    pos.colname = "pos",
                                    highlight.colname = "color",
                                    highlight.col = highlight_colormap,
                                    signif = c(5e-08))
# Plot
snp.pl = manhattan_plot(x = no_high,
                        color.by.highlight = TRUE,
                        rescale = TRUE,
                        label.colname = "label") +
  theme(axis.text.x = element_text(size = 9,
                                 face = "bold",
                                 angle = 90),
        axis.text.y = element_text(size = 9,
                                   face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11,
                                    face = "bold"))
 
# saveRDS(snp.pl,
#        "/home/vagrant/Manuscript_1/FIGURE2/Fig2C_ManhattanPlot.rds")
snp.pl = readRDS("/home/vagrant/Manuscript_1/FIGURE2/Fig2C_ManhattanPlot.rds")

# Place the plot
plotGG(plot = snp.pl, 
       x = 4.5, 
       y = 2.25, 
       width = 3.75, 
       height = 1.5,
       just = c("left", "top"), 
       default.units = "inches")

#### PANEL D - BAG3 locus #### 
# text D
plotText(label = "D", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 4.0, 
         just = "left", 
         default.units = "inches",
         fontface = "bold") 

# ATACseq signal BAG3
# path where bigwig files are
bw.path = "/home/vagrant/epifunc/"

# 1st bio replicate smNPC
bw.smNPC_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_I.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 2nd bio replicate smNPC
bw.smNPC_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_II.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 3rd bio replicate smNPC
bw.smNPC_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_III.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 1st bio replicate TH positive neurons D15
bw.TH.pos_D15_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_I.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 2nd bio replicate TH positive neurons D15
bw.TH.pos_D15_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_II.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 3rd bio replicate TH positive neurons D15
bw.TH.pos_D15_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_III.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 1st bio replicate TH positive neurons D30
bw.TH.pos_D30_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_I.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 2nd bio replicate TH positive neurons D30
bw.TH.pos_D30_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_II.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 3rd bio replicate TH positive neurons D30
bw.TH.pos_D30_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_III.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 1st bio replicate TH positive neurons D50
bw.TH.pos_D50_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 2nd bio replicate TH positive neurons D50
bw.TH.pos_D50_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# 3rd bio replicate TH positive neurons D50
bw.TH.pos_D50_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
                        chrom = "chr10",
                        chromstart = 119611305,
                        chromend = 119691505)

# Add a scale next to the bigwig files - check the maximum to choose
scale.BAG3.max = max(c(bw.smNPC_1$score,
                       bw.smNPC_2$score,
                       bw.smNPC_3$score,
                       bw.TH.pos_D15_1$score,
                       bw.TH.pos_D15_2$score,
                       bw.TH.pos_D15_3$score,
                       bw.TH.pos_D30_1$score,
                       bw.TH.pos_D30_2$score,
                       bw.TH.pos_D30_3$score,
                       bw.TH.pos_D50_1$score,
                       bw.TH.pos_D50_2$score,
                       bw.TH.pos_D50_3$score))

print(scale.BAG3.max)

# Define parameters for the regions
region.p.BAG3 = pgParams(chrom = "chr10",
                         chromstart = 119611405,
                         chromend = 119691405,
                         assembly = "hg38",
                         range = c(0, scale.BAG3.max))

# Add the genomic label
plotGenomeLabel(
  chrom = "chr10",
  chromstart = 119611405,
  chromend = 119691405,
  assembly = "hg38",
  x = 0.5,
  y = 4.125,
  length = 4.0,
  default.units = "inches",
  scale = "Mb")

# ATACseq signal smNPC_I and scale
ATAC_smNPC_I = plotSignal(
  data = bw.smNPC_1, 
  params = region.p.BAG3,
  fill = "#D3D3D3", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.375, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal smNPC_II
ATAC_smNPC_II = plotSignal(
  data = bw.smNPC_2, 
  params = region.p.BAG3,
  fill = "#C0C0C0", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.375, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal smNPC_III
ATAC_smNPC_III = plotSignal(
  data = bw.smNPC_3, 
  params = region.p.BAG3,
  fill = "#A9A9A9", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.375, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III, 
          at = c(0, scale.BAG3.max),
          fontsize = 6)

# Add the samples name
plotText(label = "smNPC", 
         fontsize = 6, 
         fontcolor = "#A9A9A9",
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 4.5, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D15 I
ATAC_TH.pos_D15_I = plotSignal(
  data = bw.TH.pos_D15_1, 
  params = region.p.BAG3,
  fill = "#B22222", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.675, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D15 II
ATAC_TH.pos_D15_II = plotSignal(
  data = bw.TH.pos_D15_2, 
  params = region.p.BAG3,
  fill = "#A52A2A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.675, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D15 III
ATAC_TH.pos_D15_III = plotSignal(
  data = bw.TH.pos_D15_3, 
  params = region.p.BAG3,
  fill = "#8B0000", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.675, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_III, 
          at = c(0, scale.BAG3.max),
          fontsize = 6)

# Annotate sample name
plotText(label = "mDAN D15", 
         fontsize = 6, 
         fontfamily = "Helvetica",
         fontcolor = "#B22222",
         x = 0.25, 
         y = 4.75, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D30 I
ATAC_TH.pos_D30_I = plotSignal(
  data = bw.TH.pos_D30_1, 
  params = region.p.BAG3,
  fill = "#FF6347", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.975, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D30 II
ATAC_TH.pos_D30_II = plotSignal(
  data = bw.TH.pos_D30_2, 
  params = region.p.BAG3,
  fill = "#FF0000", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.975, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D30 III
ATAC_TH.pos_D30_III = plotSignal(
  data = bw.TH.pos_D30_3, 
  params = region.p.BAG3,
  fill = "#DC143C", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 4.975,
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_III, 
          at = c(0, scale.BAG3.max),
          fontsize = 6)

# Annotate sample name
plotText(label = "mDAN D30", 
         fontsize = 6, 
         fontfamily = "Helvetica",
         fontcolor = "#DC143C",
         x = 0.25, 
         y = 5.0625, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D50 I
ATAC_TH.pos_D50_I = plotSignal(
  data = bw.TH.pos_D50_1, 
  params = region.p.BAG3,
  fill = "#FFA07A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 5.275, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D50 II
ATAC_TH.pos_D50_II = plotSignal(
  data = bw.TH.pos_D50_2, 
  params = region.p.BAG3,
  fill = "#FA8072", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 5.275, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D50 III
ATAC_TH.pos_D50_III = plotSignal(
  data = bw.TH.pos_D50_3, 
  params = region.p.BAG3,
  fill = "#E9967A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 5.275,
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_III, 
          at = c(0, scale.BAG3.max),
          fontsize = 6)

# Annotate sample name
plotText(label = "mDAN D50", 
         fontsize = 6, 
         fontfamily = "Helvetica",
         fontcolor = "#E9967A",
         x = 0.25, 
         y = 5.375, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# Gene track
plotGenes(chrom = "chr10",
  chromstart = 119611405,
  chromend = 119691405,
  assembly = assembly(Genome = "hg38refGene", 
                      TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", 
                      OrgDb = "org.Hs.eg.db"),
  fontcolor = "black",
  fill = "black",
  x = 0.5,
  y = 5.625,
  width = 4.0,
  height = 0.375,
  just = c("left", "top"),
  default.units = "inches")

#### PANEL E - BAG3 and LHX1 expression (RNAseq)
# text E
plotText(label = "E", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 4.75, 
         y = 4.0, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

# Load RPKM matrix obtained from Borja and rename columns
# dat.RPKM = read_delim("/home/vagrant/epifunc/RNAseq/RPKM_genename",
#                       delim = "\t", 
#                       col_names = TRUE)
# 
# dat.RPKM = dat.RPKM %>% 
#   dplyr::select(gene_id, gene_name, everything())

# RPKM data - Remove the negative sorted neurons and astrocytes
# dat.RPKM.filt = dat.RPKM %>% 
#   dplyr::select(!contains("negsort"), -starts_with("ASTRO")) %>% 
#   pivot_longer(cols = ends_with("bam"), names_to = "Sample", values_to = "RPKM") %>% 
#   mutate(Cond = case_when(str_detect(Sample, "D15_possort") ~ "Positively sorted neurons D15",
#                          str_detect(Sample, "D30_possort") ~ "Positively sorted neurons D30",
#                          str_detect(Sample, "D50_possort") ~ "Positively sorted neurons D50",
#                          TRUE ~ "smNPCs"),
#          Cond = factor(Cond, levels = c("smNPCs", "Positively sorted neurons D15",
#                                   "Positively sorted neurons D30", 
#                                   "Positively sorted neurons D50")),
#          gene_id = gsub("\\..*", "", gene_id))

# Save the matrix of expression
# write_rds(dat.RPKM.filt,
#           "/home/vagrant/Manuscript_1/FIGURE2/mat_RPKM.rds")

# Load it for making the plots
dat.RPKM.filt = read_rds("/home/vagrant/Manuscript_1/FIGURE2/mat_RPKM.rds")

# Expression for BAG3 and LHX1
BAG3.LHX1.exp = ggboxplot(dat.RPKM.filt %>% 
                            dplyr::filter(gene_name %in% c("BAG3", "LHX1")),
                          x = "Cond", 
                          y = "RPKM", 
                          add = "jitter", 
                          fill = "Cond", 
                          palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"),
                          facet.by = "gene_name",
                          xlab = "", 
                          ylab = "Reads Per Kilobase Million (RPKM)") +
  scale_x_discrete(labels = c("smNPC",
                              expression("mDAN D15"),
                              expression("mDAN D30"),
                              expression("mDAN D50"))) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
  

# Place the plot
plotGG(plot = BAG3.LHX1.exp, 
       x = 5.0, 
       y = 4.125, 
       width = 3, 
       height = 3.5,
       just = c("left", "top"), 
       default.units = "inches")

# Add the SNP position that zoom in the TFBS
SNP_param = pgParams(chrom = "chr10",
                     chromstart = 119651405,
                     chromend = 119651405,
                     assembly = "hg38")

annoHighlight(plot = ATAC_smNPC_I, 
              params = SNP_param,
              fill = "#404788FF",
              y = 4.5, 
              height = 1.5, 
              just = c("left", "top"), 
              default.units = "inches")

annoZoomLines(plot = ATAC_smNPC_I, 
              params = SNP_param,
              y0 = 6.0, 
              x1 = c(2.25, 2.75), 
              y1 = 6.125, 
              default.units = "inches")

# Add the transcription factor binding motif (TFBS) for LHX1
pfm.LHX1 = getMatrixByID(JASPAR2020, 
                         ID = "MA1518.1")
LHX1 = new.env()

LHX1$LHX1 = pfm.LHX1@profileMatrix

LHX1 = as.list(LHX1)

# Remove the nucleotide position
LHX1_TFBS = ggseqlogo(LHX1) +
  annotate("rect", 
           xmin = 4.5, 
           xmax = 5.5, 
           ymin = -0.05, 
           ymax = 2.5, 
           alpha = .1, 
           col = "black", 
           fill = "#404788FF") +
  theme(axis.text.x = element_blank())

# Place the plot
plotGG(plot = LHX1_TFBS, 
       x = 1.4925, 
       y = 5.9375, 
       width = 1.5, 
       height = 1.5,
       just = c("left", "top"), 
       default.units = "inches")

# text F
plotText(label = "F", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 7.75, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

# ATACseq signal SCARB2
# 1st bio replicate smNPC
bw.SCARB2.smNPC_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_I.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 2nd bio replicate smNPC
bw.SCARB2.smNPC_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_II.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 3rd bio replicate smNPC
bw.SCARB2.smNPC_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/smNPC_III.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 1st bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_I.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 2nd bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_II.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 3rd bio replicate TH positive neurons D15
bw.SCARB2.TH.pos_D15_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D15_POS_III.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 1st bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_I.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 2nd bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_II.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 3rd bio replicate TH positive neurons D30
bw.SCARB2.TH.pos_D30_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/D30_POS_III.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 1st bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_1 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 2nd bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_2 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# 3rd bio replicate TH positive neurons D50
bw.SCARB2.TH.pos_D50_3 = readBigwig(file = paste0(bw.path, 
                                      "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
                        chrom = "chr4",
                        chromstart = 76173617,
                        chromend = 76253817)

# Add a scale next to the bigwig files - check the maximum to choose
scale.SCARB2.max = max(c(bw.SCARB2.smNPC_1$score,
                       bw.SCARB2.smNPC_2$score,
                       bw.SCARB2.smNPC_3$score,
                       bw.SCARB2.TH.pos_D15_1$score,
                       bw.SCARB2.TH.pos_D15_2$score,
                       bw.SCARB2.TH.pos_D15_3$score,
                       bw.SCARB2.TH.pos_D30_1$score,
                       bw.SCARB2.TH.pos_D30_2$score,
                       bw.SCARB2.TH.pos_D30_3$score,
                       bw.SCARB2.TH.pos_D50_1$score,
                       bw.SCARB2.TH.pos_D50_2$score,
                       bw.SCARB2.TH.pos_D50_3$score))

print(scale.SCARB2.max)

# Define parameters for the regions
region.p.SCARB2 = pgParams(chrom = "chr4",
                         chromstart = 76173717,
                         chromend = 76253717,
                         assembly = "hg38",
                         range = c(0, scale.SCARB2.max))

# Add the genomic label
plotGenomeLabel(
  chrom = "chr4",
  chromstart = 76173717,
  chromend = 76253717,
  assembly = "hg38",
  x = 0.5,
  y = 7.875,
  length = 4.0,
  default.units = "inches",
  scale = "Mb")

# ATACseq signal smNPC_I
ATAC_smNPC_I_SCARB2 = plotSignal(
  data = bw.SCARB2.smNPC_1, 
  params = region.p.SCARB2,
  fill = "#D3D3D3", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.125, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal smNPC_II
ATAC_smNPC_II_SCARB2 = plotSignal(
  data = bw.SCARB2.smNPC_2, 
  params = region.p.SCARB2,
  fill = "#C0C0C0", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.125, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal smNPC_III
ATAC_smNPC_III_SCARB2 = plotSignal(
  data = bw.SCARB2.smNPC_3, 
  params = region.p.SCARB2,
  fill = "#A9A9A9", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.125, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III_SCARB2, 
          at = c(0, scale.SCARB2.max),
          fontsize = 6)

# Add the samples name
plotText(label = "smNPC", 
         fontsize = 6, 
         fontcolor = "#A9A9A9",
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 8.250, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D15 I
ATAC_TH.pos_D15_I_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D15_1, 
  params = region.p.SCARB2,
  fill = "#B22222", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.4375, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D15 II
ATAC_TH.pos_D15_II_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D15_2, 
  params = region.p.SCARB2,
  fill = "#A52A2A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.4375, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D15 III
ATAC_TH.pos_D15_III_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D15_3, 
  params = region.p.SCARB2,
  fill = "#8B0000", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.4375,  
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D15_III_SCARB2, 
          at = c(0, scale.SCARB2.max),
          fontsize = 6)

plotText(label = "mDAN D15", 
         fontsize = 6, 
         fontcolor = "#B22222",
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 8.5625, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D30 I
ATAC_TH.pos_D30_I_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D30_1, 
  params = region.p.SCARB2,
  fill = "#FF6347", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.75, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D30 II
ATAC_TH.pos_D30_II_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D30_2, 
  params = region.p.SCARB2,
  fill = "#FF0000", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.75, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D30 III
ATAC_TH.pos_D30_III_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D30_3, 
  params = region.p.SCARB2,
  fill = "#DC143C", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 8.75, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D30_III_SCARB2, 
          at = c(0, scale.SCARB2.max),
          fontsize = 6)

plotText(label = "mDAN D30", 
         fontsize = 6, 
         fontcolor = "#DC143C",
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 8.875, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# ATACseq signal TH positive neurons D50 I
ATAC_TH.pos_D50_I_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D50_1, 
  params = region.p.SCARB2,
  fill = "#FFA07A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 9.0625, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")


# ATACseq signal TH positive neurons D50 II
ATAC_TH.pos_D50_II_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D50_2, 
  params = region.p.SCARB2,
  fill = "#FA8072", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 9.0625, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

# ATACseq signal TH positive neurons D50 III
ATAC_TH.pos_D50_III_SCARB2 = plotSignal(
  data = bw.SCARB2.TH.pos_D50_3, 
  params = region.p.SCARB2,
  fill = "#E9967A", 
  alpha = 0.7, 
  linecolor = NA,
  x = 0.5, 
  y = 9.0625, 
  width = 4.0, 
  height = 0.25,
  just = c("left", "top"), 
  default.units = "inches")

annoYaxis(plot = ATAC_TH.pos_D50_III_SCARB2, 
          at = c(0, scale.SCARB2.max),
          fontsize = 6)

plotText(label = "mDAN D50", 
         fontsize = 6, 
         fontcolor = "#E9967A",
         fontfamily = "Helvetica",
         x = 0.25, 
         y = 9.125, 
         just = c("left", "top"),
         default.units = "inches",
         fontface = "bold")

# Gene track
plotGenes(
  chrom = "chr4",
  chromstart = 76173717,
  chromend = 76253717,
  assembly = assembly(Genome = "hg38refGene", 
                      TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", 
                      OrgDb = "org.Hs.eg.db"),
  fontcolor = "black",
  fill = "black",
  x = 0.5,
  y = 9.5,
  width = 4.0,
  height = 0.375,
  just = c("left", "top"),
  default.units = "inches")

# Add the SNP position that zoom in the TFBS
SNP_param_SCARB2 = pgParams(chrom = "chr4",
                     chromstart = 76213717,
                     chromend = 76213717,
                     assembly = "hg38")

annoHighlight(plot = ATAC_smNPC_I_SCARB2, 
              params = SNP_param_SCARB2,
              fill = "#404788FF",
              y = 8.25, 
              height = 1.75, 
              just = c("left", "top"), 
              default.units = "inches")

annoZoomLines(plot = ATAC_smNPC_I_SCARB2, 
              params = SNP_param_SCARB2,
              y0 = 10.0, 
              x1 = c(2.25, 2.75), 
              y1 = 10.125, 
              default.units = "inches")

# text G
plotText(label = "G", 
         fontsize = 12,
         fontfamily = "Helvetica",
         x = 4.75, 
         y = 7.75, 
         just = "left", 
         default.units = "inches",
         fontface = "bold")

SCARB2.NR2C2.exp = ggboxplot(dat.RPKM.filt %>% 
                            dplyr::filter(gene_name %in% c("SCARB2", "NR2C2", "FAM47E")) %>% 
                              mutate(gene_name = factor(gene_name,
                                                        levels = c("SCARB2", "NR2C2", "FAM47E"))),
                          x = "Cond", 
                          y = "RPKM", 
                          add = "jitter", 
                          fill = "Cond", 
                          palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"),
                          facet.by = "gene_name",
                          xlab = "", 
                          ylab = "Reads Per Kilobase Million (RPKM)") +
  scale_x_discrete(labels = c("smNPC",
                              "mDAN D15",
                              "mDAN D30",
                              "mDAN D50")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
  

# Place the plot
plotGG(plot = SCARB2.NR2C2.exp, 
       x = 5.0, 
       y = 8.0, 
       width = 3, 
       height = 3.5,
       just = c("left", "top"), 
       default.units = "inches")

# Add the transcription factor binding motif (TFBS) for NR2C2
pfm.NR2C2 = getMatrixByID(JASPAR2020, 
                         ID = "MA1536.1")

# Reverse to position weight matrix as the SNP is on the reverse strand
pfm.NR2C2 = reverseComplement(pfm.NR2C2)

NR2C2 = new.env()

NR2C2$NR2C2 = pfm.NR2C2@profileMatrix

NR2C2 = as.list(NR2C2)

# Remove nucleotide
NR2C2_TFBS = ggseqlogo(NR2C2) +
  annotate("rect", 
           xmin = 3.5, 
           xmax = 4.5, 
           ymin = -0.05, 
           ymax = 2.0, 
           alpha = .1, 
           col = "black", 
           fill = "#404788FF") +
  theme(axis.text.x = element_blank()) 

# Place the plot
plotGG(plot = NR2C2_TFBS, 
       x = 1.5975, 
       y = 10.0, 
       width = 1.5, 
       height = 1.5,
       just = c("left", "top"), 
       default.units = "inches")

pageGuideHide()
dev.off()
```
