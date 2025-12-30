## Allelic Imbalance of PD SNPs


``` r
PD_rSNPs_254 <- read_tsv("../SNPs/254rSNPs_GWAS_pval_rsID.bed", 
                         col_names = c("chr", "start", "end", "SNP_position", "name"))


# Make a Granges object of the 254 PD rSNPs (width has to be one)
searchArea_254 <- makeGRangesFromDataFrame(PD_rSNPs_254, 
                                           starts.in.df.are.0based = TRUE,
                                           keep.extra.columns = TRUE)

reads_254 <- impBamGAL(UserDir = "ATAC_BAMs", 
                       files = c("D15_pos.sorted.bam", "D30_pos.sorted.bam", "D50_pos.sorted.bam" , "smNPC.sorted.bam"),
                       searchArea = searchArea_254, 
                       verbose = TRUE)

# Count the number of reads for the 2 different alleles at that position
countList_ase_254 <- ASEsetFromCountList(searchArea_254, 
                                         getAlleleCounts(BamList = reads_254, 
                                                         GRvariants = searchArea_254, 
                                                         verbose = TRUE))

# Now check how many are truly heterozygous thanks to the whole genome sequencing data (WGS)
# Use table browser to retrieve according to this link http://genome.ucsc.edu/FAQ/FAQreleases.
# Load now the coordinates that are in hg19
PD_rSNPs_254_hg19 <- read_tsv("../SNPs/rsIDs_254PDrSNPs_hg19_coord.txt",
                               show_col_types = FALSE)

# Only keep SNPs that are on canonical chromosomes
PD_rSNPs_254_hg19.GR <- PD_rSNPs_254_hg19 |> 
  dplyr::filter(!str_detect(`#chrom`, "_")) |> 
  dplyr::select(`#chrom`:name) |> 
  dplyr::rename(chr = `#chrom`) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqnames.field = "chr",
                           start.field = "chromStart",
                           end.field = "chromEnd",
                           starts.in.df.are.0based = TRUE)

# Load WGS data of TH REP1 mCHERRY cell line (hg19 genome version)
dat.SNPs.only.GR <- read_rds("../SNPs/E19D017a78.merge.qc.SNP_INDEL.hg19.annotated.rds") |> 
  dplyr::select(CHR:POS, REF:FILTER, E19D017a78) |> 
  mutate(START = POS - 1) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqnames.field = "CHR",
                           start.field = "START",
                           end.field = "POS",
                           starts.in.df.are.0based = TRUE)

GenomeInfoDb::seqlevelsStyle(dat.SNPs.only.GR) <- "UCSC"

# Overlap
hits <- findOverlaps(PD_rSNPs_254_hg19.GR, dat.SNPs.only.GR)
rSNP.254_in_E19D017a78 <- PD_rSNPs_254_hg19.GR[queryHits(hits)]
E19D017a78_w_rSNP.254 <- dat.SNPs.only.GR[subjectHits(hits)]

mcols(rSNP.254_in_E19D017a78) <- cbind(mcols(rSNP.254_in_E19D017a78),
                                       mcols(E19D017a78_w_rSNP.254))

rSNP.254_in_E19D017a78 <- as.data.frame(rSNP.254_in_E19D017a78) |> 
  as_tibble()
rSNP.254_in_E19D017a78
```

```
# A tibble: 242 × 10
   seqnames    start      end width strand name       REF   ALT   FILTER E19D017a78                                                         
   <fct>       <int>    <int> <int> <fct>  <chr>      <chr> <chr> <chr>  <chr>                                                              
 1 chr17    43930238 43930238     1 *      rs79589869 C     A     PASS   ./.:.:.:.:.:.:.:.:.:.:.                                            
 2 chr17    44345063 44345063     1 *      rs2732651  C     T     PASS   0/1:43,25:0.368:68:17,16:26,9:450,0,450:99:18,25,12,13:887,0,1384:…
 3 chr17    43503000 43503000     1 *      rs62064651 G     A     PASS   0/1:11,11:0.5:22:7,5:4,6:319,0,392:99:4,7,5,6:319,0,392:4,7,7,4    
 4 chr17    43503284 43503284     1 *      rs76344126 G     A     PASS   0|1:14,13:0.481:27:4,5:10,8:450,0,450:99:7,7,7,6:504,0,1016:435032…
 5 chr17    43503294 43503294     1 *      rs7209501  A     C     PASS   0|1:14,13:0.481:27:4,5:10,8:450,0,450:99:7,7,7,6:504,0,1016:435032…
 6 chr17    44019643 44019643     1 *      rs62062769 G     A     PASS   ./.:.:.:.:.:.:.:.:.:.:.                                            
 7 chr17    44019680 44019680     1 *      rs62062770 T     C     PASS   ./.:.:.:.:.:.:.:.:.:.:.                                            
 8 chr17    43828935 43828935     1 *      rs17426106 G     C     PASS   ./.:.:.:.:.:.:.:.:.:.:.                                            
 9 chr17    43829353 43829353     1 *      rs62054442 A     G     PASS   ./.:.:.:.:.:.:.:.:.:.:.                                            
10 chr17    43825339 43825339     1 *      rs62054435 C     G     PASS   ./.:.:.:.:.:.:.:.:.:.:. 
```

Out of 254 PD-rSNPs, we have information about 242

``` r

# Check how many are "unknown" (./.:.)
rSNP.254_in_E19D017a78 |> 
  mutate(Genotype = gsub("\\:.*", "\\1", E19D017a78)) |> 
  summarise(GT.num = n(), .by = Genotype) 
```

```
# A tibble: 5 × 2
  Genotype GT.num
  <chr>     <int>
1 ./.         164
2 0/1          45
3 0|1           9
4 1/1          21
5 1|1           3
```

54 PD-rSNPs are heterozygous in TH REP1 mCHERRY cell line, 24 are homozygous for the alternative allele and 164 are unknown

- Check with a chi-squared test how many of the 54 heterozygous PD-rSNPs lead to chromatin allelic imbalance in at least one of the 4 samples (smNPC, mDAN-D15, mDAN-D30, mDAN-D50)

``` r
# Check if the 54 PD-rSNPs that are heterozygous in TH REP1 mCHERRY cell line also lead to chromatin allelic imbalance
# rsIDs of the 54 heterozygous PD-rSNPs 
rSNP.254_in_E19D017a78_0.1 <- rSNP.254_in_E19D017a78 |> 
  mutate(Genotype = gsub("\\:.*", "\\1", E19D017a78),
         Genotype2 = gsub("\\|", "\\/", Genotype)) |> 
  dplyr::filter(str_detect(Genotype2, "0/1")) |>
  dplyr::select(name)

# Get their coordinates
hetero_54_rSNP_in_E19D017a78.hg38 <- searchArea_254 |> 
  as_tibble() |> 
  dplyr::select(-start) |> # Not the real start when converting a GRanges to a tibble 
  mutate(start = end - 1) |> 
  dplyr::select(seqnames, start, end, name) |> 
  dplyr::filter(name %in% rSNP.254_in_E19D017a78_0.1$name) |> 
  arrange(seqnames, start) |> 
  unite(col = "coordToInves", seqnames, end, sep = "_")

# Extract the 54 heterozygous PD-rSNPs from the ASet object
PDrSNPs_54_het <- countList_ase_254[names(countList_ase_254) %in% hetero_54_rSNP_in_E19D017a78.hg38$coordToInves, ]


searchArea.0.1 <- searchArea_254 |> 
  as_tibble() |> 
  mutate(location = paste0(seqnames, '_', start)) |> 
  dplyr::filter(name %in% hetero_54_rSNP_in_E19D017a78.hg38$name) |> 
  arrange(factor(location, levels = rownames(PDrSNPs_54_het))) |> 
  makeGRangesFromDataFrame(keep.extra.columns = FALSE,
                           seqnames.field = "seqnames",
                           start.field = "start",
                           end.field = "end",
                           starts.in.df.are.0based = FALSE)


allele_counts <- getAlleleCounts(BamList = reads_254, 
                                 GRvariants = searchArea_254, 
                                 verbose = TRUE)

allele_counts_54 <- allele_counts[names(allele_counts) %in% rownames(PDrSNPs_54_het)]


# And make it a ASet object back and proceed with a chi-squared test
PDrSNPs_54_het.ASE <- ASEsetFromCountList(searchArea.0.1, allele_counts_54, verbose = TRUE)

PDrSNPs_54_het.ASE.chisq <- chisq.test(PDrSNPs_54_het.ASE) |> 
  as_tibble() |> 
  mutate(Sample = colnames(PDrSNPs_54_het.ASE), .before = 1)

# Consider that a PD-rSNP lead to allelic imbalance if at least one sample (smNPC, mDAN-D15, mDAN-D30, mDAN-D50) shows allelic imbalance
PDrSNPs_54_het.ASE.chisq |> 
  select_if(\(x) any(x < 0.05))
# 46 among 54 PD-rSNPs show allelic imbalance
write_tsv(PDrSNPs_54_het.ASE.chisq, "Gerard_et_al_2026/FIGURE2/PDrSNPs_54_het.ASE.chisq.tsv")
```

## Manhattan plot

Produce as a PNG to be used in the figure panel.

``` r
# Associated then the significant GWAS pvalue
snp_2_high <- read_tsv("sneep_254_gwas_p_rsid.tsv", show_col_types = FALSE) |> 
  filter(!is.na(name),
                gwas_pvalue <= 5e-08) |> 
  select(chr,start, snp = name, gwas_pvalue)
  


snp_to_pl <- read_tsv("nalls_allSNPs_hg38.tsv.gz", show_col_types = FALSE) |>
  select(chr = seqnames, start, gwas_pvalue = p) |> 
  bind_rows(snp_2_high) |>
  mutate(p = -log10(gwas_pvalue),
         chr = factor(chr,
                        levels = paste0("chr",
                                        1:22)),
         color = if_else(is.na(snp), "no_high", "high"),
         label = if_else(snp %in% c("rs1465922", "rs144814361"), snp, NA_character_))

highlight_colormap <- c("no_high" = adjustcolor( "grey", 
                                                 alpha.f = 0.2), 
                        "high" = "#6600FF")

no_high <- manhattan_data_preprocess(snp_to_pl,
                                     pval.colname = "gwas_pvalue",
                                     chr.colname = "chr",
                                     pos.colname = "start",
                                     highlight.colname = "color",
                                     highlight.col = highlight_colormap,
                                     signif = 5e-08)
# Plot
snp.pl <- manhattan_plot(x = no_high,
                         color.by.highlight = TRUE,
                         rescale = TRUE,
                         label.font.size = 2.5,
                         label.colname = "label") +
  theme(axis.text.x = element_text(size = 7, hjust = 1, vjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))

ggsave("manhattan_plot.png", plot = snp.pl, dpi = "print", height = 2, width = 3.75)
```