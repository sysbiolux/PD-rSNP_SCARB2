
## Creating contact matrices with HicExplorer (h5 format)

Using https://hicexplorer.readthedocs.io/en/latest/

Work on HPC after installing hicexplorer in a micromamba environment

### Prepare restriction sites file

``` bash
micromamba activate hicexp # v3.7.6
seqtk  seq -l 120 juicer/references/GRCh38.genome.fa.gz > GRCh38.genome.fa
hicFindRestSite -f GRCh38.genome.fa -p GATC -o hg38_MboI.bed
```

### Align reads and build contact matrices

``` bash
# for each replicate
module load env/deprecated
module load bio/SAMtools
module load bio/BWA

micromamba activate hicexp

for FOLDER in LowC_mDAN.D30_0? LowC_smNPC_0?
  do echo $FOLDER
  bwa mem -A 1 -B 4 -E 50 -L 0 -t 12 /work/projects/lowc_mdan/GRCh38.genome.fa ${FOLDER}/*R1*fastq.gz | \
    samtools view -Shb - > ${FOLDER}/R1.bam
  bwa mem -A 1 -B 4 -E 50 -L 0 -t 12 /work/projects/lowc_mdan/GRCh38.genome.fa ${FOLDER}/*R2*fastq.gz | \
    samtools view -Shb - > ${FOLDER}/R2.bam

  # build matrix with HicExplorer, MQ 15 is the default
  hicBuildMatrix --samFiles ${FOLDER}/R1.bam ${FOLDER}/R2.bam --outBam ${FOLDER}/hicexp.bam \
  --outFileName ${FOLDER}/hic_10kb.h5 --QCfolder ${FOLDER}/HiC_file_10kb_QC --binSize 10000 \
  --restrictionSequence GATC --danglingSequence GATC \
  --restrictionCutFile /work/projects/lowc_mdan/hg38_MboI.bed --threads 8 --inputBufferSize 400000 \
  --genomeAssembly hg38 --minMappingQuality 15
done

```

Send over the `.h5` files to Jafar Sharif for downstream analysis.


### RNA-seq signal 

500kb bins needed for correlation with Hi-C data

- Convert to `saf` format

``` bash
module load bio/R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
awk 'BEGIN{OFS="\t"}{print $1"_"$2"_"$3, $1, $2, $3, "."}' 500kb_bins.bed > 500kb_bins.saf
```

Count reads in each bin with `Rsubread`

``` bash
library(Rsubread)

bam <- c("20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/smNPC_20190312_reads.Aligned.sortedByCoord.out.q30.bam",
         "20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/D30_possort_20190228_reads.Aligned.sortedByCoord.out.q30.bam",
         "20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/smNPC_20190319_reads.Aligned.sortedByCoord.out.q30.bam",
         "20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/D30_possort_20190328_reads.Aligned.sortedByCoord.out.q30.bam",
         "20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/smNPC_20190322_reads.Aligned.sortedByCoord.out.q30.bam",
         "20190712_RNA_seq_Samples_NESC_D15_D30_D50_ASTRO/STAR/q30/D30_possort_20190321_reads.Aligned.sortedByCoord.out.q30.bam")


message("BAMS: ", head(bam), length(bam))

ref <- "500kb_bins.saf"

fc <- featureCounts(
  as.character(bam), isPairedEnd = FALSE,
  annot.ext = ref,
  isGTFAnnotationFile = FALSE,
  minMQS = 15,
  strandSpecific = 2,
  nthreads = 8,
  useMetaFeatures = FALSE, allowMultiOverlap = FALSE)

colnames(fc$counts) <- sub("_reads.Aligned.sortedByCoord.out.q30.bam", "",  colnames(fc$counts))

head(fc$stat)

saveRDS(fc, "rna_seq_500k_bins_fc.rds")

```


## Compartment analysis and TAD calling by Jafar Sharif

### Sum matrices

``` bash
# hicexplorer, version 3.7.2, bioconda python3.8
$ hicSumMatrices -m replicate1_10kb.h5 replicate2_10kb.h5 replicate3_10kb.h5  -o merged_10kb.h5
```


### merge matrix bins

``` bash
# hicexplorer, version 3.7.2, bioconda python3.8
$ hicMergeMatrixBins -m merged_10kb.h5 --numBins 5 -o merged_50kb.h5
$ hicMergeMatrixBins -m merged_10kb.h5 --numBins 50 -o merged_500kb.h5
```

### correct matrices

``` bash
# hicexplorer, version 3.7.2, bioconda python3.8
$ hicCorrectMatrix diagnostic_plot -m merged_10kb.h5 --plotName diagnostic_plot_merged_10kb.png --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
$ hicCorrectMatrix correct -m merged_10kb.h5 --filterThreshold "lower limit" "upper limit" --perchr -o corrected_merged_10kb.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
$ hicCorrectMatrix correct -m merged_50kb.h5 --filterThreshold "lower limit" "upper limit" --perchr -o corrected_merged_50kb.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
$ hicCorrectMatrix correct -m merged_500kb.h5 --filterThreshold "lower limit" "upper limit" --perchr -o corrected_merged_500kb.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
```

### convert format

``` bash
# hicexplorer, version 3.7.2, bioconda python3.8
$ hicConvertFormat -m corrected_merged_10kb.h5 --inputFormat h5 --outputFormat cool -o corrected_merged_10kb.cool
$ hicConvertFormat -m corrected_merged_50kb.h5 --inputFormat h5 --outputFormat cool -o corrected_merged_50kb.cool
$ hicConvertFormat -m corrected_merged_500kb.h5 --inputFormat h5 --outputFormat cool -o corrected_merged_500kb.cool
```

### find TADs

``` bash
# hicexplorer, version 3.7.2, bioconda python3.8
$ hicFindTADs -m corrected_merged_10kb.h5 --minDepth 30000 --maxDepth 60000 --numberOfProcessors 16 --outPrefix matrix_10kb --minBoundaryDistance 200000 --correctForMultipleTesting fdr --threshold 0.05
```

### calculate PC1

``` bash
# fanc, version 0.9.23b, bioconda python3.7
# calculated PC1 values were manually checked to correct the direction of the eigenvector (+ or -) for each chromosome
$ fanc compartments corrected_merged_500kb.cool corrected_merged_500kb.ab
$ fanc compartments -v corrected_merged_500kb.ev.txt corrected_merged_500kb.cool corrected_merged_500kb.ab

# fanc plot, for whole chromosomes or selected regions
# fanc, version 0.9.23b, bioconda python3.7
# O/E matrix
$ fancplot -o region_plot.png chrXX:start-end -p triangular -vmin -2 -vmax 2 -c bwr -e corrected_merged_50kb.cool
$ fancplot -o whole_chromosome_plot.png chrXX -p square -c bwr -e corrected_merged_50kb.cool -vmin -2 -vmax 2
```


### calculation of intra- or inter- compartment strengths

```
# pentad, https://github.com/magnitov/pentad, python3.7
# cis calculation (intra-compartment interactions)
$ python src/get_pentad_cis.py corrected_merged_500kb.cool corrected_merged_500kb_PC1.txt --out_pref corrected_merged_500kb_cis
$ python src/plot_pentad.py  corrected_merged_500kb_cis.json --title corrected_merged_500kb_cis --out_pref corrected_merged_500kb_cis_plot --vmin "lower limit" --vmax "upper limit"
# distance calculation (inter-compartment interactions)
$ python src/get_pentad_distance.py corrected_merged_500kb.cool corrected_merged_500kb_PC1.txt --out_pref corrected_merged_500kb_distance --distance 10 25
$ python src/plot_pentad.py corrected_merged_500kb_distance.json --title corrected_merged_500kb_distance --out_pref corrected_merged_500kb_distance_plot --vmin "lower limit" --vmax "upper limit"
```

### PC1 dotplot

RStudio, R version 4.2.1

``` r
PC1_group_dotplot <- read_csv("PC1_group_dotplot.csv")
sp <- ggplot(PC1_group_dotplot, aes(x = PC1_NPC, y = PC1_DAN, color = group)) +
  geom_point(aes(color = Group), alpha = 0.6) +
  scale_color_manual(values = c('Blue', 'gray', 'Red')) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    panel.border = element_rect(colour = "black", fill = NA)
  )
sp
```


### pentads heatmap plot

RStudio, R version 4.2.1

``` r
library(gplots)
intra_comp_matrix <- read_csv("R_heatmap_DAN_01_intra_A.csv")

data_matrix <- data.matrix(intra_comp_matrix, rownames.force = NA)
heatmap.2(
  data_matrix,
  dendrogram = 'none',
  Rowv = FALSE,
  Colv = FALSE,
  trace = 'none',
  col = colorRampPalette(c('royalblue1', 'white', 'firebrick2'))(1000),
  breaks = seq(0.7, 1.2, length.out = 1001),
  density.info = 'none'
)
heatmap.2(
  data_matrix,
  dendrogram = 'none',
  Rowv = FALSE,
  Colv = FALSE,
  trace = 'none',
  col = colorRampPalette(c('springgreen4', 'white', 'lightpink3'))(1000),
  breaks = seq(-0.2, 0.2, length.out = 1001),
  density.info = 'none'
)
```


