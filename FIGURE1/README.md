
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