
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