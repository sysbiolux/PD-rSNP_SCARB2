## FIGURE 2

author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "31 January, 2025"

### Low-C data processing for STARE analysis



Low-C technique have been applied to 110K TH REP1 mCHERRY smNPC (3 biological replicates), 110K TH REP1 mCHERRY mDANs differentiated for 30 days (1 biological replicate) using [Reinhardt differentiation protocol](https://doi.org/10.1371/journal.pone.0059252).

A Python environment has been set up to run the [FAN-C software](https://github.com/vaquerizaslab/fanc).\
FAN-C requires the installation of HDF5


``` bash
# In a specific directory
cd $HOME/Tools
mkdir hdf5-build
cd hdf5-build

# Download version 1.10.5
wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz

# Unpack
tar xzf hdf5-1.10.5.tar.gz
rm hdf5-1.10.5.tar.gz
cd hdf5-1.10.5

# Install 
./configure --prefix=$HOME
make
make install
```

Create the virtual Python environment


``` bash
# Load modules to use Python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Create the environment in my home directory with ther other environments
python3 -m venv ~/Environment/FANC

# Activate the newly created environment
source ~/Environment/FANC/bin/activate

# And install FAN-C
pip3 install fanc

# Check that FAN-C version
fanc --version
```

There is an error "ImportError: cannot import name 'GC' from 'Bio.SeqUtils'. The version of biopython is 1.83 and the GC function is still present in the 1.75 version but I cannot find it in the 1.83 version. Remove the 1.83 version and install the 1.75.


``` bash
# Uninstall the 1.83 version
pip uninstall biopython

# And install the 1.75
pip install biopython==1.75
```

Check FAN-C now


``` bash
fanc --version
```

**Conclusion** : FAN-C version 0.9.27 has been successfully installed!

#### *1. Run fastqc on the samples to check their quality*

Run the script for the 1st biological replicate (smNPC (N1 + N2) and mDAN D30 (N1))


``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/fastqc.sh
```

Display the script


``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/fastqc.sh
```

```
#!/bin/bash -l
#SBATCH -J N1_smNPC_mDAN_fastqc
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --time=08:00:00
#SBATCH -p batch
#SBATCH --qos=normal

# Load module containing fastqc
module load bio/FastQC/0.11.9-Java-11

# Check fastqc version
fastqc --v

# Perform quality control using fastqc module.
for i in $SCRATCH/LowC_smNPC_THpos_neur/fastq/*.gz
do
	fastqc -o $SCRATCH/LowC_smNPC_THpos_neur/FASTQC_res/ $i
done
```
Since I have to use the same genome version as the one used for the ATACseq data and as it contains a lot of contigs that are irrelevant for downstream Hi-C analysis, limit the analysis to the canonical chromosomes and use the `fanc fragments` command first to generate a "map" of canonical chromosomes that will be passed to `fanc auto` using the `-g` parameter.

``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_in_silico_digestion.sh   # Takes 5 minutes booking a full node on iris
```
Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_in_silico_digestion.sh
```

```
#!/bin/bash -l
#SBATCH -J N1_smNPC_mDAN_in_silico_digestion
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --time=06:00:00

# This script is for performing in silico digestion of the genome using the restriction enzyme MboI (the same enzyme that been used for the wetlab experiment) #
# Load modlue containing python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Activate the python environment for FAN-C
source ~/Environment/FANC/bin/activate

# Check FAN-C version
fanc --version

# Run in-silico genome digestion
fanc fragments -c 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' \
$SCRATCH/bwa_index/GRCh38.genome.fa.gz \
'MboI' \
$SCRATCH/FANC_hg38.p1.cano.chr.bed

deactivate
```

#### *2. Run FAN-C in auto mode to do the mapping, normalisation and generating hic matrices*
Process each of the biological replicates separately as suggested by [the FAN-C authors](https://github.com/vaquerizaslab/fanc/issues/24).  

Run the script for the 1st biological replicate (smNPC and mDAN D30).

``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_launcher_N1.sh
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_launcher_N1.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_N2_smNPC_N1_mDAN_smNPC_Op
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=96:00:00
#SBATCH -p batch
#SBATCH --qos=long

# Load modlue containing python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Activate the python environment for FAN-C
source ~/Environment/FANC/bin/activate

# Check FAN-C version
fanc --version

# As FAN-C needs bwa for the mapping step, load the module containing bwa
module load bio/BWA/0.7.17-GCC-10.2.0

# And check version
bwa

# Set a TMPDIR variable
export TMPDIR=$SCRATCH/tempfiles/
echo $TMPDIR

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# Run FAN-C in automatic mode smNPC - Biological replicate 1 - Normalisation method is Knight-Ruiz
fanc auto $SCRATCH/LowC_smNPC_THpos_neur/fastq/N1_TH_REP1_mCHERRY_smNPC_LOW_C_LowTH_SB_S2_R1_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/fastq/N1_TH_REP1_mCHERRY_smNPC_LOW_C_LowTH_SB_S2_R2_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC \
-g $SCRATCH/FANC_hg38.p1.cano.chr.bed \
-i $SCRATCH/bwa_index/GRCh38.genome.fa.gz \
-r MboI \
-n N1_smNPC \
-t 14 \
--le-inward-cutoff 5000 \
--le-outward-cutoff 5000 \
--fanc-parallel \
--norm-method KR \
--iterative \
--split-ligation-junction \
-q 3 \
-tmp

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

deactivate
```

The run is done for N1_smNPC but not for N1_TH^+^ day 30 neurons (max 48 hours on qos batch on iris).  

Run the script for the 1st biological replicate of neurons, the 2nd biological replicate of smNPC and the smNPC sample that has been generated during the optimisation 

``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_launcher_N1_TH_N2_smNPC.sh
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_launcher_N1_TH_N2_smNPC.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_N2_smNPC_N1_mDAN_smNPC_Op
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=96:00:00
#SBATCH -p batch
#SBATCH --qos=long

# Load modlue containing python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Activate the python environment for FAN-C
source ~/Environment/FANC/bin/activate

# Check FAN-C version
fanc --version

# As FAN-C needs bwa for the mapping step, load the module containing bwa
module load bio/BWA/0.7.17-GCC-10.2.0

# And check version
bwa

# Set a TMPDIR variable
export TMPDIR=$SCRATCH/tempfiles/
echo $TMPDIR

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# Run FAN-C in automatic mode for mDAN neurons day 30 - Biological replicate 1 - Normalisation method is Knight-Ruiz
fanc auto $SCRATCH/LowC_smNPC_THpos_neur/fastq/N1_TH_REP1_mCHERRY_TH_Neur_day30_LOW_C_LowTH_SB_S1_R1_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/fastq/N1_TH_REP1_mCHERRY_TH_Neur_day30_LOW_C_LowTH_SB_S1_R2_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30 \
-g $SCRATCH/FANC_hg38.p1.cano.chr.bed \
-i $SCRATCH/bwa_index/GRCh38.genome.fa.gz \
-r MboI \
-n N1_mDAN_D30 \
-t 14 \
--le-inward-cutoff 5000 \
--le-outward-cutoff 5000 \
--fanc-parallel \
--norm-method KR \
--iterative \
--split-ligation-junction \
-q 3 \
-tmp

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# Run FAN-C in automatic mode for smNPC - Biological replicate 2 - Normalisation method is Knight-Ruiz
fanc auto $SCRATCH/LowC_smNPC_THpos_neur/fastq/N2_TH_REP1_mCHERRY_smNPC_LOW_C_LowTH_SB_S3_R1_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/fastq/N2_TH_REP1_mCHERRY_smNPC_LOW_C_LowTH_SB_S3_R2_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC \
-g $SCRATCH/FANC_hg38.p1.cano.chr.bed \
-i $SCRATCH/bwa_index/GRCh38.genome.fa.gz \
-r MboI \
-n N2_smNPC \
-t 14 \
--le-inward-cutoff 5000 \
--le-outward-cutoff 5000 \
--fanc-parallel \
--norm-method KR \
--iterative \
--split-ligation-junction \
-q 3 \
-tmp

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# Run FAN-C in automatic mode for smNPC - Sample used for the optimisation - Normalisation method is Knight-Ruiz
fanc auto $SCRATCH/LowCO/LowCO/fastq/Optimisation_LowC_TH_REP1_mCHERRY_smNPC_lowC_LowCO_SB_S2_R1_001.fastq.gz \
$SCRATCH/LowCO/LowCO/fastq/Optimisation_LowC_TH_REP1_mCHERRY_smNPC_lowC_LowCO_SB_S2_R2_001.fastq.gz \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC \
-g $SCRATCH/FANC_hg38.p1.cano.chr.bed \
-i $SCRATCH/bwa_index/GRCh38.genome.fa.gz \
-r MboI \
-n opti_smNPC \
-t 14 \
--le-inward-cutoff 5000 \
--le-outward-cutoff 5000 \
--fanc-parallel \
--norm-method KR \
--iterative \
--split-ligation-junction \
-q 3 \
-tmp

deactivate
```

Convert newly created pairs files into hic files readable by other tools (like Juicer)

``` bash
# Copy the pairs files to be able to work on parallel with them
cp $SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/pairs/N1_smNPC.pairs $SCRATCH/
cp $SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/pairs/N1_mDAN_D30.pairs $SCRATCH/
cp $SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/pairs/N2_smNPC.pairs $SCRATCH/
cp $SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/pairs/opti_smNPC.pairs $SCRATCH/

# And run
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_hic_to_juicer.sh
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_hic_to_juicer.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_hic_to_juicer
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH -p batch
#SBATCH --qos=long

# Load module containing python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Activate the python environment for FAN-C
source ~/Environment/FANC/bin/activate

# Check FAN-C version
fanc --version

# As FAN-C needs bwa for the mapping step, load the module containing bwa
module load bio/BWA/0.7.17-GCC-10.2.0

# And check version
bwa

# Remove tempfiles and set a TMPDIR variable
rm -rf $SCRATCH/tempfiles/*
export TMPDIR=$SCRATCH/tempfiles/
echo $TMPDIR

# Pairs files are in SCRATCH
cd $SCRATCH

# 5kb resolution
#parallel -j 4 "fanc to-juicer {} {.}.juicer.5kb.hic --juicer-tools-jar $HOME/juicer_tools.2.20.00.jar -tmp -r 5000" ::: *.pairs

# 10kb resolution
#parallel -j 4 "fanc to-juicer {} {.}.juicer.10kb.hic --juicer-tools-jar $HOME/juicer_tools.2.20.00.jar -tmp -r 10000" ::: *.pairs

# 25kb resolution
#parallel -j 4 "fanc to-juicer {} {.}.juicer.25kb.hic --juicer-tools-jar $HOME/juicer_tools.2.20.00.jar -tmp -r 25000" ::: *.pairs

# 50kb resolution
parallel -j 4 "fanc to-juicer {} {.}.juicer.50kb.hic --juicer-tools-jar $HOME/juicer_tools.2.20.00.jar -tmp -r 50000" ::: *.pairs

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

deactivate
```


#### Merge 3 biological replicates of smNPC into a mega map. Do the same for the mDAN D30

Recap from https://github.com/sysbiolux/PD-rSNP_SCARB2/tree/master/Juicer_pipeline:

Create a juicer apptainer image following https://github.com/aidenlab/juicer/tree/main/Docker:

``` bash
module load tools/Apptainer
singularity pull juicer.sif docker://aidenlab/juicer:v2.0.1
cd /mnt/aiongpfs/projects/lowc_mdan
apptainer exec --no-home juicer.sif ls -l
```

All those steps are in the `launcher_merge_mDAN.sh` (copy over for `smNPC`) that call for the container `run_merge_juicer.sh`

Output files were renamed before sending to Dennis as they have identical names in the their sub-folders:

``` bash
4bc93055d922c0591740db908d6127c5  mega_mDAN.D30_inter_q1.hic (4.2GB)
cd2f7547ec62ad6f66445135363e8802  mega_mDAN.D30_inter_q30.hic (3.2GB)
4c2c5d4eddf31c9fd89ac2bf6cc1e5e2  mega_smNPC_inter_q1.hic (7.8GB)
5908d20f63151318cd52ef01292a830f  mega_smNPC_inter_q30.hic (7.2GB)
```


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