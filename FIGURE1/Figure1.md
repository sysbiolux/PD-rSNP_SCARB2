---
title: "FIGURE 1"
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

### FIGURE 1


``` r
########## FIGURE 1 ########## Save as TIFF, 300 ppi
tiff("/home/vagrant/Manuscript_1/FIGURE1/FIGURE1.tiff", width = 8.27, height = 9.5,
    units = "in", res = 300, compression = "lzw")

# Create a A4 blank page
pageCreate(width = 8.27, height = 9.5, default.units = "inches", showGuides = TRUE)

#### PANEL A - text Figure 1
plotText(label = "Figure 1", fontsize = 14, fontfamily = "Helvetica", x = 0.25, y = 0.25,
    just = "left", default.units = "inches", fontface = "bold")

# text A
plotText(label = "A", fontsize = 12, fontfamily = "Helvetica", x = 0.25, y = 0.5,
    just = "left", default.units = "inches", fontface = "bold")

########################################## Figure 1A Data generation scheme
########################################## #### Figure 1A generated in
########################################## Biorender
fig1A = readPNG("/home/vagrant/Manuscript_1/FIGURE1/Fig1A_Biorender_vertical.png")

plotRaster(image = fig1A, x = 0.5, y = 2.5, default.units = "inches", width = 2.5,
    height = 3.75, just = "left", interpolate = TRUE)

########################################## Figure 1B Odd ratio
########################################## ################# text B
plotText(label = "B", fontsize = 12, fontfamily = "Helvetica", x = 3, y = 0.5, just = "left",
    default.units = "inches", fontface = "bold")

# Load data obtained from Jochen
RNAseq_odd = read_delim("/home/vagrant/Manuscript_1/FIGURE1/REFORMAT_tpm_gene_sets.gsa.txt",
    delim = "\t", col_names = TRUE)

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
    nudge_y = 0.4) + xlab("Odds ratio \n(95% Confidence Interval)") + ylab("")

# Place the plot
plotGG(plot = fig1B, x = 3.25, y = 0.5, width = 4.5, height = 4.5, just = c("left",
    "top"), default.units = "inches")
```

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

Run `fanc expected` command which calculates the expected values of the contact matrix

``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_hic_expected.sh    # Needs 5 minutes using ine full node on iris
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_hic_expected.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_hic_O/E
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=60:00:00
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

# Run fan-c expected command for every hic files and resolutions
## N1-smNPC and N1-mDAN D30
### 1 mb resolution
echo "Calculating expected values for the 1mb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_1mb" "N1-mDAN D30_1mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/1mb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_1mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_1mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/1mb_N1-smNPC_N1-mDAN-D30_OE.txt

### 2 mb resolution
echo "Calculating expected values for the 2mb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_2mb" "N1-mDAN D30_2mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/2mb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_2mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_2mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/2mb_N1-smNPC_N1-mDAN-D30_OE.txt

### 5 mb resolution
echo "Calculating expected values for the 5mb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_5mb" "N1-mDAN D30_5mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5mb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_5mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_5mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5mb_N1-smNPC_N1-mDAN-D30_OE.txt

### 5 kb resolution
echo "Calculating expected values for the 5kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_5kb" "N1-mDAN D30_5kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_5kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_5kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 10 kb resolution
echo "Calculating expected values for the 10kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_10kb" "N1-mDAN D30_10kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/10kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_10kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_10kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/10kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 25 kb resolution
echo "Calculating expected values for the 25kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_25kb" "N1-mDAN D30_25kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/25kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_25kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_25kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/25kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 50 kb resolution
echo "Calculating expected values for the 50kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_50kb" "N1-mDAN D30_50kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/50kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_50kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_50kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/50kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 100 kb resolution
echo "Calculating expected values for the 100kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_100kb" "N1-mDAN D30_100kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/100kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_100kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_100kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/100kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 250 kb resolution
echo "Calculating expected values for the 250kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_250kb" "N1-mDAN D30_250kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/250kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_250kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_250kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/250kb_N1-smNPC_N1-mDAN-D30_OE.txt

### 500 kb resolution
echo "Calculating expected values for the 500kb resolution on N1-smNPCs and N1-mDAN D30"
fanc expected -l "N1-smNPCs_500kb" "N1-mDAN D30_500kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/500kb_N1-smNPC_N1-mDAN-D30_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/binned/N1_smNPC_500kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/binned/N1_mDAN_D30_500kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/500kb_N1-smNPC_N1-mDAN-D30_OE.txt

# Run fan-c expected command for every hic files and resolutions
## N2-smNPC and opt-smNPC
### 1 mb resolution
echo "Calculating expected values for the 1mb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_1mb" "opt-smNPC_1mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/1mb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_1mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_1mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/1mb_N2-smNPC_N1-opt-smNPC_OE.txt

### 2 mb resolution
echo "Calculating expected values for the 2mb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_2mb" "opt-smNPC_2mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/2mb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_2mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_2mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/2mb_N2-smNPC_N1-opt-smNPC_OE.txt

### 5 mb resolution
echo "Calculating expected values for the 5mb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_5mb" "opt-smNPC_5mb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5mb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_5mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_5mb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5mb_N2-smNPC_N1-opt-smNPC_OE.txt

### 5 kb resolution
echo "Calculating expected values for the 5kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_5kb" "opt-smNPC_5kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_5kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_5kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/5kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 10 kb resolution
echo "Calculating expected values for the 10kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_10kb" "opt-smNPC_10kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/10kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_10kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_10kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/10kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 25 kb resolution
echo "Calculating expected values for the 25kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_25kb" "opt-smNPC_25kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/25kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_25kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_25kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/25kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 50 kb resolution
echo "Calculating expected values for the 50kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_50kb" "opt-smNPC_50kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/50kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_50kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_50kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/50kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 100 kb resolution
echo "Calculating expected values for the 100kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_100kb" "opt-smNPC_100kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/100kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_100kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_100kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/100kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 250 kb resolution
echo "Calculating expected values for the 250kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_250kb" "opt-smNPC_250kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/250kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_250kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_250kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/250kb_N2-smNPC_N1-opt-smNPC_OE.txt

### 500 kb resolution
echo "Calculating expected values for the 500kb resolution on N2-smNPC and opt-smNPC"
fanc expected -l "N2-smNPC_500kb" "opt-smNPC_500kb" \
-p $SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/500kb_N2-smNPC_opt-smNPC_OE.pdf \
-tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/binned/N2_smNPC_500kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/binned/opti_smNPC_500kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/Obs_Expected_res/500kb_N2-smNPC_N1-opt-smNPC_OE.txt

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

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

#### *3. hic files are now available per replicate. Merge them*

``` bash
# And run
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_replicate_merge_smNPC.sh
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_replicate_merge_smNPC.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_merge_replicate_smNPC
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

# Run FAN-C in automatic mode for merging hic files for smNPCs (3 biological replicates)
fanc auto $SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/N1_smNPC.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_N2_smNPC/hic/N2_smNPC.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_output_opti_smNPC/hic/opti_smNPC.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC \
-n merged_smNPC \
-t 14

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

deactivate
```

#### *4. hic files input for SNEEP*
It is better to run SNEEP on the hic files from smNPC that have been merged. The command `fanc to-juicer` works using `.pairs` files. Using the merge smNPC file, convert it to a text file using `fanc dump`
It works but Dennis needs the 5kb and 10kb matrix with bins. Rerun.

``` bash
sbatch $SCRATCH/LowC_smNPC_THpos_neur/scripts/FANC_hic_to_txt_and_bin_nolog.sh
```

Display the script

``` bash
cat ~/Desktop/Manuscript_1/FIGURE1/scripts/FANC_hic_to_txt_and_bin_nolog.sh
```

```
#!/bin/bash -l
#SBATCH -J FANC_from_hic_to_txt_and_bin_nolog
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --time=02:00:00

# This script is for converting a hic file from 3 hic files into a text file#
# Load modlue containing python3
module load lang/Python/3.8.6-GCCcore-10.2.0

# Activate the python environment for FAN-C
source ~/Environment/FANC/bin/activate

# Check FAN-C version
fanc --version

# Set a TMPDIR variable
export TMPDIR=$SCRATCH/tempfiles/
echo $TMPDIR

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# Run the fan dump command using the 5kb and 10 kb bins for Dennis
# 5kb
fanc dump -tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_5kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_5kb_matrix_noOE_nolog.txt \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_5kb_bins_noOE_nolog.txt

echo "Remove first temporary files previously written"
rm -rf $SCRATCH/tempfiles/*

# 10kb
fanc dump -tmp \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_10kb.hic \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_10kb_matrix_noOE_nolog.txt \
$SCRATCH/LowC_smNPC_THpos_neur/FANC_merged_smNPC/hic/binned/merged_smNPC_10kb_bins_noOE_nolog.txt

deactivate
```


``` r
#################################### Figure 1C LowC - Chr4 ######### text C
plotText(label = "C", fontsize = 12, fontfamily = "Helvetica", x = 0.25, y = 5.25,
    just = "left", default.units = "inches", fontface = "bold")

# Load TH REP1 mCHERRY mDAN neurons D30 (N1) at 50kb resolution
mDAN.D30.50kb = readHic(file = "/home/vagrant/epifunc/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/juicer_like/N1_mDAN_D30.juicer.50kb.hic",
    chrom = "chr4", assembly = "hg38", resolution = 50000, res_scale = "BP", norm = "NONE",
    matrix = "observed")

# Load TH REP1 mCHERRY mDAN neurons D30 (N1) at 25kb resolution
mDAN.D30.25kb = readHic(file = "/home/vagrant/epifunc/LowC_smNPC_THpos_neur/FANC_output_N1_mDAN_D30/hic/juicer_like/N1_mDAN_D30.juicer.25kb.hic",
    chrom = "chr4", assembly = "hg38", resolution = 25000, res_scale = "BP", norm = "NONE",
    matrix = "observed")

# Load TH REP1 mCHERRY smNPCs (N1) at 50kb resolution
smNPC.50kb = readHic(file = "/home/vagrant/epifunc/LowC_smNPC_THpos_neur/FANC_output_N1_smNPC/hic/juicer_like/N1_smNPC.juicer.50kb.hic",
    chrom = "chr4", assembly = "hg38", resolution = 50000, res_scale = "BP", norm = "NONE",
    matrix = "observed")

# Full chr4 to display
region.p.chr4 = pgParams(chrom = "chr4", assembly = assembly(Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), just = c("left",
    "top"), width = 3.5, fontcolor = "black", fill = "black")

# Plot full chr4 - smNPC on top of the diagonale and mDAN D30 at the bottom of
# the diagonale
LowC_smNPC.top = plotHicSquare(data = smNPC.50kb, params = region.p.chr4, zrange = c(0,
    5), resolution = "auto", half = "top", x = 0.5, y = 5.5, height = 3.5)

annoHeatmapLegend(plot = LowC_smNPC.top, fontsize = 7, x = 0.375, y = 6, width = 0.07,
    height = 0.5, just = c("left", "top"), default.units = "inches")

# smNPC
plotText(label = "smNPC", fontsize = 12, fontfamily = "Helvetica", x = 0.5625, y = 5.625,
    just = "left", default.units = "inches", fontface = "bold")

LowC_mDAND30.bottom = plotHicSquare(data = mDAN.D30.50kb, params = region.p.chr4,
    zrange = c(0, 5), resolution = "auto", half = "bottom", x = 0.5, y = 5.5, height = 3.5)

annoHeatmapLegend(plot = LowC_mDAND30.bottom, fontsize = 7, x = 4.0625, y = 6, width = 0.07,
    height = 0.5, just = c("left", "top"), default.units = "inches")

# mDAN D30
plotText(label = "mDAN D30", fontsize = 12, fontfamily = "Helvetica", x = 3, y = 8.9375,
    just = "left", default.units = "inches", fontface = "bold")

# Add genome label of chr4
annoGenomeLabel(plot = LowC_mDAND30.bottom, scale = "Mb", axis = "x", x = 0.5, y = 9.0625,
    just = c("left", "top"))

# Annotate SNCA domain that will be zoomed in fid1D
SNCA.TAD = GRanges("chr4", ranges = IRanges(start = 8.9e+07, end = 9.1e+07))

domainAnno = annoDomains(plot = LowC_mDAND30.bottom, data = SNCA.TAD, half = "bottom",
    linecolor = "red")

################################### Figure 1D SNCA locus ######### Zoom on SNCA
################################### an add ATACseq data of the locus Define
################################### parameters for a small part of chr4 where
################################### SNCA is Add ATACseq track of smNPC, mDAN
################################### D15, mDAN D30, mDAN D50 text D
plotText(label = "D", fontsize = 12, fontfamily = "Helvetica", x = 4.125, y = 5.25,
    just = "left", default.units = "inches", fontface = "bold")

bw.path = "/home/vagrant/epifunc/"

# smNPC N1
bw.SNCA.smNPC.1 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_I.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# smNPC N2
bw.SNCA.smNPC.2 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_II.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# smNPC N3
bw.SNCA.smNPC.3 = readBigwig(file = paste0(bw.path, "BIGWIG/smNPC_III.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D15 N1
bw.SNCA.mDAN.D15.1 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_I.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D15 N2
bw.SNCA.mDAN.D15.2 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_II.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D15 N3
bw.SNCA.mDAN.D15.3 = readBigwig(file = paste0(bw.path, "BIGWIG/D15_POS_III.bw"),
    chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D30 N1
bw.SNCA.mDAN.D30.1 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_I.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D30 N2
bw.SNCA.mDAN.D30.2 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_II.bw"), chrom = "chr4",
    chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D30 N3
bw.SNCA.mDAN.D30.3 = readBigwig(file = paste0(bw.path, "BIGWIG/D30_POS_III.bw"),
    chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D50 N1
bw.SNCA.mDAN.D50.1 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S1.bw"),
    chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D50 N2
bw.SNCA.mDAN.D50.2 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S3.bw"),
    chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07)

# mDAN D50 N3
bw.SNCA.mDAN.D50.3 = readBigwig(file = paste0(bw.path, "BIGWIG/HFFTHmCherry_D50_possort_S4.bw"),
    chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07)

# Add a scale next to the bigwig files - check the maximum to choose
scale.SNCA.max = max(c(bw.SNCA.smNPC.1$score, bw.SNCA.smNPC.2$score, bw.SNCA.smNPC.3$score,
    bw.SNCA.mDAN.D15.1$score, bw.SNCA.mDAN.D15.2$score, bw.SNCA.mDAN.D15.3$score,
    bw.SNCA.mDAN.D30.1$score, bw.SNCA.mDAN.D30.2$score, bw.SNCA.mDAN.D30.3$score,
    bw.SNCA.mDAN.D50.1$score, bw.SNCA.mDAN.D50.2$score, bw.SNCA.mDAN.D50.3$score))

print(scale.SNCA.max)

# Define a small region of chr4
region.p.chr4.small = pgParams(chrom = "chr4", chromstart = 8.9e+07, chromend = 9.1e+07,
    assembly = assembly(Genome = "hg38refGene", TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
        OrgDb = "org.Hs.eg.db"), just = c("left", "top"), width = 3.5, fontcolor = "black",
    fill = "black", range = c(0, scale.SNCA.max))


# Plot Low-C data SNCA in mDAN D30 (10kb resolution)
chr4.SNCA.sq = plotHicTriangle(data = mDAN.D30.25kb, params = region.p.chr4.small,
    height = 1.75, resolution = 25000, x = 4.5, y = 5.5, zrange = c(0, 5), just = c("left",
        "top"), default.units = "inches", palette = colorRampPalette(brewer.pal(n = 9,
        "YlGnBu")))


# ATACseq signal smNPC_I and scale
ATAC_smNPC_I = plotSignal(data = bw.SNCA.smNPC.1, params = region.p.chr4.small, fill = "#313695",
    alpha = 0.7, linecolor = NA, x = 4.5, y = 7.25, height = 0.25, just = c("left",
        "top"), default.units = "inches")

# ATACseq signal smNPC_II and scale
ATAC_smNPC_II = plotSignal(data = bw.SNCA.smNPC.2, params = region.p.chr4.small,
    fill = "#313695", alpha = 0.6, linecolor = NA, x = 4.5, y = 7.25, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal smNPC_III and scale
ATAC_smNPC_III = plotSignal(data = bw.SNCA.smNPC.3, params = region.p.chr4.small,
    fill = "#313695", alpha = 0.5, linecolor = NA, x = 4.5, y = 7.25, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_smNPC_III, at = c(0, scale.SNCA.max), fontsize = 6)

plotText(label = "smNPC", fontsize = 6, fontcolor = "#313695", fontfamily = "Helvetica",
    x = 4.125, y = 7.3125, just = c("left", "top"), default.units = "inches", fontface = "bold")

# ATACseq signal mDAN D15 I and scale
ATAC_mDAN.D15_I = plotSignal(data = bw.SNCA.mDAN.D15.1, params = region.p.chr4.small,
    fill = "#053061", alpha = 0.7, linecolor = NA, x = 4.5, y = 7.625, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D15 II and scale
ATAC_mDAN.D15_II = plotSignal(data = bw.SNCA.mDAN.D15.2, params = region.p.chr4.small,
    fill = "#053061", alpha = 0.6, linecolor = NA, x = 4.5, y = 7.625, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D15 III and scale
ATAC_mDAN.D15_III = plotSignal(data = bw.SNCA.mDAN.D15.3, params = region.p.chr4.small,
    fill = "#053061", alpha = 0.5, linecolor = NA, x = 4.5, y = 7.625, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_mDAN.D15_III, at = c(0, scale.SNCA.max), fontsize = 6)

plotText(label = "mDAN D15", fontsize = 6, fontcolor = "#053061", fontfamily = "Helvetica",
    x = 4.125, y = 7.6875, just = c("left", "top"), default.units = "inches", fontface = "bold")

# ATACseq signal mDAN D30 I and scale
ATAC_mDAN.D30_I = plotSignal(data = bw.SNCA.mDAN.D30.1, params = region.p.chr4.small,
    fill = "#2D004B", alpha = 0.7, linecolor = NA, x = 4.5, y = 8, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D30 II and scale
ATAC_mDAN.D30_II = plotSignal(data = bw.SNCA.mDAN.D30.2, params = region.p.chr4.small,
    fill = "#2D004B", alpha = 0.6, linecolor = NA, x = 4.5, y = 8, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D30 III and scale
ATAC_mDAN.D30_III = plotSignal(data = bw.SNCA.mDAN.D30.3, params = region.p.chr4.small,
    fill = "#2D004B", alpha = 0.5, linecolor = NA, x = 4.5, y = 8, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_mDAN.D30_III, at = c(0, scale.SNCA.max), fontsize = 6)

plotText(label = "mDAN D30", fontsize = 6, fontcolor = "#2D004B", fontfamily = "Helvetica",
    x = 4.125, y = 8.125, just = c("left", "top"), default.units = "inches", fontface = "bold")

# ATACseq signal mDAN D50 I and scale
ATAC_mDAN.D50_I = plotSignal(data = bw.SNCA.mDAN.D50.1, params = region.p.chr4.small,
    fill = "#003C30", alpha = 0.7, linecolor = NA, x = 4.5, y = 87.375, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D50 II and scale
ATAC_mDAN.D50_II = plotSignal(data = bw.SNCA.mDAN.D50.2, params = region.p.chr4.small,
    fill = "#003C30", alpha = 0.6, linecolor = NA, x = 4.5, y = 8.375, height = 0.25,
    just = c("left", "top"), default.units = "inches")

# ATACseq signal mDAN D50 III and scale
ATAC_mDAN.D50_III = plotSignal(data = bw.SNCA.mDAN.D50.3, params = region.p.chr4.small,
    fill = "#003C30", alpha = 0.5, linecolor = NA, x = 4.5, y = 8.375, height = 0.25,
    just = c("left", "top"), default.units = "inches")

annoYaxis(plot = ATAC_mDAN.D50_III, at = c(0, scale.SNCA.max), fontsize = 6)

plotText(label = "mDAN D50", fontsize = 6, fontcolor = "#003C30", fontfamily = "Helvetica",
    x = 4.125, y = 8.5, just = c("left", "top"), default.units = "inches", fontface = "bold")

# Add gene name and genome labels for the chr4 Add gene name
SNCA_g = plotGenes(params = region.p.chr4.small, x = 4.5, y = 8.8125, height = 0.5,
    geneHighlights = data.frame(gene = c("SNCA"), color = "red"), assembly = assembly(Genome = "hg38refGene",
        TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), geneBackground = "black")


# Add genome label
annoGenomeLabel(plot = SNCA_g, scale = "Mb", axis = "x", x = 4.5, y = 8.6875, just = c("left",
    "top"))

pageGuideHide()
dev.off()
```
