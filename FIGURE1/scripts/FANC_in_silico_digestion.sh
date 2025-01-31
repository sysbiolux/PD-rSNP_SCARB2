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
