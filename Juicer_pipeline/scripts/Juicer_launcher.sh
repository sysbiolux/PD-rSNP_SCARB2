#!/bin/bash -l
#SBATCH -J Juicer_pipeline_fromBAM
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH -c 28
#SBATCH --time=144:00:00
#SBATCH -p bigmem
#SBATCH --qos=iris-bigmem-long

# Activate the python environment for Juicer
source $HOME/Environment/JUICER/bin/activate

# Load needed modules
module load env/legacy/2020b
module load bio/SAMtools/1.12-GCC-10.2.0
module load bio/BWA/0.7.17-GCC-10.2.0

# Go to SCRATCH folder for the run
cd $SCRATCH

# And run
# smNPC N1
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_smNPC_01/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28

# smNPC N2
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_smNPC_02/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28

# smNPC N3
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_smNPC_03/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28

# mDAN D30 N1
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_mDAN.D30_01/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28

# mDAN D30 N2
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_mDAN.D30_02/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28

# mDAN D30 N3
$SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
-d $SCRATCH/LowC_mDAN.D30_03/ -p $SCRATCH/hg38.chrom.sizes \
-s MboI \
-y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
-z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
-t 28



deactivate