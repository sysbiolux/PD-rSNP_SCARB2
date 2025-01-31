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
