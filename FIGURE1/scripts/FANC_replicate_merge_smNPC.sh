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
