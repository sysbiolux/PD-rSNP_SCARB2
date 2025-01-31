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
