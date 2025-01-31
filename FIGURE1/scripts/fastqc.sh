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
