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
