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
