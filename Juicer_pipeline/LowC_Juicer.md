---
title: "LowC_Juicer_AllRep"
author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "08 July, 2025"
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



#### *1. Install Juicer (CPU version) on HPC (Aion)*

``` bash
# Connect to HPC 
ssh aion-cluster

# Book some ressources for the installation
si
```


``` bash
# Load modules to use Python3
module load env/development/2024a
module load lang/Python/3.12.3-GCCcore-13.3.0

# Create the environment in my home directory with ther other environments
python3 -m venv ~/Environment/JUICER

# Activate the newly created environment
#source $HOME/Environment/JUICER/bin/activate

# And install Juicer according to Github
cd $SCRATCH
git clone https://github.com/theaidenlab/juicer.git
cd $SCRATCH/juicer

ln -s $SCRATCH/juicer/CPU scripts
cd scripts/common
#wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
#ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
cd ../..
mkdir references
cd references
mv $SCRATCH/bwa_index/ references/
mv  -v $SCRATCH/juicer/references/bwa_index/* references/

# Generate site position of the MboI enzyme on hg38
# Book some ressources
salloc -p interactive -N 1 -c 4 --qos debug -C batch -t 01:45:00

# Unzip (and then rezip) the hg38 pastch 1 fasta file and create the restriction file of the enzyme of interest
gunzip $SCRATCH/juicer/references/GRCh38.genome.fa.gz
python $SCRATCH/juicer/misc/generate_site_positions.py MboI hg38 $SCRATCH/juicer/references/GRCh38.genome.fa
gzip $SCRATCH/juicer/references/GRCh38.genome.fa

cd $SCRATCH/juicer
# Juicer must have a specific folder structure. Create restriction sites folder
mkdir restriction_sites 

# Move the restriction file to the right folder
mv $SCRATCH/hg38_MboI.txt $SCRATCH/juicer/restriction_sites

# Make chromosome size from the restriction file (as advised by juicer)
awk 'BEGIN{OFS="\t"}{print $1, $NF}' $SCRATCH/juicer/restriction_sites/hg38_MboI.txt > $SCRATCH/hg38.chrom.sizes

# Move the juicer tool jar previsouly downloaded to the script directory (under common)
mv $HOME/juicer_tools.2.20.00.jar $SCRATCH/juicer/CPU
mv $SCRATCH/juicer/scripts/juicer_tools.2.20.00.jar $SCRATCH/juicer/scripts/common
ln -s $SCRATCH/juicer/scripts/common/juicer_tools.2.20.00.jar  juicer_tools.jar

# this is optional, only needed for fragment-delimited files (this is my case)
ln -s $SCRATCH/juicer/restriction_sites restriction_sites

# The fastq folder needed by juicer has to be top directory. Move the necessary fastq files there
cd $SCRATCH
mkdir fastq

# Move smNPC samples first
mv /scratch/users/dgerard/LowCO/LowCO/fastq/*smNPC* $SCRATCH/fastq/
mv /scratch/users/dgerard/LowC_smNPC_THpos_neur/fastq/*smNPC* $SCRATCH/fastq/

# Move mDAN D30 now
mv /scratch/users/dgerard/LowC_smNPC_THpos_neur/fastq/*day30* $SCRATCH/fastq/
mv /scratch/users/dgerard/LowC2/fastq/*D30* $SCRATCH/fastq/

# Change the permissions
chmod 750 N4*
chmod 750 N3*
```

#### *2. Run Juicer*
Run the batch script for juicer

``` bash
# Transfer to HPC
rsync -avzPhu /Volumes/deborah.gerard/Documents/epifunc/LowC_Juicer/scripts/Juicer_launcher.sh iris-cluster:/scratch/users/dgerard/

# Run
sbatch $SCRATCH/Juicer_launcher.sh

# Take the slurm log
rsync -avzPhu iris-cluster:/scratch/users/dgerard/slurm-4176155.out ~/Desktop/
```

Every samples were merged and I had a out of memory error. Make directory of biological replicates per samples and rerun Juicer per sample

``` bash
# on Iris in $SCRATCH
salloc -p interactive -N 1 -c 4 --qos debug -C batch -t 00:45:00

cd $SCRATCH
# For 3 biological replicate of smNPCs
mkdir -p LowC_smNPC_01/fastq
mkdir -p LowC_smNPC_02/fastq
mkdir -p LowC_smNPC_03/fastq

# For 3 biological replicate of mDAN D30
mkdir -p LowC_mDAN.D30_01/fastq
mkdir -p LowC_mDAN.D30_02/fastq
mkdir -p LowC_mDAN.D30_03/fastq

# Move fastq files to respective folder
mv ./fastq/N1_TH_REP1_mCHERRY_smNPC_LOW_C* $SCRATCH/LowC_smNPC_01/fastq/
mv ./fastq/N2_TH_REP1_mCHERRY_smNPC_LOW_C* $SCRATCH/LowC_smNPC_02/fastq/
mv ./fastq/Optimisation* $SCRATCH/LowC_smNPC_03/fastq/

mv ./fastq/N1_TH_REP1_mCHERRY_TH_Neur_day30* $SCRATCH/LowC_mDAN.D30_01/fastq/
mv ./fastq/N3_TH_REP1_mCHERRY_THpos_NEUR_D30* $SCRATCH/LowC_mDAN.D30_02/fastq/
mv ./fastq/N4_TH_REP1_mCHERRY_THpos_NEUR_D30* $SCRATCH/LowC_mDAN.D30_03/fastq/

# Transfer the fastq files for the 1st biological replicate of smNPC
rsync -avzPhu /Volumes/deborah.gerard/Documents/epifunc/LowC_smNPC_THpos_neur/fastq/N1_TH_REP1_mCHERRY_smNPC_LOW_C* iris-cluster:/scratch/users/dgerard/LowC_smNPC_01/fastq/

# Juicer launch script transfer
rsync -avzPhu /Volumes/deborah.gerard/Documents/epifunc/LowC_Juicer/scripts/Juicer_launcher.sh iris-cluster:/scratch/users/dgerard/

# Launch
sbatch $SCRATCH/Juicer_launcher.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/epifunc/LowC_Juicer/scripts/Juicer_launcher.sh
```

```
## #!/bin/bash -l
## #SBATCH -J Juicer_pipeline_fromBAM
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH -c 28
## #SBATCH --time=144:00:00
## #SBATCH -p bigmem
## #SBATCH --qos=iris-bigmem-long
## 
## # Activate the python environment for Juicer
## source $HOME/Environment/JUICER/bin/activate
## 
## # Load needed modules
## module load env/legacy/2020b
## module load bio/SAMtools/1.12-GCC-10.2.0
## module load bio/BWA/0.7.17-GCC-10.2.0
## 
## # Go to SCRATCH folder for the run
## cd $SCRATCH
## 
## # And run
## # smNPC N1
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_smNPC_01/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## # smNPC N2
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_smNPC_02/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## # smNPC N3
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_smNPC_03/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## # mDAN D30 N1
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_mDAN.D30_01/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## # mDAN D30 N2
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_mDAN.D30_02/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## # mDAN D30 N3
## $SCRATCH/juicer/scripts/juicer.sh -D $SCRATCH/juicer/ \
## -d $SCRATCH/LowC_mDAN.D30_03/ -p $SCRATCH/hg38.chrom.sizes \
## -s MboI \
## -y $SCRATCH/juicer/restriction_sites/hg38_MboI.txt \
## -z $SCRATCH/juicer/references/GRCh38.genome.fa.gz \
## -t 28
## 
## 
## 
## deactivate
```

Download the slurm output

``` bash
rsync -avzPhu iris-cluster:/scratch/users/dgerard/slurm-4187610.out ~/Desktop/LowC_Juicer/
```

Copy the results in /work/projects/lowc_mdan

``` bash
# Book some ressources
si -t 01:05:00
# For smNPC N1
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_smNPC_01/aligned $SCRATCH/LowC_smNPC_01/splits /work/projects/lowc_mdan/LowC_smNPC_01/

# For smNPC N2
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_smNPC_02/aligned $SCRATCH/LowC_smNPC_02/splits /work/projects/lowc_mdan/LowC_smNPC_02/

# For smNPC N3
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_smNPC_03/aligned $SCRATCH/LowC_smNPC_03/splits /work/projects/lowc_mdan/LowC_smNPC_03/

# For mDAN D30 N1
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_mDAN.D30_01/aligned $SCRATCH/LowC_mDAN.D30_01/splits /work/projects/lowc_mdan/LowC_mDAN.D30_01/

# For mDAN D30 N2
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_mDAN.D30_02/aligned $SCRATCH/LowC_mDAN.D30_02/splits /work/projects/lowc_mdan/LowC_mDAN.D30_02/

# For mDAN D30 N3
rsync -avzP --no-p --no-g --chmod=ug=rwX $SCRATCH/LowC_mDAN.D30_03/aligned $SCRATCH/LowC_mDAN.D30_03/splits /work/projects/lowc_mdan/LowC_mDAN.D30_03/
```

#### *3. Merge 3 biological replicates of smNPC into a mega map. Do the same for the mDAN D30*

Create a juicer apptainer image following https://github.com/aidenlab/juicer/tree/main/Docker:

don't mount home but bind the current project folder in the container

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
