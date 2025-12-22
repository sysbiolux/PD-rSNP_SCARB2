#!/bin/bash -l 
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=aurelien.ginolhac@uni.lu
#SBATCH -J juicer_mDAN
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 20
#SBATCH --mem=170GB
#SBATCH --time=0-03:00:00
#SBATCH -p batch

module load env/release
module load tools/Apptainer

# Apptainer `--no-home` does not bind home but the current directory.
# However it **cannot resolve links** so replace with real path
cd /mnt/aiongpfs/projects/lowc_mdan

# Done manually for mDNA.D30
#outputdir="mega_mDAN.D30"
outputdir="mega_smNPC"
mkdir -p $outputdir
echo "Working in $outputdir"

# quality 30
merged_names30=$(find LowC_smNPC_* | grep merged30.txt | tr '\n' ' ')

# above sort from mega.sh assumes that individual merged are sorted (-m) which is not true
# applying solution from https://groups.google.com/g/3d-genomics/c/axgvyaDbcLg/m/HpZ1RHyIDAAJ
# sorting individually first
for MERGE in $merged_names30
  do echo $MERGE
  OUT=$(echo $MERGE | awk -F "/" '{print $1}')
  sort -T /dev/shm -k2,2d -k6,6d -k3,3n -k7,7n $MERGE > ${outputdir}/${OUT}_merged30_sort.txt
done

sort  -m -k2,2d -k6,6d -k3,3n -k7,7n ${outputdir}/*_merged30_sort.txt > ${outputdir}/merged30.txt

apptainer exec --no-home juicer.sif ./run_merge_juicer.sh ${outputdir} 30

# quality 1

merged_names=$(find LowC_smNPC_* | grep merged1.txt | tr '\n' ' ')

# above sort from mega.sh assumes that individual merged are sorted (-m) which is not true
# applying solution from https://groups.google.com/g/3d-genomics/c/axgvyaDbcLg/m/HpZ1RHyIDAAJ
# sorting individually first
for MERGE in $merged_names
  do echo $MERGE
  OUT=$(echo $MERGE | awk -F "/" '{print $1}')
  sort -T /dev/shm -k2,2d -k6,6d -k3,3n -k7,7n $MERGE > ${outputdir}/${OUT}_merged1_sort.txt
done

sort  -m -k2,2d -k6,6d -k3,3n -k7,7n ${outputdir}/*_merged1_sort.txt > ${outputdir}/merged1.txt

apptainer exec --no-home juicer.sif ./run_merge_juicer.sh ${outputdir} 1
#
