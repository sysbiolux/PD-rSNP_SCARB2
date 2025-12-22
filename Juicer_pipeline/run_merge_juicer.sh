#!/bin/bash

set -eu

if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters, folder_name, 1 or 30" >&2
    exit 2
fi
 

outputdir="$1"
qual="$2" # 1 or 30

echo "Output $outputdir for quality $qual"

juicer_tools statistics juicer/restriction_sites/hg38_MboI.txt \
  $outputdir/inter_${qual}.txt $outputdir/merged${qual}.txt hg38

# indexing for multi-threading, probably useless as not working
time index_by_chr.awk ${outputdir}/merged${qual}.txt 500000 > ${outputdir}/merged${qual}_index.txt
# multi-threaded pre is NOT working, use one thread only
time juicer_tools pre -n -f juicer/restriction_sites/hg38_MboI.txt -s ${outputdir}/inter_${qual}.txt \
  -g ${outputdir}/inter_${qual}_hists.m -q ${qual} -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 \
  --threads 1 -i ${outputdir}/merged${qual}_index.txt -t /dev/shm \
  ${outputdir}/merged${qual}.txt ${outputdir}/inter_${qual}.hic hg38

time juicer_tools addNorm --threads 10 ${outputdir}/inter_${qual}.hic 

