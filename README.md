# Identification of Parkinson’s disease-associated regulatory variants in human dopaminergic neurons reveals modulators of SCARB2 and BAG3 expression
Please cite the following reference: 

All the scripts have been run into a Docker container based on [`r2u`](https://eddelbuettel.github.io/r2u/).   
This repository contains the code and data to reproduce the figures of the manuscript.
Every figure has its own folder with the related data to it.  

## Docker container

- To build the image from the folder containing the Dockerfile: `docker buildx build -t 250603:24.04 .`

Of note, ensure that `renv` is not invoked if you use it: `mv .Rprofile Rprofile`

- To start a container from the image, binding the local directory: `docker run -u $(id -u):$(id -g) -ti -v `pwd`:/mnt 250603:24.04`

Once in the container, for Figure3 for example: `cd /mnt/FIGURE3 ;  Rscript --vanilla Figure3.R`

## Data availability

Sequencing of epigenetics data can be found at [EGA-archive](https://ega-archive.org) under controlled access.

- [EGAD00001009288](https://ega-archive.org/datasets/EGAD00001009288), associated [publication](https://europepmc.org/abstract/MED/38177910).
    + ATAC/ChIP
    + RNA-seq
- [EGAD50000002258](https://ega-archive.org/datasets/EGAD50000002258) which contains new sequencing data:
    + 3D chromatin contact using LowC, 6 samples. 3 replicates for each condition: smNPC (controls) cells and differentiated neurons after 30 days. Resulting paired-end are provided as BAM files.
    + In order to confirm the predicted impact of the PD-associated allele of rs144814361 on BAG3 promoter, we proceeded with genome editing of the TH-REP1 cell line by using prime editing to insert the “T” allele at the position chr10:119651405 in the BAG3 promoter. Chromatin accessibility using ATAC-seq, 9 samples of derived cell line TH-REP1. 3 replicates for each condition: 1) iPSC Wild Type 2) SNP-BAG3 variant in iPSC 3) SNP-BAG3 variant in smNPC. BAG3 variant: heterozygous for the SNP variant: rs144814361 (chr10 119651405). Resulting sequencing are provided as paired-end FASTQ files.

## Acknowledgements

- The sequencing was performed at the Genomics Platform of Luxembourg Centre for System Biomedicine.
- The computational analysis presented in this paper were carried out using the HPC facilities of the University of Luxembourg https://hpc.uni.lu. 
- Figures created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
- Figure 1A was created with [BioRender](https://www.biorender.com/)  using license Catillon, M. (2025) https://BioRender.com/5kzs36o.

## References


- R Core Team (2025). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing,
  Vienna, Austria. https://www.R-project.org.
- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen
  TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H
  (2019). "Welcome to the tidyverse." _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 https://doi.org/10.21105/joss.01686.
- Nicole E Kramer, Eric S Davis, Craig D Wenger, Erika M Deoudes, Sarah M Parker, Michael I Love, Douglas H Phanstiel, Plotgardener: cultivating precise multi-panel figures in R, _Bioinformatics_, 2022.
- Jesper R Gådin, Ferdinand van't Hooft, Per Eriksson and Lasse Folkersen (2015): AllelicImbalance: an R/bioconductor package
  for detecting, managing, and visualizing allele expression imbalance data from RNA sequencing. BMC _Bioinformatics_.
- Kassambara A (2025). _ggpubr: 'ggplot2' Based Publication Ready Plots_. R package version 0.6.2,
  https://rpkgs.datanovia.com/ggpubr.
