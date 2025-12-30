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

Sequencing of epigenetics data can be found at [EGA-archive](https://ega-archive.org/datasets/EGAD00001009288) under controlled access.

## Acknowledgements

- Figures created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
- Figure 1A was created with [BioRender](https://www.biorender.com/).

## References


- R Core Team (2025). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing,
  Vienna, Austria. https://www.R-project.org.
- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen
  TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H
  (2019). "Welcome to the tidyverse." _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 https://doi.org/10.21105/joss.01686.
- Nicole E Kramer, Eric S Davis, Craig D Wenger, Erika M Deoudes, Sarah M Parker, Michael I Love, Douglas H Phanstiel, Plotgardener: cultivating precise multi-panel figures in R, _Bioinformatics_, 2022.
- Jesper R Gådin, Ferdinand van't Hooft, Per Eriksson and Lasse Folkersen (2015): AllelicImbalance: an R/bioconductor package
  for detecting, managing, and visualizing allele expression imbalance data from RNA sequencing. BMC _Bioinformatics_.