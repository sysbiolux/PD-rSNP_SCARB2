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

## Acjnowledgements

- Figures created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
- Figure 1A and Figure 2E were created with [BioRender](https://www.biorender.com/).
