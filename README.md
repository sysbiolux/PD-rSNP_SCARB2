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

The data can be found at XXX

## Figure 1
Figure 1 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
Figure 1A has been created with [BioRender](https://www.biorender.com/)

## Figure 2
Figure 2 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
Figure 2A has been created with [BioRender](https://www.biorender.com/)

## Figure 3
Figure 3 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 
Figure 3A has been created with [BioRender](https://www.biorender.com/)

## Figure 4
Figure 4 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 

## Supplementary Figure 1
Supplementary Figure 1 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 

## Supplementary Figure 2
Supplementary Figure 2 has been created with the [plotgardener R package](https://phanstiellab.github.io/plotgardener/index.html). 

