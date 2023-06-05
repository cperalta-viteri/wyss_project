# Introduction

This repository contains R-based functions to help in accessing and analyzing big
brain -omics datasets. The repository is organized this way:

- `Rpackages/`: contains R packages 
- `Rpackages/ABAData`: is a R-package with functions to access Allen Brain Atlas reference datasets.
- `Rpackages/scml`: is a R-package with functions to train deep-learning models on single-cell datasets.
- `apps/`: contains shiny apps
- `apps/geo_download`: is a shiny app to search GEO and download GEO datasets.
- `apps/standardize_metadata`: is a shiny app to standardize metadata of single-cell experiments.
- `container/`: contains a description of a docker container to work on single-cell -omics.


See pakages vignettes for more documentations on the packages.


# R-packages installation
```
R CMD INSTALL Rpackages/scml
R CMD INSTALL Rpackages/ABAData
```


# TODO

- Add ISH data to ABAData
- Create `data/` with ABA datasets
- Build `.tar.gz` files of the packages
- `standardize_metadata` shiny app
- Example code to create a SingleCellExperiment from 10x or CSV



# Quick Start

  1. Search and download a dataset from GEO using the shiny app GUI `apps/geo_download`:
```
shiny::runApp("apps/geo_download/")
```

The dataset you have downloaded with the app will be made available into directory `apps/geo_downloads/downloads/`. 
At this step, depending on the format of the data deposited on GEO, you have to implement a custom R-script to produce `SingleCellExperiment` objects from the downloaded datasets. To help you in this step, you can reuse example scripts available into `apps/geo_downloads/data-examples/`.





```
shiny::runApp("apps/standardize_metadata/")
```

Run

the download shiny app

<!--
1. Install docker for your OS <https://www.docker.com>

2. In a terminal


```

docker run -it --rm --gpus all -p 8888:8888 -v $PWD:/home/user/workdir unigebsp/neurosuite:cuda11.6_focal_R4.2_BioC3.16_Py3.10_latest

```

3. In a web browser open the URL http://localhost:8888/


# Build the working environment

1.  Clone the repository and launch the build command
```{bash}
git clone git@github.com:Christ14n97/wyss_project.git
docker-compose run --rm build
```


```bash
# STILL TO TEST: how to mount a config file to save customization on RS IDE.
docker run -it --rm -p 8787:8787 -v $PWD/.config:/user/rstudio/.config -v $PWD:/home/user/workdir dockerenvi:latest 
```


# Run working environment on Baobab 

## Launch the working environment using Singularity
```{bash}
ssh prados@baobab2.hpc.unige.ch

# Convert the docker container to SIF image for singularity
ml GCC/9.3.0 Singularity/3.7.3-Go-1.14
singularity run -e PASSWORD=bioc -p 8787:8787 docker://bioconductor/bioconductor_docker
singularity pull docker://unigebsp/neurosuite

# Launch the server 
srun -n1 --cpus-per-task=8 --time 4:00:00 --partition=shared-bigmem --mem=64G singularity run neurosuite_latest.sif

# Create a tunnel between your workstation and baobab wrker node
ssh -L localhost:8888:cpu013.baobab:8888 prados@baobab2.hpc.unige.ch

# Open your web-browser to http://localhost:8888/

```

## Launch R terminal with GUI forwarding
```{bash}
ssh -Y prados@baobab2.hpc.unige.ch
ml GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1
salloc -n1 --cpus-per-task=8 --time 4:00:00 --partition=shared-bigmem --mem=64G --x11
R
```
-->
