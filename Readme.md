# Quick Start

1. Install docker for your OS <https://www.docker.com>

2. In a terminal
```
docker run --rm -p8888:8888 unigebsp/neurosuite
```

3. In a web browser open the URL http://localhost:8888/


# Build the working environment

1.  Clone the repository and launch the build command
```{bash}
git clone git@github.com:Christ14n97/wyss_project.git
docker-compose run --rm build
```

# Rocket Project

Currently, it only has RStudio as IDE. Python3 is available to use, though.

Build the image:

```bash
docker build --pull --rm -f "container/Dockerfile" -t neuro_suite:latest "container"
```

Run container

```bash
# this is how I mount my directory on my container using simple repo2docker command
docker run -it --rm -p 8787:8787 -v $PWD:/home/rstudio/workdir neuro_suite:latest 
```

```bash
# STILL TO TEST: how to mount a config file to save customization on RS IDE.
docker run -it --rm -p 8787:8787 -v $PWD/.config:/user/rstudio/.config -v $PWD:/home/rstudio/workdir dockerenvi:latest 
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
