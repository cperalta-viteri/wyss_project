# take this image to create a non-cuda environment
#FROM bioconductor/bioconductor_docker:RELEASE_3_16

FROM nvidia/cuda:12.0.0-base-ubuntu20.04

RUN TZ=Europe/Madrid && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && apt-get install -y \
    software-properties-common aptitude cuda-toolkit-11-6

# adding R key
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# adding R repo
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# adding python repo to install python
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update &&\
    apt-get -y install python3.10 python3-pip

# installing R 4.2
RUN aptitude -y install r-base r-base-dev build-essential

#creating user for autamate login into shell 
### testing fase
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su &&\
    sed 's/required/sufficient/g' -i /etc/pam.d/chsh &&\
    useradd -ms /usr/bin/bash -p "psswd" user

ENV HOME=/home/user

# installing rstudio server
RUN apt-get -y install gdebi-core wget &&\
    wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2022.12.0-353-amd64.deb -P $HOME/Downloads/ &&\
    gdebi -n $HOME/Downloads/rstudio-server-2022.12.0-353-amd64.deb

COPY jupyter_lab_config.py /home/user/.jupyter/jupyter_lab_config.py
COPY cudnn-local-repo-ubuntu2004-8.7.0.84_1.0-1_amd64.deb $HOME/Downloads
RUN gdebi -n $HOME/Downloads/cudnn-local-repo-ubuntu2004-8.7.0.84_1.0-1_amd64.deb


RUN chown -R user $HOME && chgrp -R user $HOME


#installing packages as user, not root
USER user
RUN pip install jupyterlab jupyter-rsession-proxy
RUN echo 'export PATH="$HOME/.local/bin:$PATH"' >> $HOME/.bashrc
RUN echo 'export PATH=/usr/local/cuda-11.6/bin${PATH:+:${PATH}}' >> $HOME/.bashrc
RUN pip install charset-normalizer==2.1.1 \
    jsonpointer fqdn isoduration webcolors uri-template



USER root
# correcting locales
RUN aptitude install -y locales locales-all
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

#installing missing dependencies
RUN aptitude -y install \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    cmake \
    libssl-dev \
    libxml2-dev \
    libmagick++-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libpoppler-cpp-dev
    
RUN R -q -e " \
    install.packages(c('Matrix','curl','XML','rjson','jsonlite','rhdf5','igraph','reshape2','Rcpp')); \
    install.packages(c('LowRankQP','bmrm','limma')); \
    install.packages(c('scales','RColorBrewer','wesanderson')); \
    "

RUN R -q -e " \
    install.packages(c('ggplot2','patchwork','ggrepel','ggraph','cowplot','ggnewscale'));"


RUN R -q -e "\
    install.packages(c('jpeg','png','ggpubr')); \
    install.packages(c('irlba', 'umap', 'Rtsne')); \
    "

# setting environmental variables (add before, at line 62)
ENV PATH="/usr/local/cuda-11.6/bin${PATH:+:${PATH}}"
ENV PATH="$HOME/.local/bin:$PATH"

RUN R -q -e " \
    install.packages(c('torch','luz'));torch::install_torch();"


RUN R -q -e "\
    install.packages(c('tidyverse','tidygraph')); \
    install.packages(c('reticulate','Seurat')); \
    "


RUN R -q -e " \
    # bioconductor packages \
    install.packages(c('BiocManager'));\
    BiocManager::install(version='3.16');\
    BiocManager::install(c('ShortReads','Rsamtools','rtracklayer')); \
    BiocManager::install(c('S4Vectors','IRanges','GenomicRanges','GenomicFeatures','GenomicAlignments','SummarizedExperiment')); \
    BiocManager::install(c('Rhdf5lib','rhdf5','HDF5Array')); \
    BiocManager::install(c('SingleCellExperiment','SpatialExperiment','LoomExperiment')); \
    BiocManager::install(c('DropletUtils','scuttle','scran','scater','uwot','AUCell','SingleR')); \
    BiocManager::install(c('biomaRt')); \
    BiocManager::install(c('edgeR','GEOquery','zellkonverter','sceasy')); \
    BiocManager::install(c('splatter'))"

RUN R -q -e "BiocManager::install(c('demuxmix'))"

RUN R -q -e "\
    install.packages(c('devtools')); \
    devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always'); \
    devtools::install_github('cellgeni/sceasy')\
    "


USER user
COPY requirements.txt /tmp
RUN pip install -r /tmp/requirements.txt


USER root
EXPOSE 8888


# removing cache files
RUN rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
# removing temporary folders
RUN rm -rf /tmp/downloaded_packages && \
    rm -rf /tmp/RtmpIj1cu8/downloaded_packages && \
    rm -rf /tmp/RtmpmxTsSH/downloaded_packages

USER user

WORKDIR $HOME
CMD ["/home/user/.local/bin/jupyter", "lab"]


