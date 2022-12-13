FROM bioconductor/bioconductor_docker:RELEASE_3_16

COPY requirements.txt /tmp


RUN R -q -e "install.packages(c('ggplot2', 'tidyverse', 'scales', \ 
    'rcolorbrewer', 'wesanderson', 'jpeg', 'rjson', 'torch', 'luz', 'cowplot', \
    'Matrix', 'LowRankQP', 'umap', 'irlba', 'Seurat', 'ggnewscale', 'ggpubr', \
    'ggraph', 'tidygraph')); \
    # finishing torch installation
    torch::install_torch(); \
    # bioconductor packages
    BiocManager::install(c('S4Vectors', 'DropletUtils', 'HDF5Array','rhdf5','Rhdf5lib', 'scuttle', 'IRanges', 'SummarizedExperiment'));\
    # github packages
    devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')"

# removing temporary folders
RUN rm -rf /tmp/downloaded_packages && rm -rf /tmp/RtmpIj1cu8/downloaded_packages

# update-upgrade | installing docker and CUDA
RUN apt update && apt -y upgrade && apt -y install docker.io nvidia-cuda-toolkit

RUN pip install -r /tmp/requirements.txt

ENV DISABLE_AUTH=true
