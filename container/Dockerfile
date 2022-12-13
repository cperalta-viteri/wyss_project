FROM bioconductor/bioconductor_docker:RELEASE_3_16

# this set the WD to run commands. 
#WORKDIR /rstudio

COPY requirements.txt /tmp


RUN R -q -e "install.packages(c('ggplot2', 'tidyverse', 'scales', \ 
    'rcolorbrewer', 'wesanderson', 'jpeg', 'rjson', 'torch', 'luz', 'cowplot', \
    'Matrix', 'LowRankQP', 'umap', 'irlba')); \
    torch::install_torch()"

RUN R -q -e "install.packages(c('Seurat', 'ggnewscale'))"

RUN R -q -e "BiocManager::install(c('S4Vectors', 'DropletUtils', 'HDF5array', \
    'rhdf5','Rhdf5lib', 'scuttle', 'IRanges', 'SummarizedExperiment'))"

RUN R -q -e "devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')"

RUN rm -rf /tmp/downloaded_packages && rm -rf /tmp/RtmpIj1cu8/downloaded_packages

RUN apt update && apt -y upgrade && apt -y install docker.io

RUN pip install -r /tmp/requirements.txt

RUN R -q -e "install.packages(c('ggpubr', 'ggraph', 'tidygraph'))"

ENV DISABLE_AUTH=true

## BioConductor
# cellkonverter edgeR 