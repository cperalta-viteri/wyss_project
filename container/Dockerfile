FROM bioconductor/bioconductor_docker:RELEASE_3_16

# update-upgrade | installing docker and CUDA
RUN apt-get update && apt-get install -y \
    nvidia-cuda-toolkit \
    docker.io && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

RUN R -q -e " \
    install.packages(c('Matrix','curl','XML','rjson','jsonlite','rhdf5','igraph','reshape2','Rcpp')); \
    install.packages(c('LowRankQP','bmrm','limma')); \
    install.packages(c('scales','RColorBrewer','wesanderson')); \
"
RUN R -q -e " \
    install.packages(c('ggplot2','patchwork','ggrepel','ggraph','cowplot','ggnewscale','ggpubr')); \
    install.packages(c('jpeg','png')); \
    install.packages(c('irlba', 'umap', 'Rtsne')); \
"
RUN R -q -e " \
    install.packages(c('torch','luz'));torch::install_torch(); \
    install.packages(c('tidyverse','tidygraph')); \
    install.packages(c('reticulate','Seurat')); \
"
RUN R -q -e " \
    # bioconductor packages \
    BiocManager::install(c('ShortReads','Rsamtools','rtracklayer')); \
    BiocManager::install(c('S4Vectors','IRanges','GenomicRanges','GenomicFeatures','GenomicAlignments','SummarizedExperiment')); \
    BiocManager::install(c('Rhdf5lib','rhdf5','HDF5Array')); \
    BiocManager::install(c('SingleCellExperiment','SpatialExperiment','LoomExperiment')); \
    BiocManager::install(c('DropletUtils','scuttle','scran','scater','uwot','AUCell','SingleR')); \
    BiocManager::install(c('edgeR','GEOquery','zellkonverter','sceasy')); \
\
    # github packages \
    devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always') \
"

# removing temporary folders
#RUN rm -rf /tmp/downloaded_packages && rm -rf /tmp/RtmpIj1cu8/downloaded_packages


#COPY requirements.txt /tmp
#RUN pip install -r /tmp/requirements.txt

ENV DISABLE_AUTH=true
