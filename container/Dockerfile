FROM bioconductor/bioconductor_docker:RELEASE_3_16


WORKDIR /usr

ENV DISABLE_AUTH=true

COPY requirements.txt /tmp


RUN R -q -e "install.packages(c('ggplot2', 'tidyverse', 'scales', 'seurat', \ 
    'rcolorbrewer', 'wesanderson', 'jpeg', 'rjson', 'torch', 'luz', 'cowplot', \
    'Matrix', 'LowRankQP', 'umap', 'irlba')); \
    torch::install_torch()"

RUN rm -rf /tmp/downloaded_packages

RUN apt -y install docker.io

RUN pip install --no-cache-dir -r /tmp/requirements.txt