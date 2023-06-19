
FROM rocker/ml:4.2.2

RUN apt-get update \
    && apt-get install -y \
      wget \
      bzip2 \
      ca-certificates \
      sudo \
      locales \
      fonts-liberation \
      run-one \
      pigz \
      docker.io \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN usermod -aG docker rstudio

# installing packages as user, not root
USER rstudio
RUN pip install --upgrade pip \
   && pip install --upgrade jupyterlab \
      jupyterlab-spreadsheet-editor \
      jupyterlab-github \
      jupyterlab-system-monitor \
      jupyter-rsession-proxy 

COPY install.R /tmp
#RUN R -q -f "/tmp/install.R"

EXPOSE 8888
CMD ["/home/rstudio/.local/bin/jupyter", "lab", "--port=8888", "--ip=0.0.0.0", "--allow-root"]


