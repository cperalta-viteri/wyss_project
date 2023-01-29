#!/bin/bash

# Download Allen data for Mouse whole cortex % hippocampus from GEO
mkdir data/GSE185862
wget -P data/GSE185862 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/*
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.json' > data/GSE185862/dend.json
