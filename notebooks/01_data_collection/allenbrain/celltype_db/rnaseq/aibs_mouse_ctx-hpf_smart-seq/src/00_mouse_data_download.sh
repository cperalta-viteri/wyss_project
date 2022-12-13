#!/bin/zsh

# Download Allen 2019 - Mouse whole cortex % hippocampus SMART-SEQ
# with 10X-SMART-SEQ TAXONOMY (2021)
# https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x

mkdir ../data

curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/expression_matrix.hdf5' | pigz > ../data/matrix.hdf5 &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/expression_matrix_introns.csv' | pigz > ../data/intron.csv.gz &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/expression_matrix_exons.csv' | pigz > ../data/exon.csv.gz &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/trimmed_means.csv' | pigz > ../data/trimmed_means.csv.gz &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/medians.csv' | pigz > ../data/medians.csv.gz &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/metadata.csv' | pigz > ../data/metadata.csv.gz &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/dend.json' > ../data/dend.json &
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/tsne.csv' | pigz > ../data/tsne.csv.gz &
#curl 'https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip' > data/taxonomy.zip
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/5f/8f/5f8fb8ff-1133-4f3b-bdfa-43a0144bd054/readme_mouse_smart-seq.txt' > data/readme_human_ctx.txt
#curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip' > data/sample-exp_component_mapping_human_smart-seq_oct_2019.zip

