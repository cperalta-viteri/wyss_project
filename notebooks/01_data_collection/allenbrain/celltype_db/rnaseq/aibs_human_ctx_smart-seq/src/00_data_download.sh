

# Download Allen 2019 Human Multiple Cortical Areas SMART-see
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv' | pigz > data/matrix.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/intron.csv' | pigz > data/intron.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/exon.csv' | pigz > data/exon.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/trimmed_means.csv' | pigz > data/trimmed_means.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/medians.csv' | pigz > data/medians.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv' | pigz > data/metadata.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/dend.json' > data/dend.json
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv' | pigz > data/tsne.csv.gz
curl 'https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip' > data/taxonomy.zip
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/30/80/30803016-9515-45fd-a6af-8595902c176f/readme_human_ctx.txt' > data/readme_human_ctx.txt
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip' > data/sample-exp_component_mapping_human_smart-seq_oct_2019.zip

