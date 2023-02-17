

# Download Allen 2019 Human Multiple Cortical Areas SMART-see
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv' | pigz > data/human_ctx_sseq/matrix.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/intron.csv' | pigz > data/human_ctx_sseq/intron.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/exon.csv' | pigz > data/human_ctx_sseq/exon.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/trimmed_means.csv' | pigz > data/human_ctx_sseq/trimmed_means.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/medians.csv' | pigz > data/human_ctx_sseq/medians.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv' | pigz > data/human_ctx_sseq/metadata.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/dend.json' > data/human_ctx_sseq/dend.json
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv' | pigz > data/human_ctx_sseq/tsne.csv.gz
curl 'https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip' > data/human_ctx_sseq/taxonomy.zip
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/30/80/30803016-9515-45fd-a6af-8595902c176f/readme_human_ctx.txt' > data/human_ctx_sseq/readme_human_ctx.txt
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip' > data/human_ctx_sseq/sample-exp_component_mapping_human_smart-seq_oct_2019.zip



# Download Allen 2020 10x human dataset
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv' | pigz > data/human_ctx_sseq/matrix.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv' | pigz > data/human_ctx_sseq/metadata.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json' | pigz > data/human_ctx_sseq/dend.json.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv' | pigz > data/human_ctx_sseq/tsne.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv' | pigz > data/human_ctx_sseq/medians.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv' | pigz > data/human_ctx_sseq/trimmed_means.csv.gz
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds' > data/human_ctx_sseq/human_dendrogram.rds
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip' > data/human_ctx_sseq/sample-exp_component_mapping_human_10x_apr2020.zip
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt' > data/human_ctx_sseq/readme-m1-10.txt
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/38/05/3805d61e-5c4d-446d-9331-cb61eb181198/taxonomy.txt' > data/human_ctx_sseq/taxonomy.txt

