

# Download Allen 2020 10x human dataset
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv' | pigz > data/matrix.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv' | pigz > data/metadata.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json' | pigz > data/dend.json.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv' | pigz > data/tsne.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv' | pigz > data/medians.csv.gz
curl 'https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv' | pigz > data/trimmed_means.csv.gz
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds' > data/human_dendrogram.rds
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip' > data/sample-exp_component_mapping_human_10x_apr2020.zip
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt' > data/readme-m1-10.txt
curl 'https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/38/05/3805d61e-5c4d-446d-9331-cb61eb181198/taxonomy.txt' > data/taxonomy.txt


