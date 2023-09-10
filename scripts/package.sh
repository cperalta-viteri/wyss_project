

tar -C ../ \
  -cvf 2023-30-06_wyss_project_single_cell_transcriptomics.tar \
  --exclude _doc/ \
  --exclude .Rproj.user/ \
  --exclude apps/geo_download/downloads/ \
  --exclude .DS_Store \
  --exclude .Rhistory \
  --exclude .DS_Store \
  --exclude .git/ \
  --exclude 2023-30-06_wyss_project_single_cell_transcriptomics.tar \
  wyss_project/




scp Rpackages/scml_1.2.0.tar.gz prados@baobab2.hpc.unige.ch:~/scratch/
scp scripts/baobab_train_aba.R prados@baobab2.hpc.unige.ch:~/scratch/

ssh prados@baobab2.hpc.unige.ch
ml GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1 CUDA/11.7.0
R CMD INSTALL ~/scratch/scml_1.2.0.tar.gz


ml GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1 CUDA/11.7.0
srun -p public-longrun-cpu --time=4-00:00:00 --mem=32G -n1 --cpus-per-task=4 Rscript ~/scratch/baobab_train_aba.R &> ~/scratch/baobab_train_aba2.log &

salloc -p shared-cpu --time=3:00:00 --mem=64G -n1 --cpus-per-task=4 
R
Rscript ~/scratch/baobab_train_aba.R &> ~/scratch/baobab_train_aba.log &



