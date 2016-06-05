#!/bin/bash
#$ -S /bin/bash
#cd /home/rprunier/structure
#$ -cwd
#$ -M rachel.prunier@gmail.com
#$ -m besa
#$ -N f_stats_shuffled


module load JAGS/4.0.0

Rscript f-statistics-cluster.R
Rscript extract-thetas.R
Rscript detect-outlier.R