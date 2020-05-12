#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N project-network-gan
#$ -cwd
#$ -l h_rt=200:00:00
#$ -M fabio.cassano@ed.ac.uk
#$ -pe smp 8

source ~fcassano/.bashrc
conda activate pygna

snakemake --snakefile Snakefile_paper --configfile config_paper_multi.yaml --cores 1

source deactivate
module unload anaconda
