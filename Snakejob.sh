#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N workflow-pygna
#$ -cwd
#$ -l h_rt=200:00:00
#$ -M fabio.cassano@ed.ac.uk
#$ -pe smp 8

source ~fcassano/.bashrc
conda activate pygna38

snakemake --snakefile Snakefile --cores 16

source deactivate
module unload anaconda
