# PyGNA Snakemake workflow:

Current version:  1.0.1-dev

A simple Snakemake workflow to perform network analyses using the Python Geneset Network Analysis ([PyGNA](https://github.com/stracquadaniolab/pygna)) package.

## Authors

* Viola Fanfani, v.fanfani@sms.ed.ac.uk (lead developer)
* Giovanni Stracquadanio, giovanni.stracquadanio@ed.ac.uk

## Overview

![dag.png](dag.png)

## Usage

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/stracquadaniolab/workflow-pygna/releases).

In any case, if you use this workflow in a paper, please cite our PyGNA as follows:


#### Step 2: Configure workflow

Configure the workflow according to your needs, via editing the file `config.yaml`.
Now a template config file is `config_temp.yaml`.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stablve/executable.html) for further details.

#### Step 4: Check results

Results are stored in the `results` folder. 

## Benchmarking pipeline

### Stochastic Block Model 

We provide a full pipeline for the benchmark of GNT and GNA methods through the SBM generated data.

First check the `config_sbm.yaml` configfile and use the desired parameters, then the whole pipeline can be run with:

    snakemake --snakefile Snakefile_sbm --configfile config_sbm.yaml --cores 1

### High Degree Nodes

For the HDN simulations please refere to 
the previous paper snakefile `Snakefile_paper_old`

## Paper Results, TCGA

The entire pipeline for the processing of cancer data can be run as follows 

**Please note:** the `scripts/tcga-download.R` script downloads and processes the RNAseq dataset. This step is time and memory consuming, but, most importantly, we have noticed that it is difficult to be able to replicate the exact environment/TCGA version being installed (there are many issues of this kind raised on the TCGAbiolinks github repo ). For reproducibility, we provide the differential expression tables already preprocessed. We have generated this data with R=4 and TCGAbiolinks v2.14.0. 

  use the `--use-conda` flag to install the environment with the same parameters


Following the next steps you should be able to run the paper pipelines:

1. First download the data folder from [Zenodo](https://zenodo.org/record/3574027#.XfObSZP7RTY)   

2. The pipeline uses the data folder path, you can:  
    A. add the data folder inside the workflow-pygna folder 
    B. change the relative path in the config file  

3. Check the `config_paper.yaml` configuration files. They include number of permutations and cores parameters, tweak them as needed (for the moment we have set 3 and 1
so that you can quickly check if all results are generated. ).

4. As we provide intermediate files whose generation can be very time consuming      (differential expression results and SP/RWR matrices), run 

        snakemake --snakefile Snakefile_paper  --configfile config_paper.yaml -t 

    to update all files time tag and avoid recreating them (snakemake would try to run the pipeline rule again if the files have been recently modified).


## Paper Results, previous version

We provide a Snakefile to replicate the results of the paper:  

    snakemake --snakefile Snakefile_paper --configfile config_paper --use-conda --cores $N
- *single_geneset*: with this suffix we refer to all the results of for the high-throughput    experiments taken from TCGA biolinks (Fig. 2).   
  **Please note:** the `scripts/tcga_rnaseq.R` script downloads and processes the BLCA RNAseq dataset. This step is time and memory consuming, but, most importantly, we have noticed that it is difficult to be able to replicate the exact environment/TCGA version being installed (there are many issues of this kind raised on the TCGAbiolinks github repo ). For reproducibility, we provide the blca_diffexp.csv file that contains the full differential expression results. We have generated this data with R=3.6 and TCGAbiolinks v2.14.0. 
  use the `--use-conda` flag to install the environment with the same parameters
- *multi*: refers to the analysis of multiple geneset from the Bailey et al. paper (Fig. 4).

Following the next steps you should be able to run the paper pipelines:

1. First download the data folder from [Zenodo](https://zenodo.org/record/3574027#.XfObSZP7RTY)   

2. The pipeline uses the data folder path, you can:  
    A. add the data folder inside the workflow-pygna folder 
    B. change the relative path in the config file  

3. Check the `config_paper_single.yaml` and `config_paper_multi.yaml` configuration files. They include number of permutations and cores parameters, tweak them as needed (for the moment we have set 3 and 1
so that you can quickly check if all results are generated. ).

4. As we provide intermediate files whose generation can be very time consuming      (differential expression results and SP/RWR matrices), run 

        snakemake --snakefile Snakefile_paper <analysis>_all --configfile config_paper_<analysis>.yaml -t 

   with `<analisys>` being one of single or multi, to update all files time tag and avoid recreating them (snakemake would try to run the pipeline rule again if the files have been recently modified).

5. To obtain all the results for the single geneset (avoid the first step to have the full regeneration of all files):

        snakemake snakemake --snakefile Snakefile_paper single_all --configfile config_paper_single.yaml -t 
        
        snakemake --snakefile Snakefile_paper single_all --configfile config_paper_single.yaml --use-conda

6. To obtain the results for the multi geneset

        snakemake snakemake --snakefile Snakefile_paper multi_all --configfile config_paper_multi.yaml -t 
        
        snakemake --snakefile Snakefile_paper multi_all --configfile config_paper_multi.yaml

7. To obtain the results for the hdn simulations  

        snakemake --snakefile Snakefile_paper hdn_all --configfile config_paper_hdn.yaml

8. To obtain the results for the multi geneset  

        snakemake --snakefile Snakefile_paper sbm_all --configfile config_paper_sbm.yaml

