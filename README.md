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


## Paper Results

We provide a Snakefile to replicate the results of the paper:  

- *single_geneset*: with this suffix we refer to all the results of for the high-throughput    experiments taken from TCGA biolinks (Fig. 2).   
  **Please note:** the `scripts/tcga_rnaseq.R` script downloads and processes the BLCA RNAseq dataset. This step is time and memory consuming, but, most importantly, we have noticed that it is dvifficult to be able to replicate the exact environment/TCGA version being installed (there are many issues of this kind raised on the TCGAbiolinks github repo ). For reproducibility, inside `data/GDCdata` we provide the BLCA expression data on which our analysis is done and the blca_diffexp.csv file contains the full differential expression results. We have generated this data with TCGAbiolinks v2.15.2. 

- *multi*: refers to the analysis of multiple geneset from the Bailey et al. paper (Fig. 4).



Following the next steps you should be able to run the paper pipelines:

1. First download the data folder from [data repo](https://add_our_data)   

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
        
        snakemake --snakefile Snakefile_paper single_all --configfile config_paper_single.yaml

6. To obtain the results for the multi geneset

        snakemake snakemake --snakefile Snakefile_paper multi_all --configfile config_paper_multi.yaml -t 
        
        snakemake --snakefile Snakefile_paper multi_all --configfile config_paper_multi.yaml

7. To obtain the results for the hdn simulations  

        snakemake --snakefile Snakefile_paper hdn_all --configfile config_paper_hdn.yaml

8. To obtain the results for the multi geneset  

        snakemake --snakefile Snakefile_paper sbm_all --configfile config_paper_sbm.yaml

