# PyGNA Snakemake workflow:

Current version:  0.1.1-dev

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

- *single_geneset*: with this suffix we refer to all the results of for the high-throughput experiments taken from TCGA biolinks (Fig. 2).    
- *multi*: refers to the analysis of multiple geneset from the Bailey et al. paper (Fig. 4).

Following the next steps you should be able to run the paper pipelines:

1. First download the data folder from [data repo](https://add_our_data)   

2. The pipeline uses a relative path for the data repo, you can:  
    A. add the data folder in the same location you have the workflow-pygna repo  
    B. change the relative path in the config file  

3. Check the `config_paper_single.yaml` and `config_paper_multi.yaml` configuration files. They include number of permutations and cores parameters, tweak them as needed (for the moment we have set 3 and 1
so that you can quickly check if all results are generated. ).

4. To obtain all the results for the single geneset

    snakemake --snakefile Snakefile_paper single_all --configfile config_paper_single.yaml

5. To obtain the results for the multi geneset  

    snakemake --snakefile Snakefile_paper multi_all --configfile config_paper_multi.yaml

6. To obtain the results for the hdn simulations  

    snakemake --snakefile Snakefile_paper hdn_all --configfile config_paper_hdn.yaml

7. To obtain the results for the multi geneset  

    snakemake --snakefile Snakefile_paper sbm_all --configfile config_paper_sbm.yaml