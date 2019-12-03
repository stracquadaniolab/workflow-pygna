'''
    Rules for generating HDN and SBM simulations
    and to run the analysis
'''


# SBM parameters
n_nodes= config["sbm_parameters"]["n_nodes"]
T=config["sbm_parameters"]["theta0"].split(',')
P=config["sbm_parameters"]["percentage"].split(',')
D=config["sbm_parameters"]["density"].split(',')
S=range(config["sbm_parameters"]["number_simulations"])

NOP = config["topology"]["number_of_permutations"]
ANALYSIS = config["topology"]["analyse"]

rule generate_sbm:
   params:
       folder=OUTPATH+"sbm_sim/",
       nodes=n_nodes,
       theta0=config["sbm_parameters"]["theta0"],
       percentage=config["sbm_parameters"]["percentage"],
       density=config["sbm_parameters"]["density"],
       n_s=config["sbm_parameters"]["number_simulations"]
   output:
       network=OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_network.tsv",
       geneset=OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_genes.gmt"
   shell:
       "pygna generate-sbm2-network {params.folder} --n-nodes {params.nodes} -t {params.theta0} --percentage {params.percentage} -d {params.density} --n-simulations {params.n_s}"



'''
    Creates SP and RWR matrices
'''
rule create_SP_matrix_sbm:
    input:
        OUTPATH+"sbm_sim/{suffix}_network.tsv"
    output:
        OUTPATH+"sbm_sim/matrices/{suffix}_SP.hdf5"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule create_RWR_matrix_sbm:
    input:
        OUTPATH+"sbm_sim/{suffix}_network.tsv"
    output:
        OUTPATH+"sbm_sim/matrices/{suffix}_RWR.hdf5"
    shell:
        "pygna build-rwr-diffusion {input} --output-file {output}"


rule analyse_module_sbm:
    input:
        network=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_genes.gmt")
    params:
        cores=config["topology"]["cores"],
        nop=NOP
    output:
        OUTPATH+"results/sbm_t_{t}_p_{p}_d_{d}_s_{s}_table_module.csv"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {output}  --number-of-permutations {params.nop} --size-cut 1 --cores {params.cores}"


rule analyse_sp_sbm:
    input:
        network=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_genes.gmt"),
        matrix=ancient(OUTPATH+"sbm_sim/matrices/sbm_t_{t}_p_{p}_d_{d}_s_{s}_SP.hdf5")
    params:
        nop=NOP,
        cores=config["topology"]["cores"],
    output:
        OUTPATH+"results/sbm_t_{t}_p_{p}_d_{d}_s_{s}_table_sp.csv"
    shell:
        "pygna test-topology-sp {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule analyse_RWR_sbm:
    input:
        network=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"sbm_sim/sbm_t_{t}_p_{p}_d_{d}_s_{s}_genes.gmt"),
        matrix=ancient(OUTPATH+"sbm_sim/matrices/sbm_t_{t}_p_{p}_d_{d}_s_{s}_RWR.hdf5")
    params:
        nop=NOP,
        cores=config["topology"]["cores"],
    output:
        OUTPATH+"results/sbm_t_{t}_p_{p}_d_{d}_s_{s}_table_rwr.csv"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule final_sbm:
    input:
        expand(OUTPATH+"results/sbm_t_{t}_p_{p}_d_{d}_s_{s}_table_{b}.csv", t=T,p=P, d=D, s=S , b=ANALYSIS)
    output:
        OUTPATH+"final_sbm.csv"
    shell:
        """head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done """