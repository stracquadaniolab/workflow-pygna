'''
    Rules for generating HDN and SBM simulations
    and to run the analysis
'''


# HDN parameters
n_nodes= config["hdn_parameters"]["n_nodes"]
NP=[str(i) for i in config["hdn_parameters"]["network_probability"]]
VP=[str(i) for i in config["hdn_parameters"]["hdn_probability"]]
P=[str(i) for i in config["hdn_parameters"]["hdn_percentage"]]
S=range(config["hdn_parameters"]["number_simulations"])
NOP = config["topology"]["number_of_permutations"]
ANALYSIS = config["topology"]["analyse"]


'''
    Creates SP and RWR matrices
'''
rule create_SP_matrix_sim:
    input:
        OUTPATH+"hdn_sim/{suffix}_network.tsv"
    output:
        OUTPATH+"hdn_sim/matrices/{suffix}_SP.hdf5"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule create_RWR_matrix_sim:
    input:
        OUTPATH+"hdn_sim/{suffix}_network.tsv"
    output:
        OUTPATH+"hdn_sim/matrices/{suffix}_RWR.hdf5"
    shell:
        "pygna build-rwr-diffusion {input} --output-file {output}"


rule generate_hdn:
   params:
       folder=OUTPATH+"hdn_sim/",
       prefix="sim_hdn_np_{np}_vp_{vp}_p_{p}",
       nodes=n_nodes,
       np="{np}",
       vp="{vp}",
       p="{p}",
       s=config["hdn_parameters"]["number_simulations"]
   output:
       network=OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv",
       geneset=OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"
   shell:
       "pygna generate-hdn-network {params.folder} {params.prefix} --n-nodes {params.nodes} --network-prob {params.np} --vip-prob {params.vp} --vip-percentage {params.p} --number-of-simulations {params.s}"



rule analyse_module_hdn:
    input:
        network=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt")
    params:
        cores=config["topology"]["cores"],
        nop=NOP
    output:
        OUTPATH+"results/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_table_module.csv"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {output}  --number-of-permutations {params.nop} --size-cut 1 --cores {params.cores}"


rule analyse_sp_hdn:
    input:
        network=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"),
        matrix=ancient(OUTPATH+"hdn_sim/matrices/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_SP.hdf5")
    params:
        nop=NOP,
        cores=config["topology"]["cores"],
    output:
        OUTPATH+"results/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_table_sp.csv"
    shell:
        "pygna test-topology-sp {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule analyse_RWR_hdn:
    input:
        network=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(OUTPATH+"hdn_sim/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"),
        matrix=ancient(OUTPATH+"hdn_sim/matrices/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_RWR.hdf5")
    params:
        nop=NOP,
        cores=config["topology"]["cores"],
    output:
        OUTPATH+"results/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_table_rwr.csv"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule final_hdn:
    input:
        expand(OUTPATH+"results/sim_hdn_np_{np}_vp_{vp}_p_{p}_s_{s}_table_{b}.csv", np=NP,vp=VP, p=P, s=S , b=ANALYSIS)
    output:
        OUTPATH+"final_hdn.csv"
    shell:
        """head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done """