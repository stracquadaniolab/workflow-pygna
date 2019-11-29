'''
    Rules for generating HDN and SBM simulations
    and to run the analysis
'''


# VIP parameters
n_nodes= config["vip_parameters"]["n_nodes"]
NP=[str(i) for i in network_probability]
VP=[str(i) for i in vip_probability]
P=[str(i) for i in vip_percentage]
S=range(number_of_simulations)


rule generate_vip:
   params:
       folder=hdn_data_folder+"",
       prefix="sim_vip_np_{np}_vp_{vp}_p_{p}",
       nodes=n_nodes,
       np="{np}",
       vp="{vp}",
       p="{p}",
       s=number_of_simulations
   output:
       network=hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv",
       geneset=hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"
   shell:
       "pygna generate-vip-network {params.folder} {params.prefix} --n-nodes {params.nodes} --network-prob {params.np} --vip-prob {params.vp} --vip-percentage {params.p} --number-of-simulations {params.s}"




rule analyse_module:
    input:
        network=ancient(hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(hdn_data_folder+"new_genesets/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt")
    params:
        nop=NOP,
        cores=config["topology"]["cores"]
    output:
        sim_path+"results/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_table_module.csv"
    shell:
        "pygna analyse-module {input.network} {input.geneset} {params.out} {params.s} --number-of-permutations {params.nop} --size-cut 1"

rule analyse_internal_degree:
    input:
        network=ancient(hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(hdn_data_folder+"new_genesets/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt")
    params:
        s="sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}",
        out=OUTPATH_HDN
        nop=NOP
    output:
        sim_path+"results/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_table_internal_degree.csv"
    shell:
        "pygna analyse-internal-degree {input.network} {input.geneset} {params.out} {params.s} --number-of-permutations {params.nop} --size-cut 1"

rule analyse_location:
    input:
        network=ancient(hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(hdn_data_folder+"new_genesets/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"),
        matrix=ancient(sim_path+"matrices/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_SP.hdf5")
    params:
        s="sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}",
        out=OUTPATH_HDN
        nop=NOP
    output:
        sim_path+"results/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_table_SP.csv"
    shell:
        "pygna analyse-location {input.network} {input.matrix} {input.geneset} {params.out} {params.s} --number-of-permutations {params.nop} --cores 2"

rule analyse_RWR_vip:
    input:
        network=ancient(hdn_data_folder+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_network.tsv"),
        geneset=ancient(hdn_data_folder+"new_genesets/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_genes.gmt"),
        matrix=ancient(sim_path+"matrices/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_RWR.hdf5")
    params:
        s="sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}",
        out=OUTPATH_HDN,
        nop=NOP
    output:
        sim_path+"results/sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_table_RW.csv"
    shell:
        "pygna analyse-RW {input.network} {input.geneset} {input.matrix} {params.out} {params.s} --number-of-permutations {params.nop} --cores 2"

rule final_vip:
    input:
        expand(OUTPATH_HDN+"sim_vip_np_{np}_vp_{vp}_p_{p}_s_{s}_table_{b}.csv", np=NP,vp=VP, p=P, s=S , b=ANALYSIS)
    output:
        OUTPATH_HDN+"final_vip_module.csv"
    shell:
        "head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done "
