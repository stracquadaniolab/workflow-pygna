
rule create_SP_mix_matrix:
    input:
        sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"
    output:
        sim_path+"matrices/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_sp.hdf5"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule create_rwr_mix_matrix:
    input:
        sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"
    output:
        sim_path+"matrices/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_rwr.hdf5"
    shell:
        "pygna build-rwr-diffusion {input} --output-file {output}"




rule generate_gnt_sbm:
   params:
       nodes=n_nodes,
       np="{np}",
       nc="{nc}",
       nn="{nn}",
       s="{s}",
   output:
       network=sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv",
       geneset=sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt"
   shell:
       """
       pygna generate-gnt-sbm {output.network} {output.geneset}\
       --N {params.nodes} -f {params.nc} --d {params.np} --block-size {params.nn}"""

rule analyse_mix_module:
    input:
        network=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt")
    params:
        s="gnt_sbm_np_{np}_vp_{nc}_p_{nn}_s_{s}",
        out=OUTPATH,
        nop=NOP
    output:
        OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_module.csv"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --size-cut 1"

rule analyse_mix_total_degree:
    input:
        network=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_total_degree.csv"
    shell:
        "pygna test-topology-total-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --size-cut 1"

rule analyse_mix_internal_degree:
    input:
        network=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_internal_degree.csv"
    shell:
        "pygna test-topology-internal-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --size-cut 1"

rule analyse_mix_location:
    input:
        network=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt"),
        matrix=ancient(sim_path+"matrices/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_sp.hdf5")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_sp.csv"
    shell:
        "pygna test-topology-sp {input.network}  {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores 2 --size-cut 1"

rule analyse_mix_rwr_vip:
    input:
        network=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_genes.gmt"),
        matrix=ancient(sim_path+"matrices/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_rwr.hdf5")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_rwr.csv"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores 2 --size-cut 1"


rule final_mix_single_analysis:
    input:
        expand(OUTPATH+"pygna_results/gnt_sbm_np_{np}_fc_{nc}_nn_{nn}_s_{s}_table_{{b}}.csv", np=NPm,nc=NCm, nn=NNm, s=Sm)
    output:
        expand("{sim_path}final_gnt_{{b}}.csv", sim_path=OUTPATH)
    shell:
        "head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done "



rule mix_final:
    input:
        expand(OUTPATH+"final_gnt_{b}.csv",b=gnt_analysis)
    output:
        OUTPATH + "final_gnt.csv"
    shell:
        "head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done "
