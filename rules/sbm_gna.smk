
rule create_SP_ct_matrix:
    input:
        sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_network.tsv"
    output:
        sim_path+"matrices/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_sp.hdf5"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule create_rwr_ct_matrix:
    input:
        sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_network.tsv"
    output:
        sim_path+"matrices/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_rwr.hdf5"
    shell:
        "pygna build-rwr-diffusion {input} --output-file {output}"


rule generate_ct_sbm:
   params:
       nodes = n_nodes,
       np    = "{np}",
       cis   = "{cis}",
       trans = "{trans}",
       nn    = "{nn}",
       s     = "{s}",
   output:
       network=sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_network.tsv",
       geneset=sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes.gmt",
       geneset2 = sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes_B.gmt"
   shell:
       """
       pygna generate-gna-sbm {output.network} {output.geneset}\
        --output-gmt2 {output.geneset2}\
       --N {params.nodes} --d {params.np} --fc-trans {params.trans} --fc-cis {params.cis} --block-size {params.nn}\
       """

rule analyse_gna_location:
    input:
        network=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes.gmt"),
        genesetB=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes_B.gmt"),
        matrix=ancient(sim_path+"matrices/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_sp.hdf5")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_table_association_sp.csv"
    shell:
        """
        pygna test-association-sp {input.network}  {input.geneset} {input.matrix} {output}\
        --file-geneset-b {input.genesetB} --number-of-permutations {params.nop} --cores 1 --size-cut 1
        """

rule analyse_gna_rwr:
    input:
        network=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_network.tsv"),
        geneset=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes.gmt"),
        genesetB=ancient(sim_path+"gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_genes_B.gmt"),
        matrix=ancient(sim_path+"matrices/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_rwr.hdf5")
    params:
        nop=NOP
    output:
        OUTPATH+"pygna_results/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_table_association_rwr.csv"
    shell:
        """
        pygna test-association-rwr {input.network} {input.geneset} {input.matrix} {output}\
        --file-geneset-b {input.genesetB} --number-of-permutations {params.nop} --cores 1 --size-cut 1
        """


rule final_ct_single_analysis:
    input:
        expand(OUTPATH+"pygna_results/gna_sbm_np_{np}_cis_{cis}_trans_{trans}_nn_{nn}_s_{s}_table_association_{{b}}.csv", np=NPct,cis=cis, trans=trans, nn=NNct , s=Sct)
    output:
        expand("{sim_path}final_gna_{{b}}.csv", sim_path=OUTPATH)
    shell:
        "head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done "



rule ct_final:
    input:
        expand(OUTPATH+"final_gna_{b}.csv",b=gna_analysis)
    output:
        OUTPATH + "final_gna.csv"
    shell:
        "head -1 {input[0]} > {output}; for filename in {input}; do sed 1d $filename >> {output}; done "
