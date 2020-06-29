'''
    All rules for pygna
'''

rule generate_data:
    output:
        OUTPATH+"{n}/{n}.csv",
        OUTPATH+"{n}/log.csv",
        OUTPATH+"{n}/{n}_datafilt.csv"
    conda: "../envs/tcgabiolinks.yaml"
    params:
        name= "{n}",
        folder= OUTPATH+"{n}",
        raw_data = OUTPATH+"raw_data/"
    script:
        "../scripts/tcga-download.R"

rule generate_gmt:
    input:
         OUTPATH+"{n}/{n}.csv"
    output:
        OUTPATH+"{n}/{n}.gmt"
    conda: "../envs/pygna.yaml"
    params:
        dataset="{n}"
    shell:
        "pygna geneset-from-table {input} {params.dataset} --output-gmt {output} -f significant -d significant -n genes.Entrezid -t 0.5 -a greater"

rule merge_gmt:
    input:
        expand(OUTPATH+"{n}/{n}.gmt", n=GENESET),
    output:
        OUTPATH+"merged.gmt"
    shell:
        "cat {input} >> {output}"

rule generate_matrix_sp:
    input:
        NETWORK
    output:
        protected(SP_MATRIX)
    conda: "../envs/pygna.yaml"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule generate_matrix_rwr:
    input:
        NETWORK
    output:
        protected(RWR_MATRIX)
    conda: "../envs/pygna.yaml"
    shell:
        "pygna build-rwr-diffusion {input} --output-file {output}"

rule topology_module:
    input:
        network=NETWORK,
        geneset=OUTPATH+"merged.gmt"
    output:
        OUTPATH+"table_topology_module.csv"
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"],
        diagnostic_folder = config["parameters"]["diagnostic_folder"]
    conda:
        "../envs/pygna.yaml"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores}    --diagnostic-null-folder {params.diagnostic_folder}"

rule topology_internal_degree:
    input:
        network=NETWORK,
        geneset=OUTPATH+"merged.gmt"
    output:
        OUTPATH+"table_topology_internal_degree.csv"
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"],
        diagnostic_folder = config["parameters"]["diagnostic_folder"]
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-internal-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores}  --diagnostic-null-folder {params.diagnostic_folder}"

rule topology_total_degree:
    input:
        network=NETWORK,
        geneset=OUTPATH+"merged.gmt"
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"],
        diagnostic_folder = config["parameters"]["diagnostic_folder"]
    output:
        OUTPATH+"table_topology_total_degree.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-total-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores} --diagnostic-null-folder {params.diagnostic_folder}"

rule topology_sp:
    input:
        network=NETWORK,
        geneset=OUTPATH+"merged.gmt",
        matrix=SP_MATRIX
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"],
        diagnostic_folder = config["parameters"]["diagnostic_folder"]
    output:
        OUTPATH+"table_topology_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-sp {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores} --diagnostic-null-folder {params.diagnostic_folder}"

rule topology_rwr:
    input:
        network=NETWORK,
        geneset=OUTPATH+"merged.gmt",
        matrix=RWR_MATRIX
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"],
        diagnostic_folder = config["parameters"]["diagnostic_folder"]
    output:
	    OUTPATH+"table_topology_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores} --diagnostic-null-folder {params.diagnostic_folder}"

rule association_RW:
    input:
        network=NETWORK,
        A=OUTPATH+"{n}/{n}.gmt",
        B=GENESETB,
        matrix=RWR_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"]
    output:
	    OUTPATH+"{n}/table_association_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {output} --file-geneset-b {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"



rule association_SP:
    input:
        network=NETWORK,
        A=OUTPATH+"{n}/{n}.gmt",
        B=GENESETB,
        matrix=SP_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"],
        diagnostic_folder = config["parameters"]["outpath"]
    output:
	    OUTPATH+"{n}/table_association_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} --file-geneset-b {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"



# SINGLE GENESET COMPARISON
rule within_comparison_RW:
    input:
        network=NETWORK,
        A=OUTPATH+"merged.gmt",
        matrix=RWR_MATRIX
    params:
        nop=config["within_comparison"]["number_of_permutations"],
        cores=config["within_comparison"]["cores"]
    output:
	    OUTPATH+"table_within_comparison_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"


rule within_comparison_SP:
    input:
        network=NETWORK,
        A=OUTPATH+"merged.gmt",
        matrix=SP_MATRIX
    params:
        nop=config["within_comparison"]["number_of_permutations"],
        cores=config["within_comparison"]["cores"]
    output:
	    OUTPATH+"table_within_comparison_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"



rule plot_topology_module:
    input:
        OUTPATH+"table_topology_module.csv"
    output:
        OUTFIGURES+"barplot_module.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_internal_degree:
    input:
        OUTPATH+"table_topology_internal_degree.csv",
    output:
        OUTFIGURES+"barplot_internal_degree.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_total_degree:
    input:
        OUTPATH+"table_topology_total_degree.csv",
    output:
        OUTFIGURES+"barplot_total_degree.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_sp:
    input:
        OUTPATH+"table_topology_sp.csv",
    output:
        OUTFIGURES+"barplot_sp.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_rwr:
    input:
        OUTPATH+"table_topology_rwr.csv",
    output:
        OUTFIGURES+"barplot_rwr.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_within_comparison_rwr:
    input:
        OUTPATH+"table_within_comparison_rwr.csv",
    output:
        OUTFIGURES+"heatmap_within_comparison_rwr.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --rwr --single-geneset --annotate"

rule plot_within_comparison_sp:
    input:
        OUTPATH+"table_within_comparison_sp.csv",
    output:
        OUTFIGURES+"heatmap_within_comparison_sp.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --single-geneset --annotate"

rule plot_association_rwr:
    input:
        OUTPATH+"{n}/table_association_rwr.csv",
    output:
        OUTFIGURES+"{n}_volcano_association_rwr.{e}"
    shell:
        "pygna paint-volcano-plot {input} {output} --rwr --annotate --threshold-y 1"

rule plot_association_sp:
    input:
        OUTPATH+"{n}/table_association_sp.csv",
    output:
        OUTFIGURES+"{n}_volcano_association_sp.{e}"
    shell:
        "pygna paint-volcano-plot {input} {output} --annotate --threshold-y 1"
