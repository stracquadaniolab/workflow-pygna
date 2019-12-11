'''
    All rules for pygna
'''

rule generate_data:
    input:
        OUTPATH
    output:
        GENESET_CSV
    script:
        "../scripts/TCGA_vignette.R"

rule generate_gmt:
    input:
        GENESET_CSV
    output:
        GENESET
    shell:
        "pygna geneset-from-table {input} tcga_biolink_brca --output-gmt {output} -f FDR -d significant -n genes.Entrezid"


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
        geneset=GENESET
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
        geneset=GENESET
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
        geneset=GENESET
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
        geneset=GENESET,
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
        geneset=GENESET,
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
        A=GENESET,
        B=config["association"]["geneset_B"],
        matrix=RWR_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"]
    output:
	    OUTPATH+"table_association_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {output} --file-geneset-b {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"



rule association_SP:
    input:
        network=NETWORK,
        A=GENESET,
        B=config["association"]["geneset_B"],
        matrix=SP_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"],
        diagnostic_folder = config["parameters"]["outpath"]
    output:
	    OUTPATH+"table_association_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} --file-geneset-b {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"


# DIFFUSION HOTNET
rule test_diffusion_hotnet:
    input:
        network=NETWORK,
        A=GENESET_CSV,
        matrix=RWR_MATRIX,
    params:
        nop=config["within_comparison"]["number_of_permutations"],
        cores=config["within_comparison"]["cores"]
    output:
	    OUTPATH+"table_diffusion.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-diffusion-hotnet {input.network} {input.A} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores} --name-column genes.Entrezid --weight-column logFC --filter-column PValue --normalise"


# SINGLE GENESET COMPARISON
rule within_comparison_RW:
    input:
        network=NETWORK,
        A=GENESET,
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
        A=GENESET,
        matrix=SP_MATRIX
    params:
        nop=config["within_comparison"]["number_of_permutations"],
        cores=config["within_comparison"]["cores"]
    output:
	    OUTPATH+"table_within_comparison_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"
