'''
    All rules for pygna
'''


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
        OUTPATH+"_topology_module.csv"
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"]
    conda: 
        "../envs/pygna.yaml"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_internal_degree:
    input:
        network=NETWORK,
        geneset=GENESET
    output:
        OUTPATH+"_topology_internal_degree.csv"
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"]
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-internal-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_total_degree:
    input:
        network=NETWORK,
        geneset=GENESET
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"]
    output:
        OUTPATH+"_topology_total_degree.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-total-degree {input.network} {input.geneset} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_sp:
    input:
        network=NETWORK,
        geneset=GENESET,
        matrix=SP_MATRIX
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"]
    output:
        OUTPATH+"_topology_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-sp {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_rwr:
    input:
        network=NETWORK,
        geneset=GENESET,
        matrix=RWR_MATRIX
    params:
        nop=config["topology"]["number_of_permutations"],
        cores=config["topology"]["cores"]
    output:
	    OUTPATH+"_topology_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"

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
	    OUTPATH+"_association_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {output} -B {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"


rule association_SP:
    input:
        network=NETWORK,
        A=GENESET,
        B=config["association"]["geneset_B"],
        matrix=SP_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"]
    output:
	    OUTPATH+"_association_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} -B {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"

# COMPARISON
rule comparison_RW:
    input:
        network=NETWORK,
        A=GENESET,
        B=config["association"]["geneset_B"],
        matrix=RWR_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"]
    output:
	    OUTPATH+"_comparison_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {output} -B {input.B} --number-of-permutations {params.nop} --cores {params.cores}"


rule comparison_SP:
    input:
        network=NETWORK,
        A=GENESET,
        B=config["association"]["geneset_B"],
        matrix=RWR_MATRIX
    params:
        nop=config["association"]["number_of_permutations"],
        cores=config["association"]["cores"]
    output:
	    OUTPATH+"_comparison_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} -B {input.B} --number-of-permutations {params.nop} --cores {params.cores}"


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
	    OUTPATH+"_within_comparison_rwr.csv"
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
	    OUTPATH+"_within_comparison_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {output} --number-of-permutations {params.nop} --cores {params.cores}"