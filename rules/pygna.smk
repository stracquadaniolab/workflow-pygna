'''
    Analysis of the genesesets with the merged file
'''

rule generate_matrix_sp:
    input:
        NETWORK
    output:
        SP_MATRIX
    conda: "../envs/pygna.yaml"
    shell:
        "pygna build-distance-matrix {input} {output}"

rule generate_matrix_rwr:
    input:
        NETWORK
    output:
        RWR_MATRIX
    conda: "../envs/pygna.yaml"
    shell:
        "pygna build-RWR-diffusion {input} --output-file {output}"

rule topology_module:
    input:
        network=NETWORK,
        geneset=GENESET
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
        OUTPUT_FOLDER+PREFIX+"_table_module.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-module {input.network} {input.geneset} {params.out} {params.prefix} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_internal_degree:
    input:
        network=NETWORK,
        geneset=GENESET
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
        OUTPUT_FOLDER+PREFIX+"_table_internal_degree.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-internal-degree {input.network} {input.geneset} {params.out} {params.prefix} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_total_degree:
    input:
        network=NETWORK,
        geneset=GENESET
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
        OUTPUT_FOLDER+PREFIX+"_table_total_degree.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-total-degree {input.network} {input.geneset} {params.out} {params.prefix} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_sp:
    input:
        network=NETWORK,
        geneset=GENESET,
        matrix=SP_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
        OUTPUT_FOLDER+PREFIX+"_table_SP.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-sp {input.network} {input.geneset}  {input.matrix} {params.out} {params.prefix} --number-of-permutations {params.nop} --cores {params.cores}"

rule topology_rwr:
    input:
        network=NETWORK,
        geneset=GENESET,
        matrix=RWR_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
	    OUTPUT_FOLDER+PREFIX+"_table_RW.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-topology-rwr {input.network} {input.geneset} {input.matrix} {params.out} {params.prefix} --number-of-permutations {params.nop} --cores {params.cores}"

rule association_RW:
    input:
        network=NETWORK,
        A=GENESET,
        B=GENESET_B,
        matrix=RWR_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
	    OUTPUT_FOLDER+PREFIX+"_table_association_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {params.out} {params.prefix} -B {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"


rule association_SP:
    input:
        network=NETWORK,
        A=GENESET,
        B=GENESET_B,
        matrix=SP_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
	    OUTPUT_FOLDER+PREFIX+"_table_association_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {params.out} {params.prefix} -B {input.B} --keep --number-of-permutations {params.nop} --cores {params.cores}"


rule comparison_RW:
    input:
        network=NETWORK,
        A=GENESET,
        B=GENESET_B,
        matrix=RWR_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
	    OUTPUT_FOLDER+PREFIX+"_table_comparison_rwr.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-rwr {input.network} {input.A} {input.matrix} {params.out} {params.prefix} -B {input.B} --number-of-permutations {params.nop} --cores {params.cores}"


rule comparison_SP:
    input:
        network=NETWORK,
        A=GENESET,
        B=GENESET_B,
        matrix=RWR_MATRIX
    params:
        prefix=PREFIX,
        out=OUTPUT_FOLDER,
	    nop=NOP,
        cores=CORES
    output:
	    OUTPUT_FOLDER+PREFIX+"_table_comparison_sp.csv"
    conda: "../envs/pygna.yaml"
    shell:
        "pygna test-association-sp {input.network} {input.A} {input.matrix} {params.out} {params.prefix} -B {input.B} --number-of-permutations {params.nop} --cores {params.cores}"
