'''
    All rules for pygna visualisation
'''



rule plot_topology_module:
    input:
        OUTPATH+"table_topology_module.csv",
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
        OUTPATH+"table_association_rwr.csv",
    output:
        OUTFIGURES+"heatmap_association_rwr.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --rwr --annotate"

rule plot_association_sp:
    input:
        OUTPATH+"table_association_sp.csv",
    output:
        OUTFIGURES+"heatmap_association_sp.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --annotate"