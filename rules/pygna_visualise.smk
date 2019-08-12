'''
    All rules for pygna visualisation
'''



rule plot_topology_module:
    input:
        OUTPATH+"_topology_module.csv",
    output:
        OUTFIGURES+"_barplot_module.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_internal_degree:
    input:
        OUTPATH+"_topology_internal_degree.csv",
    output:
        OUTFIGURES+"_barplot_internal_degree.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_total_degree:
    input:
        OUTPATH+"_topology_total_degree.csv",
    output:
        OUTFIGURES+"_barplot_total_degree.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_sp:
    input:
        OUTPATH+"_topology_sp.csv",
    output:
        OUTFIGURES+"_barplot_sp.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_topology_rwr:
    input:
        OUTPATH+"_topology_rwr.csv",
    output:
        OUTFIGURES+"_barplot_rwr.{e}"
    shell:
        "pygna paint-datasets-stats {input} {output}"

rule plot_within_comparison_rwr:
    input:
        OUTPATH+"_within_comparison_rwr.csv",
    output:
        OUTFIGURES+"_heatmap_within_comparison_rwr.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --rwr --single-geneset --annotate"

rule plot_within_comparison_sp:
    input:
        OUTPATH+"_within_comparison_sp.csv",
    output:
        OUTFIGURES+"_heatmap_within_comparison_sp.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --single-geneset --annotate"

rule plot_association_rwr:
    input:
        OUTPATH+"_association_rwr.csv",
    output:
        OUTFIGURES+"_heatmap_association_rwr.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --rwr --annotate"

rule plot_association_sp:
    input:
        OUTPATH+"_association_sp.csv",
    output:
        OUTFIGURES+"_heatmap_association_sp.{e}"
    shell:
        "pygna paint-comparison-matrix {input} {output} --annotate"