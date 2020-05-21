configfile: "config_temp.yaml"
import os

# Network Parameters
NETWORK=config["parameters"]["network_file"]
OUTPATH=config["parameters"]["outpath"]
GENESET=config["parameters"]["tcga_dataset"]
GENESETB = config["association"]["geneset_B"]

TOPOLOGY = config["topology"]["analyse"]
ASSOCIATION = config["association"]["analyse"]
COMPARISON = config["within_comparison"]["analyse"]

EXTENSIONS=config['figures']['extension']
OUTFIGURES=OUTPATH + config['figures']['outpath']

os.mkdir(config["parameters"]["diagnostic_folder"])
# If the matrices for the full analysis are not specified, they are generated by the rule
if len(config["parameters"]["sp_matrix"])>0:
    SP_MATRIX=config["parameters"]["sp_matrix"]
else:
    SP_MATRIX=OUTPATH+"_sp_matrix.hdf5"

if len(config["parameters"]["rwr_matrix"])>0:
    RWR_MATRIX=config["parameters"]["rwr_matrix"]
else:
    RWR_MATRIX=OUTPATH+"_rwr_matrix.hdf5"

#include: "rules/pygna.smk"
include: "rules/pygna_visualise.smk"
include: "rules/paper_analysis.smk"

rule all:
    input:
        expand(OUTPATH+"{n}/{n}.csv", n=GENESET),
        expand(OUTPATH+"{n}/table_topology_{t}.csv", t=TOPOLOGY, n=GENESET),
        expand(OUTFIGURES+"{n}/barplot_{t}.{e}", t=TOPOLOGY, e=EXTENSIONS,n=GENESET),
        expand(OUTPATH+"{n}/{m}/table_association_{t}.csv", t=ASSOCIATION, n=GENESET, m=GENESETB),
        expand(OUTFIGURES+"{n}/{m}/heatmap_association_{t}.{e}", t=ASSOCIATION, e=EXTENSIONS, n=GENESET, m=GENESETB),
        expand(OUTPATH+"{n}/table_within_comparison_{t}.csv", t=COMPARISON, n=GENESET),
        expand(OUTFIGURES+"{n}/heatmap_within_comparison_{t}.{e}", t=COMPARISON, e=EXTENSIONS,n=GENESET),

rule GNT_all:
    input:
        expand(OUTPATH+"{n}/table_topology_{t}.csv", t=TOPOLOGY,n=GENESET),
        expand(OUTFIGURES+"{n}/barplot_{t}.{e}", t=TOPOLOGY, e=EXTENSIONS,n=GENESET),

rule GNA_association_all:
    input:
        expand(OUTPATH+"{n}/table_association_{t}.csv", t=ASSOCIATION,n=GENESET),
        expand(OUTFIGURES+"{n}/heatmap_association_{t}.{e}", t=ASSOCIATION, e=EXTENSIONS,n=GENESET),

rule single_geneset:
    input:
        expand(OUTPATH+"{n}/table_topology_{t}.csv", t=TOPOLOGY,n=GENESET),
        expand(OUTFIGURES+"{n}/barplot_{t}.{e}", t=TOPOLOGY, e=EXTENSIONS,n=GENESET),
        expand(OUTPATH+"{n}/table_association_{t}.csv", t=ASSOCIATION,n=GENESET),
        expand(OUTFIGURES+"{n}/heatmap_association_{t}.{e}", t=ASSOCIATION, e=EXTENSIONS,n=GENESET),
