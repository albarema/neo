#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""


import pandas as pd

configfile: "config.yaml"

## --------------------------------------------------------------------------------
##### Modules #####
include: "rules/pcascan.smk"
include: "rules/polyadapt.smk"

## --------------------------------------------------------------------------------
##### Wildcards #####
wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+"

PANNAME=["euras10","euras50"]
## --------------------------------------------------------------------------------
##### Target rules #####
def run_all_rules(_):
    inputs = []

    # this will trigger rules get_plink_euras, run_pca and plot_pca
    #inputs.append("pca/plots/neo.impute-euras.pcadapt-k3.pca.plot.pdf")

    # annotate top genes
    inputs.extend(expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-ALLk.manhattan.plot.pdf", panel=PANNAME))

    # this will run one of the phenotypes
    #"UKBiobank/data/gwasfreqs-pops-102_irnt.tsv.gz",

    # this will run all the phenotypes
    for pheno in pd.read_table('phenoname.txt')['phenoname'].tolist():
        tsv = checkpoints.polyAdapt_freqs.get(pheno=pheno, level='pops').output.candi

        with open(tsv) as fin:
            if len(fin.readlines()) > 3:
                inputs.append("UKBiobank/selection_UKBV2/Genscores-pops-{pheno}.txt".format(pheno=pheno))
    #qx
    #expand("UKBiobank/selection_UKBV2/Genscores-pops-{pheno}.txt",pheno=pd.read_table('phenoname.txt')['phenoname'].tolist())

    return inputs


## targets
rule all:
    input:
        run_all_rules

        # expand("annotate_genes/neo.impute-{panel}-pcadapt-ALLk.manhattan.plot.pdf", panel=PANNAMES)
        # this will trigger rules get_plink_euras, run_pca and plot_pca
        #"annotate_genes/all-k.manhattan.plot.pdf"

        # annotate top genes
        #"annotate_genes/neo.impute-euras.pcadapt.manhattan_annotated.png",
        # this will run one of the phenotypes
        #"UKBiobank/data/gwasfreqs-pops-102_irnt.tsv.gz",

        # this will run all the phenotypes
        #expand("UKBiobank/data/gwasfreqs-pops-{pheno}.tsv.gz",pheno=pd.read_table('phenoname.txt')['phenoname'].tolist())
        #qx
        #expand("UKBiobank/selection_UKBV2/Genscores-pops-{pheno}.txt",pheno=pd.read_table('phenoname.txt')['phenoname'].tolist())
