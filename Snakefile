#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""


import pandas as pd
import csv

configfile: "config.yaml"

## --------------------------------------------------------------------------------
##### Modules #####
include: "rules/pcascan.smk"
include: "rules/polyadapt.smk"
include: "rules/freqs_fst.smk"
# include: "rules/vcf2acf.smk"
## --------------------------------------------------------------------------------
##### Wildcards #####
wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+"

PANNAME=['eurasGP08pass']
## --------------------------------------------------------------------------------
##### Target rules #####
def run_all_rules(_):
    inputs = []

    # this will trigger rules get_plink_euras, run_pca and plot_pca
    # inputs.append("pca/plots/neo.impute-euras.pcadapt-k3.pca.plot.pdf")

    # annotate top genes
    # inputs.extend(expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-ALLk.manhattan.plot.pdf", panel=PANNAME))

    # this will run one of the phenotypes
    #"UKBiobank/data/gwasfreqs-pops-102_irnt.tsv.gz",

    # this will run all the phenotypes
    with open("phenotypes_qx_eurasGP08pass.txt", "w") as fout:
        for pheno in pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist():
            tsv = checkpoints.polyAdapt_freqs.get(pheno=pheno, level='pops', panel='eurasGP08pass').output.candi
            # tsv = "UKBiobank/data/euras/gwasfreqs_candidates-pops-" + str(pheno) + ".tsv"
            with open(tsv) as fin:
                if len(fin.readlines()) > 10:
                    inputs.append("UKBiobank/selection_UKBV2/eurasGP08pass/Genscores-pops-{pheno}.txt".format(pheno=pheno))
                    inputs.append("UKBiobank/selection_UKBV2/eurasGP08pass/QX_fm_report-pops-{pheno}.txt".format(pheno=pheno))
                    fout.write(str(pheno)+ '\n')
    # get qx plots
    # inputs.append("UKBiobank/Qx-pvalue/Qx_allvals_sigtraits.txt")

    # get polygenic scores plots
    # inputs.append("UKBiobank/polyscores_sigtraits_5e-8.txt")
    return inputs



## targets
rule allrules:
    input:
        # expand("UKBiobank/selection_UKBV2/Genscores-inds-{pheno}.txt", pheno=pd.read_table('phenoname_inds.txt')['phenoname'].tolist()),
        run_all_rules

        # expand("annotate_genes/neo.impute-{panel}-pcadapt-ALLk.manhattan.plot.pdf", panel=PANNAMES)
        # this will trigger rules get_plink_euras, run_pca and plot_pca
        # "annotate_genes/all-k.manhattan.plot.pdf"

        # annotate top genes
        # "annotate_genes/neo.impute-euras.pcadapt.manhattan_annotated.png",
        # this will run one of the phenotypes
        # "UKBiobank/data/gwasfreqs-pops-102_irnt.tsv.gz",

        # this will run all the phenotypes. Run this before running rull_all_rules (to speed it up)
        # expand("UKBiobank/data/gwasfreqs-pops-{pheno}.tsv.gz",pheno=pd.read_table('phenoname.txt')['phenoname'].tolist())
        # qx
        # expand("UKBiobank/selection_UKBV2/Genscores-pops-{pheno}.txt",pheno=pd.read_table('phenoname.txt')['phenoname'].tolist())
