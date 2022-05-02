#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""

import pandas as pd
import csv

configfile: "config.yaml"

# ALL COMBINATIONS
#PHENOS=pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist()
PANAMES=['euras','europe', 'eurasGP08pass', 'europeGP08pass']
LEVELS=['pops']
# Phenos pass > 10 snps
PHENOS=pd.read_table('phenotypes_qx_euras.txt')['phenotype'].tolist()


## rules
rule betas:
    input:
        # candidates
        # expand("UKBiobank/data/europe/gwasfreqs_candidates-pops-{pheno}.tsv", level=LEVELS, pheno=PHENOS),
        # Frequency matched qx
        expand("UKBiobank/selection_UKBV2/{panel}/QX_fm_report-{level}-{pheno}.txt",panel=['euras','europe'], level='pops', pheno=PHENOS),
        # Genscores
        # expand("UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}.txt",panel= 'europe', level='pops', pheno=PHENOS),
        # "UKBiobank/europe/polyscores_sigtraits_5e-8.txt"


## global parameters

# rule polar_beta:
#     input:
#         os.path.join(config['uk_dir'],"{pheno}.gwas.imputed_v3.both_sexes.tsv.gz")
#     output:
#         ungz=temp(os.path.join(config['uk_dir'], "{pheno}.flipped.byP")),
#         gz=os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz")
#     shell:
#         """
#         python scripts/FlipMinorUKB.py -i {input} -o {output.ungz}
#         cat <(head -1 {output.ungz}) <(tail -n+2 {output.ungz}| sort -k2,2 -k3,3g) | bgzip -c > {output.gz}
#         tabix -s 2 -b 3 -e 3 {output.gz}
#         """


checkpoint polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/{panel}-{level}_Filtered.acf.gz",
        lbd=config['lbd']
    output:
        #freqs="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP_Alba.py  -a {input.popfile} -g {input.infile} -o UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        cat <(head -1 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv) <(tail -n+2 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv| sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        python scripts/partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p 5e-08
        python scripts/extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p 0.00001
        """

rule polyAdapt_qx:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
    output:
        qx="UKBiobank/selection_UKBV2/{panel}/QX_report-{level}-{pheno}.txt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """

rule polyAdapt_gbr:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        gbr="paneldir/gbr.tsv.gz"
    output:
        qxfm="UKBiobank/selection_UKBV2/{panel}/QX_fm_report-{level}-{pheno}.txt",
    shell:
        """
        Rscript scripts/CalcQX_GBR-matched_Alba.R -w {input.candi} -e {input.neut} -a {input.gbr} -n 1000 -m {output.qxfm} -j 1000
        """

# For plotting run polyadapt_plotting.smk





