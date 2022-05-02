#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""

import pandas as pd
import os
configfile: "config.yaml"

##
wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+",
    pcs="[^-]+"

## global parameters
# PCA_K = range(2,11)
PCA_K=range(2,7)
PANNAMES=["eurasModern", "europeModern"] 
TEST=['pcadapt'] # pcadapt

rule plot:
    input:
        # Get manhattan and pca plots
        expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.vals.tsv", panel=PANNAMES, k=PCA_K)



        # Get pca coloured by polygenic scores
        
rule get_plink_euras:
    input:
        bed="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.bed",
        bim="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.bim",
        fam="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.fam",
        inds="{panel}-inds.txt"
    output:
        "plink/{prefix}-{panel}.bed",
        "plink/{prefix}-{panel}.bim",
        "plink/{prefix}-{panel}.fam"
    params:
        miss=[0.5]
    shell:
        "plink --bfile /science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{wildcards.prefix}.1000g --keep-fam {input.inds} --maf 0.05 "
        "   --geno {params.miss} --make-bed --out plink/{wildcards.prefix}-{wildcards.panel}"


rule run_pca:
    input:
        "plink/{prefix}-{panel}.bed"        
    output:
        expand("pca/{panel}/{prefix}-{panel}-k{k}.{test}", k=PCA_K, allow_missing=True)
    params:
        min_k=min(PCA_K),
        max_k=max(PCA_K),
        meth=['mahalanobis'] # mahalanobis
    shell:
        """
        Rscript scripts/pca_k.R {input} {params.min_k} {params.max_k} {wildcards.panel} {params.meth} {wildcards.test}
        """

rule plot_pca:
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.{test}"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv",
    params:
        max_k=max(PCA_K)
    shell:
        """
        Rscript scripts/pcadapt_plot.R {wildcards.prefix}-{wildcards.panel} {input} {wildcards.test} {wildcards.panel} {params.max_k}
        """

rule distr_pcadapt:
    input:
        expand("pca/{panel}/{prefix}-{panel}-k{k}.{test}", prefix=['neo.impute'], allow_missing=True),
    output:
        out1="plots/{panel}/{panel}-k{k}-{test}.png",
        out2="plots/{panel}/chi.distr-{panel}-k{k}-{test}.png",
    shell:
        """
        Rscript scripts/pcadapt_distr_chiq.R {wildcards.k}
        """

# rule chiq2pdf:
#     input:
#         expand("plots/{panel}/{prefix}-{panel}-{test}-k{k}-chisq_chi2stat.png", k=PCA_K, allow_missing=True)
#     output:
#         "plots/{panel}/{prefix}-{panel}-{test}-ALLk.chisq.plot.pdf"
#     shell:
#         """
#         Rscript scripts/png_2_pdf.R {output} {input}
#         """

rule plot_man:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.manhattan.plot.png"
    shell:
        """
        Rscript scripts/manplot_pcadapt.R {input} {wildcards.prefix}-{wildcards.panel}-{wildcards.test} {wildcards.k} {wildcards.panel} {wildcards.test}
        """

rule png2pdf:
    input:
        expand("annotate_genes/{panel}/{prefix}-{panel}-pcadapt-k{k}.manhattan.plot.png", k=PCA_K, allow_missing=True)
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-pcadapt-ALLk.manhattan.plot.pdf"
    shell:
        """
        Rscript scripts/png_2_pdf.R {output} {input}
        """