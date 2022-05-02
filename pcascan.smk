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
PANNAMES=["eurasGP08pass", "europeGP08pass"] # TEST DIFF MISSINIGNESS COEFFICIENTS
TEST=['pcadapt'] # pcadapt

rule plot:
    input:
        #
        #expand("plink/neo.impute-{panel}GP08pass.bed", panel=['euras', 'europe']),
        # Get manhattan and pca plots
        # expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.vals.tsv", panel=PANNAMES, k=PCA_K)
        expand("annotate_genes/{panel}/neo.impute-{panel}-{test}-ALLk.manhattan.plot.pdf", panel=PANNAMES, test=TEST),
        # pcadapt values distribution
        #expand("plots/{panel}/neo.impute-{panel}-{test}-ALLk.chisq.plot.pdf", panel=PANNAMES, test=TEST),
        # pcadapt get PC2 pvalues which separates hg
        #expand("annotate_genes/{panel}/neo.impute-{panel}-{test}-k3-PC2.vals.tsv",panel=PANNAMES, test=TEST)
        # Get correlation pca loadings and degas loadings
        # "degas/corr_pcs.txt"

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

# B. SKIP
# # Extra filters to test: it removes all sites that don't passed that cutoff (sites >10% individuals with GP < 0.8)
# rule rm_inds_GP:
#     input:
#         bed="plink/{prefix}-{panel}.bed",
#         bim="plink/{prefix}-{panel}.bim",
#         fam="plink/{prefix}-{panel}.fam",
#         rmid="{panel}.rsid.GP09.txt" #gp=08/09
#     output:
#         "plink/{prefix}-{panel}GP09pass.bed",
#         "plink/{prefix}-{panel}GP09pass.bim",
#         "plink/{prefix}-{panel}GP09pass.fam"
#     shell:
#         "plink --bfile plink/{wildcards.prefix}-{wildcards.panel} --extract {input.rmid} "
#          " --make-bed --out plink/{wildcards.prefix}-{wildcards.panel}GP09pass"

rule run_pca:
    input:
        "plink/{prefix}-{panel}.bed"        
    output:
        expand("pca/{panel}/{prefix}-{panel}-k{k}.{test}", k=PCA_K, allow_missing=True),
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
         "pca/{panel}/{prefix}-{panel}-k{k}.{test}" # TODO do you mean to hard code k=8 here?
    output:
        # "pca/plots/{panel}/{prefix}-{panel}-{test}-k10.pca.plot.pdf",
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv",
    params:
        max_k=max(PCA_K)
    shell:
        """
        Rscript scripts/pcadapt_plot.R {wildcards.prefix}-{wildcards.panel} {input} {wildcards.test} {wildcards.panel} {params.max_k}
        """

rule pca_pvals: # ony for pcadapt, for the component wise we already ave the pvalues
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.{test}"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pcs}.vals.tsv",
    params:
        pc_k="{pcs}"
    shell:
        """
        Rscript scripts/GetPvalues_pcadapt.R {input} {params.pc_k} {wildcards.panel} {output}
        """

rule distr_pcadapt:
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.{test}", prefix=['neo.impute'] , allow_missing=True),
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

rule annotate_genes: # OLD ANNOTATION
    input:
        posval="annotate_genes/{prefix}-{panel}-{test}.vals.tsv",
        hg19=config['hg19_file']
    output:
        posval="annotate_genes/{prefix}-{panel}-{test}.annotated.vals.tsv",
        top="annotate_genes/{prefix}-{panel}-{test}.annotated.top.vals.tsv",
        merged="annotate_genes/{prefix}-{panel}-{test}.merged.top.vals.tsv",
        plot="annotate_genes/{prefix}-{panel}-{test}.manhattan_annotated.png"
    shell:
        """
        python scripts/IDprocess_human.py {input.posval} {input.hg19} {output.posval}
        Rscript scripts/get.top20.R {output.posval} {output.top}
        python scripts/merge_for_manhattan_Alba.py {output.top} {input.posval} {output.merged}
        Rscript scripts/manhattan_with_genes.R {output.merged} {output.plot}
        """
