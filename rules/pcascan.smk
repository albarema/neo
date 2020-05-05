#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""

import os
## rules

## global parameters
PCA_K = range(2,12)
PANNAMES=["euras10", "euras50"]


rule get_plink_euras:
    input:
        bed="plink/{prefix}.bed",
        bim="plink/{prefix}.bim",
        fam="plink/{prefix}.fam",
        inds="euras.inds.txt"
    output:
        "plink/{prefix}-{panel}.bed",
        "plink/{prefix}-{panel}.bim",
        "plink/{prefix}-{panel}.fam"
    params:
        miss=[0.1]
    shell:
        "plink --bfile plink/{wildcards.prefix} --keep-fam {input.inds} --maf 0.05 --geno {params.miss} --make-bed "
        "      --out plink/{wildcards.prefix}-{wildcards.panel}"

rule run_pca:
    input:
        "plink/{prefix}-{panel}.bed"
    output:
        expand("pca/{panel}/{prefix}-{panel}-k{k}.pcadapt", k=PCA_K, allow_missing=True),
        "pca/plots/{panel}/{prefix}-{panel}.pcadapt.screeplot.pdf"
    params:
        min_k=min(PCA_K),
        max_k=max(PCA_K)
    shell:
        """
        Rscript pca_k.R {input} {params.min_k} {params.max_k} {wildcards.panel}
        """

rule plot_pca:
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.pcadapt"  # TODO do you mean to hard code k=8 here?
    output:
        "pca/plots/{panel}/{prefix}-{panel}-{test}-k{k}.pca.plot.pdf",
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv",
    shell:
        """
        Rscript pcadapt_plot.R {wildcards.prefix}-{wildcards.panel} {input} {wildcards.test} {wildcards.panel}
        """
rule plot_man:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.manhattan.plot.png"
    shell:
        """
        Rscript manplot_pcadapt.R {input} {wildcards.prefix}-{wildcards.panel}-{wildcards.test} {wildcards.k} {wildcards.panel}
        """

rule png2pdf:
    input:
        expand("annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.manhattan.plot.png", k=PCA_K, allow_missing=True)
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-ALLk.manhattan.plot.pdf"
    shell:
        """
        Rscript png_2_pdf.R {output} {input}
        """
rule annotate_genes:
    input:
        posval="annotate_genes/{prefix}-{panel}.{test}.vals.tsv",
        hg19=config['hg19_file']
    output:
        posval="annotate_genes/{prefix}-{panel}.{test}.annotated.vals.tsv",
        top="annotate_genes/{prefix}-{panel}.{test}.annotated.top.vals.tsv",
        merged="annotate_genes/{prefix}-{panel}.{test}.merged.top.vals.tsv",
        plot="annotate_genes/{prefix}-{panel}.{test}.manhattan_annotated.png"
    shell:
        """
        python IDprocess_human.py {input.posval} {input.hg19} {output.posval}
        Rscript get.top20.R {output.posval} {output.top}
        python merge_for_manhattan_Alba.py {output.top} {input.posval} {output.merged}
        Rscript manhattan_with_genes.R {output.merged} {output.plot}
        """
