#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""
import pandas as pd
import csv
import os

configfile: "config.yaml"

PCA_K = range(1,21)
SET=['eur']# ,'gbr'
STEPS=[5e6]

rule all:
    input:
        expand("1000genomes/pca/plots/{set}.pcadapt.screeplot.pdf", set='all')


rule get_plink:
    input:
        vcf="1000genomes/1000g.annot.vcf.gz",
        # ind="1000genomes/{set}.txt"
    output:
        expand("1000genomes/plink/{set}.{format}", format=['bed','bim','fam'],allow_missing=True)
    shell:
        "plink --vcf {input.vcf} --maf 0.05 --make-bed --out 1000genomes/plink/{wildcards.set}" # --keep {input.ind}


rule run_pcadapt:
    input:
        "1000genomes/plink/{set}.bed"
    output:
        # expand("1000genomes/pca/{set}-k{k}.pcadapt", k=PCA_K, allow_missing=True),
        "1000genomes/pca/{set}-k20.pcadapt",
        "1000genomes/pca/plots/{set}.pcadapt.screeplot.pdf"
    params:
        min_k=min(PCA_K),
        max_k=max(PCA_K)
    shell:
        """
        Rscript scripts/pca_k_1000genomes.R {input} {params.min_k} {params.max_k} {wildcards.set}
        """

#rule corr_degas:
# todo :
# 1. Rscripts scripts/corr_pearson_degas.R
# 2. Rscripts scripts/degas_rand_corr.R	
# 3. Rscript scripts/pcadapt_plot_1000g.R : plot PCA coloured by popid
# 4.
rule randcor:
    input:
        expand("degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv", k=6, allow_missing=True)
    output:
        expand("degas/{panel}/randpval.LD{step}.loadings.{load}.tsv",load=[1,2], allow_missing=True)
    shell:
        "Rscript scripts/degas_rand_corr.R {input} {wildcards.panel} {wildcards.step}"

rule sigtab:
    input:
        p="degas/{panel}/randpval.LD{step}.loadings.{load}.tsv",
        cor=expand("degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv",k=6, allow_missing=True)
    output:
        "degas/{panel}/sig.loadings.{load}.txt"
    shell:
        "Rscript scripts/SigDegas.R "
        " -c {input.cor}"
        " -p {input.p}"
        " -o {output}"


