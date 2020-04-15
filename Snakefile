#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020

@author: gsd818
"""

## --------------------------------------------------------------------------------
## global parameters

PREFIX="neo.impute"
PANEL="euras.neo.impute"

PCA_K = range(8,12)
ADMIXTURE_K = range(7,12)
CHROMS = range(1, 22)

configfile: "phenoname.yaml"


# dirs

GLACTDIR="/willerslev/software/glactools"
VCFDIR="/willerslev/users-shared/science-snm-willerslev-dsw670/projects/neo/data/20200221-impute.filter"
UKDIR="/willerslev/datasets/UKBiobank/NealeV2"

## files
POPFILE="/willerslev/users-shared/science-snm-willerslev-gsd818/neo/impute/UKBiobank/data/popfile_pass.acf"
LBD="/willerslev/datasets/ldetect-data/EUR/fourier_ls-all_nochr.bed"
PANELFILE="euras.clusters.panelfile.txt"



## --------------------------------------------------------------------------------
## output file sets

plink_imp = expand("plink/" + PREFIX + ".{ext}", ext=['bed', 'bim', 'fam'])

plink_euras = expand("plink/" + PANEL + ".{ext}", ext=['bed', 'bim', 'fam'])

pcadapt_k = expand("pca/" + PANEL + ".k{k}.pcadapt", k=PCA_K)


## --------------------------------------------------------------------------------
## targets

rule plink_all:
    input:
        plink_euras

rule pca_all:
    input:
        "pca/plots/" + PANEL + ".pca.plot.pdf"
        
rule all:
    input:
        "UKBiobank/data/gwasfreqs_{pheno}.tsv.gz"
        
## --------------------------------------------------------------------------------
## rules

rule get_plink_euras:
    input:
        bed="plink/" + PREFIX + ".bed",
        bim="plink/" + PREFIX + ".bim",
        fam="plink/" + PREFIX + ".fam",
        inds="euras.inds.txt"
    output:
        plink_euras
    shell:
        """
        plink --bfile plink/{PREFIX} --keep-fam {input.inds} --maf 0.01 --geno 0.01 --make-bed --out plink/{PANEL}
        """
        
rule run_pca:
    input:
        bed="plink/" + PANEL + ".bed",
    output:
        pcadapt_k,
        "pca/plots/" + PANEL + ".pca.screeplot.pdf"
    shell:
        """
        Rscript pca_k.R {input.bed} 
        """
        
rule plot_pca:
    input: 
        pcak="pca/" + PANEL + ".k8.pcadapt"
    output:
        "pca/plots/" + PANEL + ".pca.plot.pdf"
    shell:
        """
        Rscript pcadapt_plot.R {PANEL} {input.pcak}
        """        


rule polyAdapt_freqs_1:
    input:
    	lambda wildcards: config["phenoname"][wildcards.pheno],
        infile=os.path.join(UKDIR,"{pheno}.flipped.byP.gz"),
        popfile=POPFILE,
        lbd=LBD
    output:
    	"UKBiobank/data/gwasfreqs_{pheno}.tsv.gz"
    shell:
        """
        python2 acf2ukbfreq_byP.py  -a {input.popfile} -g {input.infile} -o {output}
        """              