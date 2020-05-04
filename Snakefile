#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:01:16 2020
@author: gsd818
"""

import os
import pandas as pd

configfile: "config.yaml"

## --------------------------------------------------------------------------------
## global parameters
PCA_K = range(2,12)
PANNAMES=["euras10", "euras50"]

wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+"

## --------------------------------------------------------------------------------
## targets
rule all:
    input:
        expand("annotate_genes/neo.impute-{panel}-pcadapt-ALLk.manhattan.plot.pdf", panel=PANNAMES)
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

## --------------------------------------------------------------------------------
## rules

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
        expand("pca/{prefix}-{panel}-k{k}.pcadapt", k=PCA_K, allow_missing=True),
        "pca/plots/{prefix}-{panel}.pcadapt.screeplot.pdf"
    params:
        min_k=min(PCA_K),
        max_k=max(PCA_K)
    shell:
        """
        Rscript pca_k.R {input} {params.min_k} {params.max_k} {wildcards.panel}
        """

rule plot_pca:
    input:
        "pca/{prefix}-{panel}-k{k}.pcadapt"  # TODO do you mean to hard code k=8 here?
    output:
        "pca/plots/{prefix}-{panel}-{test}-k{k}.pca.plot.pdf",
        "annotate_genes/{prefix}-{panel}-{test}-k{k}.vals.tsv",
    shell:
        """
        Rscript pcadapt_plot.R {wildcards.prefix}-{wildcards.panel} {input} {wildcards.test}
        """

rule plot_man:
    input:
        "annotate_genes/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{prefix}-{panel}-{test}-k{k}.manhattan.plot.png"
    shell:
        """
        Rscript manplot_pcadapt.R {input} {wildcards.prefix}-{wildcards.panel}-{wildcards.test} {wildcards.k}
        """

rule png2pdf:
    input:
        expand("annotate_genes/{prefix}-{panel}-{test}-k{k}.manhattan.plot.png", k=PCA_K, allow_missing=True)
    output:
        "annotate_genes/{prefix}-{panel}-{test}-ALLk.manhattan.plot.pdf"
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

rule polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz"),
        popfile="paneldir/{level}-euras.clusters.acf.gz",
        lbd=config['lbd']
    output:
        freqs="UKBiobank/data/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python acf2ukbfreq_byP.py  -a {input.popfile} -g {input.infile} -o {output.freqs}
        cat <(head -1 {output.freqs}) <(tail -n+2 {output.freqs} | sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        python partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p 5e-08
        python extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p 0.00001
        """

checkpoint phenos_ids:
    input:
        "phenoname.txt"
    output:
        "phenoname_polyscores.txt"
    shell:
        """
        cd UKBiobank/data
        cat <(echo 'phenoname') <(find gwasfreqs_candidates-pops* -type f | xargs wc -l | awk '$1 > 3' | sed 's/gwasfreqs_candidates-pops-//' | sed 's/.tsv//'| awk '{print $2}') > {output}
        """

rule polyAdapt_qx:
    input:
        pd.read_table('phenoname_polyscores.txt')['phenoname'].tolist()
        neut="UKBiobank/data/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/gwasfreqs_candidates-{level}-{pheno}.tsv",
        gbr="paneldir/gbr.tsv.gz"
    output:
        qx="UKBiobank/selection_UKBV2/QX_report-{level}-{pheno}.txt",
        qxfm="UKBiobank/selection_UKBV2/QX_fm_report-{level}-{pheno}.txt",
        scores="UKBiobank/selection_UKBV2/Genscores-{level}-{pheno}.txt"
    shell:
        """
        Rscript CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        Rscript CalcQX_GBR-matched_Alba.R -w {input.candi} -e {input.neut} -a {input.gbr} -n 1000 -m {output.qxfm} -j 1000
        """
