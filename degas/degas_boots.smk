import pandas as pd
import os
configfile: "config.yaml"

#PANNAMES=['eurasGP08pass','europeGP08pass', 'eur', 'gbr']
PANNAMES=['europe', 'euras']
STEPS=[5e6]
PSEUDO=5000

rule all:
    input:
        # Bootstrapping
        #expand("degas/{panel}/sig.boots.loadings.{load}.{panel}.txt",panel=PANNAMES, load=[1,2]),
        expand("degas/{panel}/sig.boots.s{pseudo}.LD{step}.loadings.{load}.{panel}.txt",panel=['europe'], load=[1,2], step=STEPS, pseudo=PSEUDO)


rule boots:
    input:
        expand("degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv", k=6, allow_missing=True)
    output:
        expand("degas/{panel}/bootsrandpval.loadings.{load}.tsv",load=[1,2], allow_missing=True)
    shell:
        "Rscript scripts/boostrap_degas.R {input} {wildcards.panel}"

rule boots_blocks:
    input:
        expand("degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv", k=6, allow_missing=True)
    output:
        expand("degas/{panel}/bootsrandpval.s{pseudo}.LD{step}.loadings.{load}.tsv",load=[1,2], allow_missing=True)
    shell:
        "Rscript scripts/boostrap_degas_blocks.R {input} {wildcards.panel} {wildcards.step} {wildcards.pseudo}"

rule sigtab_block:
    input:
        p="degas/{panel}/bootsrandpval.s{pseudo}.LD{step}.loadings.{load}.tsv",
        cor="degas/{panel}/topcorr.loadings{load}.tsv"
    output:
        "degas/{panel}/sig.boots.s{pseudo}.LD{step}.loadings.{load}.{panel}.txt"
    shell:
        "Rscript scripts/SigDegas.R "
        " -c {input.cor}"
        " -p {input.p}"
        " -o {output}"

rule sigtab:
    input:
        p="degas/{panel}/bootsrandpval.loadings.{load}.tsv",
        cor="degas/{panel}/topcorr.loadings{load}.tsv"
    output:
        "degas/{panel}/sig.boots.loadings.{load}.{panel}.txt"
    shell:
        "Rscript scripts/SigDegas.R "
        " -c {input.cor}"
        " -p {input.p}"
        " -o {output}"




