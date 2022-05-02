import pandas as pd
import csv

configfile: "config.yaml"

# ALL COMBINATIONS
PANAMES=['europe', 'euras','europeGP08pass', 'eurasGP08pass']

## rules
rule all_plots:
    input:
        # qx stats plots
        expand("UKBiobank/Qx-pvalue/{panel}/Qx_allvals_sigtraits.txt", panel=PANAMES),
        # PRS plots
        expand("UKBiobank/{panel}/polyscores_sigtraits_5e-8.txt", panel=PANAMES)

rule qx_distr:
    input:
        "phenotypes_qx_{panel}.txt"
    output:
        allqx="UKBiobank/Qx-pvalue/{panel}/Qx_allvals_alltraits.txt",
        sigqx="UKBiobank/Qx-pvalue/{panel}/Qx_allvals_sigtraits.txt"
    shell:
        """
        Rscript scripts/QX_distr_plots.R {input} {output.allqx} {output.sigqx} {wildcards.panel}
        """

rule pol_scores:
    input:
        cols="neo.impute.filter.clusterInfo.tsv",
        qx="UKBiobank/Qx-pvalue/{panel}/Qx_allvals_sigtraits.txt",
        categ="traits_descript_categories.tsv"
    output:
        "UKBiobank/{panel}/polyscores_sigtraits_5e-8.txt"
    shell:
        """
        Rscript scripts/scoresPlot.R {input.qx} {input.cols} {input.categ} {wildcards.panel}
        """


