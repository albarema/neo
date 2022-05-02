import pandas as pd
import csv

configfile: "config.yaml"

# ALL COMBINATIONS
#PHENOS=pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist()
PANAMES=['euras']
LEVELS=['pops']

# 17-08-20
PHENOS='G35'


## rules
rule candis:
    input:
        expand("UKBiobank/selection_UKBV2/{panel}/Genscores-pops-nogwt-{pheno}-CI.txt", pheno=PHENOS, panel=PANAMES),
        # expand("pca/plots/{panel}/neo.impute-{panel}-k6-pca.plot_by{pheno}_5e-8.nolabels.png", pheno=PHENOS, panel=PANAMES)


rule polyAdapt_qx:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-pops-G35.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
    output:
        qx="UKBiobank/selection_UKBV2/{panel}/QX_report-{level}-{pheno}-CI.txt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}-CI.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba_CI.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """

rule pca_byPolscore:
    input:
        pca="pca/{panel}/{prefix}-{panel}-k6.pcadapt",
    output:
        "pca/plots/{panel}/{prefix}-{panel}-k6-pca.plot_by{pheno}_5e-8.nolabels.png",
    params:
        pheno="{pheno}"
    shell:
        """
        Rscript scripts/pcadapt_byPolscores_plot.R {wildcards.prefix}-{wildcards.panel} {input.pca} {params.pheno} 1 {wildcards.panel}
        """
