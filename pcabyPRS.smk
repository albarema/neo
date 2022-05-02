
import pandas as pd
configfile: "config.yaml"
#
## --------------------------------------------------------------------------------
##### Wildcards #####
wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+",
    pcs="[^-]+"

## --------------------------------------------------------------------------------
##### global parameters #####
PANNAMES=[ "eurasGP08pass"] # TEST DIFF MISSINIGNESS COEFFICIENTS# "eurasGP09pass", "europeGP09pass","euras","europeGP08pass",
PHENOS=pd.read_table('phenotypes_qx_eurasGP08pass_sig_BC.txt')['phenotype'].tolist()

## --------------------------------------------------------------------------------
#####
rule all:
    input:
#        expand("UKBiobank/data/{panel}/weightedFreqs-inds-50_irnt.tsv.gz", panel='europeGP08pass'),
        # get individual weigthed PRS
        expand("UKBiobank/selection_UKBV2/{panel}/Genscores-inds-{pheno}.txt", panel=PANNAMES, pheno=PHENOS),
        # get pca plots by PRS
        expand("pca/plots/{panel}/neo.impute-{panel}-k2-pca.plot_by{pheno}_5e-8.nolabels.png",panel=PANNAMES, pheno=PHENOS),
        # get map by PRS
        expand("plots/{panel}/neo.impute-{panel}.by-{pheno}.map.plot.pdf", panel=PANNAMES, pheno=PHENOS)


rule weighted_effect:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/{panel}-pops_Filtered.acf.gz",
        vcf="vcf/{panel}-qx.vcf.gz",
    output:
        tsv=temp("UKBiobank/data/{panel}/weightedFreqs-{level}-{pheno}.tsv"),
        gzip="UKBiobank/data/{panel}/weightedFreqs-{level}-{pheno}.tsv.gz",
    shell:
        """
        python scripts/weigthedES_byGP.py -a {input.popfile} -v {input.vcf} -g {input.infile} -o {output.tsv}
        cat <(head -1 {output.tsv}) <(tail -n+2 {output.tsv} | sort -k1,1 -k2,2g) | bgzip -c > {output.gzip}
        tabix -s 1 -b 2 -e 2 {output.gzip}
        """

rule get_can:
    input:
        wes="UKBiobank/data/{panel}/weightedFreqs-{level}-{pheno}.tsv.gz",
        lbd=config['lbd']
    output:
        can="UKBiobank/data/{panel}/weightedFreqs_candidates-{level}-{pheno}.tsv",
    shell:
        "python scripts/partitionUKB_inds_byP.py"
        " -i {input.wes}"
        " -b {input.lbd}"
        " -o {output.can}"
        " -p 5e-08"

rule get_PRS:
    input:
        can="UKBiobank/data/{panel}/weightedFreqs_candidates-{level}-{pheno}.tsv",
    output:
        prs="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}.txt"
    shell:
        "Rscript scripts/CalcPRS.R"
        " -w {input.can} "
        " -s {output.prs}"

rule pca_byPolscore:
    input:
        pca="pca/{panel}/{prefix}-{panel}-k{k}.pcadapt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-inds-{pheno}.txt"
    output:
        "pca/plots/{panel}/{prefix}-{panel}-k{k}-pca.plot_by{pheno}_5e-8.nolabels.png",
    params:
        pheno="{pheno}"
    shell:
        """
        Rscript scripts/pcadapt_byPolscores_plot.R {wildcards.prefix}-{wildcards.panel} {input.pca} {params.pheno} 1 {wildcards.panel} {input.scores}
        """

rule map_byPolscore:
    input:
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-inds-{pheno}.txt",
    output:
        "plots/{panel}/{prefix}-{panel}.by-{pheno}.map.plot.pdf"
    shell:
         "Rscript scripts/map_byPolscores_plot.R"
         " -i {wildcards.prefix}-{wildcards.panel}"
         " -p {wildcards.panel}"
         " -s {input.scores}"
         " -t {wildcards.pheno}"
         " -o {output}"


