import pandas as pd
import os
configfile: "config.yaml"

#PANNAMES=['eurasGP08pass','europeGP08pass', 'eur', 'gbr']
PANNAMES=['europe', 'europeGP09pass']
PCA_K=range(2,7)
LOADS=range(1,3)
#STEPS=[1,1e6,1e2, 1e3, 1e4, 1e5]
STEPS=[1e6,5e6, 1e5]

rule all:
    input:
        # expand("pca/plots/{panel}-maf01/neo.impute-{panel}-maf01.pcadapt.screeplot.pdf", panel=PANNAMES),
        # expand("plots/degas/degas.bar.{panel}.k{k}.png", k=6, panel=PANNAMES),
        #expand("degas/{panel}/randpval.LD.loadings.1.tsv", panel=PANNAMES, load=[1,2])
        # significant correlations info
        expand("degas/{panel}/sig.LD{step}.loadings.{load}.{panel}.txt", panel=PANNAMES, load=[1,2], step=STEPS),
        # Bootstrapping
        #expand("degas/{panel}/bootsrandpval.loadings.{load}.tsv",panel=PANNAMES, load=[1,2])
        expand("degas/{panel}/bootsrandpval.LD{step}.loadings.{load}.tsv",panel=PANNAMES, load=[1,2], step=STEPS)



# pcascan.smk needs to be run first
rule get_plink_euras:
    input:
        bed="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.bed",
        bim="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.bim",
        fam="/science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{prefix}.1000g.fam",
        inds=expand("{pan}-inds.txt", pan='europe', allow_missing=True)
    output:
        "plink/{prefix}-{panel}.bed",
        "plink/{prefix}-{panel}.bim",
        "plink/{prefix}-{panel}.fam"
    params:
        miss=[0.5]
    shell:
        "plink --bfile /science/willerslev/scratch/neo/sg_freeze_20200615/1000g/plink/{wildcards.prefix}.1000g --keep-fam {input.inds} " #--maf 0.01
        "  --geno {params.miss} --make-bed --out plink/{wildcards.prefix}-{wildcards.panel}"

# rule run_pca:
#     input:
#         "plink/{prefix}-{panel}.bed"
#     output:
#         expand("pca/{panel}/{prefix}-{panel}-k{k}.pcadapt", k=PCA_K, allow_missing=True),
#         "pca/plots/{panel}/{prefix}-{panel}.pcadapt.screeplot.pdf"
#     params:
#         min_k=min(PCA_K),
#         max_k=max(PCA_K)
#     shell:
#         """
#         Rscript scripts/pca_k.R {input} {params.min_k} {params.max_k} {wildcards.panel}
#         """
#
# rule degas: #DO NOT RUN, maybe append to corr_pcs based on k
#     input:
#         pca="pca/{panel}/neo.impute-{panel}-k{k}.pcadapt",
#         degas="DeGAs/contrib_var.csv",
#         bim="plink/neo.impute-{panel}.bim",
#         fam="plink/neo.impute-{panel}.fam"
#     output:
#         plot="plots/degas/degas.bar.{panel}.k{k}.png",
#         tsv="degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv",
#     shell:
#         "Rscript scripts/corr_pearson_degas.R"
#          " -p {input.pca}"
#          " -d {input.degas}"
#          " -i {input.bim}"
#          " -f {input.fam}"
#          " -s {wildcards.panel}"
#          " -o {output.plot}"

#Rscript scripts/corr_pearson_degas.R -p pca/europe/neo.impute-europe-k6.pcadapt -d DeGAs/contrib_var.csv -i plink/neo.impute-europe.bim -f plink/neo.impute-europe.fam -s europe -o plots/degas/degas.bar.europe.k6.png

rule randcor_blocks:
    input:
        expand("degas/{panel}-pcadapt-vs-degas-loadings-k{k}.tsv", k=6, allow_missing=True)
    output:
        expand("degas/{panel}/randpval.LD{step}.loadings.{load}.tsv",load=[1,2], allow_missing=True)
    shell:
        "Rscript scripts/degas_rand_corr.R {input} {wildcards.panel} {wildcards.step}"

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
        expand("degas/{panel}/bootsrandpval.LD{step}.loadings.{load}.tsv",load=[1,2], allow_missing=True)
    shell:
        "Rscript scripts/boostrap_degas_blocks.R {input} {wildcards.panel} {wildcards.step}"

rule sigtab:
    input:
        p="degas/{panel}/randpval.LD{step}.loadings.{load}.tsv",
        cor="degas/{panel}/topcorr.loadings{load}.tsv"
    output:
        "degas/{panel}/sig.LD{step}.loadings.{load}.{panel}.txt"
    shell:
        "Rscript scripts/SigDegas.R "
        " -c {input.cor}"
        " -p {input.p}"
        " -o {output}"


