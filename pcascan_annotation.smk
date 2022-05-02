
# PCADAPT FOR COMPONENT WISE ANALYSES #
import pandas as pd
configfile: "config.yaml"

##
wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
    level="[^-]+",
    k="[^-]+",
    pc="[^-]+",


wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",

PANNAME=['eurasGP08pass', 'europeGP08pass']
pops=pd.read_table('popdir/neopops_list.txt')['files'].tolist()
GWASCAT_PVALUE = "5e-8"
TEST=['pcadapt.cw']
# PCA_K = range(2,11)
PCA_K=range(2,7)
#K=max(PCA_K)
K=3
PC=range(1,K+1)

rule all_freqs:
    input:
        # manplots
        expand("plots/{panel}/pcadapt_k{k}/neo.impute-{panel}-{test}-pc{pc}-listTopchr.txt",k=K, pc=PC,panel=PANNAME,test=TEST),
        expand("annotate_genes/{panel}/neo.impute-{panel}-{test}-k{k}-pc{pc}.manhattan.plot.png", k=K, pc=PC,panel=PANNAME,test=TEST),
        expand("annotate_genes/{panel}/neo.impute-{panel}-{test}-k{k}-pc{pc}.AnnotatedGenes.tsv", panel=PANNAME, k=K, pc= PC, test=TEST) # run only 1-2 jobs at a time (or won't annotate the regions)
        # Get pops freqs for top regions and plot
        #expand("freqs/{panel}/neo.impute-{panel}-qx-popFreqs.tsv", panel=PANNAME),
        #expand("plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-listTopchr.txt", panel=PANNAME, k=3, prefix='neo.impute'),

        # Annotate top region candidates
        #expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.topRegions.tsv", k=PCA_K, panel=PANNAME),
        #expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.AnnotatedGenes.tsv", panel=PANNAME, k=3) # run only 1-2 jobs at a time (or won't annotate the regions)


rule run_pca:
    input:
        "plink/{prefix}-{panel}.bed"
    output:
        expand("pca/{panel}/{prefix}-{panel}-k{k}.{test}", k=PCA_K, allow_missing=True),
    params:
        min_k=min(PCA_K),
        max_k=max(PCA_K),
        meth=['componentwise']
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

rule distr_pcadapt:
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.{test}",
    output:
        out1="plots/{panel}/{panel}-k{k}-{test}.png",
        out2="plots/{panel}/chi.distr-{panel}-k{k}-{test}.png",
    shell:
        """
        Rscript scripts/pcadapt_distr_chiq.R {wildcards.k}
        """

rule plot_man_cw:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.manhattan.plot.png" #pc ? range(1,k+1)
    shell:
        """
        Rscript scripts/manplot_pcadapt.R {input} {wildcards.prefix}-{wildcards.panel}-{wildcards.test} {wildcards.k} {wildcards.panel} {wildcards.test}
        """

rule get_topRegion:
    input:
        pca="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv",
        #expand("annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv", pc=pc, allow_missing=True)
    shell:
        """
        Rscript scripts/TopRegions_pcadapt.R {input.pca} {wildcards.panel} {wildcards.test} {wildcards.k}
        """

rule freqs_plot:
    input:
        freqs="freqs/{panel}/{prefix}-{panel}-qx-popFreqs.tsv",
        top="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv"
    params:
        k="{k}"
    output:
        "plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-{test}-pc{pc}-listTopchr.txt",
#        "plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-{test}-qx-popFreqs-chr{chrom}.png" # don't uncomment, don't know which chrom are top
    shell:
        """
        Rscript scripts/Plot_TopRegionsFreqs.R {input.freqs} {input.top} {params.k} {wildcards.panel} {wildcards.test} {output}
        """

rule plot_man_locus:
    input:
        vals="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv",
        top="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv"
    output:
        manp=""
    shell:
        """"
        
        """"


rule annotate_genes_top:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.AnnotatedGenes.tsv"
    resources:
        ensembl_api=1
    shell:
        "Rscript scripts/annotate_genes.R {input} {output}"



