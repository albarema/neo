PCA_K = range(2,12)
PANNAME=['europeGP08pass', 'eurasGP08pass']#

rule all_freqs:
    input:
#        expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k3-PC2.AnnotatedGenes.tsv", panel='eurasGP08pass'),
        expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k4-PC{pc}.AnnotatedGenes.tsv", pc=[1,2], panel=PANNAME)

rule get_pvals_pcs:
    input:
        "pca/{panel}/{prefix}-{panel}-k{k}.pcadapt"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.vals.tsv"
    shell:
        "Rscript scripts/GetPvalues_pcadapt.R {input} {wildcards.pc} {wildcards.panel} {output}"

rule get_topRegion:
    input:
        pca="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.topRegions.tsv"
    shell:
        """
        Rscript scripts/TopRegions_pcadapt.R {input.pca} {output}
        """

#TODO rule freqs_plots
rule freqs_plot: # don't run this as part of the pipeline, we don't know the {chrom}
    input:
        freqs="freqs/{panel}/{prefix}-{panel}-qx-popFreqs.tsv",
        top="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.topRegions.tsv"
    params:
        k="{k}"
    output:
        "plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-qx-popFreqs-chr{chrom}.png"
    shell:
        """
        Rscript scripts/Plot_TopRegionsFreqs.R {input.freqs} {input.top} {params.k} {wildcards.panel} {output}
        """

rule annotate_genes_top:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.topRegions.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-PC{pc}.AnnotatedGenes.tsv"
    shell:
        "Rscript scripts/annotate_genes.R {input} {output}"


