PCA_K = range(2,12)

rule all_freqs:
    input:
        expand("annotate_genes/euras50/neo.impute-euras50-pcadapt-k{k}.topRegions.tsv", k=PCA_K),
        "Fst_GBR_neoclusters.png"

rule gbr_freqs:
    input:
        "paneldir/gbr.tsv.gz"
    output:
        "gbr-qx-popFreqs.tsv"
    shell:
        """
        Rscript Get_GBRFreqs.R {input} {output}
        """

rule clusters_freqs:
    input:
        "paneldir/pops-euras.clusters.acf.gz",
    output:
        "neo.impute-euras-qx-popFreqs.tsv"
    shell:
        """
        Rscript Get_ClustersFreqs.R {input} {output}
        """

rule freqs_plot:
    input:
        freqs="neo.impute-euras-qx-popFreqs.tsv",
        pca="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    params:
        k="{k}"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    shell:
        """
        Rscript TopRegions_Freqs_Plot.R {input.freqs} {input.pca} {params.k} {output}
        """

rule get_fst:
    input:
        freqs="neo.impute-euras-qx-popFreqs.tsv",
        gbr="gbr-qx-popFreqs.tsv",
        pops="euras.pops.panelfile.txt"
    output:
          "Fst_GBR_neoclusters.tsv"
    shell:
         "Rscript CalcFst.R {input.freqs} {input.gbr} {input.pops} {output}"


rule fst_plot:
    input:
        "Fst_GBR_neoclusters.tsv"
    output:
        "Fst_GBR_neoclusters.png"
    shell:
        "Rscript manplo_FST.R {input} {output}"
