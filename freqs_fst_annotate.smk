# Run annotate_pvalues_PC2 to get annotations for a specific PC
PCA_K = range(2,7)

import itertools
import pandas as pd

wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",

PANNAME=['eurasGP08pass', 'europeGP08pass']
pops=pd.read_table('popdir/neopops_list.txt')['files'].tolist()
GWASCAT_PVALUE = "5e-8"
TEST=['pcadapt', 'pcadapt.cw']

rule all_freqs:
    input:
        expand("annotate_genes/{panel}/neo.impute-{panel}-{test}-k{k}.AnnotatedGenes.tsv", panel=PANNAME, k=3, test=TEST) # run only 1-2 jobs at a time (or won't annotate the regions)
        # Get pops freqs for top regions and plot
        #expand("freqs/{panel}/neo.impute-{panel}-qx-popFreqs.tsv", panel=PANNAME),
        #expand("plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-listTopchr.txt", panel=PANNAME, k=3, prefix='neo.impute'),

        # Annotate top region candidates
        #expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.topRegions.tsv", k=PCA_K, panel=PANNAME),
        #expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.AnnotatedGenes.tsv", panel=PANNAME, k=3) # run only 1-2 jobs at a time (or won't annotate the regions)

        # Fst
        #"fst/euras/all-pops.tsv",
        #["fst/euras/{}_vs_{}.log".format(*pair) for pair in itertools.combinations(pops, 2)],

rule gbr_freqs:
    input:
        "paneldir/gbr.tsv.gz"
    output:
        "freqs/gbr-qx-popFreqs.tsv"
    shell:
        """
        Rscript scripts/Get_GBRFreqs.R {input} {output}
        """

rule clusters_freqs:
    input:
        "paneldir/{panel}-pops_Filtered.acf.gz", # euras-inds_Filtered.acf.gz
    output:
        "freqs/{panel}/{prefix}-{panel}-qx-popFreqs.tsv"
    shell:
        """
        Rscript scripts/Get_ClustersFreqs.R {input} {output}
        """

rule get_topRegion:
    input:
        pca="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.vals.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv",
        #expand("annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}-pc{pc}.topRegions.tsv", pc=pc, allow_missing=True)
    shell:
        """
        Rscript scripts/TopRegions_pcadapt.R {input.pca} {output}
        """

rule freqs_plot:
    input:
        freqs="freqs/{panel}/{prefix}-{panel}-qx-popFreqs.tsv",
        top="annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    params:
        k="{k}"
    output:
        "plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-{test}-listTopchr.txt",
#        "plots/{panel}/pcadapt_k{k}/{prefix}-{panel}-{test}-qx-popFreqs-chr{chrom}.png" # don't uncomment, don't know which chrom are top
    shell:
        """
        Rscript scripts/Plot_TopRegionsFreqs.R {input.pca} {wildcards.panel} {wildcards.test} {wildcards.k} {output}
        """

rule annotate_genes_top:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.AnnotatedGenes.tsv"
    resources:
        ensembl_api=1
    shell:
        "Rscript scripts/annotate_genes.R {input} {output}"


