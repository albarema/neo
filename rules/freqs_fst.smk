PCA_K = range(2,12)

import itertools
import pandas as pd

wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
 
pops=pd.read_table('popdir/neopops_list.txt')['files'].tolist()

rule fst_all_pops:
    input:
        ["Fst/{}_vs_{}.log".format(*pair) for pair in itertools.combinations(pops, 2)]
  

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

rule annotate_genes:
    input:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.AnnotatedGenes.tsv"
    shell:
        "Rscript annotate_genes.R {input} {output}" #TODO: run in candy

rule fst_plink:
    input:
        bfile=expand("plink/{prefix}-{panel}.{ext}", ext=['bed', 'bim', 'fam'], allow_missing=True), # TODO change FID from .fam
        panel="Fst/euras.pops.panelfile.fst.txt",
        rmid="euras.rmid.panelfile.txt"
    output:
        "Fst/{prefix}-{panel}-{test}.{ext}"
    shell:
        """
        #plink --bfile plink/{wildcards.prefix}-{wildcards.panel} --keep-fam {input.rmid} --make-bed --out Fst/{wildcards.prefix}-{wildcards.panel}
        plink --bfile Fst_plink/{wildcards.prefix}-{wildcards.panel} --fst --within {input.panel} --out Fst_plink/{wildcards.prefix}-{wildcards.panel}-{wildcards.test}
        """

rule get_fst:
    input:
        freqs="neo.impute-euras-qx-popFreqs.tsv",
        gbr="gbr-qx-popFreqs.tsv",
        pops="euras.pops.panelfile.txt"
    output:
        fstgbr="Fst/Fst_GBR_neoclusters.tsv",
        fstneo="Fst/Fst_neoclusters.tsv"
    shell:
        "Rscript CalcFst.R {input.freqs} {input.gbr} {input.pops} {output.fstgbr} {output.fstneo}"
        

rule fst_between_pops:
    input:
        "vcf/euras.annot.vcf"
    output:
        "Fst/{pop1}_vs_{pop2}.log"
    shell:
        "/willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcftools "
        " --vcf {input}"
        " --weir-fst-pop popdir/{wildcards.pop1}.txt "
        " --weir-fst-pop popdir/{wildcards.pop2}.txt "
        " --out Fst/{wildcards.pop1}_vs_{wildcards.pop2} "
