PCA_K = range(2,7)

import itertools
import pandas as pd

wildcard_constraints:
    prefix="[^-]+",
    panel="[^-]+",
    test="[^-]+",
 
pops=pd.read_table('popdir/neopops_list.txt')['files'].tolist()
PANAME=["eurasGP08pass", "europeGP08pass"] #

rule all_freqs:
    input:
        # Annotate jobs
        #expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.topRegions.tsv", k=PCA_K, panel=PANAME),
        expand("annotate_genes/{panel}/neo.impute-{panel}-pcadapt-k{k}.AnnotatedGenes.tsv", k=PCA_K, panel=PANAME), # run only 1-2 jobs at a time (or won't annotate the regions)

        # FST jobs
        #"fst/euras/all-pops.tsv",
        #["fst/euras/{}_vs_{}.log".format(*pair) for pair in itertools.combinations(pops, 2)],

        # Get pops frequencies

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
        "paneldir/{panel}-pops.Filtered.acf.gz", # euras-inds_Filtered.acf.gz
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
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    shell:
        """
        Rscript scripts/TopRegions_pcadapt.R {input.pca} {output}
        """

rule freqs_plot: # DON?T RUN this as part of the pipeline. TO CHOOSE THE "K" BEFORE
    input:
        freqs="freqs/{panel}/{prefix}-{panel}-qx-popFreqs.tsv",
        top="annotate_genes/{panel}/{prefix}-{panel}-pcadapt-k{k}.topRegions.tsv"
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
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.topRegions.tsv"
    output:
        "annotate_genes/{panel}/{prefix}-{panel}-{test}-k{k}.AnnotatedGenes.tsv"
    resources:
        ensembl_api=1
    shell:
        "Rscript scripts/annotate_genes.R {input} {output}" #TODO: run in candy

# rule get_fst:
#     input:
#         freqs="neo.impute-euras-qx-popFreqs.tsv",
#         gbr="gbr-qx-popFreqs.tsv",
#         pops="euras.pops.panelfile.txt"
#     output:
#         fstgbr="Fst/Fst_GBR_neoclusters.tsv",
#         fstneo="Fst/Fst_neoclusters.tsv"
#     shell:
#         "Rscript scripts/CalcFst.R {input.freqs} {input.gbr} {input.pops} {output.fstgbr} {output.fstneo}"
        

rule fst_between_pops:
    input:
        "vcf/euras-fst_filtered.vcf.gz"
    output:
        weir=temp("fst/euras/{pop1}_vs_{pop2}.weir.fst")
    log:
        protected("fst/euras/{pop1}_vs_{pop2}.log")
    shell:
        "/willerslev/users-shared/science-snm-willerslev-                                                                                                                                                             gsd818/Applications/vcftools/bin/vcftools "
        " --gzvcf {input}"
        " --weir-fst-pop popdir/{wildcards.pop1}.txt "
        " --weir-fst-pop popdir/{wildcards.pop2}.txt "
        " --out fst/euras/{wildcards.pop1}_vs_{wildcards.pop2} &> {log}"


rule fst_all_vals:
    input:
        "fst/euras/{pop1}_vs_{pop2}.log" # pop1 do gbr
    output:
        temp("fst/euras/{pop1}_vs_{pop2}.tsv")
    shell:
        "python scripts/get_fst_vals.py {input} {output} {wildcards.pop1} {wildcards.pop2}"


rule fst_all_pops:
    input:
        ["fst/euras/{}_vs_{}.tsv".format(*pair) for pair in itertools.combinations(pops, 2)]
    output:
        "fst/euras/all-pops.tsv"
    shell:
        "echo 'pop1\tpop2\tmean\tweighted\tnum_kept\tnum_snps' > {output};"
        "cat {input} >> {output}"    
