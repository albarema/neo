## --------------------------------------------------------------------------------
## Packages
import itertools
import pandas as pd


## --------------------------------------------------------------------------------

CHR = range(1,23)

pops=pd.read_table('popdir/1000gpops_list.txt')['files'].tolist()
euras=pd.read_table('popdir/neopops_list.txt')['files'].tolist()

## --------------------------------------------------------------------------------
##### Target rules #####
rule all:
    input:
        expand("fst/GBR_euras/GBR_vs_{pop2}.tsv", pop2=euras),
        #expand("1000genomes/chr{chrom}.recode.vcf", chrom=CHR),
        #["fst/1000genomes/{}_vs_{}.log".format(*pair) for pair in itertools.combinations(pops, 2)],
        "fst/GBR_euras/all-pops.tsv"

## --------------------------------------------------------------------------------
## Rules
rule get_filtered_vcf:
    input:
        vcf="/willerslev/datasets/1000_Genomes_phase3_v5a/individual_chromosomes/chr{chrom}.1kg.phase3.v5a.vcf.gz",
        rmid="plink/neo.impute-euras.bim.txt"
    output:
        "1000genomes/chr{chrom}.recode.vcf"
    shell:
        "/willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcftools "
        " --gzvcf {input.vcf} "
        "--positions {input.rmid} " 
        "--recode --out 1000genomes/chr{wildcards.chrom}"

rule merge_vcf:
    input:
        expand("1000genomes/chr{chrom}.recode.vcf", chrom=CHR)
    output:
        protected("1000genomes/1000g.annot.vcf.gz")
    shell:
        """
        export PERL5LIB=/willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/src/perl/
        /willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcf-concat {input} | bgzip -c > {output}
        """

# rule ancient_modern:
#     input:
#         mod="1000genomes/1000g.annot.vcf.gz",
#         anc="vcf/euras.annot.vcf.bgz"
#     output:
#         "euras_100g_all.annot.vcf.gz"
#     shell:
#          """
#          export PERL5LIB=/willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/src/perl/
#          /willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcf-merge {input.anc} {input.mod} | bgzip -c > {output}
#          """


rule fst_between_1000gpops:
    input:
        "1000genomes/{set}.annot.vcf.gz"
    output:
        weir=temp("fst/{set}/{pop1}_vs_{pop2}.weir.fst")
    log:
        log=protected("fst/{set}/{pop1}_vs_{pop2}.log"),
    shell:
        "/willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcftools "
        " --gzvcf {input}"
        " --weir-fst-pop popdir/{wildcards.pop1}.txt "
        " --weir-fst-pop popdir/{wildcards.pop2}.txt "
        " --out fst/{wildcards.set}/{wildcards.pop1}_vs_{wildcards.pop2} &> {log} "

rule fst_all_vals:
    input:
        "fst/{set}/{pop1}_vs_{pop2}.log"
    output:
        temp("fst/{set}/{pop1}_vs_{pop2}.tsv")
    shell:
        "python scripts/get_fst_vals.py {input} {output} {wildcards.pop1} {wildcards.pop2}"


rule fst_all_pops:
    input:
        #["fst/1000genoes/{}_vs_{}.tsv".format(*pair) for pair in itertools.combinations(pops, 2)]
        expand("fst/GBR_euras/GBR_vs_{pop2}.tsv", pop2=euras)
    output:
        "fst/{set}/all-pops.tsv"
    shell:
        "echo 'pop1\tpop2\tmean\tweighted\tnum_kept\tnum_snps' > {output};"
        "cat {input} >> {output}"  
