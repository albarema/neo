# https://www.nature.com/articles/ejhg2015269
GWASCAT_PVALUE = "5e-8"

rule all:
    input:
        "gwascat/gwas_catalog_targets.tsv"

rule gwascat_fetch:
    output:
        tsv="gwascat/gwas_catalog_v1.0.2-associations_e100_r2020-06-04.tsv",
    shell:
        "wget --quiet -O {output} https://www.ebi.ac.uk/gwas/api/search/downloads/alternative/"


rule gwascat_genome_wide_significant:
    input:
        "gwascat/gwas_catalog_v1.0.2-associations_e100_r2020-06-04.tsv",
    output:
        temp("gwascat/gwas_catalog_significant.tsv"),
    params:
        col=lambda wildcards, input: open(input[0]).readline().split("\t").index("P-VALUE") + 1,
    shell:
        "awk -F'\\t' 'NR==1 || ${params.col} < {GWASCAT_PVALUE} {{print $0}}' {input} > {output}"


rule gwascat_ontology:
    input:
        "gwascat/gwas_catalog_significant.tsv",
    output:
        "gwascat/gwas_catalog_targets.tsv",
    shell:
        "python scripts/gwascat_ontology.py"
        " --tsv {input}"
        " --output {output}"
