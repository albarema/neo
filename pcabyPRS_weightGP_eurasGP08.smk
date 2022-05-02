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
#PHENOS=pd.read_table('phenotypes_qx_eurasGP08pass_sig_BC.txt')['phenotype'].tolist()
PHENOS='G35'
PANNAME=['euras', 'europe']
## --------------------------------------------------------------------------------
#####
rule all:
    input:
        #expand("UKBiobank/data/allind/weightedFreqs-inds-{pheno}.tsv.gz", pheno=PHENOS)
        expand("UKBiobank/data/allind/weightedFreqs_candidates-inds-{pheno}-2e1.tsv", pheno=PHENOS),
        # genscores
        expand("UKBiobank/selection_UKBV2/{panel}/Genscores-inds-{pheno}.txt", pheno=PHENOS, panel=PANNAME),
        # pcaplot
        expand("pca/plots/{panel}/neo.impute-{panel}-k{k}-pca.plot_by{pheno}_2e-1.nolabels.png", pheno=PHENOS, k=2, panel=PANNAME),
        # maps
        expand("plots/{panel}/neo.impute-{panel}.by-{pheno}.map.plot.pdf", pheno=PHENOS, k=2, panel=PANNAME)
#
--------------------------------------------------------------------------------
rule weighted_effect:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/euras-pops.acf.gz",
        vcf="vcf/euras-qx.vcf.gz",
    output:
        tsv=temp("UKBiobank/data/allind/weightedFreqs-{level}-{pheno}.tsv"),
        gzip="UKBiobank/data/allind/weightedFreqs-{level}-{pheno}.tsv.gz",
    shell:
        """
        python scripts/weigthedES_byGP.py -a {input.popfile} -v {input.vcf} -g {input.infile} -o {output.tsv}
        cat <(head -1 {output.tsv}) <(tail -n+2 {output.tsv} | sort -k1,1 -k2,2g) | bgzip -c > {output.gzip}
        tabix -s 1 -b 2 -e 2 {output.gzip}
        """

rule get_can:
    input:
        wes="UKBiobank/data/allind/weightedFreqs-{level}-{pheno}.tsv.gz",
        lbd=config['lbd']
    output:
        can="UKBiobank/data/allind/weightedFreqs_candidates-{level}-{pheno}-2e1.tsv",
    shell:
        "python scripts/partitionUKB_inds_byP.py"
        " -i {input.wes}"
        " -b {input.lbd}"
        " -o {output.can}"
        " -p 2e1"

rule get_PRS:
    input:
        can="UKBiobank/data/allind/weightedFreqs_candidates-{level}-{pheno}-2e1.tsv",
        inds="euras.inds.panelfile.txt",
        sites="paneldir/{panel}-pops_Filtered.acf.gz"
    output:
        prs="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}-2e1.txt"
    shell:
        "Rscript scripts/CalcPRS_allind.R"
        " -p {input.inds}"
        " -w {input.can} "
        " -r {input.sites}"
        " -s {output.prs}"

rule pca_byPolscore:
    """
    Need to change the output file inside the script
    """
    input:
        pca="pca/{panel}/{prefix}-{panel}-k{k}.pcadapt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-inds-{pheno}-2e1.txt"
    output:
        "pca/plots/{panel}/{prefix}-{panel}-k{k}-pca.plot_by{pheno}_2e-1.nolabels.png",
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
