import pandas as pd
configfile: "config.yaml"

## global parameters
CHR = range(1,23)
PANNAME=['euras','europe']
PHENOS=pd.read_table('phenotypes_qx_eurasGP08pass_sig_BC.txt')['phenotype'].tolist()

## ALL
rule all:
    input:
        #expand("UKBiobank/111121/allind/weightedFreqs-inds-{pheno}.tsv.gz", pheno=PHENOS)
        #expand("UKBiobank/111121/allind/weightedFreqs_candidates-inds-{pheno}.tsv", pheno=PHENOS),
        expand("UKBiobank/111121/{panel}/Genscores-inds-{pheno}.txt", pheno=PHENOS, panel=PANNAME),
	
## rules
rule get_filtered_vcf:
    input:
        vcf=os.path.join(config['vcf_dir'], "neo.impute.1000g.vcf.gz"),
        rmid="{panel}-inds.txt",
    output:
        "vcf/{panel}-qx.vcf.gz"
    shell:
        """
        bcftools view -q 0.05 -Q 0.95 -S {input.rmid} {input.vcf} | bgzip -c > {output}
        """

rule weighted_effect:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/euras-pops.acf.gz",
        vcf="vcf/neo.impute.1000g.vcf.gz",
    output:
        tsv="UKBiobank/111121/allind/weightedFreqs-{level}-{pheno}.tsv",
    shell:
        """
        python scripts/weigthedES_byGP.py -a {input.popfile} -v {input.vcf} -g {input.infile} -o {output.tsv}
        """

rule gzip:
    input:
        tsv="UKBiobank/111121/allind/weightedFreqs-{level}-{pheno}.tsv",
    output:
        gzip="UKBiobank/111121/allind/weightedFreqs-{level}-{pheno}.tsv.gz",   
    shell:
        """
        cat <(head -1 {input.tsv}) <(tail -n+2 {input.tsv} | sort -k1,1 -k2,2g) | bgzip -c > {output.gzip}
        tabix -s 1 -b 2 -e 2 {output.gzip}
        """ 

rule get_can:
    input:
        wes="UKBiobank/111121/allind/weightedFreqs-{level}-{pheno}.tsv.gz",
        lbd=config['lbd']
    output:
        can="UKBiobank/111121/allind/weightedFreqs_candidates-{level}-{pheno}.tsv",
    shell:
        "python scripts/partitionUKB_inds_byP.py"
        " -i {input.wes}"
        " -b {input.lbd}"
        " -o {output.can}"
        " -p 5e-8"

rule get_PRS:
    input:
        can="UKBiobank/111121/allind/weightedFreqs_candidates-{level}-{pheno}.tsv",
        inds="{panel}.inds.panelfile.txt",
        sites="paneldir/{panel}-pops_Filtered.acf.gz"
    output:
        prs="UKBiobank/111121/{panel}/Genscores-{level}-{pheno}.txt"
    shell:
        "Rscript scripts/CalcPRS_allind.R"
        " -p {input.inds}"
        " -w {input.can} "
        " -r {input.sites}"
        " -s {output.prs}"

rule pca_byPolscore:
    input:
        pca="pca/{panel}/{prefix}-{panel}-k{k}.pcadapt",
        scores="UKBiobank/111121/{panel}/Genscores-{level}-{pheno}.txt"
    output:
        "pca/plots/{panel}/{prefix}-{panel}-k{k}-pca.plot_by{pheno}_5e-8.nolabels.png",
    params:
        pheno="{pheno}"
    shell:
        """
        Rscript scripts/pcadapt_byPolscores_plot.R {wildcards.prefix}-{wildcards.panel} {input.pca} {params.pheno} 1 {wildcards.panel} {input.scores}
        """

rule pca_byPolscore_anders:
    input:
        pca="pca/{panel}/{prefix}-{panel}-k{k}.pcadapt",
        scores="UKBiobank/111121/{panel}/Genscores-{level}-{pheno}.txt"
    output:
        "anders_fisher/plots/{panel}-sigTraits.png",
    params:
        pheno="{pheno}"
    shell:
        """
        Rscript check_PRS_inds.R {wildcards.prefix}-{wildcards.panel} {input.pca} {params.pheno} 1 {wildcards.panel} {input.scores}
        """
