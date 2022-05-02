configfile: "config.yaml"
import pandas as pd

PANNAMES=['europe'] # ,'euras', 'europe'
# ALLPHENOS=pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist()
# Phenos euras and europe:
# FPHENOS=pd.read_table('phenotypes_qx_euras.txt')['phenotype'].tolist()
# Phenos GP08:
# FPHENOS=pd.read_table('phenotypes_qx_eurasGP08pass.txt')['phenotype'].tolist()

# EURAS, no genotypes set to missing! but sites are removed if more than 10% of the individuals have GP < 0.8
## ALL
## --------------------------------------------------------------------------------
##### Target rules #####
rule all_pops:
    input:
        # Get candidates
        #expand("UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",level='inds', panel=PANNAMES, pheno=ALLPHENOS),
        # Get Genscores
        # expand("UKBiobank/selection_UKBV2/{panel}/QX_fm_report-{level}-{pheno}.txt", level='inds', panel=PANNAMES, pheno=FPHENOS),
        expand("UKBiobank/selection_UKBV2/{panel}/QX_report-{level}-{pheno}-{pval}.txt", level='inds', panel=PANNAMES, pheno='G35', pval=['1e5', '1e3', '0.1'])


rule polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/{panel}-qx-{level}.acf.gz",
        lbd=config['lbd']
    output:
        #freqs="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP_Alba.py -a {input.popfile} -g {input.infile} -o UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        cat <(head -1 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv) <(tail -n+2 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv| sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        python scripts/partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p 5e-08
        python scripts/extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p 0.00001
        """

# which phenos: python scripts/PassPhenos_forQx.py --panel {panel} --outfile phenotypes_qx_{panel}.txt

rule polyAdapt_qx:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}-{pval}.tsv",
    output:
        qx="UKBiobank/selection_UKBV2/{panel}/QX_report-{level}-{pheno}-{pval}.txt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}-{pval}.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """

rule polyAdapt_gbr:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        gbr="paneldir/gbr.tsv.gz"
    output:
        qxfm="UKBiobank/selection_UKBV2/{panel}/QX_fm_report-{level}-{pheno}.txt",
    shell:
        """
        Rscript scripts/CalcQX_GBR-matched_Alba.R -w {input.candi} -e {input.neut} -a {input.gbr} -n 1000 -m {output.qxfm} -j 1000
        """

