configfile: "config.yaml"
import pandas as pd

LEVELS=['pops','inds']
PHENOS=pd.read_table('phenotypes.filtered.both_sexes.txt')['phenotype'].tolist()

PANNAME=['europeGP08pass','europeGP09pass']

## ALL
rule all_pops:
    input:
#        "vcf/europeGP09pass-qx.vcf.gz",
#        expand("paneldir/europeGP09pass-{level}_Filtered.acf.gz", level='pops'),
        expand("UKBiobank/data/{panel}/gwasfreqs_candidates-pops-{pheno}.tsv", pheno='50_irnt', panel=PANNAME),
        expand("UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-50_irnt.txt", level='pops',panel=PANNAME)

## rules
# rule get_pass_sites:
#     input:
#         vcf="vcf/europe-GP08_filtered_sites.vcf.gz", # euras-qx-infopass.vcf.gz, can be used since rmid and inds are being filtered
#         rmid="vcf/{panel}-sites.txt", #inlcude --positions {input.rmid}
#         inds="europe.rmid.panelfile.txt"
#     output:
#         passed="vcf/{panel}-qx.vcf.gz"
#     shell:
#         """
#         /willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcftools --gzvcf {input.vcf} --positions {input.rmid} --keep {input.inds} --recode --stdout | bgzip -c > {output.passed}
#         tabix -s 1 -b 2 -e 2 {output.passed}
#         """

rule get_panel_galact_v2:
    input:
        vcf="vcf/{panel}-qx.vcf.gz",
        epofile="/willerslev/datasets/EPO/all.epo.gz",
        faifile="/willerslev/datasets/reference_genomes/hs37d5/hs37d5.fa.fai",
        lbdfile="/willerslev/datasets/ldetect-data/EUR/fourier_ls-all_nochr.bed",
        panelfile="europe.{level}.panelfile.txt" #TODO pops
    output:
        tmp1=temp("tmp_vcf.{panel}.{level}.acf"),
        tmp2=temp("tmp2v2_vcf_allchr.{panel}.{level}.acf"),
        popfile="paneldir/{panel}-{level}_Filtered.acf.gz" #TODO pops
    shell:
        """
        /willerslev/software/glactools/glactools vcfm2acf --epo {input.epofile} --fai {input.faifile} <(bcftools view {input.vcf}) > {output.tmp1}
        /willerslev/software/glactools/glactools meld -f {input.panelfile} {output.tmp1} > {output.tmp2}
        cat <(/willerslev/software/glactools/glactools view -h {output.tmp2} | head -1) <(/willerslev/software/glactools/glactools view {output.tmp2}| sort -k1,1 -k2,2g) | bgzip -c > {output.popfile}
        tabix -s 1 -b 2 -e 2 {output.popfile}
        """

checkpoint polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.flipped.byP.gz"),
        popfile="paneldir/{panel}-{level}_Filtered.acf.gz",
        lbd=config['lbd']
    output:
        #freqs="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP.py -a {input.popfile} -g {input.infile} -o UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        cat <(head -1 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv) <(tail -n+2 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv| sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
        rm UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        python scripts/partitionUKB_byP.py -i {output.outfile} -b {input.lbd} -o {output.candi} -p 5e-08
        python scripts/extractneutral_byP.py -i {output.outfile} -o {output.neut} -s 20 -p 0.00001
        """

rule polyAdapt_qx:
    input:
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
    output:
        qx="UKBiobank/selection_UKBV2/{panel}/QX_report-{level}-{pheno}.txt",
        scores="UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-{pheno}.txt"
    shell:
        """
        Rscript scripts/CalcQX_edit4parallel_Alba.R -w {input.candi} -e {input.neut} -o {output.qx} -s {output.scores} -n 1000 -j 1000
        """

