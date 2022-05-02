configfile: "config.yaml"

## global parameters
CHR = range(1,23)
PANNAME=['euras','europe']
LEVEL=['pops', 'inds']

## ALL
rule all:
    input:
        #expand("paneldir/{panel}-{level}.acf.gz", level=LEVEL, panel='euras.oldbatch.nomiss'),
        expand("UKBiobank/selection_UKBV2/{panel}/Genscores-{level}-50_irnt.txt", level=LEVEL, panel='euras.oldbatch.nomiss')

# euras-fst_filtered.vcf.gz: used for euras.oldbatch results
# euras-GP09_filtered_sites.vcf.gz: masked version of euras-fst_filtered.vcf.gz.Results in: euras.oldbatchGP09
## rules
rule get_filtered_vcf:
    input:
        vcf="vcf/euras-qx.vcf.gz",
        inds="/willerslev/users-shared/science-snm-willerslev-gsd818/neo/20200309-impute/euras.rmid.panelfile_V2.txt", # 11 samples now not clustered are removed
        rmid="../20200309-impute/rsid.oldbatch.txt"
    output:
        passed="vcf/{panel}.qx.vcf.gz"
    shell:
        """
        /willerslev/users-shared/science-snm-willerslev-gsd818/Applications/vcftools/bin/vcftools --gzvcf {input.vcf} --keep {input.inds} --positions {input.rmid} --recode --stdout | bgzip -c > {output.passed}
        tabix -s 1 -b 2 -e 2 {output.passed} 
        """

rule get_panel_galact:
    input:
        vcf="vcf/{panel}.qx.vcf.gz",
        epofile="/willerslev/datasets/EPO/all.epo.gz",
        faifile="/willerslev/datasets/reference_genomes/hs37d5/hs37d5.fa.fai",
        lbdfile="/willerslev/datasets/ldetect-data/EUR/fourier_ls-all_nochr.bed",
        panelfile="/willerslev/users-shared/science-snm-willerslev-gsd818/neo/20200309-impute/euras.{level}.panelfile_V2.txt" #TODO inds
    output:
        tmp1=temp("temp_vcf-{panel}-{level}.acf"),
        tmp2=temp("temp2_vcf_allchr-{panel}-{level}.acf")
    shell:
         """
         /willerslev/software/glactools/glactools vcfm2acf --epo {input.epofile} --fai {input.faifile} <(bcftools view {input.vcf}) > {output.tmp1}
         /willerslev/software/glactools/glactools meld -f {input.panelfile} {output.tmp1} > {output.tmp2}
         """

rule merge_galact:
    input:
        acftp="temp2_vcf_allchr-{panel}-{level}.acf"
    output:
        popfile="paneldir/{panel}-{level}.acf.gz"
    shell:
         """
        cat <(/willerslev/software/glactools/glactools view -h {input.acftp} | head -1) <(/willerslev/software/glactools/glactools view {input.acftp}| sort -k1,1 -k2,2g) | bgzip -c > {output.popfile}
        tabix -s 1 -b 2 -e 2 {output.popfile}
        #rm {input.acftp}
        """

checkpoint polyAdapt_freqs:
    input:
        infile=os.path.join(config['uk_dir'], "{pheno}.gwas.imputed_v3.both_sexes.tsv.gz"),
        popfile="paneldir/{panel}-{level}.acf.gz",
        lbd=config['lbd']
    output:
        #freqs="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv",
        outfile="UKBiobank/data/{panel}/gwasfreqs-{level}-{pheno}.tsv.gz",
        candi="UKBiobank/data/{panel}/gwasfreqs_candidates-{level}-{pheno}.tsv",
        neut="UKBiobank/data/{panel}/gwasfreqs_neutral-{level}-{pheno}.tsv"
    shell:
        """
        python scripts/acf2ukbfreq_byP_Alba.py  -a {input.popfile} -g {input.infile} -o UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv
        cat <(head -1 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv) <(tail -n+2 UKBiobank/data/{wildcards.panel}/gwasfreqs-{wildcards.level}-{wildcards.pheno}.tsv| sort -k1,1 -k2,2g) | bgzip -c > {output.outfile}
        tabix -s 1 -b 2 -e 2 {output.outfile}
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


