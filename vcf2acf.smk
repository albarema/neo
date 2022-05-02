configfile: "config.yaml"

## global parameters
CHR = range(1,23)
PANNAME=['euras','europe']

## ALL
rule all:
    input:
        #"vcf/euras-qx.vcf"
        #expand("temp2_vcf_chr{chrom}.acf", chrom=CHR),
        "paneldir/euras-inds.acf.gz"
	
## rules
rule get_filtered_vcf:
    input:
        vcf=os.path.join(config['vcf_dir'], "neo.impute.1000g.vcf.gz"),
        rmid="{panel}.rmid.panelfile.txt",
    output:
        "vcf/{panel}-qx.vcf.gz"
    shell:
        """
        bcftools view -S {input.rmid} {input.vcf} | bgzip -c > {output}
        """

rule get_panel_galact:
    input:
        vcf="vcf/euras-qx.vcf.gz",
        epofile="/willerslev/datasets/EPO/all.epo.gz",
        faifile="/willerslev/datasets/reference_genomes/hs37d5/hs37d5.fa.fai",
        lbdfile="/willerslev/datasets/ldetect-data/EUR/fourier_ls-all_nochr.bed",
        panelfile="euras.inds.panelfile.txt" #TODO pops
    output:
        tmp1=temp("temp_vcf.acf"),
        tmp2="temp2_vcf_allchr.acf"
    shell:
         """
         /willerslev/software/glactools/glactools vcfm2acf --epo {input.epofile} --fai {input.faifile} <(bcftools view {input.vcf}) > {output.tmp1}
         /willerslev/software/glactools/glactools meld -f {input.panelfile} {output.tmp1} > {output.tmp2}
         """

rule merge_galact:
    input:
        acftp="temp2_vcf_allchr.acf"
    output:
        popfile=protected("paneldir/{panel}-inds.acf.gz") #TODO pops
    shell:
         """
        cat <(/willerslev/software/glactools/glactools view -h {input.acftp} | head -1) <(/willerslev/software/glactools/glactools view {input.acftp}| sort -k1,1 -k2,2g) | bgzip -c > {output.popfile}
        tabix -s 1 -b 2 -e 2 {output.popfile}
        #rm {input.acftp}
        """
