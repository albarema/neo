CHR = range(1,23)

rule all:
    input:
        expand("temp2_vcf_chr{chrom}.acf", chrom=CHR),
        "paneldir/inds-euras.clusters2.acf.gz"
	
	
rule get_filtered_vcf:
    input:
        vcf="vcf/chr{chrom}.maf1pm.bcf",
        rmid="euras.rmid.panelfile.txt"
    output:
          "vcf/chr{chrom}.maf1pm-euras-qx.vcf"
    shell:
        """
        bcftools view -S {input.rmid} {input.vcf} > {output}
        """

rule get_panel_galact:
    input:
        vcf="vcf/chr{chrom}.maf1pm-euras-qx.vcf",
        epofile="/willerslev/datasets/EPO/all.epo.gz",
        faifile="/willerslev/datasets/reference_genomes/hs37d5/hs37d5.fa.fai",
        lbdfile="/willerslev/datasets/ldetect-data/EUR/fourier_ls-all_nochr.bed",
        panelfile="euras.inds.panelfile.txt" #TODO pops
    output:
        tmp1=temp("temp_vcf_chr{chrom}.acf"),
        tmp2="temp2_vcf_chr{chrom}.acf"
    shell:
         """
         /willerslev/software/glactools/glactools vcfm2acf --epo {input.epofile} --fai {input.faifile} <(bcftools view {input.vcf}) > {output.tmp1}
         /willerslev/software/glactools/glactools meld -f {input.panelfile} {output.tmp1} > {output.tmp2}
         """

rule merge_galact:
    input:
        tmp2=expand("temp2_vcf_chr{chrom}.acf", chrom=CHR,allow_missing=True)
    output:
        acftp=temp("temp2_vcf_allchr.acf"),
        popfile=protected("paneldir/inds-euras.clusters2.acf.gz") #TODO pops
    shell:
         """
        /willerslev/software/glactools/glactools cat {input.tmp2} > {output.acftp}
        cat <(/willerslev/software/glactools/glactools view -h {output.acftp} | head -1) <(/willerslev/software/glactools/glactools view {output.acftp}| sort -k1,1 -k2,2g) | bgzip -c > {output.popfile}
        tabix -s 1 -b 2 -e 2 {output.popfile}
        rm temp2_vcf_chr*.acf
        """
