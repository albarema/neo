
All pipelines are written in snakemake. I'd recommend to have a look at this repository for a more detailed explanation: https://github.com/albarema/GWAS_choice. 


- POPULATION-LEVEL PRS
STEP 1: VCF2ACF.smk: Get allele counts from the vcf and populations of interest
STEP 2: polyadapt_qx.smk. Get Qx statistics and polygenic scores for those phenotypes  in which we have enough "associated-SNPs" that passed whatever significant threshold we want to use (i.e.: 5e-8)

- INDIVIDUAL-LEVEL PRS
pcabyPRS_weightedGP.smk: pipeline to get the individual polygenic scores. We use the genotype probabilities to compute a weigthed version of the scores. This way, we take into account the uncertainty of the true genotypes


