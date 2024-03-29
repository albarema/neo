
# Identifying candidates for positive selection using patterns of ancient population differentiation

## Workflow to perform selection scans for a given K's

We use *pcadapt* software that performs Principal component analyses (PCA) on genomic data. More info on the software is here: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html.  

The aim is to detect positive selection by looking for loci in the genome that show strong differentiation in the allele frequencies between populations. Since we work with ancient data, those differences might also be representing allele frequency changes over time. 

When we performed PCA analyses it was very important to perform pre-filtering in our data. We applied the following filters:
- SNPs with MAF > 5%
- Genotype missingness rate < 50%
- variants where <10% of individuals had post-imputation genotype probability (GP) <= 0.8

Scree plots showing the percentage variance explained by each component can help us choose the number of K's we should use for our selection analyses. We used those components that explained at least 1% of the total variance in allele frequencies. 

Analyses performed:
- All Eurasians: which captures the west- vs. east- Eurasian differences.
- West-Eurasians: restricting our analyses to West-Eurasia. We removed those groups that were separated in the previous analyses by the first component (Siberian samples)
- West Eurasian Hunter-gatherer (Hg) vs. farmer scan: restricted to the loadings from the components separating HG and Neolithic farming peoples.

Candidate regions under selection were chosen by looking at those SNPs with the lowest p-values and concatenating those regions with significant variants. We labelled these regions with the HGNC protein-coding genes falling within the genomic coordinated (retrieved using the R package: biomaRt: https://bioconductor.org/packages/release/bioc/html/biomaRt.html).

Manhattan plot - All Eurasians. 

![manhattan](man.png)

