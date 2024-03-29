# Correlation between components of variation in population structure and components of variation in SNP-trait association
We used the Tanigawa et al., 2019 method - DeGAs - to look at components of genetic association in the UK Biobank and compare them with the variation observed in ancient population structure. The aim is to determine which trait-associated variants may have been key in human evolution thus, contributing to the genetic makeup of Europeans today.

## Principal component analyses - PCA
Firstly, we need to perform PCA on our set of ancient samples to look at the population structure for which we used pcadapt: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html.  

## Correlation between population structrure and SNP_trait associations capture in modern humans
We are mainly interested in the correlation between the first two population structure components (PC1 and PC2) and the SNP-trait association from de DeGAs analyses. The first two components separate western- vs east-Eurasia, and hunter-gatherers (HG) groups vs. farmers (FA) respectively. 

## Randomization scheme and boostrapping

We finally, randomised the direction of the correlation to test if the observed correlation was significantly different from those observed in a null model with no association between SNP-trait and population structure.

Top 4 DeGaS components associated with the HG-FA. These components of SNP-trait associated represented mainly lifestyle traits and body differences. 

![Component1](https://github.com/albarema/neo/blob/master/PS1.png)

