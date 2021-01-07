# DBSLMM
Deterministic Bayesian Sparse Linear Mixed Model

## Update log
### v0.2 
* validates the flods of heritability: 0.8, 1 and 1.2
* fits the LMM, when the chromosome without any large effect SNPs
* fits the external validation all by R code. 

## Tutorial for external test only using summary statistics
### make SNP correlation matrix
We should use reference panel to construct SNP correlation matrix. You can treat [1000 Genome Project](https://www.internationalgenome.org/data#download) as reference panel. We also need the [block information]()
