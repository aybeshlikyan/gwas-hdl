#!/bin/bash

## PLINK Linear Regression
plink --bfile 'snpdata' --linear --pheno 'phenotype.txt'

## Recoding files for R
plink --bfile 'project_3_files/snpdata' --recode A --recode-allele 'project_3_files/snpdata.bim'

## Get files from UK Biobank and uncompress
wget https://www.dropbox.com/s/65jisgxwbbdrkaw/30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

gzcat 30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz > 30760_irnt.gwas.imputed_v3.both_sexes.tsv
gzcat variants.tsv.bgz > variants.tsv

## R Linear Regression, Comparison with PLINK, and PLINK Manhattan Plot, and 
## Comparing minimum PLINK p-value from each chromosome with their Biobank SNP p-value in R
Rscript r_script.r