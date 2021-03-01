
## Loading phenotype data into R
pheno <- read.table(file="phenotype.txt", header=FALSE)
names(pheno)=c("FID","IID","Phenotype")

## Loading recoded PLINK file into R and performing linear regression
plink.traw <- file("plink.traw", "r")
line <- readLines(plink.traw, n=1)
pVals <- rep(NA, 828325)

for (i in (1:828325)){
  line <- readLines(plink.traw, n=1)
  geno.snps <- as.numeric(unlist(strsplit(line, "\t"))[2510-(2503:0)])
  fit <- lm(pheno$Phenotype~geno.snps)
  currPVal <- summary(fit)$coefficients[2,4]
  pVals[i] = currPVal
}

close(plink.traw)
rm(plink.traw)
rm(line)
rm(geno.snps)

## Loading PLINK p-values into R
library(data.table)
plink.lr <- fread("plink.assoc.linear")

## Plotting frequencies of differences between PLINK p-values and R p-values
# plot((1:828325), (pVals-plink.lr$P), 
#      main="Difference in P-Values from R and PLINK",
#      xlab="SNP", ylab="Difference")

hist(pVals-plink.lr$P,
     main="Difference in P-Values from R and PLINK",
     xlab="Difference")
summary(pVals-plink.lr$P)

## Loading function and generating Manhattan plot for PLINK p-values
source("qqman.r")
manhattan(plink.lr)

## Finding SNP with minimum p-value in each chromosome
chrs <- c()
snps <- c()
minPVals <- c()

for(i in (1:22)){
	chrRows <- plink.lr[plink.lr$CHR==i]
	minRows <-chrRows[chrRows$P==min(chrRows$P)]
	chrs <- c(chrs, i)
	snps <- c(snps, minRows$SNP)
	minPVals <-c(minPVals, minRows$P)
}

data.PVals <- data.frame(chrs, snps, minPVals)

## Loading Biobank files into R and extracted relevant rows
variants <- fread("variants.tsv", select=c("variant", "chr", "rsid"))
relevant.variants <- variants[variants$rsid %in% data.PVals$snps]
rm(variants)
biobank.pvals <- fread("30760_irnt.gwas.imputed_v3.both_sexes.tsv", select=c("variant", "pval"))
relevant.pvals <- biobank.pvals[biobank.pvals$variant %in% relevant.variants$variant]
rm(biobank.pvals)

## Adding corresponding Biobank p-values to original dataframe
data.PVals$biobankPVal <- NaN
for (i in data.PVals$snps) {
  data.PVals[data.PVals$snps == i,]$biobankPVal = relevant.pvals[relevant.pvals$variant == relevant.variants[relevant.variants$rsid == i]$variant]$pval
}

data.PVals