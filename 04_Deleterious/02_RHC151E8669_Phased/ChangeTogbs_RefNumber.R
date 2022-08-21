# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
#https://rdrr.io/bioc/SNPRelate/f/vignettes/SNPRelate.Rmd
args<-commandArgs(TRUE)
chr=args[1]

 geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/02_RHC151E8669_gvcf.367vcf/01_VCF/C151_RH_E8669_2hap_chr",chr,"_merge.g")
#geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/02_Pop367Diploid_Combined/Diploid367DP4Qaulity_PlusNewDP4_Combine_chr",chr,"_SubstrLandraceMR0.5")
vcf.fn=paste0(geno_prefix,".vcf")
gbs_file=paste(geno_prefix,"_RefNumber.gds",sep='')
snpgdsVCF2GDS(vcf.fn, gbs_file,method="copy.num.of.ref")

geno <- snpgdsOpen(gbs_file)
# Take out genotype data for the first 3 samples and the first 5 SNPs
#g <- read.gdsn(index.gdsn(geno, "genotype"), start=c(1,1), count=c(5,3))
#Or take out genotype data with sample and SNP IDs, and four possible values are returned 0, 1, 2 and NA (3 is replaced by NA):
#g <- snpgdsGetGeno(geno, sample.id=..., snp.id=...)
X_all <- snpgdsGetGeno(geno)
X_all[1:5,1:5]

Geno_0_LineNum=apply(X_all,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_all,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum<Geno_2_LineNum))

