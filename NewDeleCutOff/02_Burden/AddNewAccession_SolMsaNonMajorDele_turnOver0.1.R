# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24

library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
#https://rdrr.io/bioc/SNPRelate/f/vignettes/SNPRelate.Rmd
#args<-commandArgs(TRUE)
#chr=args[1]
#ConservationCutoff=as.numeric(args[2])
all.BurdenInfor_file="PreGeno_NonSolMajor_turnOver0.1Dele"
ConservationCutoff=2.63
allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
all.BurdenInfor=NULL
for (chr in allchrs){
geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/DerivedAlleleDeleterious/03_Burden_367Panel_PlusC10.20.RH.C151.E8669/01_VCF/chr",chr,"_merge_RefMissing_MR0.7")
bed.fn <- paste(geno_prefix,".bed",sep='')
fam.fn <- paste(geno_prefix,".fam",sep='')
bim.fn <- paste(geno_prefix,".bim",sep='')
gbs_file=paste(geno_prefix,".gds",sep='')
vcf_freq_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/DerivedAlleleDeleterious/03_Burden_367Panel_PlusC10.20.RH.C151.E8669/01_VCF/chr",chr,"_merge_RefMissing_FreqBiAllele")
GERP_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/DerivedAlleleDeleterious/01_DerivedAllele/ConstraintedSites_alleleInfor_chr",chr,".txt")
print(geno_prefix)
print(GERP_file)
print(bed.fn)
# Convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gbs_file)
#vcf.fn <- "C:/your_folder/your_vcf_file.vcf"
#snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")
geno <- snpgdsOpen(gbs_file)
# Take out genotype data for the first 3 samples and the first 5 SNPs
#g <- read.gdsn(index.gdsn(geno, "genotype"), start=c(1,1), count=c(5,3))
#Or take out genotype data with sample and SNP IDs, and four possible values are returned 0, 1, 2 and NA (3 is replaced by NA):
#g <- snpgdsGetGeno(geno, sample.id=..., snp.id=...)

##major : 0; heter 1,alt:2 ; missing NA
snp_infor_combine=snpgdsSNPList(geno, sample.id=NULL)
#snp_infor_combine$pos_allele=paste(snp_infor_combine$position,snp_infor_combine$allele,sep='_')
#write.table(snp_infor_combine,paste(geno_prefix,"SnpList.txt",sep=""),quote=F,row.name=F)
index_SNP=(!grepl ('0',snp_infor_combine$allele))
print (paste("Check Indel", sum(!index_SNP)))
#write.table(snp_infor_combine,paste(geno_prefix,"SnpList.txt",sep=""),quote=F,row.name=F)
sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
#write.table(sample.id,paste(geno_prefix,"AccessionList.txt",sep=""),quote=F,row.name=F)


#snp_infor=snp_infor_combine[(!is.na(index_367) & index_SNP) ,]
#write.table(snp_infor,paste(geno_prefix,"SnpList_367SNPs_filterPos.txt",sep=""),quote=F,row.name=F)

X_all_combine <- snpgdsGetGeno(geno)
X_all_combine[1:5,1:5]
#X_all=X_all_combine[,!is.na(index_367)]
##get the burden weoght
GERP_OneChr=read.table(GERP_file,header=T)
names(GERP_OneChr)[1:6]=c("Chrom","Pos_start","Pos","GERP","Neutral","Depth")
##########################################################################333333#####
# heter, home 2,deleter1, weight by GERP_Z,calculated by depth >=20, conserved GERP score:
###############3#########
GERP_FilterConserved=subset(GERP_OneChr,GERP>=ConservationCutoff)
index_Conserved=match(snp_infor_combine$position,GERP_FilterConserved$Pos)
w_all=GERP_FilterConserved$GERP[index_Conserved]
X_Conserved_Minor2=X_all_combine[,!is.na(index_Conserved)]
w_Conserved=w_all[!is.na(index_Conserved)]
snp_infor_Constrainted=snp_infor_combine[!is.na(index_Conserved),]

SnpNum_Combine=nrow(snp_infor_combine)
Conserved_SnpNnum=sum(!is.na(index_Conserved))

Geno_0_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum<Geno_2_LineNum))

vcf_freq=read.table(vcf_freq_file,skip=1)
names(vcf_freq)[1:8]=c("Chrom","position","N_ALLELES","N_ch","DM_Allele","RefFreq","Alt_Allele","AltFreq")

 index_RefMajor=which(vcf_freq$RefFreq>=vcf_freq$AltFreq)
 vcf_freq$MajorAllele=NA
 vcf_freq$MajorAllele[index_RefMajor]=vcf_freq$DM_Allele[index_RefMajor]
 index_Refminor=which(vcf_freq$RefFreq<vcf_freq$AltFreq)
 vcf_freq$MajorAllele[index_Refminor]=vcf_freq$Alt_Allele[index_Refminor]
 vcf_freq$MinorAllele=NA
vcf_freq$MinorAllele[which(vcf_freq$RefFreq>=vcf_freq$AltFreq)]=vcf_freq$Alt_Allele[which(vcf_freq$RefFreq>=vcf_freq$AltFreq)]
vcf_freq$MinorAllele[which(vcf_freq$RefFreq<vcf_freq$AltFreq)]=vcf_freq$DM_Allele[which(vcf_freq$RefFreq<vcf_freq$AltFreq)]

index_ConservedSNPFreq=match(snp_infor_Constrainted$position,vcf_freq$position)
sum(is.na(index_ConservedSNPFreq))

GERP_FilterConserved$MinorAllele_367PlusC10.20=vcf_freq$MinorAllele[match(GERP_FilterConserved$Pos,vcf_freq$position)]
GERP_FilterConserved$MajorAllele_367PlusC10.20=vcf_freq$MajorAllele[match(GERP_FilterConserved$Pos,vcf_freq$position)]

GERPConserved_minor_Sol100Major=GERP_FilterConserved[GERP_FilterConserved$MinorAllele_367PlusC10.20==GERP_FilterConserved$Sol100_MajorAllel & GERP_FilterConserved$MinorAllele_367PlusC10.20>0.1 ,]
index_missAsign=match(GERPConserved_minor_Sol100Major$Pos,snp_infor_Constrainted$position)
missAsignNum=length(index_missAsign)
print(paste("missAsignNum",missAsignNum))


 HomeDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==2,na.rm=T))
 HeterDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==1,na.rm=T))

X_Conserved=X_Conserved_Minor2
X_Conserved[,index_missAsign]=(X_Conserved[,index_missAsign]-2)*(-1)
sum(X_Conserved[,-index_missAsign]==X_Conserved_Minor2[,-index_missAsign],na.rm=T)

HomeDeleNum=apply(X_Conserved,1,function(x) sum(x==2,na.rm=T))
HeterDeleNum=apply(X_Conserved,1,function(x) sum(x==1,na.rm=T))
##gds from bed,bim,fam: major 0; minor:2 ; heter 1, missing NA
#heter1, home dele 2
Dele_matrix_weight_Conserved <- t(w_Conserved * t(X_Conserved))
Burden_index_EachiLine_Conserved=apply(Dele_matrix_weight_Conserved,1,function(x) sum(x,na.rm=T))
X_Conserved_heter=X_Conserved
X_Conserved_heter[X_Conserved==2]=0

Burden_index_heterMatrix_Conserved=t(w_Conserved * t(X_Conserved_heter))
Burden_index_Heter_Conserved=apply(Burden_index_heterMatrix_Conserved,1,function(x) sum(x,na.rm=T))
Burden_index_Homo_Conserved=(Burden_index_EachiLine_Conserved-Burden_index_Heter_Conserved)/2
Chrom=paste("chr",chr,sep="")
BurdenInfor=data.frame(ConservationCutoff,Chrom,SnpNum_Combine,Conserved_SnpNnum,sample.id,HeterDeleNum,HomeDeleNum,Burden_index_EachiLine_Conserved,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved)
#write.table(BurdenInfor,paste("DMV6_chr",chr,"_GATKFilterd.snp_367MR0.5_ConservationCutoff",ConservationCutoff,"_SolMsaNonMajor_BurdenInfo_filter367SNP.txt",sep=""),quote=F,row.name=F,sep='\t')
 all.BurdenInfor=rbind(all.BurdenInfor,BurdenInfor)
}


library(tidyr)
library(ggplot2)
library(dplyr)
wholeGenome=data.frame(all.BurdenInfor %>%
                         group_by(sample.id,ConservationCutoff) %>%
                         summarize(SnpNum=sum(SnpNum),
                                   Conserved_SnpNnum = sum(Conserved_SnpNnum),
                                   HeterDeleNum  = sum(HeterDeleNum),
                                   HomeDeleNum=sum(HomeDeleNum),
                                   HomeDeleNum_pre=sum(HomeDeleNum_pre),
                                   HeterDeleNum_pre=sum(HeterDeleNum_pre),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   Burden_index_Heter_Conserved = sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved  = sum(Burden_index_Homo_Conserved),
                         ))
wholeGenome$Chrom="AllChrs"
write.table(wholeGenome,paste0(all.BurdenInfor_file,"_WholeGenome.txt"),quote=F,sep='\t',row.name=F)
