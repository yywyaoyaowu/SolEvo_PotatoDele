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
ConservationCutoff=as.numeric(args[2])
#geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/05_AddNewAcession_MissRef/01_VCF/chr",chr,"_merge_RefMissing_MR0.7")
geno_prefix=paste0("DMV6_AddNewAccession_",chr,"_merge_RefMissing_MR0.7_182Landrace")
bed.fn <- paste(geno_prefix,".bed",sep='')
fam.fn <- paste(geno_prefix,".fam",sep='')
bim.fn <- paste(geno_prefix,".bim",sep='')
gbs_file=paste(geno_prefix,".gds",sep='')

GERP_file=("DMV6_AllChrs_GATKFilterd.snp_182LandraceMR0.5Mac3_Freq.bed_GERP.DerivedAlleleInfor.txt")
print(geno_prefix)
print(GERP_file)
# Convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gbs_file)
geno <- snpgdsOpen(gbs_file)
##major : 0; heter 1,alt:2 ; missing NA
snp_infor_combine=snpgdsSNPList(geno, sample.id=NULL)
index_SNP=(!grepl ('0',snp_infor_combine$allele))
print (paste("Check Indel", sum(!index_SNP)))
sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))

##########################################################################333333#####
# heter, home 2,deleter1, weight by GERP_Z,calculated by depth >=20, conserved GERP score:
###############3#########

##get the burden weight
GERP_AllChrs=read.table(GERP_file,header=T)
GERP_OneChr=GERP_AllChrs[GERP_AllChrs$Chr==paste0("chr",chr),]
index_Landrace=match(snp_infor_combine$position,GERP_OneChr$Pos)
snp_infor=snp_infor_combine[!is.na(index_Landrace),]
write.table(snp_infor,paste(geno_prefix,paste0("SnpList_182SNPs_filterPos_chr",chr,".txt")),quote=F,row.name=F)

X_all_combine <- snpgdsGetGeno(geno)
X_all_combine[1:5,1:5]
X_all=X_all_combine[,!is.na(index_Landrace)]

GERP_FilterConserved=GERP_OneChr[(GERP_OneChr$GERP>=ConservationCutoff),]
index_Conserved=match(snp_infor$position,GERP_FilterConserved$Pos)
w_all=GERP_FilterConserved$GERP[index_Conserved]

X_Conserved_Minor2=X_all[,!is.na(index_Conserved)]
snp_infor_Constrainted=snp_infor[!is.na(index_Conserved),]
w_Conserved=w_all[!is.na(index_Conserved)]

SnpNum_Combine=nrow(snp_infor_combine)
SnpNum=nrow(snp_infor)
Conserved_SnpNnum=sum(!is.na(index_Conserved))
Geno_0_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum<Geno_2_LineNum))
##gds from bed,bim,fam: major 0; minor:2 ; heter 1, missing NA
#heter1, home dele 2

GERPConserved_minor_Sol100Major=GERP_FilterConserved[GERP_FilterConserved$MinorAllele==GERP_FilterConserved$Sol100_MajorAllel,]
index_missAsign=match(GERPConserved_minor_Sol100Major$Pos,snp_infor_Constrainted$position)
missAsignNum=length(index_missAsign)
print(paste("missAsignNum",missAsignNum))
X_Conserved=X_Conserved_Minor2
X_Conserved[,index_missAsign]=(X_Conserved[,index_missAsign]-2)*(-1)
sum(X_Conserved[,-index_missAsign]==X_Conserved_Minor2[,-index_missAsign],na.rm=T)

HomeDeleNum=apply(X_Conserved,1,function(x) sum(x==2,na.rm=T))
HeterDeleNum=apply(X_Conserved,1,function(x) sum(x==1,na.rm=T))
HomeDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==2,na.rm=T))
HeterDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==1,na.rm=T))

Dele_matrix_weight_Conserved <- t(w_Conserved * t(X_Conserved))
Burden_index_EachiLine_Conserved=apply(Dele_matrix_weight_Conserved,1,function(x) sum(x,na.rm=T))
X_Conserved_heter=X_Conserved
X_Conserved_heter[X_Conserved==2]=0

Burden_index_heterMatrix_Conserved=t(w_Conserved * t(X_Conserved_heter))
Burden_index_Heter_Conserved=apply(Burden_index_heterMatrix_Conserved,1,function(x) sum(x,na.rm=T))
Burden_index_Homo_Conserved=(Burden_index_EachiLine_Conserved-Burden_index_Heter_Conserved)/2
Chrom=paste("chr",chr,sep="")
BurdenInfor=data.frame(ConservationCutoff,Chrom,SnpNum_Combine,SnpNum,Conserved_SnpNnum,sample.id,HeterDeleNum,HomeDeleNum,Burden_index_EachiLine_Conserved,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved)
write.table(BurdenInfor,paste("DMV6_chr",chr,"_GATKFilterd.snp_Landrace182MR0.5_ConservationCutoff",ConservationCutoff,"_SolMsaNonMajor_Dele_BurdenInfo_AddAccession.txt",sep=""),quote=F,row.name=F,sep='\t')
