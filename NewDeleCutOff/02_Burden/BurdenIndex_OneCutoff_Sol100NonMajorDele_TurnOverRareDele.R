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
geno_prefix=paste("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/01_Landrace_VCF_2/DMV6_chr",chr,"_GATKFilterd.snp_182LandraceMR0.5Mac3",sep="")
bed.fn <- paste(geno_prefix,".bed",sep='')
fam.fn <- paste(geno_prefix,".fam",sep='')
bim.fn <- paste(geno_prefix,".bim",sep='')
gbs_file=paste(geno_prefix,".gds",sep='')

GERP_file=("DMV6_AllChrs_GATKFilterd.snp_182LandraceMR0.5Mac3_Freq.bed_GERP.DerivedAlleleInfor.txt")
print(geno_prefix)
print(GERP_file)
#print(bed.fn)

# Convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gbs_file)
#vcf.fn <- "C:/your_folder/your_vcf_file.vcf"
#snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")
geno <- snpgdsOpen(gbs_file)
# Take out genotype data for the first 3 samples and the first 5 SNPs
#g <- read.gdsn(index.gdsn(geno, "genotype"), start=c(1,1), count=c(5,3))
#Or take out genotype data with sample and SNP IDs, and four possible values are returned 0, 1, 2 and NA (3 is replaced by NA):
#g <- snpgdsGetGeno(geno, sample.id=..., snp.id=...)
X_all <- snpgdsGetGeno(geno)
X_all[1:5,1:5]
##major : 0; heter 1,alt:2 ; missing NA
snp_infor=snpgdsSNPList(geno, sample.id=NULL)
snp_infor$genome_pos=paste(snp_infor$chromosome,snp_infor$position,sep="_")
#write.table(snp_infor,paste(geno_prefix,"SnpList.txt",sep=""),quote=F,row.name=F)
sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
#write.table(sample.id,paste(geno_prefix,"AccessionList.txt",sep=""),quote=F,row.name=F)

##get the burden weight
GERP_AllChrs=read.table(GERP_file,header=T)
GERP_OneChr=GERP_AllChrs[GERP_AllChrs$Chr==paste0("chr",chr),]
##########################################################################333333#####
# heter, home 2,deleter1, weight by GERP_Z,calculated by depth >=20, conserved GERP score:
###############3#########
GERP_FilterConserved=GERP_OneChr[(GERP_OneChr$GERP>=ConservationCutoff),]

index_Conserved=match(snp_infor$position,GERP_FilterConserved$Pos)
w_all=GERP_FilterConserved$GERP[index_Conserved]

X_Conserved_Minor2=X_all[,!is.na(index_Conserved)]
snp_infor_Constrainted=snp_infor[!is.na(index_Conserved),]
w_Conserved=w_all[!is.na(index_Conserved)]

SnpNum=nrow(snp_infor)
Conserved_SnpNnum=sum(!is.na(index_Conserved))
Geno_0_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum<Geno_2_LineNum))

#####Heter Burden
X_Conserved_heter=X_Conserved_Minor2
X_Conserved_heter[X_Conserved_Minor2==2]=0

Burden_index_heterMatrix_Conserved=t(w_Conserved * t(X_Conserved_heter))
Burden_index_Heter_Conserved=apply(Burden_index_heterMatrix_Conserved,1,function(x) sum(x,na.rm=T))
HeterDeleNum=apply(X_Conserved_heter,1,function(x) sum(x==1,na.rm=T))

GERPConserved_minor_Sol100Major=GERP_FilterConserved[GERP_FilterConserved$MinorAllele==GERP_FilterConserved$Sol100_MajorAllel & GERP_FilterConserved$MAF >0.1,]
index_missAsign=match(GERPConserved_minor_Sol100Major$Pos,snp_infor_Constrainted$position)
missAsignNum=length(index_missAsign)
print(paste("missAsignNum",missAsignNum))
X_Conserved=X_Conserved_Minor2
X_Conserved[,index_missAsign]=(X_Conserved[,index_missAsign]-2)*(-1)
sum(X_Conserved[,-index_missAsign]==X_Conserved_Minor2[,-index_missAsign],na.rm=T)

HomeDeleNum=apply(X_Conserved,1,function(x) sum(x==2,na.rm=T))
HomeDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==2,na.rm=T))
HeterDeleNum_pre=apply(X_Conserved_Minor2,1,function(x) sum(x==1,na.rm=T))

Dele_matrix_weight_Conserved <- t(w_Conserved * t(X_Conserved))
Burden_index_EachiLine_Conserved=apply(Dele_matrix_weight_Conserved,1,function(x) sum(x,na.rm=T))

print(X_Conserved_heter[1:3,c(which(X_Conserved[3,]==1)[1:5],which(X_Conserved[3,]==2)[1:5])])

Burden_index_Homo_Conserved=(Burden_index_EachiLine_Conserved-Burden_index_Heter_Conserved)/2
Chrom=paste("chr",chr,sep="")
BurdenInfor=data.frame(ConservationCutoff,Chrom,SnpNum,Conserved_SnpNnum,sample.id,HeterDeleNum,HomeDeleNum,Burden_index_EachiLine_Conserved,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved)
write.table(BurdenInfor,paste("DeleBurden/DMV6_chr",chr,"_GATKFilterd.snp_Landrace182MR0.5_ConservationCutoff",ConservationCutoff,"_SolMsaNonMajor_Dele_BurdenInfo_TurnOverRareDele.txt",sep=""),quote=F,row.name=F,sep='\t')

Dele_geno_data=data.frame(t(X_Conserved))
names(Dele_geno_data)=sample.id
snp_infor_Constrainted$GERP=w_Conserved
Dele_geno=cbind(snp_infor_Constrainted,Dele_geno_data)
write.table(Dele_geno,paste("DeleBurden/DMV6_chr",chr,"_GATKFilterd.snp_Landrace182MR0.5_ConservationCutoff",ConservationCutoff,"_SolMsaNonMajor_DeleGenotype_TurnOverRareDele.txt",sep=""),quote=F,row.name=F,sep='\t')

