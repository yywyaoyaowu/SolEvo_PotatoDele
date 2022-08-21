# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
library(dplyr)

all.BurdenInfor_file="AllAccession_DeleBurden_BothMinor_BurdenInfo_LandracePlusNewPopMR0.5Dele"
GERP_file="Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05.bed_GERP2.SolMsaInfo"
GERP_AllChrs=read.table(GERP_file,header=T)

AllCutOffs=c(2.75,2.63,3,3.5,2,2.5)
 for (ConservationCutoff in AllCutOffs){
   all.BurdenInfor=NULL

   GERP_AllChrs_Dele=GERP_AllChrs[GERP_AllChrs$GERP>=ConservationCutoff & GERP_AllChrs$MajorAllelePop==GERP_AllChrs$Sol100_MajorAllel,]

   allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
for (chr in allchrs){
  geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/02_Pop367Diploid_Combined/Diploid367DP4Qaulity_PlusNewDP4_Combine_chr",chr,"_SubstrLandraceMR0.5")
  gbs_file=paste(geno_prefix,"_MinorIs2.gds",sep='')
  vcf.fn <- paste0(geno_prefix,".vcf")
  geno <- snpgdsOpen(gbs_file)
##major : 0; heter 1,alt:2 ; missing NA
  snp_infor=snpgdsSNPList(geno, sample.id=NULL)
  index_SNP=(!grepl ('0',snp_infor$allele))
  print (paste("Check Indel", sum(!index_SNP)))
  sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))

   X_all <- snpgdsGetGeno(geno)
   X_all[1:5,1:5]

   GERP_FilterConserved=GERP_AllChrs_Dele[GERP_AllChrs_Dele$Chrom==paste0("chr",chr),]
   index_Conserved=match(snp_infor$position,GERP_FilterConserved$Pos)
   w_all=GERP_FilterConserved$GERP[index_Conserved]
   w_Conserved=w_all[!is.na(index_Conserved)]

   X_Conserved_Minor2=X_all[,!is.na(index_Conserved)]
   snp_infor_Constrainted=snp_infor[!is.na(index_Conserved),]
   GERP_FilterConserved_Sites=GERP_FilterConserved[!is.na(match(GERP_FilterConserved$Pos,snp_infor_Constrainted$position)),]

SnpNum=nrow(snp_infor)
Conserved_SnpNnum=sum(!is.na(index_Conserved))
Geno_0_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==0,na.rm=T))
Geno_2_LineNum=apply(X_Conserved_Minor2,2,function(x) sum(x==2,na.rm=T))
print("Check minor is 2")
print(sum(Geno_0_LineNum<Geno_2_LineNum))
##gds from bed,bim,fam: major 0; minor:2 ; heter 1, missing NA
#heter1, home dele 2

###Heter Burden
X_Conserved_heter=X_Conserved_Minor2
X_Conserved_heter[X_Conserved_Minor2==2]=0
Burden_index_heterMatrix_Conserved=t(w_Conserved * t(X_Conserved_heter))
Burden_index_Heter_Conserved=apply(Burden_index_heterMatrix_Conserved,1,function(x) sum(x,na.rm=T))
HeterDeleNum=apply(X_Conserved_Minor2,1,function(x) sum(x==1,na.rm=T))

X_Conserved=X_Conserved_Minor2
HomeDeleNum=apply(X_Conserved,1,function(x) sum(x==2,na.rm=T))
Dele_matrix_weight_Conserved <- t(w_Conserved * t(X_Conserved))
Burden_index_EachiLine_Conserved=apply(Dele_matrix_weight_Conserved,1,function(x) sum(x,na.rm=T))
Burden_index_Homo_Conserved=(Burden_index_EachiLine_Conserved-Burden_index_Heter_Conserved)/2
Chrom=paste("chr",chr,sep="")
BurdenInfor=data.frame(ConservationCutoff,Chrom,SnpNum,Conserved_SnpNnum,sample.id,HeterDeleNum,HomeDeleNum,Burden_index_EachiLine_Conserved,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved)
    all.BurdenInfor=rbind(all.BurdenInfor,BurdenInfor)
    #write.table(BurdenInfor,paste("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/03_AddNewAccession/DMV6_chr",chr,"_ConservationCutoff",ConservationCutoff,"_Landrace182_mac3DP4QualityDeleSites_SolMsaNonMajor_TurnOverRareDele_BurdenInfo_AddNewAccession.txt",sep=""),quote=F,row.name=F,sep='\t')
  # Dele_geno_data=data.frame(t(X_Conserved))
   #names(Dele_geno_data)=sample.id
   #snp_infor_Constrainted$GERP=w_Conserved
   #Dele_geno=cbind(snp_infor_Constrainted,GERP_FilterConserved_Sites[,c(-c(2:3))],Dele_geno_data)
  # write.table(Dele_geno,paste("DMV6_chr",chr,"_ConstraintedConservationCutoff",ConservationCutoff,"_PlusNovelSites_SolMsaNonMajor_TurnOverRareDele_DeleGenotype_AddNewAccession.txt",sep=""),quote=F,row.name=F,sep='\t')
  }
write.table(all.BurdenInfor,paste0(all.BurdenInfor_file,"GERP",ConservationCutoff),quote=F,row.name=F,sep='\t')
#all.BurdenInfor=read.table("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/03_AddNewAccession/AddNewAccession_ConstraintedDelePlusNovelSNPs_SolMsaNonMajor_TurnOverRareDele_DiffConservationCutoff_BurdenInfo.txt",header=T)

wholeGenome=data.frame(all.BurdenInfor %>%
                         group_by(sample.id,ConservationCutoff) %>%
                         summarize(SnpNum=sum(SnpNum),
                                   Conserved_SnpNnum = sum(Conserved_SnpNnum),
                                   HeterDeleNum  = sum(HeterDeleNum),
                                   HomeDeleNum=sum(HomeDeleNum),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   Burden_index_Heter_Conserved = sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved  = sum(Burden_index_Homo_Conserved),
                         ))
wholeGenome$Chrom="AllChrs"
   write.table(wholeGenome,paste0(all.BurdenInfor_file,"_WholeGenome_GERP",ConservationCutoff),quote=F,row.name=F,sep='\t')
 }