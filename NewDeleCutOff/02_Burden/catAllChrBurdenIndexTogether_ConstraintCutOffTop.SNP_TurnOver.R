# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/3/17
Cutoff=read.table("../SelectedCutOff.txt",header=T)
allBurdenIndex=NULL
for (n in 1:nrow(Cutoff)){
oneCutOff=Cutoff$GERP_Cutoff[n]
oneCutoff_id=Cutoff$Type[n]
allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
  for (ch in allchrs){
    Burden_file=paste0("DMV6_chr",ch,"_GATKFilterd.snp_Landrace182MR0.5_ConservationCutoff",oneCutOff,"_SolMsaNonMajor_Dele_BurdenInfo_TurnOverRareDele.txt")
    oneChrBurdenIndex=read.table(Burden_file,header=T)
    oneChrBurdenIndex$Cutoff_Id=oneCutoff_id
  print (nrow(oneChrBurdenIndex))
    allBurdenIndex=rbind(allBurdenIndex,oneChrBurdenIndex)
}
}

library(tidyr)
library(ggplot2)
library(dplyr)
wholeGenome=data.frame(allBurdenIndex %>%
                         group_by(sample.id,Cutoff_Id,ConservationCutoff) %>%
                         summarize(SnpNum=sum(SnpNum),
                                   Conserved_SnpNnum = sum(Conserved_SnpNnum),
                                   HeterDeleNum  = sum(HeterDeleNum ),
                                   HomeDeleNum=sum(HomeDeleNum),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   Burden_index_Heter_Conserved = sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved  = sum(Burden_index_Homo_Conserved),
                         ))
wholeGenome$Chrom="AllChrs"
index=match(names(allBurdenIndex),names(wholeGenome))
allBurdenIndex.final=rbind(allBurdenIndex,wholeGenome[,index])
write.table(allBurdenIndex.final,"DiploidLandrace182MR0.5_AllChrs_DiffConservationCutoff_BurdenInfo_SolMsaNonMajorDele_TurnOverRareDele.txt",quote=F,sep='\t',row.name=F)


####
##add new accession
Cutoff=read.table("../SelectedCutOff.txt",header=T)
allBurdenIndex=NULL
for (n in 1:nrow(Cutoff)){
  oneCutOff=Cutoff$GERP_Cutoff[n]
  oneCutoff_id=Cutoff$Type[n]
  allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
  for (ch in allchrs){
    Burden_file=paste0("DMV6_chr",ch,"_GATKFilterd.snp_Landrace182MR0.5_ConservationCutoff",oneCutOff,"_SolMsaNonMajor_Dele_BurdenInfo_AddAccession_TurnOverRareDele.txt")
    oneChrBurdenIndex=read.table(Burden_file,header=T)
    oneChrBurdenIndex$Cutoff_Id=oneCutoff_id
    print (nrow(oneChrBurdenIndex))
    allBurdenIndex=rbind(allBurdenIndex,oneChrBurdenIndex)
  }
}

wholeGenome=data.frame(allBurdenIndex %>%
                         group_by(sample.id,Cutoff_Id,ConservationCutoff) %>%
                         summarize(SnpNum_Combine=sum(SnpNum_Combine),
                                    SnpNum=sum(SnpNum),
                                   Conserved_SnpNnum = sum(Conserved_SnpNnum),
                                   HeterDeleNum  = sum(HeterDeleNum ),
                                   HomeDeleNum=sum(HomeDeleNum),
                                   HomeDeleNum_pre=sum(HomeDeleNum_pre),
                                   HeterDeleNum_pre=sum(HeterDeleNum_pre),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   Burden_index_Heter_Conserved = sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved  = sum(Burden_index_Homo_Conserved),
                         ))
wholeGenome$Chrom="AllChrs"
index=match(names(allBurdenIndex),names(wholeGenome))
allBurdenIndex.final=rbind(allBurdenIndex,wholeGenome[,index])
write.table(allBurdenIndex.final,"DiploidLandrace182MR0.5_AllChrs_DiffConservationCutoff_BurdenInfo_SolMsaNonMajorDele_TurnOverRareDele_AddNewAccession.txt",quote=F,sep='\t',row.name=F)

NewAccession=allBurdenIndex.final[!grepl("PG",allBurdenIndex.final$sample.id),]
write.table(NewAccession,"DiploidLandrace367MR0.5_Substr182AddNewAccession_DiffConservationCutoff_BurdenInfo_TurnOverRareDele.txt",quote=F,sep='\t',row.name=F)
