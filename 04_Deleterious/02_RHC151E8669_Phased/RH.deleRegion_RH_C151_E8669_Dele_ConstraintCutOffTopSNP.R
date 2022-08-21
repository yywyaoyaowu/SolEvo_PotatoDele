# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
#https://rdrr.io/bioc/SNPRelate/f/vignettes/SNPRelate.Rmd
#input
HeterRegion.all=read.table("RH_RH1015_heter_stat_50k_0.02_250k_merge.bed")
names(HeterRegion.all)=c("DM_V6_chrom","DM_V6_start","DM_V6_end")
HeterRegion.all$Loci=paste(HeterRegion.all$DM_V6_chrom,HeterRegion.all$DM_V6_start,HeterRegion.all$DM_V6_end,sep="_")
HeterRegion=HeterRegion.all[HeterRegion.all$DM_V6_start>0,]
chrs=as.character(unique(HeterRegion$DM_V6_chrom))

AllPhasedDeleStat=NULL
PhasedAccessionList=c("RH","C151","E8669")
AllPhasedDeleStat=NULL
AllCutOffs=c(2.75,3.5,2)

for (Chrom in chrs){
  HeterRegion_one.chr=HeterRegion[HeterRegion$DM_V6_chrom==Chrom,]
  for (phasedAccession in PhasedAccessionList){
    print(phasedAccession)
    allInfo_allphased=read.table(paste0("01_DeleInfor/RHC151E8669_",Chrom,"_GERP2.75_DeleGenotype.txt"),header=T)
    index_phased=grep(phasedAccession,names(allInfo_allphased))
    allInfo=allInfo_allphased[,c(1:11,index_phased,ncol(allInfo_allphased))]
    Polymorphism_Num=nrow(allInfo)
    print( names(allInfo)[c(12:13)])
    names(allInfo)[c(12:13)]=c("Hap1","Hap2")
    sample_id=phasedAccession

    for (ConservationCutoff in AllCutOffs){
      for (n in 1:nrow(HeterRegion_one.chr)){
        Chrom=HeterRegion_one.chr$DM_V6_chrom[n]
        HeterLoci=HeterRegion_one.chr$Loci[n]
        DMV6_Start=HeterRegion_one.chr$DM_V6_start[n]
        DMV6_End=HeterRegion_one.chr$DM_V6_end[n]
        all.ReginDeleInfo=allInfo[(allInfo$position >=HeterRegion_one.chr$DM_V6_start[n] & allInfo$position <=HeterRegion_one.chr$DM_V6_end[n]),]
        Polymorphism_Num=sum(!is.na(all.ReginDeleInfo$Hap1) & !is.na(all.ReginDeleInfo$Hap2))
        allInfo_Conserved=all.ReginDeleInfo[(!is.na(all.ReginDeleInfo$Hap1) & !is.na(all.ReginDeleInfo$Hap2)  & all.ReginDeleInfo$GERP>=ConservationCutoff),]
        Conserved_SnpNnum=nrow(allInfo_Conserved)
        index_heter=((allInfo_Conserved$Hap1==0 & allInfo_Conserved$Hap2==2) |(allInfo_Conserved$Hap1==2 & allInfo_Conserved$Hap2==0))
        HeterDeleNum=sum(index_heter,na.rm=T)
        HomoDeleNum=sum(allInfo_Conserved$Hap1==2 & allInfo_Conserved$Hap2==2,na.rm=T)
        Burden_index_Heter_Conserved=sum(allInfo_Conserved$GERP[index_heter],na.rm=T)
        Burden_index_Homo_Conserved=sum(allInfo_Conserved$GERP[allInfo_Conserved$Hap1==2 & allInfo_Conserved$Hap2==2],na.rm=T)
        Burden_index_EachiLine_Conserved=2*Burden_index_Homo_Conserved+Burden_index_Heter_Conserved

        Hap1_Dele=sum(allInfo_Conserved$Hap1==2,na.rm=T)
        Hap2_Dele=sum(allInfo_Conserved$Hap2==2,na.rm=T)

        allInfo_Conserved_heter=allInfo_Conserved[index_heter,]
        Hap1Dele_DownNeiHapDele=which(allInfo_Conserved_heter$Hap1[-nrow(allInfo_Conserved_heter)]==2 & allInfo_Conserved_heter$Hap2[-1]==2)
        Hap1Dele_UpNeiHapDele=which(allInfo_Conserved_heter$Hap1[-1]==2 & allInfo_Conserved_heter$Hap2[-nrow(allInfo_Conserved_heter)]==2)+1
        Hap1Dele_NeiHapDele=unique(sort(c(Hap1Dele_DownNeiHapDele,Hap1Dele_UpNeiHapDele)))
        Hap1Dele_NeiHapDele_Num=length(Hap1Dele_NeiHapDele)
        a=allInfo_Conserved_heter[unique(sort(c(Hap1Dele_DownNeiHapDele,Hap1Dele_DownNeiHapDele+1,Hap1Dele_UpNeiHapDele,Hap1Dele_UpNeiHapDele-1))),]
        RecombinationBasedHap1=length(Hap1Dele_DownNeiHapDele)+length(Hap1Dele_UpNeiHapDele)

        Hap2Dele_DownNeiHapDele=which(allInfo_Conserved_heter$Hap2[-nrow(allInfo_Conserved_heter)]==2 & allInfo_Conserved_heter$Hap1[-1]==2 )
        Hap2Dele_UpNeiHapDele=which(allInfo_Conserved_heter$Hap2[-1]==2 & allInfo_Conserved_heter$Hap1[-nrow(allInfo_Conserved_heter)]==2)+1
      Hap2Dele_NeiHapDele=unique(sort(c(Hap2Dele_DownNeiHapDele,Hap2Dele_UpNeiHapDele)))
      Hap2Dele_NeiHapDele_Num=length(Hap2Dele_NeiHapDele)
      RecombinationBasedHap2=length(Hap2Dele_DownNeiHapDele)+length(Hap2Dele_UpNeiHapDele)
      b=allInfo_Conserved_heter[unique(sort(c(Hap2Dele_DownNeiHapDele,Hap2Dele_DownNeiHapDele+1,Hap2Dele_UpNeiHapDele,Hap2Dele_UpNeiHapDele-1))),]
      print(paste (phasedAccession,Chrom,"GERP",ConservationCutoff,"Check"))
      if(nrow(a)>0 | nrow(b)>0){
        print(all(a==b),na.rm=T)
      }
      Nearby_Dele_Phased_Num=nrow(a)
        Duplicated_index=which(a$Hap1[-1]==a$Hap1[-nrow(a)] & a$Hap2[-1]==a$Hap2[-nrow(a)] )+1
        a_rmNeiSame=a[-Duplicated_index,]
        Nearby_Dele_Phased_Num_rmDump=nrow(a_rmNeiSame)
      PhasedDeleStat=data.frame(phasedAccession,HeterLoci,Polymorphism_Num,ConservationCutoff,Conserved_SnpNnum,HeterDeleNum,HomoDeleNum,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved,Burden_index_EachiLine_Conserved,RecombinationBasedHap1,RecombinationBasedHap2,Nearby_Dele_Phased_Num,Nearby_Dele_Phased_Num_rmDump,Hap1_Dele,Hap1Dele_NeiHapDele_Num,Hap2_Dele,Hap2Dele_NeiHapDele_Num)
      AllPhasedDeleStat=rbind(AllPhasedDeleStat,PhasedDeleStat)
    }
  }
}
}
write.table(AllPhasedDeleStat,"RH.HeterRegion_RHC151E8669_PhasedDele_Stat_DiffCutOff_AllChrs_Recombination_SNP.TopConstraint.txt",quote=F,sep='\t',row.name=F)
#AllPhasedDeleStat=read.table("RH.deleRegion_RHC151E8669_PhasedDele_Stat_DiffCutOff_AllChrs_Recombination_SNP.TopConstraint.txt",header=T)
library(tidyr)
library(ggplot2)
library(dplyr)
allHeterDele=data.frame(AllPhasedDeleStat %>%
                         group_by(phasedAccession,ConservationCutoff) %>%
                         summarize(Polymorphism_Num=sum(Polymorphism_Num),
                                   Conserved_SnpNnum=sum(Conserved_SnpNnum),
                                   HeterDeleNum=sum(HeterDeleNum),
                                   HomoDeleNum=sum(HomoDeleNum),
                                   Burden_index_Heter_Conserved=sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved=sum(Burden_index_Homo_Conserved),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   RecombinationBasedHap1 = sum(RecombinationBasedHap1),
                                   RecombinationBasedHap2  = sum(RecombinationBasedHap2 ),
                                   Nearby_Dele_Phased_Num=sum(Nearby_Dele_Phased_Num),
                                   Nearby_Dele_Phased_Num_rmDump=sum(Nearby_Dele_Phased_Num_rmDump),
                                   Hap1_Dele = sum(Hap1_Dele),
                                   Hap1Dele_NeiHapDele_Num  = sum(Hap1Dele_NeiHapDele_Num ),
                                   Hap2_Dele = sum(Hap2_Dele),
                                   Hap2Dele_NeiHapDele_Num  = sum(Hap2Dele_NeiHapDele_Num)))
allHeterDele$HeterLoci="AllHeter"
index=match(names(AllPhasedDeleStat),names(allHeterDele))
AllPhasedDeleStat.final=rbind(AllPhasedDeleStat,allHeterDele[,index])
write.table(AllPhasedDeleStat.final,"RH.HeterRegion_RHC151E8669_PhasedDele_Stat_DiffCutOff_AllChrs_Recombination_SNP.TopConstraint.txt",quote=F,sep='\t',row.name=F)

