# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
library(dplyr)
#https://rdrr.io/bioc/SNPRelate/f/vignettes/SNPRelate.Rmd
#HeterRegion.all=read.table("RH1015_heter_stat_50k_0.5_merge.bed2")
##input
HeterRegion_file="RH_RH1015_heter_stat_50k_0.02_250k_merge.bed"
PhasedDeleStat_allGenome_file=paste("RHC151E8669_PhasedDele_Stat_DiffCutOff_eachChrs_Recombination.final.txt_WholeGenome.txt")
RH_Heter_PhaseDele_file="RH.HeterRegion_RHC151E8669_PhasedDele_Stat_DiffCutOff_AllChrs_Recombination_SNP.TopConstraint.txt"
Genome_Len=731287687
##output
heterDele_Enrichment_file="RH_F4_HeterRegion_PhasedDele_Stat_DiffCutOff_Enrichment_ComparedHomoRegion_ConstraintCutOffTop.SNP.txt"
HomeRegionStat_file="RH_F4_HomoRegion_PhasedDele_Stat_DiffCutOff_ConstraintCutOffTop.SNP.txt"

##
HeterRegion.all=read.table(HeterRegion_file)
names(HeterRegion.all)=c("DM_V6_chrom","DM_V6_start","DM_V6_end")
HeterRegion.all$Loci=paste(HeterRegion.all$DM_V6_chrom,HeterRegion.all$DM_V6_start,HeterRegion.all$DM_V6_end,sep="_")
HeterRegion=HeterRegion.all[HeterRegion.all$DM_V6_start>0,]
HeterRegion$Length=HeterRegion$DM_V6_end-HeterRegion$DM_V6_start+1

#library(tidyverse)
PhasedDeleStat_allGenome=read.table(PhasedDeleStat_allGenome_file,header=T)
PhasedDeleStat_all=subset(PhasedDeleStat_allGenome,Chrom=="AllChrs")
PhasedDeleStat_all$DeleNum=PhasedDeleStat_all$HeterDeleNum+PhasedDeleStat_all$HomoDeleNum
RH_Heter_PhaseDele=read.table(RH_Heter_PhaseDele_file,header=T)
RH_Heter_PhaseDele$DeleNum=RH_Heter_PhaseDele$HeterDeleNum	+RH_Heter_PhaseDele$HomoDeleNum
index=match(RH_Heter_PhaseDele$HeterLoci,HeterRegion$Loci)

RH_Heter_PhaseDele$DM_V6_chrom=HeterRegion$DM_V6_chrom[index]
RH_Heter_PhaseDele$DM_V6_start=HeterRegion$DM_V6_start[index]
RH_Heter_PhaseDele$DM_V6_end=HeterRegion$DM_V6_end[index]
RH_Heter_PhaseDele$RegionLength=HeterRegion$DM_V6_end[index]-HeterRegion$DM_V6_start[index]+1
RH_Heter_PhaseDele$Cutoff_phasedAccession=paste(RH_Heter_PhaseDele$ConservationCutoff,RH_Heter_PhaseDele$phasedAccession,sep='_')
RH_Heter_PhaseDele$RegionLength[RH_Heter_PhaseDele$HeterLoci=="AllHeter"]=sum(HeterRegion$Length)
RH_Heter_PhaseDele2=RH_Heter_PhaseDele

RH_Heter_PhaseDele.all=subset(RH_Heter_PhaseDele,HeterLoci=="AllHeter")
#RH_Heter_PhaseDele.all$RegionLength=sum(HeterRegion$Length)
index.stat=match(c("RegionLength","DeleNum","Nearby_Dele_Phased_Num_rmDump"),names(RH_Heter_PhaseDele.all))
HomoRegionInfo=RH_Heter_PhaseDele.all[,c(1:5,index.stat)]

PhasedDeleStat_all$phasedAccession==HomoRegionInfo$phasedAccession
HomoRegionInfo$RegionLength=Genome_Len-RH_Heter_PhaseDele.all$RegionLength
HomoRegionInfo$DeleNum=PhasedDeleStat_all$DeleNum-RH_Heter_PhaseDele.all$DeleNum

HomoRegionInfo$Nearby_Dele_Phased_Num_rmDump=PhasedDeleStat_all$Nearby_Dele_Phased_Num_rmDump-RH_Heter_PhaseDele.all$Nearby_Dele_Phased_Num_rmDump
HomoRegionInfo$Cutoff_phasedAccession=paste(HomoRegionInfo$ConservationCutoff,HomoRegionInfo$phasedAccession,sep='_')

for (n in 1:nrow(RH_Heter_PhaseDele2)){
   i=match(RH_Heter_PhaseDele2$Cutoff_phasedAccession[n],HomoRegionInfo$Cutoff_phasedAccession)
  print(paste(n,i))
  RH_Heter_PhaseDele2$Cutoff_phasedAccession[n]
  ##per length
  if(RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n]>0){
    RH_Heter_PhaseDele2$DeleNum_EnrichFold[n]=(RH_Heter_PhaseDele2$DeleNum[n]/RH_Heter_PhaseDele2$RegionLength[n])/(HomoRegionInfo$DeleNum[i]/HomoRegionInfo$RegionLength[i])
    RH_Heter_PhaseDele2$DeleNum_EnrichP.value[n]=fisher.test(matrix(c(HomoRegionInfo$DeleNum[i],RH_Heter_PhaseDele2$DeleNum[n],HomoRegionInfo$RegionLength[i],RH_Heter_PhaseDele2$RegionLength[n]),2))$p.value

    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFold[n]=(RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n]/RH_Heter_PhaseDele2$RegionLength[n])/(HomoRegionInfo$Nearby_Dele_Phased_Num_rmDump[i]/HomoRegionInfo$RegionLength[i])
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFoldP.value[n]=fisher.test(matrix(c(HomoRegionInfo$Nearby_Dele_Phased_Num_rmDump[i],RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n],HomoRegionInfo$RegionLength[i],RH_Heter_PhaseDele2$RegionLength[n]),2))$p.value
  }else{
    RH_Heter_PhaseDele2$DeleNum_EnrichFold[n]=NA
    RH_Heter_PhaseDele2$DeleNum_EnrichP.value[n]=NA
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFold[n]=NA
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFoldP.value[n]=NA
  }
  ###Per Polymophsm
  if(RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n]>0){
    RH_Heter_PhaseDele2$DeleNum_EnrichFold_PerPolymophsm[n]=(RH_Heter_PhaseDele2$DeleNum[n]/RH_Heter_PhaseDele2$Polymorphism_Num[n])/(HomoRegionInfo$DeleNum[i]/HomoRegionInfo$Polymorphism_Num[i])
    RH_Heter_PhaseDele2$DeleNum_EnrichP.value_PerPolymophsm[n]=fisher.test(matrix(c(HomoRegionInfo$DeleNum[i],RH_Heter_PhaseDele2$DeleNum[n],HomoRegionInfo$Polymorphism_Num[i],RH_Heter_PhaseDele2$Polymorphism_Num[n]),2))$p.value
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFold_PerPolymophsm[n]=(RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n]/RH_Heter_PhaseDele2$Polymorphism_Num[n])/(HomoRegionInfo$Nearby_Dele_Phased_Num_rmDump[i]/HomoRegionInfo$Polymorphism_Num[i])
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFoldP.value_PerPolymophsm[n]=fisher.test(matrix(c(HomoRegionInfo$Nearby_Dele_Phased_Num_rmDump[i],RH_Heter_PhaseDele2$Nearby_Dele_Phased_Num_rmDump[n],HomoRegionInfo$Polymorphism_Num[i],RH_Heter_PhaseDele2$Polymorphism_Num[n]),2))$p.value
  }else{
    RH_Heter_PhaseDele2$DeleNum_EnrichFold_PerPolymophsm[n]=NA
    RH_Heter_PhaseDele2$DeleNum_EnrichP.value_PerPolymophsm[n]=NA
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFold_PerPolymophsm[n]=NA
    RH_Heter_PhaseDele2$NearbyDele_Phased_Num_EnrichFoldP.value_PerPolymophsm[n]=NA
  }
}
write.table(RH_Heter_PhaseDele2,heterDele_Enrichment_file,quote=F,sep='\t',row.name=F)
write.table(HomoRegionInfo,HomeRegionStat_file,quote=F,sep='\t',row.name=F)
