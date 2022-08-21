# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/2/24
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
GERP_file="../02_Pop367Diploid_Combined/Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05.bed_GERP2.SolMsaInfo"
GERP_AllChrs=read.table(GERP_file,header=T)
PhasedAccessionList=c("RH","C151","E8669")

AllPhasedDeleStat=NULL
hap_dele_stat=paste0("RHC151E8669_PhasedDele_Stat_DiffCutOff_eachChrs_Recombination.final.txt")
 allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
  for (chr in allchrs){
    Chrom=paste0("chr",chr)
    #input
    geno_prefix=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/02_RHC151E8669_gvcf.367vcf/01_VCF/C151_RH_E8669_2hap_chr",chr,"_merge.g")
   gbs_file=paste0(geno_prefix,"_RefNumber.gds")
#####
print(geno_prefix)
geno <- snpgdsOpen(gbs_file)
X_all <- snpgdsGetGeno(geno)
X_all[1:5,1:5]
##major : 0; heter 1,alt:2 ; missing NA
snp_infor=snpgdsSNPList(geno, sample.id=NULL)
sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
hap6_nameId=c("C151_H1","C151_H2","E8669_H1","E8669_H2","RH_H1","RH_H2")
index_6hap=match(hap6_nameId,sample.id)

    GERP_oneChr=GERP_AllChrs[GERP_AllChrs$Chrom==Chrom,]
    index=match(snp_infor$position,GERP_oneChr$Pos)
    geno_data=data.frame(t(X_all))
    names(geno_data)=sample.id
    #geno=cbind(snp_infor,GERP_oneChr[index,c(1,5,14,10,25,17)],geno_data,GERP_oneChr[index,c(14,15,18,23,24)])
    write.table(Dele_geno,paste0("01_DeleInfor/RHC151E8669_Chr",chr,"_GERP",ConservationCutoff,"_DeleGenotype.txt",sep=""),quote=F,row.name=F,sep='\t')
    AllCutOffs=c(2.75,3.5,2)
    for (ConservationCutoff in AllCutOffs){
      all.BurdenInfor=NULL
      GERP_FilterConserved=GERP_oneChr[GERP_oneChr$GERP>=ConservationCutoff & GERP_oneChr$MajorAllelePop==GERP_oneChr$Sol100_MajorAllel,]

      index_Conserved=match(snp_infor$position,GERP_FilterConserved$Pos)
      w_all=GERP_FilterConserved$GERP[index_Conserved]
      w_Conserved=w_all[!is.na(index_Conserved)]
      X_Conserved_Ref2=X_all[,!is.na(index_Conserved)]
      snp_infor_Constrainted=snp_infor[!is.na(index_Conserved),]
      GERP_FilterConserved_Sites=GERP_FilterConserved[match(snp_infor_Constrainted$position,GERP_FilterConserved$Pos),]
      print (sum(GERP_FilterConserved_Sites$Pos!=snp_infor_Constrainted$position))
      print (sum(GERP_FilterConserved_Sites$Pos==snp_infor_Constrainted$position))

      Dele_geno_Refdata=data.frame(t(X_Conserved_Ref2))
      Dele_geno_RefBased=cbind(snp_infor_Constrainted,GERP_FilterConserved_Sites[,c(1,5,14,10,25,17)],Dele_geno_Refdata,GERP_FilterConserved_Sites[,c(14,15,18,23,24)])
      write.table(Dele_geno_RefBased,paste0("01_DeleInfor/RHC151E8669_Chr",chr,"_GERP",ConservationCutoff,"_Dele_RefIs2Genotype.txt",sep=""),quote=F,row.name=F,sep='\t')

      GERPConserved_minor_Sol100Major=GERP_FilterConserved_Sites[GERP_FilterConserved_Sites$Ref==GERP_FilterConserved_Sites$Sol100_MajorAllel,]
      index_missAsign=match(GERP_FilterConserved_Sites$Pos,snp_infor_Constrainted$position)
      missAsignNum=length(index_missAsign)
      print(paste("missAsignNum",missAsignNum))
      X_Conserved=X_Conserved_Ref2
      X_Conserved[,index_missAsign]=(X_Conserved[,index_missAsign]-2)*(-1)
      sum(X_Conserved[,-index_missAsign]==X_Conserved_Ref2[,-index_missAsign],na.rm=T)
      Dele_geno_data=data.frame(t(X_Conserved))
      names(Dele_geno_data)=sample.id
      Dele_geno=cbind(snp_infor_Constrainted,GERP_FilterConserved_Sites[,c(1,5,14,10,25,17)],Dele_geno_data,GERP_FilterConserved_Sites[,c(14,15,18,23,24)])
      #write.table(Dele_geno,paste0("01_DeleInfor/RHC151E8669_Chr",chr,"_GERP",ConservationCutoff,"_DeleGenotype.txt",sep=""),quote=F,row.name=F,sep='\t')
      ##stat the dele infor
for (phasedAccession in PhasedAccessionList){
  print(phasedAccession)
  SelectedSamples_Index=match(c(paste(phasedAccession,"_H1",sep=""),paste(phasedAccession,"_H2",sep="")),names(Dele_geno))
  allInfo_Conserved=Dele_geno[,c(1:11,SelectedSamples_Index)]
  names(allInfo_Conserved)[c(12:13)]=c("Hap1","Hap2")
  allInfo=allInfo_Conserved[!is.na(allInfo_Conserved$Hap1) & !is.na(allInfo_Conserved$Hap2) & allInfo_Conserved$Hap1!=1 & allInfo_Conserved$Hap2!=1 ,]
  write.table(allInfo,paste0("01_DeleInfor/PhasedHap_",phasedAccession,"_chr",chr,"_GERP",ConservationCutoff,".txt"),quote=F,sep='\t',row.name=F)
  allInfo_Conserved=allInfo
    sample_id=phasedAccession
    Conserved_SnpNnum=nrow(allInfo)
    index_heter=((allInfo_Conserved$Hap1==0 & allInfo_Conserved$Hap2==2) |(allInfo_Conserved$Hap1==2 & allInfo_Conserved$Hap2==0))
    HeterDeleNum=sum(index_heter,na.rm=T)
   index_homo=allInfo_Conserved$Hap1==2 & allInfo_Conserved$Hap2==2
    HomoDeleNum=sum(index_homo,na.rm=T)
    Burden_index_Heter_Conserved=sum(allInfo_Conserved$GERP[index_heter],na.rm=T)
    Burden_index_Homo_Conserved=sum(allInfo_Conserved$GERP[index_homo],na.rm=T)
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

    c=allInfo_Conserved_heter[unique(sort(c(Hap1Dele_DownNeiHapDele,Hap1Dele_DownNeiHapDele+1,Hap1Dele_UpNeiHapDele,Hap1Dele_UpNeiHapDele-1,
                                Hap2Dele_DownNeiHapDele,Hap2Dele_DownNeiHapDele+1,Hap2Dele_UpNeiHapDele,Hap2Dele_UpNeiHapDele-1))),]

  print(paste (phasedAccession,Chrom,"GERP",ConservationCutoff,"Check"))
    if(nrow(a)>0 | nrow(b)>0){
      print ("phased dele same")
      print(all(a==b),na.rm=T)
    }

    Nearby_Dele_Phased_Num=nrow(a)
  Nearby_Dele_Phased_Num_all=nrow(c)
    write.table(a,paste0("01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_",phasedAccession,"_chr",chr,"_GERP",ConservationCutoff,".txt"),quote=F,sep='\t',row.name=F)
    write.table(c,paste0("01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_",phasedAccession,"_chr",chr,"_GERP",ConservationCutoff,"_upDown.txt"),quote=F,sep='\t',row.name=F)
    Duplicated_index=which(a$Hap1[-1]==a$Hap1[-nrow(a)] & a$Hap2[-1]==a$Hap2[-nrow(a)] )+1
    a_rmNeiSame=a[-Duplicated_index,]
   Nearby_Dele_Phased_Num_rmDump=nrow(a_rmNeiSame)
    write.table(a_rmNeiSame,paste0("01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_",phasedAccession,"_chr",chr,"_GERP",ConservationCutoff,"_rmNeiDump.txt"),quote=F,sep='\t',row.name=F)
    PhasedDeleStat=data.frame(phasedAccession,Chrom,ConservationCutoff,Conserved_SnpNnum,HeterDeleNum,HomoDeleNum,Burden_index_Heter_Conserved,Burden_index_Homo_Conserved,Burden_index_EachiLine_Conserved,RecombinationBasedHap1,RecombinationBasedHap2,Nearby_Dele_Phased_Num_all,Nearby_Dele_Phased_Num,Nearby_Dele_Phased_Num_rmDump,Hap1_Dele,Hap1Dele_NeiHapDele_Num,Hap2_Dele,Hap2Dele_NeiHapDele_Num)
    AllPhasedDeleStat=rbind(AllPhasedDeleStat,PhasedDeleStat)
  }
}
  }
write.table(AllPhasedDeleStat,hap_dele_stat,quote=F,sep='\t',row.name=F)

library(tidyr)
library(ggplot2)
library(dplyr)

wholeGenome=data.frame(AllPhasedDeleStat %>%
                         group_by(phasedAccession,ConservationCutoff) %>%
                         summarize(Conserved_SnpNnum=sum(Conserved_SnpNnum),
                                   HeterDeleNum=sum(HeterDeleNum),
                                   HomoDeleNum=sum(HomoDeleNum),
                                   Burden_index_Heter_Conserved=sum(Burden_index_Heter_Conserved),
                                   Burden_index_Homo_Conserved=sum(Burden_index_Homo_Conserved),
                                   Burden_index_EachiLine_Conserved=sum(Burden_index_EachiLine_Conserved),
                                   RecombinationBasedHap1 = sum(RecombinationBasedHap1),
                                   RecombinationBasedHap2  = sum(RecombinationBasedHap2 ),
                                   Nearby_Dele_Phased_Num_all=sum(Nearby_Dele_Phased_Num_all),
                                   Nearby_Dele_Phased_Num=sum(Nearby_Dele_Phased_Num),
                                   Nearby_Dele_Phased_Num_rmDump=sum(Nearby_Dele_Phased_Num_rmDump),
                                   Hap1_Dele = sum(Hap1_Dele),
                                   Hap1Dele_NeiHapDele_Num  = sum(Hap1Dele_NeiHapDele_Num ),
                                   Hap2_Dele = sum(Hap2_Dele),
                                   Hap2Dele_NeiHapDele_Num  = sum(Hap2Dele_NeiHapDele_Num)))

wholeGenome$Chrom="AllChrs"
write.table(wholeGenome,paste0(hap_dele_stat,"_WholeGenome.txt"),quote=F,sep='\t',row.name=F)
