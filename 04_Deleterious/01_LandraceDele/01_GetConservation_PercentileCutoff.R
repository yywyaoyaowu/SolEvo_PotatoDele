# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/3/17
allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
ConstraintedInfor=NULL
for (chr in allchrs){
  Constrainted_Chr=read.table(paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/02_Pop367Diploid_Combined/04_GenomeBurden/Sol_msa_ConstraintedGERP2_withDepth_withMajorAllele_chr",chr,".txt"),header=T)
  ConstraintedInfor=rbind(ConstraintedInfor,Constrainted_Chr)
}
names(ConstraintedInfor)[1:6]=c("Chrom","Pos_start","Pos","GERP","Neutral","Depth")
ConstraintedInfor$Position=paste(ConstraintedInfor[,1],ConstraintedInfor[,3],sep='_')

SNPs.infor_Freq=read.table("Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs.bed")
SNPs.infor_Freq$Position=paste(SNPs.infor_Freq[,1],SNPs.infor_Freq[,3],sep='_')
names(SNPs.infor_Freq)[1:8]=c("Chrom","Pos_start","Pos","AlleleNums","Ref","RefFreq","Alt","AltFreq")
print (table(SNPs.infor_Freq$Chrom))
SNPs.infor_Freq$MAF=apply(SNPs.infor_Freq[,c(6,8)],1,min)
SNPs.infor_Freq$RefIsMajorPop=(SNPs.infor_Freq$RefFreq >= SNPs.infor_Freq$AltFreq)
SNPs.infor_Freq$MajorAllelePop=NA
SNPs.infor_Freq$MajorAllelePop[SNPs.infor_Freq$RefIsMajorPop]=SNPs.infor_Freq$Ref[SNPs.infor_Freq$RefIsMajorPop]
SNPs.infor_Freq$MajorAllelePop[!SNPs.infor_Freq$RefIsMajorPop]=SNPs.infor_Freq$Alt[!SNPs.infor_Freq$RefIsMajorPop]
SNPs.infor_Freq$MinorAllelePop[!SNPs.infor_Freq$RefIsMajorPop]=SNPs.infor_Freq$Ref[!SNPs.infor_Freq$RefIsMajorPop]
SNPs.infor_Freq$MinorAllelePop[SNPs.infor_Freq$RefIsMajorPop]=SNPs.infor_Freq$Alt[SNPs.infor_Freq$RefIsMajorPop]

SNPs.infor=SNPs.infor_Freq
index=match(SNPs.infor$Position,ConstraintedInfor$Position)
SNPs.infor2=cbind(SNPs.infor,ConstraintedInfor[index,-c(1:3,15)])
SNPs.infor_GERP2=subset(SNPs.infor2,GERP>=2)
print(table(SNPs.infor_GERP2$Chrom))

##derived allele
SNPs.infor_GERP2$RefIsDerived=(SNPs.infor_GERP2$Ref!=SNPs.infor_GERP2$Sol100_MajorAllel)
SNPs.infor_GERP2$AltIsDerived=(SNPs.infor_GERP2$Alt!=SNPs.infor_GERP2$Sol100_MajorAllel)
SNPs.infor_GERP2$MinorIsDerived=(SNPs.infor_GERP2$MinorAllele!=SNPs.infor_GERP2$Sol100_MajorAllel)
SNPs.infor_GERP2$DerivedAlleleFreqInLandrace=NA
SNPs.infor_GERP2$DerivedAlleleFreqInLandrace[SNPs.infor_GERP2$RefIsDerived]=SNPs.infor_GERP2$RefFreq[SNPs.infor_GERP2$RefIsDerived]
SNPs.infor_GERP2$DerivedAlleleFreqInLandrace[SNPs.infor_GERP2$AltIsDerived]=SNPs.infor_GERP2$AltFreq[SNPs.infor_GERP2$AltIsDerived]
write.table(SNPs.infor_GERP2,"Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05.bed_GERP2.SolMsaInfo",quote=F,sep='\t',row.name=F)


Percentiles=c(1,0.5,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001)
all.GERP_cutInfor=NULL
polymophsm_GERP_ScoreSort=sort(SNPs.infor2$GERP,decreasing=T)
AllSites=nrow(SNPs.infor2)
Percentiles=c(1,0.5,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001)
GERP_Cutoff=c(polymophsm_GERP_ScoreSort[ceiling(AllSites*Percentiles)])
GERP_cutInfor=data.frame(Percentiles,GERP_Cutoff)
GERP_cutInfor$Type="LandraceDp4QualityPolymorphism"
all.GERP_cutInfor=rbind(all.GERP_cutInfor,GERP_cutInfor)
all.GERP_cutInfor
write.table(all.GERP_cutInfor,"Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05DP4_GERP.CutOff.txt",quote=F,sep='\t',row.name=F)


####at least two accession
library(ggplot2)

AllCutOffs=c(2,2.75,2.63,3,3.5,2.5)
for (ConservationCutoff in AllCutOffs){
  chr_vcf_freq=SNPs.infor2
chr_vcf_freq$Type=NA
chr_vcf_freq$Type="Non-constrainted"
chr_vcf_freq$Type[chr_vcf_freq$GERP>=2]="constrainted"
table(chr_vcf_freq$Type)
file_name=paste0("Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05_Conserved.nonConstrainted_MAF.hist_GERP",ConservationCutoff,".pdf")
p=ggplot(chr_vcf_freq, aes(x=MAF,colour=Type,fill=Type)) + geom_histogram(aes(y= ..density..),position="identity", alpha=0.5)+ theme_bw()
#+ scale_fill_manual(values=cbbPalette[11])
p.axis=p + theme(axis.text.x = element_text(size = 15))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("Minor Allele Frequency")+ theme(axis.title.x = element_text(size = 15))+ theme(axis.title.y = element_text(size = 15))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15),legend.position=c(.75,.85))
ggsave(p.lenged,filename =file_name,height =7,width =7)
}


#AllCutOffs=c(2,2.75,2.63,3,3.5,2.5)
AllCutOffs=c(2,2.75,3.5)
all.chr_vcf_freq_one_Constrainted=NULL
for (ConservationCutoff in AllCutOffs){
  chr_vcf_freq_one_Constrainted=subset(SNPs.infor_GERP2,GERP>=ConservationCutoff)
  chr_vcf_freq_one_Constrainted$CutOff=ConservationCutoff
  all.chr_vcf_freq_one_Constrainted=rbind(all.chr_vcf_freq_one_Constrainted,chr_vcf_freq_one_Constrainted)

  file_name=paste0("Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05_Conserved_DerivedAlleleFrequency.hist_GERP",ConservationCutoff,".pdf")
  p=ggplot(chr_vcf_freq_one_Constrainted, aes(x=DerivedAlleleFreqInLandrace)) + geom_histogram(fill="white",colour="black",position="identity")+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
  #+ scale_fill_manual(values=cbbPalette[11])
  p.axis=p + theme(axis.text.x = element_text(size = 15))+ theme(axis.text.y = element_text(size = 15))
  ##change xlab and ylab
  p.lab=p.axis + xlab("Derived Allele Frequency")  + ylab("Frequency")+ theme(axis.title.x = element_text(size = 15))+ theme(axis.title.y = element_text(size = 15))
  p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
  ggsave(p.lenged,filename =file_name,height =5,width =7)
}

##Derived allele different cutoff
all.chr_vcf_freq_one_Constrainted$CutOffID=NA
all.chr_vcf_freq_one_Constrainted$CutOffID[all.chr_vcf_freq_one_Constrainted$CutOff==2]="Minor Threshold"
all.chr_vcf_freq_one_Constrainted$CutOffID[all.chr_vcf_freq_one_Constrainted$CutOff==2.75]="Moderate Threshold"
all.chr_vcf_freq_one_Constrainted$CutOffID[all.chr_vcf_freq_one_Constrainted$CutOff==3.5]="Strong Threshold"
all.chr_vcf_freq_one_Constrainted$CutOffID=factor(all.chr_vcf_freq_one_Constrainted$CutOffID,levels=c("Minor Threshold","Moderate Threshold","Strong Threshold"))
file_name="Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05_Conserved.nonConstrainted_MAF.hist_GERPDiffCutoff.pdf"
p=ggplot(all.chr_vcf_freq_one_Constrainted, aes(x=DerivedAlleleFreqInLandrace)) + geom_histogram(aes(y= ..density..),position="identity", alpha=0.5)+ theme_bw()
#+ scale_fill_manual(values=cbbPalette[11])
p.axis=p + theme(axis.text.x = element_text(size = 15))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("Minor Allele Frequency")+ theme(axis.title.x = element_text(size = 15))+ theme(axis.title.y = element_text(size = 15))
p.lenged=p.lab+facet_wrap(~CutOffID,scales="free")+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15),legend.position=c(.75,.85))
ggsave(p.lenged,filename =file_name,height =7,width =7)


  file_name=paste0("Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05_MinorAlleleFreq.pdf")
  p=ggplot(SNPs.infor2, aes(x=MAF)) + geom_histogram(fill="white",colour="black",position="identity")+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
  #+ scale_fill_manual(values=cbbPalette[11])
  p.axis=p + theme(axis.text.x = element_text(size = 15))+ theme(axis.text.y = element_text(size = 15))
  ##change xlab and ylab
  p.lab=p.axis + xlab("Minor Allele Frequency")  + ylab("Frequency")+ theme(axis.title.x = element_text(size = 15))+ theme(axis.title.y = element_text(size = 15))
  p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
  ggsave(p.lenged,filename =file_name,height =5,width =7)



all.infor=NULL
mafs=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
all.CutOffs=c(2,2.63,2.75,3.03,3.5,3.91)
chr_vcf_freq=SNPs.infor2
for (cutoff in all.CutOffs){
  for (i in (1:(length(mafs)-1))){
    min_edge=mafs[i]
    OneMAF=mafs[i+1]
    AllSites=sum(chr_vcf_freq$MAF<=OneMAF &chr_vcf_freq$MAF>min_edge)
    ConservedSites=sum(chr_vcf_freq$MAF<=OneMAF &chr_vcf_freq$MAF>min_edge & chr_vcf_freq$GERP>=cutoff,na.rm=T)
    infor=data.frame(min_edge,OneMAF,AllSites,cutoff,ConservedSites)
    all.infor=rbind(all.infor,infor)
  }
}
all.infor$DeleRatio=all.infor$ConservedSites/all.infor$AllSites

write.table(all.infor,"Pop367_MR5_PlusNewAccession_SubstrLandrace_AllChrs_MR05_ConstraintedMAF_Stat.txt",quote=F,sep="\t",row.name=F)


AllCutOffs=c(2,2.75,2.63,3,3.5,2.5)
all.infor$MAFRegion=paste(all.infor$min_edge,all.infor$OneMAF,sep='-')
all.infor$MAFRegion=factor(all.infor$MAFRegion,levels=unique(all.infor$MAFRegion))
all.infor$Proportion=all.infor$DeleRatio*100
for (ConservationCutoff in AllCutOffs){
  MAF_Conserved_Ratio_stat=subset(all.infor,cutoff==ConservationCutoff)
p=ggplot(MAF_Conserved_Ratio_stat, aes(x=MAFRegion, y=Proportion)) + geom_bar(stat='identity',position="dodge",width=.5,fill='grey')+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 10, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 15))
  ##change xlab and ylab
p.lab=p.axis + xlab("Minor allele frequency") + ylab("Proportion of constrained sites (%)")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename =paste0("ConservedProportion_DiffMAF_bar_GERP",ConservationCutoff,".pdf"),height =4,width = 4)
}
