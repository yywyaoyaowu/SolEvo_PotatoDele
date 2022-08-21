# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/8/3

setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/02_Dele/Pop367AddNewAccession_Substr193_DP4")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
final_col=c("#C1CAC0","#DD8466","#FAE5BA","#96C9F7","#C9ACEC")
cbbPalette2=c("#56B4E9", "#C9ACEC")
cbbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
library(ggplot2)
GERP_stat=read.table("DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt",header=T)
GERP_stat$ConstraintedProportion=(GERP_stat$Conserved_num/GERP_stat$Site_num)*100
GERP_stat$DeleProportion=(GERP_stat$Dele_num/GERP_stat$SNPNum)*100
GERP_stat$Region[GERP_stat$Region=="Promoter_TSS1K"]="Promoter"

###Fig3c ##MAF>=0.01
all.CutOffs=c(2,2.75,3.5)
for (CutOff in all.CutOffs){
  figure_name=paste0("Figure3c_DeleMAF0.01_Distribution_NonSyn_GERP",CutOff,".pdf")
  GERP_oneCutoff=GERP_stat[GERP_stat$CutOff==CutOff,]
  GERP_oneCutoff_2=GERP_oneCutoff
  CDS=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="CDS"]
  Intron=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="Intron"]
  UTR=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="UTR"]
  UpDown5K=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="UpDown5K"]
  InterGenetic=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="InterGenetic"]

  NonSynonymous=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="NonSynonymous"]
  Synonymous=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="Synonymous"]
  Num = c(NonSynonymous,Synonymous,Intron,UTR,UpDown5K,InterGenetic)
  Type = c("NonSynonymous","Synonymous","Intron","UTR","UpDown5K","InterGenetic")
  piepercent = paste0(round(100*Num/sum(Num),digit=1), "%")
  print (100*Num/sum(Num))
  print (piepercent)
  lab=paste(Type,piepercent,sep=": ")
  cols=c("#000000","#696969","#A9A9A9","#C0C0C0","#D3D3D3","white")
  pdf(file=figure_name,width=6,height=4)
  pie(Num, labels=piepercent,col=cols,main = paste("Dele Distribution (GERP>=",CutOff,")"))
  legend("topright", Type, cex=0.8, fill=cols)
  dev.off()

  figure_name=paste0("Figure3c_DeleMAF0.01_Distribution_CDS_GERP",CutOff,".pdf")
  Num = c(CDS,Intron,UTR,UpDown5K,InterGenetic)
  Type = c("CDS","Intron","UTR","UpDown5K","InterGenetic")
  piepercent = paste0(round(100*Num/sum(Num),digit=1), "%")
  lab=paste(Type,piepercent,sep=": ")
  cols=c("#000000","#A9A9A9","#C0C0C0","#D3D3D3","white")
  pdf(file=figure_name,width=6,height=4)
  pie(Num, labels=piepercent,col=cols,main = paste("Dele Distribution (GERP>=",CutOff,")"))
  legend("topright", Type, cex=0.8, fill=cols)
  dev.off()
}

  ##figure3d
all.CutOffs=c(2,2.75,3.5)
for (CutOff in all.CutOffs){
  GERP_oneCutoff_2=GERP_stat[GERP_stat$CutOff==CutOff,]
  figure_name=paste0("Figure3d_DeleMAF0.01_NonSyn.NonCoding_GERP",CutOff,".pdf")
  CDS=GERP_oneCutoff_2$DeleProportion[GERP_oneCutoff_2$Region=="CDS"]
  NonSynonymous=GERP_oneCutoff_2$DeleProportion[GERP_oneCutoff_2$Region=="NonSynonymous"]
  Synonymous=GERP_oneCutoff_2$DeleProportion[GERP_oneCutoff_2$Region=="Synonymous"]
  NonCodingDele=GERP_oneCutoff_2$Dele_num[match(c("Intron","UTR","UpDown5K","InterGenetic"),GERP_oneCutoff_2$Region)]
  NonCodingSNP=GERP_oneCutoff_2$SNPNum[match(c("Intron","UTR","UpDown5K","InterGenetic"),GERP_oneCutoff_2$Region)]
  NonCoding=(sum(NonCodingDele)/sum(NonCodingSNP))*100
  DeleProportion=c(NonSynonymous,Synonymous,NonCoding)
  Region=c("NonSynonymous","Synonymous","NonCoding")
  data=data.frame(Region,DeleProportion)
  max=max(data$DeleProportion)*1.05
  data$Region=factor(data$Region,levels=Region)
  p=ggplot(data, aes(x=Region, y=DeleProportion)) + geom_bar(stat='identity',position="dodge",fill="grey",width=.6) +scale_y_continuous(limits=c(0,max))+theme_bw()
  #ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
  p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 15))
  ##change xlab and ylab
  p.lab=p.axis + xlab("") + ylab("Proportion of Deleterious Mutations(%)")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
  p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
  p.lenged
  ggsave(p.lenged,filename =figure_name,height =5,width =2.5)
}

all.CutOffs=c(2,2.75,3.5)
for (CutOff in all.CutOffs){
  GERP_oneCutoff_2=GERP_stat[GERP_stat$CutOff==CutOff,]
  figure_name=paste0("Figure3d_DeleMAF0.01_CodingNonCoding_GERP",CutOff,".pdf")
  CDS=GERP_oneCutoff_2$DeleProportion[GERP_oneCutoff_2$Region=="CDS"]
  NonCodingDele=GERP_oneCutoff_2$Dele_num[match(c("Intron","UTR","UpDown5K","InterGenetic"),GERP_oneCutoff_2$Region)]
  NonCodingSNP=GERP_oneCutoff_2$SNPNum[match(c("Intron","UTR","UpDown5K","InterGenetic"),GERP_oneCutoff_2$Region)]
  NonCoding=(sum(NonCodingDele)/sum(NonCodingSNP))*100
  DeleProportion=c(CDS,NonCoding)
  Region=c("CDS","NonCoding")
  data=data.frame(Region,DeleProportion)
  data$Region=factor(data$Region,levels=Region)
  max=max(data$DeleProportion)*1.05
  p=ggplot(data, aes(x=Region, y=DeleProportion)) + geom_bar(stat='identity',position="dodge",fill="grey",width=.6) +scale_y_continuous(limits=c(0,max))+theme_bw()
  #ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
  p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 15))
  ##change xlab and ylab
  p.lab=p.axis + xlab("") + ylab("Proportion of Deleterious Mutations(%)")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
  p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
  p.lenged
  ggsave(p.lenged,filename =figure_name,height =5,width =2.5)
}

####Figure S6, one figure.
Distributions_Regions=c("NonSynonymous","Synonymous","Intron","UTR","UpDown5K","InterGenetic")
all.CutOffs2=c(2,2.75,3.5)
GERP_stat_sub=GERP_stat[!is.na(match(GERP_stat$Region,Distributions_Regions)),]
GERP_stat_sub2=GERP_stat_sub[!is.na(match(GERP_stat_sub$CutOff,all.CutOffs2)),]
GERP_stat_sub2$Region=factor(GERP_stat_sub2$Region,levels=Distributions_Regions)

GERP_stat_sub2$GERPCutOffLevel=NA
GERP_stat_sub2$GERPCutOffLevel[GERP_stat_sub2$CutOff==2]="Minor Threshold"
GERP_stat_sub2$GERPCutOffLevel[GERP_stat_sub2$CutOff==2.75]="Moderate Threshold"
GERP_stat_sub2$GERPCutOffLevel[GERP_stat_sub2$CutOff==3.5]="Strong Threshold"
GERP_stat_sub2$GERPCutOffLevel=factor(GERP_stat_sub2$GERPCutOffLevel,levels=unique(GERP_stat_sub2$GERPCutOffLevel))
cbbPalette2=c("#A9A9A9","#808080","#000000")
max=max(GERP_stat_sub2$DeleProportion)*1.05
p=ggplot(GERP_stat_sub2, aes(x=Region, y=DeleProportion,fill=GERPCutOffLevel)) + geom_bar(stat='identity',position="dodge",width=.6) + scale_fill_manual(values=cbbPalette2) +scale_y_continuous(expand=c(0,0.2), limits=c(0,max))+theme_bw()
p.axis=p + theme(axis.text.x = element_text(size = 10, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 10))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Proportion of Constrainted (%)")+ theme(axis.title.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 6))+theme(legend.title = element_text(size = 6))
p.lenged
ggsave(p.lenged,filename ="FigureS6_ProportionOfDeleteriousSNPs_SolMsa_DiffGERP.Cutoff.pdf",height =3,width =4)

#######Figure S6
Regions=c("NonSynonymous","Synonymous","Intron","UTR","Promoter","UpDown5K","InterGenetic")
all.CutOffs=c(2,2.5,2.57,2.75,3,3.5,3.9)
for (CutOff in all.CutOffs){
  #Figure S6 4-fold; 0-fold
  figure_name=paste0("ProportionOfDeleteriousSNPsMAF0.01_SolMsa_xd_GERP",CutOff,".pdf")
  deglist_0=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="deglist_0"]
  deglist_4=GERP_oneCutoff_2$Dele_num[GERP_oneCutoff_2$Region=="deglist_4"]
  deglist_other=CDS-deglist_0-deglist_4
  figure_name=paste0("Output/NewDeleCutOff/02_Dele/Figure3c_DeleMAF0.01_Distribution_xfold_GERP",CutOff,".pdf")
  Num = c(deglist_0,deglist_other,deglist_4,Intron,UTR,UpDown5K,InterGenetic)
  Type = c("0-fold","x-fold","4-fold","Intron","UTR","UpDown5K","InterGenetic")
  piepercent = paste0(round(100*Num/sum(Num)), "%")
  lab=paste(Type,piepercent,sep=": ")
  cols=c("#000000","#696969","#808080","#A9A9A9","#C0C0C0","#D3D3D3","white")
  pdf(file=figure_name,width=6,height=4)
  pie(Num, labels=piepercent,col=cols,main = paste("Dele Distribution (GERP>=",CutOff,")"))
  legend("topright", Type, cex=0.8, fill=cols)
  dev.off()

  figure_name=paste0("Output/NewDeleCutOff/02_Dele/Figure3c_DeleMAF0.01_Distribution_CDS_GERP",CutOff,".pdf")
  Num = c(CDS,Intron,UTR,UpDown5K,InterGenetic)
  cols=c("#000000","#696969","#A9A9A9","#D3D3D3","white")
  Type = c("CDS","Intron","UTR","UpDown5K","InterGenetic")
  piepercent = paste0(round(100*Num/sum(Num)), "%")
  lab=paste(Type,piepercent,sep=": ")
  pdf(file=figure_name,width=6,height=4)
  pie(Num, labels=piepercent,col=cols,main = paste("Dele Distribution (GERP>=",CutOff,")"))
  legend("topright", Type, cex=0.8, fill=cols)
  dev.off()
}
