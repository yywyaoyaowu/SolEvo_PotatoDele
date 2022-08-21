# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/3/1
#https://emitanaka.org/posts/2022-02-20-color-considerations/
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
final_col=c("#C1CAC0","#DD8466","#FAE5BA","#96C9F7","#C9ACEC")
cbbPalette2=c("#56B4E9", "#C9ACEC")
library(ggplot2)
GERP_stat=read.table("Output/NewDeleCutOff/DM_V6_DeleStat_BurdenIndex_Landrace_DiffCutoff.txt",header=T)

#names(GERP_stat)=c("Region","AllSites","ConservedSites")
GERP_stat$Proportion=(GERP_stat$Conserved_num/GERP_stat$Site_num)*100
GERP_stat$Region[GERP_stat$Region=="Promoter_TSS1K"]="Promoter"
Regions=c("CDS","Intron","UTR","Promoter","UpDown5K","InterGenetic")
#Regions=c("CDS","deglist_0","deglist_4","Intron","UTR","Promoter","UpDown5K","InterGenetic")

cbbPalette2=c("#A9A9A9","#808080","#696969","#000000")
all.CutOffs=c(2,2.75,3.5)
for (CutOff in all.CutOffs){
figure_name=paste0("Output/NewDeleCutOff/01_ConstraintedSites/Figure2d_ProportionOfConstrainedSites_SolMsa_GERP",CutOff,".pdf")
GERP_oneCutoff=GERP_stat[GERP_stat$CutOff==CutOff,]
GERP_oneCutoff_2=GERP_oneCutoff[match(Regions,GERP_oneCutoff$Region),]
GERP_oneCutoff_2$Region=factor(GERP_oneCutoff_2$Region,levels=unique(GERP_oneCutoff_2$Region))

max=max(GERP_oneCutoff_2$Proportion)+max(GERP_oneCutoff_2$Proportion)*0.05
#()
p=ggplot(GERP_oneCutoff_2, aes(x=Region, y=Proportion)) + geom_bar(stat='identity',position="dodge",fill="grey",width=.4) +scale_y_continuous(expand=c(0,0.2), limits=c(0,max))+theme_bw()
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Proportion of constrained sites(%)")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename =figure_name,height =5,width =4)
}

####Figure S4, one figure.
figure_name=paste0("Output/NewDeleCutOff/01_ConstraintedSites/FigureS4_ProportionOfConstraintedSites_SolMsa_DiffGERP.Cutoff.pdf")
all.CutOffs2=c(2,2.75,3.5)
Regions=c("CDS","deglist_0","deglist_4","Intron","UTR","Promoter","UpDown5K","InterGenetic")
GERP_stat_sub=GERP_stat[!is.na(match(GERP_stat$Region,Regions)) & !is.na(match(GERP_stat$CutOff,all.CutOffs2)),]
GERP_stat_sub2=GERP_stat_sub[!is.na(match(GERP_stat_sub$CutOff,all.CutOffs2)),]
GERP_stat_sub2$Region=factor(GERP_stat_sub2$Region,levels=Regions)
GERP_stat_sub2$CutOff=factor(GERP_stat_sub2$CutOff,levels=all.CutOffs2)

cbbPalette2=c("#A9A9A9","#808080","#000000")
p=ggplot(GERP_stat_sub2, aes(x=Region, y=Proportion,fill=CutOff)) + geom_bar(stat='identity',position="dodge",width=.6) + scale_fill_manual(values=cbbPalette2)+theme_bw()
p.axis=p + theme(axis.text.x = element_text(size = 10, color = "black", vjust =1, hjust =1, angle = 40))+ theme(axis.text.y = element_text(size = 10))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Proportion of Constrained sites(%)")+ theme(axis.title.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 10))+theme(legend.title = element_text(size = 10))
p.lenged
ggsave(p.lenged,filename =figure_name,height =3,width =4)
