# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/8/14
library(tidyr)
library(ggplot2)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2=cbbPalette[c(3,2)]
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/02_Dele/03_RHC151E8669_phased/")
heter=read.table("RHC151E8669_PhasedDele_Stat_DiffCutOff_eachChrs_Recombination.final.txt_WholeGenome.txt",header=T)

heter$Hap1_HeterDele=heter$Hap1_Dele-heter$HomoDeleNum
heter$Hap2_HeterDele=heter$Hap2_Dele-heter$HomoDeleNum

heter$phasedAccession=factor(heter$phasedAccession,levels=c("E8669","C151","RH"))

heter_one=subset(heter,ConservationCutoff!=2)
heter_one$ConservationCutoff=factor(heter_one$ConservationCutoff,levels=c(2.75,3.5))

p=ggplot(heter_one, aes(x=phasedAccession, y=Nearby_Dele_Phased_Num,fill=ConservationCutoff)) + geom_bar(stat='identity',position="dodge",width=.4)+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()+ scale_fill_manual(values=cbbPalette2)
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, vjust = 0.5, hjust =0.5,color = "black", angle = 45))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Number of heterozygous deleterious")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15),legend.position=c(.2,.9))
p.lenged
ggsave(p.lenged,filename =paste0("replulsionPhasedHeterDele_DiffGERP.pdf"),height =6.5,width = 7)


all.CutOffs=c(2,2.75,3.5)
for (CutOff in all.CutOffs){
heter_one=subset(heter,ConservationCutoff==CutOff)
p=ggplot(heter_one, aes(x=phasedAccession, y=Nearby_Dele_Phased_Num)) + geom_bar(stat='identity',position="dodge",width=.4)+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", angle = 45))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Number of heterozygous deleterious")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename =paste0("replulsionPhasedHeterDele_GERP",CutOff,".pdf"),height =4,width = 3)


  index=match(c("Hap1_HeterDele","Hap2_HeterDele"),names(heter_one))
  heter_data=gather(heter_one,HeterType,Number,names(heter_one)[index])
  p=ggplot(heter_data, aes(x=phasedAccession, y=Number,fill=HeterType)) + geom_bar(stat='identity',position="dodge",width=.4)+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
  p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", angle = 0))+ theme(axis.text.y = element_text(size = 15))
  p.lab=p.axis + xlab("") + ylab("Number of deleterious mutation")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
  p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
  p.lenged
  ggsave(p.lenged,filename =paste0("replulsionPhasedHeterDele_GERP",CutOff,"_twoHaps.pdf"),height =4,width = 3)
}


index=match(c("Hap1_HeterDele","Hap2_HeterDele"),names(heter))
heter_data=gather(heter,HeterType,Number,names(heter)[index])
p=ggplot(heter, aes(x=phasedAccession, y=Nearby_Dele_Phased_Num)) + geom_bar(stat='identity',position="dodge",width=.4)+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", angle = 0))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Number of heterozygous deleterious")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename ="replulsionPhasedHeterDele.pdf",height =4,width = 3)



index=match(c("HeterDeleNum","Hap1_HeterDele","Hap2_HeterDele","Nearby_Dele_Phased_Num"),names(heter))
heter_data=gather(heter,HeterType,Number,names(heter)[index])

p=ggplot(heter_data, aes(x=HeterType, y=Number,fill=phasedAccession)) + geom_bar(stat='identity',position="dodge",width=.5)+scale_y_continuous(expand = expansion(mult = c(0, .05)))+ theme_bw()
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust =0.5, hjust =1, angle = 0))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Number of heterozygous deleterious")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename =".pdf",height =4,width = 4)
