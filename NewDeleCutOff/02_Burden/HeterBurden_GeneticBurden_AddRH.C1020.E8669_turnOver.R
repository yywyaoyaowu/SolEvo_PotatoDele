# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/24
library(tidyr)
library(ggplot2)
library(GGally)
library(ggpubr)
#setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper//Output/08_DeleteriousMutation/04_WholeGenomeBurden/Pre_PopMinor_Dele/01_367Panel_SNP.TopConstraint/")
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/02_Dele/Pop367AddNewAccession_Substr193")

allBurdenIndex_Constraint=NULL
AllCutOffs=c(2.75,2.63,3,3.5,2,2.5)
for (ConservationCutoff in AllCutOffs){

   BurdenIndex=read.table(paste0("AllAccession_DeleBurden_BothMinor_BurdenInfo_LandracePlusNewPopMR0.5Dele_WholeGenome_GERP",ConservationCutoff),header=T)
allBurdenIndex_Constraint=rbind(allBurdenIndex_Constraint,BurdenIndex)
}
allBurdenIndex_Constraint$Type="Landrace"
allBurdenIndex_Constraint$Type[(!grepl("PG",allBurdenIndex_Constraint$sample.id))]="Founder"
LandranceBurden=allBurdenIndex_Constraint
LandranceBurden$GeneticLoad=LandranceBurden$Burden_index_EachiLine_Conserved*0.5
LandranceBurden$HomoLoad=LandranceBurden$Burden_index_Homo_Conserved
LandranceBurden$HeterLoad=LandranceBurden$Burden_index_Heter_Conserved

LandranceBurden$RealizedLoad_h0.14=LandranceBurden$Burden_index_Homo_Conserved+0.14*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.19=LandranceBurden$Burden_index_Homo_Conserved+0.19*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.1=LandranceBurden$Burden_index_Homo_Conserved+0.1*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.2=LandranceBurden$Burden_index_Homo_Conserved+0.2*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.3=LandranceBurden$Burden_index_Homo_Conserved+0.3*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.4=LandranceBurden$Burden_index_Homo_Conserved+0.4*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.05=LandranceBurden$Burden_index_Homo_Conserved+0.05*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.15=LandranceBurden$Burden_index_Homo_Conserved+0.15*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.25=LandranceBurden$Burden_index_Homo_Conserved+0.25*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.35=LandranceBurden$Burden_index_Homo_Conserved+0.35*LandranceBurden$Burden_index_Heter_Conserved
LandranceBurden$RealizedLoad_h0.12=LandranceBurden$Burden_index_Homo_Conserved+0.12*LandranceBurden$Burden_index_Heter_Conserved

#LandranceBurden$Cutoff_Id=factor(LandranceBurden$Cutoff_Id,levels=unique(LandranceBurden$Cutoff_Id))
LandranceBurden_justL=subset(LandranceBurden,Type=="Landrace")
GERP_CutOff=as.character(unique(LandranceBurden$ConservationCutoff))
Founder=subset(LandranceBurden,Type!="Landrace")

LandranceBurden_GERP2=subset(LandranceBurden,ConservationCutoff==2)
out=boxplot.stats(LandranceBurden_GERP2$GeneticLoad)$out
out_ind=which(LandranceBurden_GERP2$GeneticLoad %in% c(out))
LandranceBurden_GERP2_rm_lan=LandranceBurden_GERP2[-out_ind,]

All.LandranceBurden_rm=NULL
for (CutOff in GERP_CutOff){
  LandranceBurden_one=LandranceBurden[(LandranceBurden$ConservationCutoff==CutOff & LandranceBurden$Type=="Landrace"),]
  #out=boxplot.stats(LandranceBurden_one$GeneticLoad)$out
  #out_ind=which(LandranceBurden_one$GeneticLoad %in% c(out))
  out_ind=c(order(LandranceBurden_one$GeneticLoad)[1:3],order(LandranceBurden_one$GeneticLoad,decreasing=T)[1:3])
  LandranceBurden_rm=LandranceBurden_one[-out_ind,]
  OutlierLength=length(out_ind)
  print(paste(CutOff,OutlierLength))
  All.LandranceBurden_rm=rbind(All.LandranceBurden_rm,LandranceBurden_rm)
}

#LandranceBurden_rm=rbind(LandranceBurden_rm,Founder)
LandranceBurden_GERP2_rm_landrace=subset(All.LandranceBurden_rm,ConservationCutoff==2.75)
founder_GERP2_infor=subset(Founder,ConservationCutoff==2.75)
founder_GERP2_infor=rbind(founder_GERP2_infor,founder_GERP2_infor,founder_GERP2_infor)
#founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="C151_NGS" | founder_GERP2_infor$sample.id=="E8669_NGS" ]="Sucess"
#founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="C10-20_NGS" | founder_GERP2_infor$sample.id=="RH_NGS" ]="Faild"
#founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="A626_NGS" | founder_GERP2_infor$sample.id=="E463_NGS" ]="Inbred"
founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="C151-NGS" | founder_GERP2_infor$sample.id=="E8669-NGS" ]="Sucess"
founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="C10-20-NGS" | founder_GERP2_infor$sample.id=="RH-NGS" ]="Faild"
founder_GERP2_infor$Type[founder_GERP2_infor$sample.id=="A626-NGS" | founder_GERP2_infor$sample.id=="E463-NGS" ]="Inbred"

LandranceBurden_GERP2_rm=rbind(LandranceBurden_GERP2_rm_landrace,founder_GERP2_infor)
Burden_ColIndex=match(c("GeneticLoad","HomoLoad","HeterLoad","RealizedLoad_h0.1"),names(LandranceBurden_GERP2_rm))
p=ggpairs(LandranceBurden_GERP2_rm,columns = Burden_ColIndex, columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"), mapping = aes(color = Type))
p.f=p+scale_colour_manual(values=c("red","green",'black',"blue"))
ggsave(p.f,filename ="Landrace_BurdenCorrelation_Sol100DerivedDele_FilterDP_PlusNovelSNPs.pdf",height =6,width =6)


p=ggscatter(LandranceBurden_GERP2_rm, x = "RealizedLoad_h0.14", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=founder_GERP2_infor,aes(x = RealizedLoad_h0.14, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+ xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP2_rmOutlier_h0.14.pdf",height =4,width =4)


Burden_index=match(c("GeneticLoad","HomoLoad","HeterLoad","RealizedLoad_h0.14"),names(LandranceBurden_GERP2_rm))
p=ggpairs(LandranceBurden_GERP2_rm_lan,columns =Burden_index, columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"))
ggsave(p,filename ="Landrace_BurdenCorrelation_LandraceDerived_justLandrace.pdf",height =5,width =5)


p=ggscatter(LandranceBurden_GERP2_rm_lan, x = "RealizedLoad_h0.14", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_final=p+xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP2_rmOutlier_h0.14.pdf",height =4,width =4)


p=ggpairs(LandranceBurden_GERP2_rm,columns = Burden_index, columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"), mapping = aes(color = Type))
p.f=p+scale_colour_manual(values=c("red",'black'))
ggsave(p.f,filename ="Landrace_BurdenCorrelation_LandraceDerived.pdf",height =6,width =6)



#####3.12
# h=0.12
LandranceBurden_GERP3.12_rm=subset(LandranceBurden_rm,ConservationCutoff==3.12)
names(LandranceBurden_GERP3.12_rm)[c(15:17,28)]
p=ggpairs(LandranceBurden_GERP3.12_rm,columns = c(15:17,28), columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"), mapping = aes(color = Type))
p.f=p+scale_colour_manual(values=c("red",'black'))
ggsave(p.f,filename ="Landrace_BurdenCorrelation_GERP.3.12_h0.12.pdf",height =6,width =6)

LandranceBurden_GERP3.12_rm_lan=subset(LandranceBurden_GERP3.12_rm,Type=="Landrace")
p=ggpairs(LandranceBurden_GERP3.12_rm_lan,columns = c(15:17,18), columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"))
ggsave(p,filename ="Landrace_BurdenCorrelation_GERP.3.12.h0.12_justLandrace.pdf",height =5,width =5)

names(LandranceBurden_rm)[c(20:22)]
LandranceBurden_diff.h=gather(LandranceBurden_rm,h,RealizedLoad,names(LandranceBurden_rm)[c(20:22)])
LandranceBurden_diff.h_Cuoff=LandranceBurden_diff.h[(LandranceBurden_diff.h$ConservationCutoff==3.12| LandranceBurden_diff.h$ConservationCutoff==2),]



LandranceBurden_diff.h_Cuoff_landrace=subset(LandranceBurden_diff.h_Cuoff,Type=="Landrace")
LandranceBurden_diff.h_Cuoff_founder=subset(LandranceBurden_diff.h_Cuoff,Type=="Founder")

landrace_GERP2=subset(LandranceBurden_diff.h_Cuoff_landrace,ConservationCutoff==2)
founder_GERP2=subset(LandranceBurden_diff.h_Cuoff_founder,ConservationCutoff==2)

p=ggscatter(landrace_GERP2, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=founder_GERP2,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_grid(~h,scales="free")+ xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.final.pdf",height =3,width =9)




landrace_GERP3.12=subset(LandranceBurden_diff.h_Cuoff_landrace,ConservationCutoff==3.12)
founder_GERP3.12=subset(LandranceBurden_diff.h_Cuoff_founder,ConservationCutoff==3.12)

p=ggscatter(landrace_GERP3.12, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=founder_GERP3.12,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_grid(~h,scales="free")+ xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP3.12_rmOutlier_diff.h0.0.25_AddC10.20.final.pdf",height =3,width =9)









LandranceBurden_rm_JustLandrace=LandranceBurden_rm[LandranceBurden_rm$Type=="Landrace",]
E8669.C1020.RH_infor=LandranceBurden_rm[LandranceBurden_rm$Type!="Landrace",]


GERP2=subset(LandranceBurden_diff.h,ConservationCutoff==2)
GERP2_rm=subset(LandranceBurden_diff.h,ConservationCutoff==2)
GERP2_rm_h0.3=GERP2_rm[(GERP2_rm$h!="RealizedLoad_h0.35" & GERP2_rm$h!="RealizedLoad_h0.4" & GERP2_rm$Type=="Landrace"), ]
E8669.C1020.RH_infor_0.0.3=GERP2_rm[(GERP2_rm$h!="RealizedLoad_h0.35" & GERP2_rm$h!="RealizedLoad_h0.4" & GERP2_rm$Type=="Selected"),]
GERP2_rm_h0.4=GERP2_rm[GERP2_rm$Type=="Landrace", ]
E8669.C1020.RH_infor_0.0.4=GERP2_rm[ GERP2_rm$Type=="Selected",]

p=ggscatter(GERP2_rm_h0.4, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "HeterLoad", y = "RealizedLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_RealizedLoad.VS.HeterLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "HeterLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.HeterLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "HomoLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.HomoLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.3, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticLoad.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.3, x = "HeterLoad", y = "RealizedLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_RealizedLoad.VS.HeterLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.3, x = "HomoLoad", y = "RealizedLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_RealizedLoad.VS.HomoLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.4, x = "RealizedLoad", y = "HomoLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_RealizedLoad.VS.HomoLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.4, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)



p=ggscatter(GERP2_rm_h0.3, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)


GERP2_rm_h0.3_2=subset(GERP2_rm_h0.3,Type!="Selected")
p=ggscatter(GERP2_rm_h0.3_2, x = "RealizedLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = RealizedLoad, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_GERP2_rmOutlier_diff.h0.0.3_AddC10.20_2.pdf",height =6,width =15)



p=ggscatter(LandranceBurden_rm, x = "RealizedLoad_h0.2", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = RealizedLoad_h0.2, y = GeneticLoad),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.2_DifferentCutOff_rmOutlier_AddC10.20.pdf",height =3,width =15)





p=ggscatter(LandranceBurden_rm, x = "RealizedLoad_h0.1", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.1_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden_rm, x = "RealizedLoad_h0.3", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.3_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden_rm, x = "RealizedLoad_h0.4", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.4_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "HeterLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.HeterBurden_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "HomoLoad", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.HomoLoad_DifferentCutOff.pdf",height =3,width =15)

p=ggscatter(LandranceBurden, x = "RealizedLoad_h0.2", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.2_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "RealizedLoad_h0.2", y = "FutureGeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_FutureGeneticLoad.VS.RealizedLoad_h0.2_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "RealizedLoad_h0.1", y = "GeneticLoad",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.RealizedLoad_h0.1_DifferentCutOff.pdf",height =3,width =15)


allBurden_GERP2Landrace=BurdenPerAccession[(BurdenPerAccession$GroupDom=="Landrance" & BurdenPerAccession$ConservationCutoff==2) ,]
lm (allBurden_GERP2Landrace$Burden_index_Heter_Conserved ~ allBurden_GERP2Landrace$Burden_index_Homo_Conserved)

lm (allBurden_GERP2Landrace$Burden_index_EachiLine_Conserved ~ allBurden_GERP2Landrace$Burden_index_Homo_Conserved)



allBurden_GERP3.12=allBurdenIndex.final[(allBurdenIndex.final$Chrom=="AllChrs" & allBurdenIndex.final$ConservationCutoff==3.12) ,]
cutoffs=unique(allBurdenIndex.final$ConservationCutoff)[-c(1:2)]

AllCutOffResult=NULL
all.MoreHeterAccessionNumStat=NULL
all.result=NULL
for (Cutoff in cutoffs){
  GERP_DEle=subset(BurdenPerAccession,ConservationCutoff==Cutoff)
  GERP_DEle.order=GERP_DEle[match(allBurden_GERP2$sample.id,GERP_DEle$sample.id),]
  print(sum(GERP_DEle.order$sample.id==allBurden_GERP2$sample.id))
  one.result=NULL
