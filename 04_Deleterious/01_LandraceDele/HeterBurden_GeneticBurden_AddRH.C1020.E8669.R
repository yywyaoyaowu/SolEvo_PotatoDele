# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/24
library(tidyr)
library(ggplot2)
library(GGally)
library(ggpubr)
#setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper//Output/08_DeleteriousMutation/04_WholeGenomeBurden/Pre_PopMinor_Dele/01_367Panel_SNP.TopConstraint/")
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/02_Dele/Pop367AddNewAccession_Substr193_DP4")

founder_sucess=c("C151-NGS","E8669-NGS")
founder_Faild=c("C10-20-NGS","RH-NGS")
inbred=c("A626-NGS","E463-NGS")

allBurdenIndex_Constraint=NULL
all.data_rm=NULL
AllCutOffs=c(2,2.75,3.5)
for (ConservationCutoff in AllCutOffs){
  LandranceBurden=read.table(paste0("AllAccession_DeleBurden_BothMinor_BurdenInfo_LandracePlusNewPopMR0.5Dele_WholeGenome_GERP",ConservationCutoff),header=T)
  # allBurdenIndex_Constraint=rbind(allBurdenIndex_Constraint,LandranceBurden)

  LandranceBurden$Type="Landrace"
  LandranceBurden$Type[(!grepl("PG",LandranceBurden$sample.id))]="Founder"
  LandranceBurden$B_Genetic=LandranceBurden$Burden_index_EachiLine_Conserved*0.5
  LandranceBurden$B_Homo=LandranceBurden$Burden_index_Homo_Conserved
  LandranceBurden$B_Heter=LandranceBurden$Burden_index_Heter_Conserved
  LandranceBurden$B_Expressed_h0.15=LandranceBurden$Burden_index_Homo_Conserved+0.15*LandranceBurden$Burden_index_Heter_Conserved
  LandranceBurden$B_Expressed_h0.1=LandranceBurden$Burden_index_Homo_Conserved+0.1*LandranceBurden$Burden_index_Heter_Conserved
  LandranceBurden$B_Expressed_h0.2=LandranceBurden$Burden_index_Homo_Conserved+0.2*LandranceBurden$Burden_index_Heter_Conserved
  LandranceBurden$B_Expressed_h0.09=LandranceBurden$Burden_index_Homo_Conserved+0.09*LandranceBurden$Burden_index_Heter_Conserved

  LandranceBurden2=LandranceBurden[is.na(match(LandranceBurden$sample.id,c(inbred,founder_Faild,founder_sucess))),]
  out_ind_rm6=c(order(LandranceBurden2$B_Genetic)[1:3],order(LandranceBurden2$B_Genetic,decreasing=T)[1:3])
  LandranceBurden_rm6=LandranceBurden2[-out_ind_rm6,]

  outlier=boxplot.stats(LandranceBurden2$B_Genetic)$out
  outlier_ind=which(LandranceBurden2$B_Genetic %in% c(outlier))
  LandranceBurden_rmOutlier=LandranceBurden2[-outlier_ind,]
  outlier_ind_num=length(outlier_ind)
  data_rm=data.frame(ConservationCutoff,outlier_ind_num)
  all.data_rm=rbind(all.data_rm,data_rm)

  Inbred=LandranceBurden[!is.na(match(LandranceBurden$sample.id,inbred)),]
  founder_S=LandranceBurden[!is.na(match(LandranceBurden$sample.id,founder_sucess)),]
  founder_F=LandranceBurden[!is.na(match(LandranceBurden$sample.id,founder_Faild)),]
  founder_S$Type="Sucess"
  founder_F$Type="Faild"
  Inbred$Type="Inbred"

  LandranceBurden_rm=LandranceBurden_rmOutlier
  LandranceBurden_rm_lan=rbind(LandranceBurden_rm,founder_S,founder_F)
  LandranceBurden_rm_lan$Type='Landrace'

  LandranceBurden_rm_all=rbind(LandranceBurden_rm_lan,founder_S,founder_F,Inbred,founder_S,founder_F,Inbred,founder_S,founder_F,Inbred)

  Burden_ColIndex=match(c("B_Heter","B_Homo","B_Genetic"),names(LandranceBurden_rm_all))
  p=ggpairs(LandranceBurden_rm_all,columns = Burden_ColIndex, columnLabels = c("B_Heter","B_Homo","B_Genetic"), mapping = aes(color = Type))
  p.f=p+scale_colour_manual(values=c("red","#009E73",'black',"#0072B2"))
  #ggsave(p.f,filename =paste0("Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rm6_test.pdf"),height =6,width =6)
  ggsave(p.f,filename =paste0("Burden_Cor/Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rmOutlier_test.pdf"),height =6,width =6)

  Mark=rbind(founder_S,founder_F,Inbred)
  Mark$Type="Founder"
  LandranceBurden_rm_all2=rbind(LandranceBurden_rm_lan,Mark)
  p=ggpairs(LandranceBurden_rm_all2,columns = Burden_ColIndex, columnLabels = c("B_Heter","B_Homo","B_Genetic"), mapping = aes(color = Type))
  p.f=p+scale_colour_manual(values=c('red','black'))
  #ggsave(p.f,filename =paste0("Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rm6.pdf"),height =4,width =4)
  ggsave(p.f,filename =paste0("Burden_Cor/Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rmOutlier.pdf"),height =6,width =6)

  LandranceBurden_rm_lan2=rbind(LandranceBurden_rm_lan,Inbred)
  LandranceBurden_rm_lan2$Type="Landrace"
  p=ggpairs(LandranceBurden_rm_lan2,columns = Burden_ColIndex, columnLabels = c("B_Heter","B_Homo","B_Genetic"), mapping = aes(color = Type))
  p.f=p+scale_colour_manual(values=c('black'))
  #ggsave(p.f,filename =paste0("Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rm6_justLandrace.pdf"),height =5,width =5)
  ggsave(p.f,filename =paste0("Burden_Cor/Diploid367AddNew_SubstrLandraceMR0.5_GeneticBurdenGERP",ConservationCutoff,"_rmOutlier_justLandrace.pdf"),height =6,width =6)

###Expressed burden
  LandranceBurden_rm_mark=rbind(LandranceBurden_rm_lan,founder_S,founder_F,Inbred)
  index_Burdens=match(c("B_Genetic","B_Homo","B_Heter"),names(LandranceBurden_rm_mark))
  LandranceBurden_Expressed=gather(LandranceBurden_rm_mark,BurdenType,Burden,names(LandranceBurden_rm_mark)[index_Burdens])
  LandranceBurden_Expressed_Landrace=subset(LandranceBurden_Expressed,Type=="Landrace")

  LandranceBurden_Expressed_Sucess=subset(LandranceBurden_Expressed,Type=="Sucess")
  LandranceBurden_Expressed_Faild=subset(LandranceBurden_Expressed,Type=="Faild")
  LandranceBurden_Expressed_Inbred=subset(LandranceBurden_Expressed,Type=="Inbred")

  p=ggscatter(LandranceBurden_Expressed_Landrace, x = "B_Expressed_h0.1", y = "Burden",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE)
  p_addpoint_S=p + geom_point(data=LandranceBurden_Expressed_Sucess,aes(x = B_Expressed_h0.1, y = Burden),size =1, color = '#0072B2')
  p_addpoint_F=p_addpoint_S + geom_point(data=LandranceBurden_Expressed_Faild,aes(x = B_Expressed_h0.1, y = Burden),size =1, color = 'red')
  p_addpoint_Inbred=p_addpoint_F + geom_point(data=LandranceBurden_Expressed_Inbred,aes(x = B_Expressed_h0.1, y = Burden),size =1, color = '#009E73')
  p_final=p_addpoint_Inbred+theme_bw()+facet_wrap(~BurdenType,scales="free")+ xlab("Expressed burden") + ylab("Burden")
  #ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rm6_GERP",ConservationCutoff,".pdf"),height =5,width =6)
  ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rmOutlier_GERP",ConservationCutoff,".pdf"),height =2,width =7)

  ###Expressed burden
  index_h=match(c("B_Expressed_h0.1","B_Expressed_h0.15",  "B_Expressed_h0.2"),names(LandranceBurden_Expressed))
  LandranceBurden_Expressed_h=gather(LandranceBurden_Expressed,h,ExpressedBurden,names(LandranceBurden_Expressed)[index_h])

  LandranceBurden_Expressed_h_Landrace=subset(LandranceBurden_Expressed_h,Type=="Landrace")
  LandranceBurden_Expressed_h_mark=subset(LandranceBurden_Expressed_h,Type!="Landrace")
  LandranceBurden_Expressed_h_Sucess=subset(LandranceBurden_Expressed_h,Type=="Sucess")
  LandranceBurden_Expressed_h_Faild=subset(LandranceBurden_Expressed_h,Type=="Faild")
  LandranceBurden_Expressed_h_Inbred=subset(LandranceBurden_Expressed_h,Type=="Inbred")

  p=ggscatter(LandranceBurden_Expressed_h_Landrace, x = "ExpressedBurden", y = "Burden",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE)

  p_addpoint_S=p + geom_point(data=LandranceBurden_Expressed_h_Sucess,aes(x = ExpressedBurden, y = Burden),size =1, color = '#0072B2')
  p_addpoint_F=p_addpoint_S + geom_point(data=LandranceBurden_Expressed_h_Faild,aes(x = ExpressedBurden, y = Burden),size =1, color = 'red')
  p_addpoint_Inbred=p_addpoint_F + geom_point(data=LandranceBurden_Expressed_h_Inbred,aes(x = ExpressedBurden, y = Burden),size =1, color = '#009E73')
  p_final=p_addpoint_Inbred+theme_bw()+facet_grid(vars(BurdenType),vars(h), scales="free")+ xlab("Expressed burden") + ylab("Burden")
  #ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rm6_GERP",ConservationCutoff,".pdf"),height =5,width =6)
  ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rmOutlier_GERP",ConservationCutoff,"_different.h2.pdf"),height =5.5,width =6)

  p=ggscatter(LandranceBurden_Expressed_h_Landrace, x = "ExpressedBurden", y = "Burden",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE)
  p_final=p+theme_bw()+facet_grid(vars(BurdenType),vars(h), scales="free")+ xlab("Expressed burden") + ylab("Burden")
  #ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rm6_GERP",ConservationCutoff,".pdf"),height =5,width =6)
  ggsave(p_final,filename =paste("Landrance_B_ExpressedVSBurden_rmOutlier_GERP",ConservationCutoff,"_different.h2_noMark.pdf"),height =5.5,width =6)

}
write.table(all.data_rm,"GeneticBurden_OutlierStat.txt",quote=F,sep='\t',row.name=F)


Burden_index=match(c("B_Genetic","B_Homo","B_Heter","B_Expressed_h0.14"),names(LandranceBurden_GERP2_rm))
p=ggpairs(LandranceBurden_GERP2_rm_lan,columns =Burden_index, columnLabels = c("Genetic Burden","Homo Burden","Heter Burden","Expressed burden"))
ggsave(p,filename ="Landrace_BurdenCorrelation_LandraceDerived_justLandrace.pdf",height =5,width =5)


p=ggscatter(LandranceBurden_GERP2_rm_lan, x = "B_Expressed_h0.14", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_final=p+xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Expressed_GERP2_rmOutlier_h0.14.pdf",height =4,width =4)


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
LandranceBurden_diff.h=gather(LandranceBurden_rm,h,B_Expressed,names(LandranceBurden_rm)[c(20:22)])
LandranceBurden_diff.h_Cuoff=LandranceBurden_diff.h[(LandranceBurden_diff.h$ConservationCutoff==3.12| LandranceBurden_diff.h$ConservationCutoff==2),]



LandranceBurden_diff.h_Cuoff_landrace=subset(LandranceBurden_diff.h_Cuoff,Type=="Landrace")
LandranceBurden_diff.h_Cuoff_founder=subset(LandranceBurden_diff.h_Cuoff,Type=="Founder")

landrace_GERP2=subset(LandranceBurden_diff.h_Cuoff_landrace,ConservationCutoff==2)
founder_GERP2=subset(LandranceBurden_diff.h_Cuoff_founder,ConservationCutoff==2)

p=ggscatter(landrace_GERP2, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=founder_GERP2,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_grid(~h,scales="free")+ xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.final.pdf",height =3,width =9)




landrace_GERP3.12=subset(LandranceBurden_diff.h_Cuoff_landrace,ConservationCutoff==3.12)
founder_GERP3.12=subset(LandranceBurden_diff.h_Cuoff_founder,ConservationCutoff==3.12)

p=ggscatter(landrace_GERP3.12, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=founder_GERP3.12,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_grid(~h,scales="free")+ xlab("Expressed burden") + ylab("Genetic burden")
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Expressed_GERP3.12_rmOutlier_diff.h0.0.25_AddC10.20.final.pdf",height =3,width =9)









LandranceBurden_rm_JustLandrace=LandranceBurden_rm[LandranceBurden_rm$Type=="Landrace",]
E8669.C1020.RH_infor=LandranceBurden_rm[LandranceBurden_rm$Type!="Landrace",]


GERP2=subset(LandranceBurden_diff.h,ConservationCutoff==2)
GERP2_rm=subset(LandranceBurden_diff.h,ConservationCutoff==2)
GERP2_rm_h0.3=GERP2_rm[(GERP2_rm$h!="B_Expressed_h0.35" & GERP2_rm$h!="B_Expressed_h0.4" & GERP2_rm$Type=="Landrace"), ]
E8669.C1020.RH_infor_0.0.3=GERP2_rm[(GERP2_rm$h!="B_Expressed_h0.35" & GERP2_rm$h!="B_Expressed_h0.4" & GERP2_rm$Type=="Selected"),]
GERP2_rm_h0.4=GERP2_rm[GERP2_rm$Type=="Landrace", ]
E8669.C1020.RH_infor_0.0.4=GERP2_rm[ GERP2_rm$Type=="Selected",]

p=ggscatter(GERP2_rm_h0.4, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "B_Heter", y = "B_Expressed",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Expressed.VS.B_Heter_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "B_Heter", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Heter_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.4, x = "B_Homo", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Homo_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.3, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Genetic.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.3, x = "B_Heter", y = "B_Expressed",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Expressed.VS.B_Heter_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)

p=ggscatter(GERP2_rm_h0.3, x = "B_Homo", y = "B_Expressed",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.3,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Expressed.VS.B_Homo_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.4, x = "B_Expressed", y = "B_Homo",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_B_Expressed.VS.B_Homo_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)


p=ggscatter(GERP2_rm_h0.4, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor_0.0.4,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.4_AddC10.20.pdf",height =6,width =15)



p=ggscatter(GERP2_rm_h0.3, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.3_AddC10.20.pdf",height =6,width =15)


GERP2_rm_h0.3_2=subset(GERP2_rm_h0.3,Type!="Selected")
p=ggscatter(GERP2_rm_h0.3_2, x = "B_Expressed", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = B_Expressed, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~h,scales="free",nrow=2)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_GERP2_rmOutlier_diff.h0.0.3_AddC10.20_2.pdf",height =6,width =15)



p=ggscatter(LandranceBurden_rm, x = "B_Expressed_h0.2", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_addpoint=p + geom_point(data=E8669.C1020.RH_infor,aes(x = B_Expressed_h0.2, y = B_Genetic),size =1, color = 'red')
p_final=p_addpoint+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.2_DifferentCutOff_rmOutlier_AddC10.20.pdf",height =3,width =15)





p=ggscatter(LandranceBurden_rm, x = "B_Expressed_h0.1", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.1_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden_rm, x = "B_Expressed_h0.3", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.3_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden_rm, x = "B_Expressed_h0.4", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.4_DifferentCutOff_rmOutlier.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "B_Heter", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.HeterBurden_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "B_Homo", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)
p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Homo_DifferentCutOff.pdf",height =3,width =15)

p=ggscatter(LandranceBurden, x = "B_Expressed_h0.2", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.2_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "B_Expressed_h0.2", y = "FutureB_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_FutureB_Genetic.VS.B_Expressed_h0.2_DifferentCutOff.pdf",height =3,width =15)


p=ggscatter(LandranceBurden, x = "B_Expressed_h0.1", y = "B_Genetic",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 5, label.sep = "\n")
)

p_final=p+theme_bw()+facet_wrap(~Cutoff_Id,scales="free",nrow=1)
ggsave(p_final,filename ="Landrance_GeneticBurden.VS.B_Expressed_h0.1_DifferentCutOff.pdf",height =3,width =15)


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
