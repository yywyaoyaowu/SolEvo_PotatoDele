# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/23
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/F2/")
#Model_h=read.table("GP_fixEffect_parameter_Stat.txt.final.txt",header=T)
#GP_r=read.table("GP_differntGERP.CutOff_cross.validation.mean.final.txt",header=T)
Model_h=read.table("GP_fixEffect_parameter_Stat_seq2.3.5_0.2.txt",header=T)
Model_h$h=Model_h$L.Coefficients_Heter_rmOutlier/Model_h$L.Coefficients_Homo_rmOutlier
library(ggplot2)

p=ggplot(Model_h, aes(x=Constraint_Cutoff, y=h)) +geom_line()+theme_bw()#x为药剂剂量，并非连续型变量
p.axis=p +facet_wrap(~trait,scales="free",nrow=1)+ theme(axis.text.x = element_text(size = 10, color = "black", vjust =1, hjust =1, angle = 0))+ theme(axis.text.y = element_text(size = 10))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("h (estimate dominance coefficient)")+ theme(axis.title.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 10))+theme(legend.title = element_text(size = 15))
ggsave(p.lenged,filename ="Model.estimated.h_rm6Outlier_GERP.2to3.5.pdf",height = 3,width = 12)

Model_h_yield=subset(Model_h,trait=="TotalYield")
p=ggplot(Model_h_yield, aes(x=Constraint_Cutoff, y=h)) +geom_line()+theme_bw()#x为药剂剂量，并非连续型变量
p.axis=p + theme(axis.text.x = element_text(size = 10, color = "black", vjust =1, hjust =1, angle = 0))+ theme(axis.text.y = element_text(size = 10))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("h (estimate dominance coefficient)")+ theme(axis.title.x = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 10))+theme(legend.title = element_text(size = 15))
ggsave(p.lenged,filename ="Model.estimated.h_rm6Outlier_GERP.2to3.5_yield.pdf",height = 3,width = 2.5)
