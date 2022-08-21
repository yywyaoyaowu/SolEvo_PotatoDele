# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/4/25
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/F2")
library(GGally)
library(ggplot2)
library(ggdensity)

library(ppcor)
##input
phenotype_file="/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/09_GenomicPrediction/F2/F2.phenotype_normalization.txt"
##Burden per Accession
Burden_all=read.table("E463_A626_F2_Burden_Eachline_GERP_AllDifferentCutOff.final.txt",header=T)
Burden_all$HomoBurden=(Burden_all$TotalAddBurden_PerAccession-Burden_all$TotalDomBurden_PerAccession)/2
Burden_all$HeterBurden=Burden_all$TotalDomBurden_PerAccession
all.phenotype=read.table(phenotype_file,header=T,sep="\t")
all.phenotype2_ForPlot=all.phenotype[,c(4,2,3,5,1,6,7)]
names(all.phenotype2_ForPlot)[1:5]=c("Yield","Plant height","Tuber number","Tuber size","Flowering time")

for (n in 1:5){
  trait_order=n
  trait=names(all.phenotype)[n]
  all.phenotype2_ForPlot[c(order(all.phenotype2_ForPlot[,n],decreasing = T)[1:4],order(all.phenotype2_ForPlot[,n])[1:4]),n]=NA
}

cor_plot <-gather(all.phenotype2_ForPlot,Trait,value,c(1:5))
##figure 4a
all.data_justpheno=all.data[all.data$BurdenType==as.character(unique(all.data$BurdenType)[1]),]
p=ggplot(cor_plot, aes(x=value))+ geom_histogram(fill="white", alpha=0.5, position="identity",colour='black') + theme_bw()
p.axis=p +facet_wrap(.~Trait,nrow=1,scales="free_x")+ theme(axis.text.x = element_text(size = 10))+ theme(axis.text.y = element_text(size = 10))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Frequency")+ theme(axis.title.x = element_text(size = 10))+ theme(axis.title.y = element_text(size = 10))
p.lenged=p.lab+theme(legend.text = element_text(size = 10))+theme(legend.title = element_text(size = 10))
ggsave(p.lenged,filename ="PhenotypeDistribution_rm6Outlier.pdf",height =2,width =9)



##figure 4b

index=match(Burden_all$samples,all.phenotype2_ForPlot$AccessionId)
Burden_all_Pheno=cbind(Burden_all,all.phenotype2_ForPlot[index,c(1:5)])
names(Burden_all_Pheno)[7:8]=c("HomoBurden","HeterBurden")
all.data_1 <-gather(Burden_all_Pheno,BurdenType,Burden,c(7:8))

all.data  <-gather(all.data_1,Trait,Pheno,c(7:11))
table(all.data$Trait)

##figure 4b
all.data$BurdenType=factor(all.data$BurdenType,levels=unique(all.data$BurdenType))
all.data$Trait=factor(all.data$Trait,levels=unique(all.data$Trait))
library(ggpubr)
AllCutOffs=c(2,2.75,3.5)
all.pcor_result=NULL
for (ConservationCutoff in AllCutOffs){
  all.data_constrainted=all.data[all.data$Constraint_Cutoff==ConservationCutoff,]
  Burden_all_Pheno_OneCutoff=subset(Burden_all_Pheno,Constraint_Cutoff==ConservationCutoff)

  p=ggscatter(all.data_constrainted, x = "Pheno", y = "Burden",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE)
  p_final=p+theme_bw()+facet_grid(BurdenType~Trait,scales="free")
  ggsave(p_final,filename =paste0("Cor_Burden_Traits_rm6Outlier_GERP",ConservationCutoff,".pdf"),height =4.5,width =9)

######Homo burden VS trait
  all.data_constrainted_Homo=all.data_constrainted[all.data_constrainted$BurdenType=="HomoBurden",]
  p=ggscatter(all.data_constrainted_Homo, x = "Pheno", y = "Burden",
            color = "black", size = .3,# Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE,
            cor.coef = TRUE)
  p_final=p+theme_bw()+facet_grid(BurdenType~Trait,scales="free")
  ggsave(p_final,filename =paste0("Cor_HomoBurden_Traits_rm6Outlier_GERP",ConservationCutoff,".pdf"),height =3,width =9)
########Heter burden VS trait
  all.data_constrainted_Heter=all.data_constrainted[all.data_constrainted$BurdenType=="HeterBurden",]
  p=ggscatter(all.data_constrainted_Heter, x = "Pheno", y = "Burden",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE,
              cor.coef = TRUE)
  p_final=p+theme_bw()+facet_grid(BurdenType~Trait,scales="free")
  ggsave(p_final,filename =paste0("Cor_HeterBurden_Traits_rm6Outlier_GERP",ConservationCutoff,".pdf"),height =2,width =9)

  ###Homo VS heter
  p=ggscatter(Burden_all_Pheno_OneCutoff, x = "HomoBurden", y = "HeterBurden",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE,
              cor.coef = TRUE)
  p_final=p+theme_bw()
  ggsave(p_final,filename =paste0("Cor_HomoHeterBurden_GERP",ConservationCutoff,".pdf"),height =4,width =4)


  ###residual
HomoBurden_ResidualByHeter_all=summary(lm(Burden_all_Pheno_OneCutoff$HomoBurden ~ Burden_all_Pheno_OneCutoff$HeterBurden))$residuals
HeterBurden_ResidualByHomo_all=summary(lm(Burden_all_Pheno_OneCutoff$HeterBurden ~Burden_all_Pheno_OneCutoff$HomoBurden))$residuals

  all.Tait_phenoBurdenResidual=NULL
for (n in 9:13){
  trait_order=n
  trait=names(Burden_all_Pheno_OneCutoff)[n]
  index=!is.na(Burden_all_Pheno_OneCutoff[,n])
  #By heter
  TraitPheno_Residual=summary(lm(Burden_all_Pheno_OneCutoff[,n] ~ Burden_all_Pheno_OneCutoff$HeterBurden))$residuals
  Burden_Residual=HomoBurden_ResidualByHeter_all[!is.na(Burden_all_Pheno_OneCutoff[,n] )]
  BurdenType="HomoBurden_Residual"
  one.traint_pheno=data.frame(trait,ConservationCutoff,TraitPheno_Residual,Burden_Residual,BurdenType)
  all.Tait_phenoBurdenResidual=rbind(all.Tait_phenoBurdenResidual,one.traint_pheno)

  #by homo
  TraitPheno_Residual=summary(lm(Burden_all_Pheno_OneCutoff[,n] ~ Burden_all_Pheno_OneCutoff$HomoBurden))$residuals
  Burden_Residual=HeterBurden_ResidualByHomo_all[!is.na(Burden_all_Pheno_OneCutoff[,n] )]
  BurdenType="HeterBurden_Residual"
  one.traint_pheno=data.frame(trait,ConservationCutoff,TraitPheno_Residual,Burden_Residual,BurdenType)
  all.Tait_phenoBurdenResidual=rbind(all.Tait_phenoBurdenResidual,one.traint_pheno)

  HomoBurden_Phenotype_cor=cor(Burden_all_Pheno_OneCutoff[index,n] , Burden_all_Pheno_OneCutoff$HomoBurden[index])
  HeterBurden_Phenotype_cor=cor(Burden_all_Pheno_OneCutoff[index,n] , Burden_all_Pheno_OneCutoff$HeterBurden[index])
  data_por=data.frame(Burden_all_Pheno_OneCutoff[!is.na(Burden_all_Pheno_OneCutoff[,n]),c(n,7,8)])
  pcor_fit=pcor(data_por)

  HomoBurden_Phenotype_estimate=pcor(data_por)$estimate[2,1]
  HeterBurden_Phenotype_estimate=pcor(data_por)$estimate[3,1]
  HomoBurden_Phenotype_p.value=pcor(data_por)$p.value[2,1]
  HeterBurden_Phenotype_p.value=pcor(data_por)$p.value[3,1]
  pcor_result=data.frame(trait,ConservationCutoff,HomoBurden_Phenotype_cor,HomoBurden_Phenotype_estimate,HomoBurden_Phenotype_p.value,HeterBurden_Phenotype_cor,HeterBurden_Phenotype_estimate,HeterBurden_Phenotype_p.value)
  all.pcor_result=rbind(all.pcor_result,pcor_result)
}

  ######Homo burden VS trait
  all.Tait_phenoBurdenResidual$BurdenType=factor(all.Tait_phenoBurdenResidual$BurdenType,levels=unique(all.Tait_phenoBurdenResidual$BurdenType))
  all.Tait_phenoBurdenResidual$trait=factor(all.Tait_phenoBurdenResidual$trait,levels=unique(all.Tait_phenoBurdenResidual$trait))
  p=ggscatter(all.Tait_phenoBurdenResidual, x = "TraitPheno_Residual", y = "Burden_Residual",
              color = "black", size = .3,# Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x = 5)
  )
  p_final=p+theme_bw()+facet_grid(BurdenType~trait,scales="free")
  print(nrow(all.Tait_phenoBurdenResidual))
  ggsave(p_final,filename =paste0("ResidualCor_Burden_Traits_rm6Outlier_GERP",ConservationCutoff,".pdf"),height =5,width =9)
}
write.table(all.pcor_result,"all.pcor_result_final.txt",quote=F,sep='\t',row.name=F)
