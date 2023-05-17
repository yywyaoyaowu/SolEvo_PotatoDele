# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/27
all.mean=NULL
all.cross=NULL
AllCutOffs=c(2,2.75,3.5)
for (ConservationCutoff in AllCutOffs){
  for (p in 1:100){
    one_median=read.table(paste0("FixEffect_GP.Result_GERP",ConservationCutoff,"/GP_differntGERP_Permutation_",p,"_cross.validation.mean.txt"),header=T)
    one_cross=read.table(paste0("FixEffect_GP.Result_GERP",ConservationCutoff,"/GP_differntGERP_Permutation_",p,"_each.cross.validation.txt"),header=T)
    all.mean=rbind(all.mean,one_median)
    all.cross=rbind(all.cross,one_cross)
  }
}
all.cross$Model="Permutation"

write.table(all.cross,"F2_500K.SNP.Block.100TimesPermutation_AddGRM_HomoHeter_DiffCutOffs_eachCrossvalidation.txt",quote=F,sep='\t',row.name=F)
write.table(all.mean,"F2_500K.SNP.Block.100TimesPermutation_AddGRM_HomoHeter_DiffCutOffs_eachPermutationMean.txt",quote=F,sep='\t',row.name=F)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette3=cbbPalette[c(3,2,4)]
GP_cutoffs=read.table("GP_differntGERP.CutOff_each.cross.validation.final_NewCutOff.txt",header=T)
GP_2.75=read.table("GP_differntGERP.CutOff_each.cross.validation.final_NewCutOff.txt_GERP2.75",header=T)
GP=rbind(GP_cutoffs,GP_2.75)
GP_baseModel=GP[,c(1,2,8)]
GP_baseModel$Model="Baseline"
GP_BurdenModel=GP[,c(1,2,10)]
GP_BurdenModel$Model="DeleteriousBurden"
PermutatedBurden=all.cross[,c(1,2,6,7)]
names(GP_baseModel)=names(GP_BurdenModel)
all.data=rbind(GP_baseModel,GP_BurdenModel,PermutatedBurden)
names(all.data)[3]="PredictionAcccuracy_r2"


library(ggplot2)
all.data$trait=factor(all.data$trait,levels=unique(all.data$trait)[c(4,2,3,5,1)])
AllCutOffs=c(2,2.75,3.5)
for (ConservationCutoff in AllCutOffs){
all.data_GERPConstrainted=subset(all.data,Constraint_Cutoff==ConservationCutoff)
p=ggplot(all.data_GERPConstrainted, aes(x=Model,y=PredictionAcccuracy_r2,fill=Model)) + geom_boxplot(width=.6,outlier.colour = NA)+facet_wrap(.~trait,nrow=1) + scale_fill_manual(values=cbbPalette3)+coord_cartesian(ylim = c(0, 0.45))+theme_bw()
#+coord_cartesian(ylim = c(0.1, 0.6))
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust =0.5, hjust =0.5, angle = 0))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("") + ylab("Prediction Accuracy (r2)")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
ggsave(p.lenged,filename =paste0("GP_BurdenFixModel_5FoldCrossValidation_100Permutations_r2_GERP",ConservationCutoff,".pdf"),height =3,width =9)
}


AllResults_AllCutoff=NULL
AllCutOffs=c(2,2.75,3.5)
alltraits=as.character(unique(all.data$trait))
allModel=as.character(unique(all.data$Model))
for (ConservationCutoff in AllCutOffs){
  AllResults=NULL
for (Trait in alltraits){
  all.data_GERPConstrainted=subset(all.data,Constraint_Cutoff==ConservationCutoff)
  r_stat_trait=subset(all.data_GERPConstrainted,trait==Trait)
  GP_baselineModel_trait=subset(r_stat_trait,Model=="Baseline")
  GP_BurdenModel_trait=subset(r_stat_trait,Model=="DeleteriousBurden")
  GP_PermutationModel_trait=subset(r_stat_trait,Model=="Permutation")

  Absolute_P.value=t.test(GP_BurdenModel_trait$PredictionAcccuracy_r2,GP_baselineModel_trait$PredictionAcccuracy_r2,alternative ="greater")$p.value
  Relative_P.value=t.test(GP_BurdenModel_trait$PredictionAcccuracy_r2,GP_PermutationModel_trait$PredictionAcccuracy_r2,alternative ="greater")$p.value

  GP_baselineModel_trait_median=median(GP_baselineModel_trait$PredictionAcccuracy_r2)
  GP_BurdenModel_trait_median=median(GP_BurdenModel_trait$PredictionAcccuracy_r2)
  GP_PermutationModel_trait_median=median(GP_PermutationModel_trait$PredictionAcccuracy_r2)

  GP_baselineModel_trait_mean=mean(GP_baselineModel_trait$PredictionAcccuracy_r2)
  GP_BurdenModel_trait_mean=mean(GP_BurdenModel_trait$PredictionAcccuracy_r2)
  GP_PermutationModel_trait_mean=mean(GP_PermutationModel_trait$PredictionAcccuracy_r2)

  Result=data.frame(Trait,ConservationCutoff,GP_baselineModel_trait_median,GP_BurdenModel_trait_median,GP_PermutationModel_trait_median,Absolute_P.value,Relative_P.value,GP_baselineModel_trait_mean,GP_BurdenModel_trait_mean,GP_PermutationModel_trait_mean)
  AllResults=rbind(AllResults,Result)
}

  AllResults=AllResults[c(4,2,3,5,1),]
  AllResults$Absolute_increase_r2=(AllResults$GP_BurdenModel_trait_median/AllResults$GP_baselineModel_trait_median-1)
  AllResults$Relative_increase_r2=(AllResults$GP_BurdenModel_trait_median/AllResults$GP_PermutationModel_trait_mean-1)
  AllResults_AllCutoff=rbind(AllResults_AllCutoff,AllResults)
}

write.table(AllResults_AllCutoff,"GP.summary_BurdenFixEffect_R.square.P.value_DiffCutOffs.txt",quote=F,sep='\t',row.name=F)

