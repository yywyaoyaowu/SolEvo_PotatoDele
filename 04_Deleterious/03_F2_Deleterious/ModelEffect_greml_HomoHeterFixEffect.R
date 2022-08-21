# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/19
library(SNPRelate)
library(regress)
library(tidyr)
library(ggplot2)
library(qgg)
require(caret)

###input
phenotype_file="/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/09_GenomicPrediction/02_F2/F2.phenotype_normalization.txt"
genotype_additiveGRM_file="/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/10_NewBin_GP/F2_bin_GRM_FromGuillaume_NewBin2603.txt"
F2_BurdenEachAccession_file_prefix="E463_A626_F2_Burden_Eachline_GERP"
#BurdenEachAccession_AllCutOff.file="E463_A626_F2_Burden_Eachline_GERP_AllDifferentCutOff.final.txt"
#BurdenEachAccession_AllCutOff.file="E463_A626_F2_Burden_Eachline_GERP_AllDifferentCutOff.final_NewCutOff.txt"
BurdenEachAccession_AllCutOff.file="E463_A626_F2_Burden_Eachline_GERP_AllDifferentCutOff.final.txt"
#weightedAddGRM_file="/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/10_NewBin_GP/01_F2_Burden/F2_AdditiveBurdenMatrix_GERP2_weightedGRM.txt"
#weightedHomoGRM_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/10_NewBin_GP/01_F2_Burden/F2_HomoBurdenMatrix_GERP2_weightedGRM.txt")
#weightedDomGRM_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/10_NewBin_GP/01_F2_Burden/F2_DomainceBurdenMatrix_GERP2_weightedGRM.txt")
#BurdenEachAccession_file="/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/10_NewBin_GP/01_F2_Burden//E463_A626_F2_Burden_Eachline_2603Bin.txt"

#output
Model_parameter_file="GP_fixEffect_parameter_Stat_NewCutOff.txt"
#all.GPresult "GP_differntGERP.CutOff_each.cross.validation.txt"
#all.GP.median "GP_differntGERP.CutOff_cross.validation.mean.txt"

##### read file
#Phenotype
all.phenotype=read.table(phenotype_file,header=T,sep="\t")
#GRM matrix
GRM_Genome=read.table(genotype_additiveGRM_file,header=T)
row.names(GRM_Genome)=names(GRM_Genome)
sample.id=names(GRM_Genome)
#ConstraitCutoff_infor=read.table(ConstraitCutoff_file,header=T)
BurdenEachAccession_AllCutOff=read.table(BurdenEachAccession_AllCutOff.file,header=T)
#allCutOffs=seq(2,4.1,by=0.1)
#allCutOffs=c(2.75)
#allCutOffs=c(2.5,2.63,2.73,2.75,3,3.03,3.5,3.91)
allCutOffs=c(3.5,2,2.75)
all.result=NULL
all.GPresult=NULL
all.GP.median=NULL
##sample with phenotype
for (n in 1:5){
 trait_order=n
trait=names(all.phenotype)[n]
##get all the samples with phenotype data
  sample_withpheno=all.phenotype$AccessionId[!is.na(all.phenotype[,n])]
  y=all.phenotype[match(sample_withpheno,all.phenotype$AccessionId),n]
  # whole genome GRM
  index=match(sample_withpheno,sample.id)
  GRM_all=as.matrix(GRM_Genome[index,index])
##Burden for each individual
#for (i in 1:length(allCutOffs)){
for (Constraint_Cutoff in allCutOffs){
  #Constraint_Cutoff=allCutOffs[i]
  BurdenEachAccession=BurdenEachAccession_AllCutOff[BurdenEachAccession_AllCutOff$Constraint_Cutoff==Constraint_Cutoff,]
  #print (paste(n,i))
  print (paste(n,Constraint_Cutoff))
##fix effect by homo and heter burden
AddBurden=BurdenEachAccession$TotalAddBurden_PerAccession[match(sample_withpheno,BurdenEachAccession$samples)]
HeterBurden=BurdenEachAccession$TotalDomBurden_PerAccession[match(sample_withpheno,BurdenEachAccession$samples)]
HomoBurden=(AddBurden-HeterBurden)*0.5
#HomoHeter_fix=cbind("(Intercept)"=1,data.frame(HomoBurden,HeterBurden))

greml_HomoHeter_AddGRM=greml(y = y, X =cbind(1, HomoBurden, HeterBurden), GRM = list(GRM=GRM_all))
V_GRM=greml_HomoHeter_AddGRM$theta[1]
V_E=greml_HomoHeter_AddGRM$theta[2]
L.Coefficients_Homo=greml_HomoHeter_AddGRM$b[2]
L.Coefficients_Heter=greml_HomoHeter_AddGRM$b[3]
  b_hat <- greml_HomoHeter_AddGRM$b   # beta estimates (coefficient)
  b_se <- sqrt(diag(greml_HomoHeter_AddGRM$vb))   # standard error of beta estimates (SE)
  z <- c(b_hat/b_se)   # z-score of beta estimates (coefficient/SE)
  pval <- pchisq(z^2, 1, lower.tail=FALSE)
  L.Coefficients_Homo_p.value=pval[2]
  L.Coefficients_Heter_p.value=pval[3]

indexRm=c(order(y,decreasing=T)[1:3],order(y)[1:3])
y_rm=y[-indexRm]
HomoBurden_rm=HomoBurden[-indexRm]
HeterBurden_rm=HeterBurden[-indexRm]
GRM_all_rm=GRM_all[-indexRm,-indexRm]
greml_HomoHeter_AddGRM_rm=greml(y = y_rm, X =cbind(1, HomoBurden_rm, HeterBurden_rm), GRM = list(GRM=GRM_all_rm))
V_GRM_rmOutlier=greml_HomoHeter_AddGRM_rm$theta[1]
V_E_rmOutlier=greml_HomoHeter_AddGRM_rm$theta[2]
L.Coefficients_Homo_rmOutlier=greml_HomoHeter_AddGRM_rm$b[2]
L.Coefficients_Heter_rmOutlier=greml_HomoHeter_AddGRM_rm$b[3]
  b_hat_rm <- greml_HomoHeter_AddGRM_rm$b   # beta estimates (coefficient)
  b_se_rm <- sqrt(diag(greml_HomoHeter_AddGRM_rm$vb))   # standard error of beta estimates (SE)
  z_rm <- c(b_hat_rm/b_se_rm)   # z-score of beta estimates (coefficient/SE)
  pval_rm <- pchisq(z_rm^2, 1, lower.tail=FALSE)
  L.Coefficients_Homo_rmOutlier_p.value=pval_rm[2]
  L.Coefficients_Heter_rmOutlier_p.value=pval_rm[3]

  Result=data.frame(trait,Constraint_Cutoff,V_GRM,V_E,L.Coefficients_Homo,L.Coefficients_Heter,L.Coefficients_Homo_p.value,L.Coefficients_Heter_p.value,V_GRM_rmOutlier,V_E_rmOutlier,L.Coefficients_Homo_rmOutlier,L.Coefficients_Heter_rmOutlier,L.Coefficients_Homo_rmOutlier_p.value,L.Coefficients_Heter_rmOutlier_p.value)
all.result=rbind(all.result,Result)

##genomic prediction and selection
  k <- 1
  set.seed(k)
  sampleSize=length(y)
  fold =5
  nvalid=5*20
  validate <- replicate(nvalid, sample(1:sampleSize, as.integer(sampleSize / fold)))
  cv_HomoHeterBurdenFix_addGRM <- greml(y = y, X =cbind(1, HomoBurden, HeterBurden), GRM = list(GRM_all), validate = validate)
  cor_HomoHeterBurdenFix_addGRM=cv_HomoHeterBurdenFix_addGRM$accuracy$Corr
  r2_HomoHeterBurdenFix_addGRM=cor_HomoHeterBurdenFix_addGRM*cor_HomoHeterBurdenFix_addGRM
  cor_HomoHeterBurdenFix_addGRM_median20crossvaliation=median(cor_HomoHeterBurdenFix_addGRM)
  r2_HomoHeterBurdenFix_addGRM_median20crossvaliation=median(r2_HomoHeterBurdenFix_addGRM)

  X_base <- matrix(rep(1),length(y))
  cv_addGRM <- greml(y = y, X =X_base, GRM = list(GRM_all), validate = validate)
  cor_addGRM=cv_addGRM$accuracy$Corr
  r2_addGRM=cor_addGRM*cor_addGRM
  cor_addGRM_median20crossvaliation=median(cor_addGRM)
  r2_addGRM_median20crossvaliation=median(r2_addGRM)

  k <- 1
  set.seed(k)
  sampleSize_rm=length(y_rm)
  fold =5
  nvalid=5*20
  validate_rm <- replicate(nvalid, sample(1:sampleSize_rm, as.integer(sampleSize_rm / fold)))
  cv_HomoHeterBurdenFix_addGRM_rmOutlier <- greml(y = y_rm, X =cbind(1, HomoBurden_rm, HeterBurden_rm), GRM = list(GRM=GRM_all_rm), validate = validate_rm)
  cor_HomoHeterBurdenFix_addGRM_rmOutlier=cv_HomoHeterBurdenFix_addGRM_rmOutlier$accuracy$Corr
  r2_HomoHeterBurdenFix_addGRM_rmOutlier=cor_HomoHeterBurdenFix_addGRM_rmOutlier*cor_HomoHeterBurdenFix_addGRM_rmOutlier
  cor_HomoHeterBurdenFix_addGRM_rmOutlier_median20crossvaliation=median(cor_HomoHeterBurdenFix_addGRM_rmOutlier)
  r2_HomoHeterBurdenFix_addGRM_rmOutlier_median20crossvaliation=median(r2_HomoHeterBurdenFix_addGRM_rmOutlier)


  X_base_rm <- matrix(rep(1),length(y_rm))
  cv_addGRM_rm <- greml(y = y_rm, X =X_base_rm, GRM = list(GRM_all_rm), validate = validate_rm)
  cor_addGRM_rm=cv_addGRM_rm$accuracy$Corr
  r2_addGRM_rm=cor_addGRM_rm*cor_addGRM_rm
  cor_addGRM_median20crossvaliation_rmOutlier=median(cor_addGRM_rm)
  r2_addGRM_median20crossvaliation_rmOutlier=median(r2_addGRM_rm)

  GPresult=data.frame(trait,Constraint_Cutoff,
                      cor_addGRM,r2_addGRM,
                      cor_HomoHeterBurdenFix_addGRM,r2_HomoHeterBurdenFix_addGRM,
                      cor_addGRM_rm,r2_addGRM_rm,
                      cor_HomoHeterBurdenFix_addGRM_rmOutlier,r2_HomoHeterBurdenFix_addGRM_rmOutlier)

  GP.median=data.frame(trait,Constraint_Cutoff,
                       cor_addGRM_median20crossvaliation,r2_addGRM_median20crossvaliation,
                       cor_HomoHeterBurdenFix_addGRM_median20crossvaliation,r2_HomoHeterBurdenFix_addGRM_median20crossvaliation,
                       cor_addGRM_median20crossvaliation_rmOutlier,r2_addGRM_median20crossvaliation_rmOutlier,
                       cor_HomoHeterBurdenFix_addGRM_rmOutlier_median20crossvaliation,r2_HomoHeterBurdenFix_addGRM_rmOutlier_median20crossvaliation)

  all.GPresult=rbind(all.GPresult,GPresult)
  all.GP.median=rbind(all.GP.median,GP.median)
}
}
#all.result$CutOff_ID=ConstraitCutoff_infor$CutOff_ID[match(all.result$Constraint_Cutoff,ConstraitCutoff_infor$CutOff_GERP.score)]
write.table(all.result,Model_parameter_file,quote=F,sep='\t',row.name=F)

write.table(all.GPresult,"GP_differntGERP.CutOff_each.cross.validation.final_NewCutOff.txt",quote=F,sep='\t',row.name=F)
write.table(all.GP.median,"GP_differntGERP.CutOff_cross.validation.mean.final_NewCutOff.txt",quote=F,sep='\t',row.name=F)
