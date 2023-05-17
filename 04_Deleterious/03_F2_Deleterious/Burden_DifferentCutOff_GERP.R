# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/4/24
library(qgg)
#input
snp_infor_homo_DeleStat_file="A626E454_homoSNP_Constraint_infor_AllChrs.txt"
Bin_geno_file="A.E-F2_genotype_RefNum_CombinBin_2603Bin.txt"
#####
#output
F2_BurdenEachAccession_file_prefix="E463_A626_F2_Burden_Eachline_GERP"
F2_bin_burden_file_prefix="E463_A626_F2_2603BinBurden_GERP"
BurdenMatrix_prefix="F2_BurdenMatrix_GERP"
#ConstraitCutoff_file="F2_SNP_GERP.cutoff.txt"
#######
#read files
Bin_geno=read.table(Bin_geno_file,header=T)
samples=names(Bin_geno)[-c(1:7)]
Bin_geno_matrix=Bin_geno[,-c(1:7)]

snp_infor_homo_AllChrs=read.table(snp_infor_homo_DeleStat_file,header=T,sep='\t')
snp_infor_homo_AllChrs$A626_AlleleDerived=(snp_infor_homo_AllChrs$A626_Allele!=snp_infor_homo_AllChrs$SolMsa_Major & !is.na(snp_infor_homo_AllChrs$SolMsa_Major))
snp_infor_homo_AllChrs$E454_AlleleDerived=(snp_infor_homo_AllChrs$E454_Allele!=snp_infor_homo_AllChrs$SolMsa_Major & !is.na(snp_infor_homo_AllChrs$SolMsa_Major))
snp_infor_homo_AllChrs_Heter=snp_infor_homo_AllChrs[(!(snp_infor_homo_AllChrs$A626_AlleleDerived & snp_infor_homo_AllChrs$E454_AlleleDerived) & !is.na(snp_infor_homo_AllChrs$SolMsa_Major)) ,]

snp_infor_homo_AllChrs_constrainted=subset(snp_infor_homo_AllChrs,Constraint>=2)
snp_infor_homo_AllChrs_constrainted.sort=snp_infor_homo_AllChrs_constrainted[order(snp_infor_homo_AllChrs_constrainted$Constraint,decreasing=T),]
#all.Cutoff=read.table(ConstraitCutoff_file,header=T)
SNPNum=nrow(snp_infor_homo_AllChrs_Heter)
ConstraintedSNPNUM=sum(snp_infor_homo_AllChrs_Heter$Constraint>=2)

#ConstraitCutoff_infor=all.Cutoff[all.Cutoff$GERP_Cutoff>=2 & !duplicated(all.Cutoff$GERP_Cutoff) ,]

all.stat_info=NULL
BurdenPerAccession=NULL
all.BurdenPerAccession=NULL
all.Bin_Dele_AllCutOff=NULL
F2_bin_burden_all=NULL

all.Bin_Dele_AllCutOff=read.table(paste0(F2_bin_burden_file_prefix,"_AllDiffCutOffs.final_NewCutOff.txt"),header=T)
all.BurdenPerAccession=read.table(paste0(F2_BurdenEachAccession_file_prefix,"_AllDifferentCutOff.final_NewCutOff.txt"),header=T)

#allCutOffs=seq(2,4.1,by=0.1)
allCutOffs=c(2.75)
#allCutOffs=c(2.5,2.63,2.73,3,3.03,3.5,3.91)
#all.Bin_Dele_AllCutOff=NULL
for (i in 1:length(allCutOffs)){

  Constraint_Cutoff=allCutOffs[i]
  all.Bin_Dele_Onecutoff=NULL
  Bin_Burden_matrix=matrix(rep(NA,nrow(Bin_geno_matrix)*ncol(Bin_geno_matrix)),ncol=ncol(Bin_geno_matrix))
  Bin_DomanianceBurden_matrix=matrix(rep(NA,nrow(Bin_geno_matrix)*ncol(Bin_geno_matrix)),ncol=ncol(Bin_geno_matrix))
  Bin_HomoBurden_matrix=matrix(rep(NA,nrow(Bin_geno_matrix)*ncol(Bin_geno_matrix)),ncol=ncol(Bin_geno_matrix))
  Bin_HeterSNP_matrix=matrix(rep(NA,nrow(Bin_geno_matrix)*ncol(Bin_geno_matrix)),ncol=ncol(Bin_geno_matrix))

for (n in 1:nrow(Bin_geno)){
  ##get the bin burden
  Chr=Bin_geno$Chr[n]
  Bin=Bin_geno$bin[n]
  Start=Bin_geno$NewBinStart[n]
  End=Bin_geno$NewBinEnd[n]

  Chrom_E464=paste("E4-63",strsplit(Chr,'_')[[1]][2],sep='_')
  snp_infor_bin=snp_infor_homo_AllChrs_Heter[(snp_infor_homo_AllChrs_Heter$position>=Start & snp_infor_homo_AllChrs_Heter$position<=End & snp_infor_homo_AllChrs_Heter$Chrom==Chrom_E464),]
  polymorphismNum=nrow(snp_infor_bin)
  snp_infor_bin_Constraint=subset(snp_infor_bin,Constraint>=Constraint_Cutoff)
  if(nrow(snp_infor_bin_Constraint)>0){
    polymorphismNum_Constrait=nrow(snp_infor_bin_Constraint)
    Burden=sum(snp_infor_bin_Constraint$Constraint)
    Burden_E454=sum(snp_infor_bin_Constraint$E454Allele_Burden,na.rm=T)
    Burden_A626=sum(snp_infor_bin_Constraint$A626Allele_Burden,na.rm=T)
  }else{
    polymorphismNum_Constrait=Burden=Burden_E454=Burden_A626=0
  }

  Bin_Dele=data.frame(Chr,Bin,Start,End,Constraint_Cutoff,polymorphismNum,polymorphismNum_Constrait,Burden,Burden_E454,Burden_A626)
  all.Bin_Dele_Onecutoff=rbind(all.Bin_Dele_Onecutoff,Bin_Dele)

  index_E454=(Bin_geno_matrix[n,]==2)
  index_H=(Bin_geno_matrix[n,]==1)
  index_A626=(Bin_geno_matrix[n,]==0)
  Bin_Burden_matrix[n,index_E454]=Burden_E454*2
  Bin_Burden_matrix[n,index_H]=Burden_E454+Burden_A626
  Bin_Burden_matrix[n,index_A626]=Burden_A626*2

  Bin_DomanianceBurden_matrix[n,index_H]=Burden_E454+Burden_A626
  Bin_DomanianceBurden_matrix[n,index_E454]=0
  Bin_DomanianceBurden_matrix[n,index_A626]=0

    Bin_HomoBurden_matrix[n,index_H]=0
    Bin_HomoBurden_matrix[n,index_E454]=Burden_E454*2
    Bin_HomoBurden_matrix[n,index_A626]=Burden_A626*2

  Bin_HeterSNP_matrix[n,index_E454]=0
  Bin_HeterSNP_matrix[n,index_H]=polymorphismNum
  Bin_HeterSNP_matrix[n,index_A626]=0

}

  #write.table(all.Bin_Dele_Onecutoff,paste0(F2_bin_burden_file_prefix,Constraint_Cutoff,".txt"),quote=F,sep='\t',row.name=F)
  all.Bin_Dele_AllCutOff=rbind(all.Bin_Dele_AllCutOff,all.Bin_Dele_Onecutoff)

  ##Burnden matrix
  Bin_chr_Burden_infor=cbind(all.Bin_Dele_Onecutoff,Bin_Burden_matrix)
 
  ##Domainace matrix
  Bin_chr_DomainaceBurden_infor=cbind(all.Bin_Dele_Onecutoff,Bin_DomanianceBurden_matrix)
 
  ##Homo matrix
  colnames(Bin_HomoBurden_matrix)=samples
  Bin_chr_HomoBurden_infor=cbind(all.Bin_Dele_Onecutoff,Bin_HomoBurden_matrix)

  ##HeterSNP matrix
  colnames(Bin_HeterSNP_matrix)=samples
  Bin_chr_HeterSNP_infor=cbind(all.Bin_Dele_Onecutoff,Bin_HeterSNP_matrix)
  
  TotalAddBurden_PerAccession=apply(Bin_Burden_matrix,2,sum)
  TotalDomBurden_PerAccession=apply(Bin_DomanianceBurden_matrix,2,sum)
  TotalHeterSNP_PerAccession=apply(Bin_HeterSNP_matrix,2,sum)

  AllSNPs=sum(all.Bin_Dele_Onecutoff$polymorphismNum)
  BurdenPerAccession=data.frame(Constraint_Cutoff,samples,TotalAddBurden_PerAccession,TotalDomBurden_PerAccession,TotalHeterSNP_PerAccession,AllSNPs)

  all.BurdenPerAccession=rbind(all.BurdenPerAccession,BurdenPerAccession)

}
write.table(all.Bin_Dele_AllCutOff,paste0(F2_bin_burden_file_prefix,"_AllDiffCutOffs.final_NewCutOff.txt"),quote=F,sep='\t',row.name=F)
write.table(all.BurdenPerAccession,paste0(F2_BurdenEachAccession_file_prefix,"_AllDifferentCutOff.final_NewCutOff.txt"),quote=F,row.name=F,sep='\t')
