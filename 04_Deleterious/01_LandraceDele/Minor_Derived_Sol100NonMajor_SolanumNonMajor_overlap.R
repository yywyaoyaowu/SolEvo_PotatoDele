# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/7/8
allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
#snp_infor_homo_DeleStat_file="/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/09_GenomicPrediction/02_F2_E454Ref_minorDele/A626E454_homoSNP_Constraint_infor_AllChrs.txt"

ConstraintedGERP2=read.table("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol_msa/WholeGenome/Sol_msa_AllChrs_GERP_withDepth.bed_Conserved2")
all.GERP.f=NULL
for (ch in allchrs[1:12]){
  chr=paste0("chr",ch)
  print(chr)
  #input
  msa_freq_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol_msa/SolMsa_100Species_",chr,".fa_BiAllelefreq")
  SolanumMsa_freq_file=paste0("/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol100_Solanum_msa/SolMsa_100Species_substrSolanum_",chr,".fa_BiAllelefreq")

  #####
  msa_freq_all=read.table(msa_freq_file,skip=1)
  names(msa_freq_all)=c("Chr","Pos","Ref","RefIsMajor","Solanum_tuberosum_stenotomumA626.IsMajor","Solanum_tuberosum_stenotomumE454.IsMajor","Allele1","Allele1Count","Allele2","Allele2Count")
  msa_freq_all$MajorAllel=msa_freq_all$Allele2
  msa_freq_all$MajorAllel[which(msa_freq_all$Allele1!="-")]=msa_freq_all$Allele1[which(msa_freq_all$Allele1!="-")]
  msa_freq_all$MajorCount=msa_freq_all$Allele2Count
  msa_freq_all$MajorCount[which(msa_freq_all$Allele1!="-")]=msa_freq_all$Allele1Count[which(msa_freq_all$Allele1!="-")]
  sum(msa_freq_all$MajorAllel!=msa_freq_all$Allele2)
  sum( msa_freq_all$MajorCount!=msa_freq_all$Allele2Count)
  msa_freq=msa_freq_all[(msa_freq_all$MajorAllel=="A"|msa_freq_all$MajorAllel=="T"|msa_freq_all$MajorAllel=="C"|msa_freq_all$MajorAllel=="G"),]

  SolanumMsa_freq_all=read.table(SolanumMsa_freq_file,skip=1)
  names(SolanumMsa_freq_all)=c("Chr","Pos","Ref","RefIsMajor","Solanum_lycopersicum","Solanum_lycopersicum.IsMajor","Solanum_melongena","Solanum_melongena.IsMajor","Allele1","Allele1Count","Allele2","Allele2Count")
  print(sum(SolanumMsa_freq_all$Allele1!="-"))
  SolanumMsa_freq_all$MajorAllel=SolanumMsa_freq_all$Allele2
  SolanumMsa_freq_all$MajorAllel[which(SolanumMsa_freq_all$Allele1!="-")]=SolanumMsa_freq_all$Allele1[which(SolanumMsa_freq_all$Allele1!="-")]
  SolanumMsa_freq_all$MajorCount=SolanumMsa_freq_all$Allele2Count
  SolanumMsa_freq_all$MajorCount[which(SolanumMsa_freq_all$Allele1!="-")]=SolanumMsa_freq_all$Allele1Count[which(SolanumMsa_freq_all$Allele1!="-")]
  sum(SolanumMsa_freq_all$MajorAllel!=SolanumMsa_freq_all$Allele2)
  sum( SolanumMsa_freq_all$MajorCount!=SolanumMsa_freq_all$Allele2Count)
  SolanumMsa_freq=SolanumMsa_freq_all[(SolanumMsa_freq_all$MajorAllel=="A"|SolanumMsa_freq_all$MajorAllel=="T"|SolanumMsa_freq_all$MajorAllel=="C"|SolanumMsa_freq_all$MajorAllel=="G"),]

  GERP=ConstraintedGERP2[ConstraintedGERP2[,1]==chr,]
  indexSol100=match(GERP[,3],msa_freq$Pos)
  GERP$Sol100_MajorAllel=msa_freq$MajorAllel[indexSol100]
  GERP$Sol100_MajorCount=msa_freq$MajorCount[indexSol100]
  GERP$Solanum_melongenaFromSol100=msa_freq$Solanum_melongena[indexSol100]
  GERP$Solanum_lycopersicumFromSol100=msa_freq$Solanum_lycopersicum[indexSol100]

  index_Solanum=match(GERP[,3],SolanumMsa_freq$Pos)
  SolanumMsa_connstrainted_infor=SolanumMsa_freq[index_Solanum,c(3,4,5,7,13,14)]
  names(SolanumMsa_connstrainted_infor)=paste0(names(SolanumMsa_connstrainted_infor),"FromSolanum")
  GERP_constrainted_final=cbind(GERP,SolanumMsa_connstrainted_infor)

  write.table(GERP_constrainted_final,paste0("Sol_msa_ConstraintedGERP2_withDepth_withMajorAllele_",chr,".txt"),quote=F,sep='\t',row.name=F)
  all.GERP.f=rbind(all.GERP.f,GERP_constrainted_final)
}
write.table(all.GERP.f,"Sol_msa_ConstraintedGERP2_withDepth_withMajorAllele_AllChrs.txt",quote=F,sep='\t',row.name=F)

