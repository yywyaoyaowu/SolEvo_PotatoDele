#!/bin/bash
#SBATCH --job-name=SolSelectedGERPStat
#SBATCH --partition=queue1,low,big,smp01
#SBATCH --qos=queue1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=52
#SBATCH --error=%J-SolSelectedGERPStat.err
#SBATCH --output=%J-SolSelectedGERPStat.out

###input
##cat Diploid367DP4Qaulity_PlusNewDP4_Combine_chr*_SubstrLandraceMR0.5.vcf.snpEff | grep -v '^#' | sed 's#ANN=#\t#g' | awk '{print $1,$2-1,$2,$4,$5,$9}' > Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs_snpEff.bed

GERP_all=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol_msa/WholeGenome/Sol_msa_AllChrs_GERP_withDepth.bed
Poly_Freq_landrace=Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs.bed
snp_infor_bed=Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs_snpEff.bed

#genome region
CDS=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_CDS.sort.merg.bed
deglist_0=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/deglist_0fold_sorted.bed_just.bed
deglist_4=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/deglist_4fold_sorted.bed_just.bed
Intron=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_intron.sort.merge.bed
UTR=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_prime_UTR.sort.merge.bed
FiveUTR=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_five_prime_UTR.sort.merge.bed
promoter=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_Promoter_TSS1K.bed_merge
GeneticRegion=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_genePlusUpDown5k.sort.merge.bed
UpDown5K=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_UpDown5K.sort.merge.bed

########################################
#output
echo CutOff_ID CutOff WholeGenome Region Site_num Conserved_num Dele_num Dele_Rare05_num Dele_Rare01_num SNPNum Dele_BurdenIndex Dele_Rare0.05BurdenIndex Dele_Rare0.01BurdenIndex SNPNumRare05 SNPNumRare01 > DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
##genome sites
 Site_num=`cat  ${GERP_all} | wc -l`
 CDS_Site_num=`cat   ${GERP_all}_CDS | wc -l`
 deglist_0_Site_num=`cat ${GERP_all}_deglist_0 | wc -l`
 deglist_4_Site_num=`cat ${GERP_all}_deglist_4 | wc -l`
 Intron_Site_num=`cat  ${GERP_all}_Intron | wc -l`
 UTR_Site_num=`cat   ${GERP_all}_UTR | wc -l`
 Promoter_TSS1K_Site_num=`cat ${GERP_all}_Promoter_TSS1K | wc -l`
 UpDown5K_Site_num=`cat  ${GERP_all}_UpDown5K | wc -l`
 InterGenetic_Site_num=`cat   ${GERP_all}_InterGenetic | wc -l`

###SNP sites
bedtools intersect -a $Poly_Freq_landrace -b $CDS > ${Poly_Freq_landrace}_CDS
bedtools intersect -a $Poly_Freq_landrace -b ${deglist_0} > ${Poly_Freq_landrace}_deglist_0
bedtools intersect -a $Poly_Freq_landrace -b ${deglist_4} > ${Poly_Freq_landrace}_deglist_4
bedtools intersect -a $Poly_Freq_landrace -b $Intron > ${Poly_Freq_landrace}_Intron
bedtools intersect -a $Poly_Freq_landrace -b $UTR > ${Poly_Freq_landrace}_UTR
bedtools intersect -a $Poly_Freq_landrace -b $promoter > ${Poly_Freq_landrace}_Promoter_TSS1K
bedtools intersect -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_Genetic
bedtools intersect -a $Poly_Freq_landrace -b $UpDown5K > ${Poly_Freq_landrace}_UpDown5K
bedtools subtract -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_InterGenetic
bedtools intersect -a $Poly_Freq_landrace -b $snp_infor_bed -wb > ${Poly_Freq_landrace}_snpEffInfo

Poly_num=`cat $Poly_Freq_landrace | wc -l`
CDS_Poly_num=`cat ${Poly_Freq_landrace}_CDS | wc -l`
deglist_0_Poly_num=`cat ${Poly_Freq_landrace}_deglist_0 | wc -l`
deglist_4_Poly_num=`cat ${Poly_Freq_landrace}_deglist_4 | wc -l`
Intron_Poly_num=`cat ${Poly_Freq_landrace}_Intron | wc -l`
UTR_Poly_num=`cat ${Poly_Freq_landrace}_UTR | wc -l`
Promoter_TSS1K_Poly_num=`cat ${Poly_Freq_landrace}_Promoter_TSS1K | wc -l`

#Genetic_Poly_num=`cat ${Poly_Freq_landrace}_Genetic | wc -l`

UpDown5K_Poly_num=`cat ${Poly_Freq_landrace}_UpDown5K | wc -l`
InterGenetic_Poly_num=`cat ${Poly_Freq_landrace}_InterGenetic | wc -l`

#geneRegion=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_gene.sort.merge.bed
#bedtools intersect -a $Poly_Freq_landrace -b $geneRegion > ${Poly_Freq_landrace}_GeneRegion
#GeneRegion_Poly_num=`cat ${Poly_Freq_landrace}_GeneRegion | wc -l`
NonSynonymous_Site_num=`cat ${Poly_Freq_landrace}_snpEffInfo | egrep 'MODERATE|HIGH' | wc -l`
Synonymous_Site_num=`grep synonymous_variant ${Poly_Freq_landrace}_snpEffInfo | wc -l`
Synonymous_Poly_num=$Synonymous_Site_num
NonSynonymous_Poly_num=$NonSynonymous_Site_num

Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}' $Poly_Freq_landrace | wc -l`
CDS_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_CDS | wc -l`
deglist_0_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_deglist_0 | wc -l`
deglist_4_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_deglist_4 | wc -l`
Intron_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_Intron | wc -l`
UTR_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_UTR | wc -l`
Promoter_TSS1K_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_Promoter_TSS1K | wc -l`
UpDown5K_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_UpDown5K | wc -l`
InterGenetic_Poly_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_InterGenetic | wc -l`
NonSynonymous_Site_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'   ${Poly_Freq_landrace}_snpEffInfo | egrep 'MODERATE|HIGH' | wc -l`
Synonymous_Site_num_Rare05=`awk '($6<0.05 || $8<0.05){print}'  ${Poly_Freq_landrace}_snpEffInfo | grep synonymous_variant | wc -l`
Synonymous_Poly_num_Rare05=$Synonymous_Site_num_Rare05
NonSynonymous_Poly_num_Rare05=$NonSynonymous_Site_num_Rare05

###########
#rare 0.01
Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}' $Poly_Freq_landrace | wc -l`
CDS_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_CDS | wc -l`
deglist_0_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_deglist_0 | wc -l`
deglist_4_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_deglist_4 | wc -l`
Intron_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_Intron | wc -l`
UTR_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_UTR | wc -l`
Promoter_TSS1K_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_Promoter_TSS1K | wc -l`
UpDown5K_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_UpDown5K | wc -l`
InterGenetic_Poly_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_InterGenetic | wc -l`
NonSynonymous_Site_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'   ${Poly_Freq_landrace}_snpEffInfo| egrep 'MODERATE|HIGH' | wc -l`
Synonymous_Site_num_Rare01=`awk '($6<0.01 || $8<0.01){print}'  ${Poly_Freq_landrace}_snpEffInfo | grep synonymous_variant | wc -l`
Synonymous_Poly_num_Rare01=$Synonymous_Site_num_Rare01
NonSynonymous_Poly_num_Rare01=$NonSynonymous_Site_num_Rare01

#####################################
#constrainted/dele sites
bedtools intersect -a ${GERP_all}_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_Conserved2
bedtools intersect -a ${GERP_all}_CDS_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_CDS_Conserved2
bedtools intersect -a ${GERP_all}_deglist_0_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_deglist_0_Conserved2
bedtools intersect -a ${GERP_all}_deglist_4_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_deglist_4_Conserved2
bedtools intersect -a ${GERP_all}_Promoter_TSS1K_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2
  bedtools intersect -a ${GERP_all}_Intron_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_Intron_Conserved2
  bedtools intersect -a ${GERP_all}_UTR_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_UTR_Conserved2
  bedtools intersect -a ${GERP_all}_UpDown5K_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_UpDown5K_Conserved2
  bedtools intersect -a ${GERP_all}_InterGenetic_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_InterGenetic_Conserved2
 bedtools intersect -a ${GERP_all}_Conserved2 -b ${Poly_Freq_landrace}_snpEffInfo -wb > ${Poly_Freq_landrace}_Conserved2_snpEffInfo
bedtools intersect -a ${GERP_all} -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_Conserved

cat SelectedCutOff.txt | awk 'NR>1' | while read CutOff CutOff_ID
do
 Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_Conserved2 | wc -l`
 Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2 | wc -l`
 Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | wc -l`
 Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}' ${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | wc -l`

 CDS_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_CDS_Conserved2 | wc -l`
 CDS_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_CDS_Conserved2 | wc -l`
 CDS_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`
 CDS_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

 deglist_0_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_deglist_0_Conserved2 | wc -l`
 deglist_0_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2 | wc -l`
 deglist_0_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`
 deglist_0_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

 deglist_4_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_deglist_4_Conserved2 | wc -l`
 deglist_4_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_deglist_4_Conserved2 | wc -l`
 deglist_4_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`
deglist_4_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

 Promoter_TSS1K_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_Promoter_TSS1K_Conserved2 | wc -l`
 Promoter_TSS1K_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2 | wc -l`
 Promoter_TSS1K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}'|wc -l`
 Promoter_TSS1K_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}'|wc -l`

 Intron_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_Intron_Conserved2 | wc -l`
 Intron_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | wc -l`
 Intron_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_Intron_Conserved2| awk '($4>="'${CutOff}'"){print}' |wc -l`
 Intron_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_Intron_Conserved2| awk '($4>="'${CutOff}'"){print}' |wc -l`

  UTR_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_UTR_Conserved2 | wc -l`
  UTR_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_UTR_Conserved2 | wc -l`
  UTR_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' |wc -l`
  UTR_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' |wc -l`

 UpDown5K_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_UpDown5K_Conserved2 | wc -l`
 UpDown5K_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_UpDown5K_Conserved2 | wc -l`
 UpDown5K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`
 UpDown5K_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

 InterGenetic_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_InterGenetic_Conserved2 | wc -l`
 InterGenetic_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_InterGenetic_Conserved2 | wc -l`
 InterGenetic_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`
 InterGenetic_Dele_Rare01_num=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

NonSynonymous_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo   | egrep 'MODERATE|HIGH' | wc -l`
NonSynonymous_Dele_Rare05_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo   | awk '($12<0.05 || $14<0.05){print}'  | egrep 'MODERATE|HIGH' | wc -l`
NonSynonymous_Dele_Rare01_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo   | awk '($12<0.01 || $14<0.01){print}'  | egrep 'MODERATE|HIGH' | wc -l`
NonSynonymous_Conserved_num=$NonSynonymous_Dele_num

Synonymous_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | grep 'synonymous_variant' | wc -l`
Synonymous_Dele_Rare05_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.05 || $14<0.05){print}'  | grep 'synonymous_variant' | wc -l`
Synonymous_Dele_Rare01_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.01 || $14<0.01){print}'  | grep 'synonymous_variant' | wc -l`
Synonymous_Conserved_num=$Synonymous_Dele_num
#Stat Burden Index
Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2 | awk '{sum+=$4};END{print sum}'`
Dele_Rare05BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
Dele_Rare01BurdenIndex=`awk '($12<0.01  || $14<0.01){print}' ${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 CDS_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_CDS_Conserved2 | awk '{sum+=$4};END{print sum}'`
 CDS_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
 CDS_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 deglist_0_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2 | awk '{sum+=$4};END{print sum}'`
 deglist_0_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
 deglist_0_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

deglist_4_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_deglist_4_Conserved2 | awk '{sum+=$4};END{print sum}'`
 deglist_4_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}'| awk '{sum+=$4};END{print sum}'`
deglist_4_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}'| awk '{sum+=$4};END{print sum}'`

 Promoter_TSS1K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2 | awk '{sum+=$4};END{print sum}'`
 Promoter_TSS1K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
 Promoter_TSS1K_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 Intron_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | awk '{sum+=$4};END{print sum}'`
 Intron_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`
 Intron_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

 UTR_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_UTR_Conserved2 | awk '{sum+=$4};END{print sum}'`
 UTR_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`
 UTR_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

 UpDown5K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved2 | awk '{sum+=$4};END{print sum}'`
 UpDown5K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
 UpDown5K_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 InterGenetic_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_InterGenetic_Conserved2 | awk '{sum+=$4};END{print sum}'`
 InterGenetic_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`
 InterGenetic_Dele_Rare01_BurdenIndex=`awk '($12<0.01 || $14<0.01){print}'  ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

NonSynonymous_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo | egrep 'MODERATE|HIGH' | awk '{sum+=$4};END{print sum}'`
NonSynonymous_Dele_Rare05_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.05 || $14<0.05){print}'  | egrep 'MODERATE|HIGH' | awk '{sum+=$4};END{print sum}'`
NonSynonymous_Dele_Rare01_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.01 || $14<0.01){print}'  | egrep 'MODERATE|HIGH' | awk '{sum+=$4};END{print sum}'`

Synonymous_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | grep 'synonymous_variant' | awk '{sum+=$4};END{print sum}'`
Synonymous_Dele_Rare05_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.05 || $14<0.05){print}'  | grep 'synonymous_variant' | awk '{sum+=$4};END{print sum}'`
Synonymous_Dele_Rare01_BurdenIndex=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2_snpEffInfo  | awk '($12<0.01 || $14<0.01){print}'  | grep 'synonymous_variant'| awk '{sum+=$4};END{print sum}'`

echo $CutOff_ID ${CutOff} WholeGenome All $Site_num $Conserved_num $Dele_num $Dele_Rare05_num $Dele_Rare01_num $Poly_num $Dele_BurdenIndex $Dele_Rare05BurdenIndex $Dele_Rare01BurdenIndex $Poly_num_Rare05 $Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome CDS $CDS_Site_num $CDS_Conserved_num $CDS_Dele_num $CDS_Dele_Rare05_num $CDS_Dele_Rare01_num $CDS_Poly_num $CDS_Dele_BurdenIndex $CDS_Dele_Rare05_BurdenIndex $CDS_Dele_Rare01_BurdenIndex $CDS_Poly_num_Rare05 $CDS_Poly_num_Rare01  >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome deglist_0 $deglist_0_Site_num $deglist_0_Conserved_num $deglist_0_Dele_num $deglist_0_Dele_Rare05_num $deglist_0_Dele_Rare01_num $deglist_0_Poly_num $deglist_0_Dele_BurdenIndex $deglist_0_Dele_Rare05_BurdenIndex $deglist_0_Dele_Rare01_BurdenIndex $deglist_0_Poly_num_Rare05 $deglist_0_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome deglist_4 $deglist_4_Site_num $deglist_4_Conserved_num $deglist_4_Dele_num $deglist_4_Dele_Rare05_num $deglist_4_Dele_Rare01_num $deglist_4_Poly_num $deglist_4_Dele_BurdenIndex $deglist_4_Dele_Rare05_BurdenIndex $deglist_4_Dele_Rare01_BurdenIndex $deglist_4_Poly_num_Rare05 $deglist_4_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome Intron $Intron_Site_num $Intron_Conserved_num  $Intron_Dele_num $Intron_Dele_Rare05_num $Intron_Dele_Rare01_num $Intron_Poly_num $Intron_Dele_BurdenIndex $Intron_Dele_Rare05_BurdenIndex $Intron_Dele_Rare01_BurdenIndex $Intron_Poly_num_Rare05 $Intron_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome UTR $UTR_Site_num $UTR_Conserved_num $UTR_Dele_num $UTR_Dele_Rare05_num $UTR_Dele_Rare01_num $UTR_Poly_num $UTR_Dele_BurdenIndex $UTR_Dele_Rare05_BurdenIndex $UTR_Dele_Rare01_BurdenIndex $UTR_Poly_num_Rare05 $UTR_Poly_num_Rare01  >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome Promoter_TSS1K $Promoter_TSS1K_Site_num $Promoter_TSS1K_Conserved_num  $Promoter_TSS1K_Dele_num $Promoter_TSS1K_Dele_Rare05_num $Promoter_TSS1K_Dele_Rare01_num $Promoter_TSS1K_Poly_num $Promoter_TSS1K_Dele_BurdenIndex $Promoter_TSS1K_Dele_Rare05_BurdenIndex $Promoter_TSS1K_Dele_Rare01_BurdenIndex $Promoter_TSS1K_Poly_num_Rare05 $Promoter_TSS1K_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome UpDown5K $UpDown5K_Site_num $UpDown5K_Conserved_num  $UpDown5K_Dele_num $UpDown5K_Dele_Rare05_num $UpDown5K_Dele_Rare01_num $UpDown5K_Poly_num $UpDown5K_Dele_BurdenIndex $UpDown5K_Dele_Rare05_BurdenIndex $UpDown5K_Dele_Rare01_BurdenIndex $UpDown5K_Poly_num_Rare05 $UpDown5K_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome InterGenetic $InterGenetic_Site_num $InterGenetic_Conserved_num $InterGenetic_Dele_num $InterGenetic_Dele_Rare05_num $InterGenetic_Dele_Rare01_num $InterGenetic_Poly_num $InterGenetic_Dele_BurdenIndex $InterGenetic_Dele_Rare05_BurdenIndex $InterGenetic_Dele_Rare01_BurdenIndex $InterGenetic_Poly_num_Rare05 $InterGenetic_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome Synonymous $Synonymous_Site_num $Synonymous_Conserved_num $Synonymous_Dele_num $Synonymous_Dele_Rare05_num $Synonymous_Dele_Rare01_num $Synonymous_Poly_num $Synonymous_Dele_BurdenIndex $Synonymous_Dele_Rare05_BurdenIndex $Synonymous_Dele_Rare01_BurdenIndex $Synonymous_Poly_num_Rare05 $Synonymous_Poly_num_Rare01  >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
echo $CutOff_ID ${CutOff} WholeGenome NonSynonymous $NonSynonymous_Site_num $NonSynonymous_Conserved_num $NonSynonymous_Dele_num $NonSynonymous_Dele_Rare05_num $NonSynonymous_Dele_Rare01_num $NonSynonymous_Poly_num $NonSynonymous_Dele_BurdenIndex $NonSynonymous_Dele_Rare05_BurdenIndex $NonSynonymous_Dele_Rare01_BurdenIndex $NonSynonymous_Poly_num_Rare05 $NonSynonymous_Poly_num_Rare01 >> DM_V6_Diploid367AddNew_SubstrLandrace_BurdenStat_Landrace_DiffCutoff.txt
done
