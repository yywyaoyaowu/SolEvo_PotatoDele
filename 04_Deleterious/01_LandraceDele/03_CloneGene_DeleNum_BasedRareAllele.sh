GERP_all=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol_msa/WholeGenome/Sol_msa_AllChrs_GERP_withDepth.bed
Poly_Freq_landrace=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/06_Landrace/02_Pop367Diploid_Combined/Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs.bed

CDS=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_CDS.sort.merg.bed
Intron=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_intron.sort.merge.bed
UTR=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_prime_UTR.sort.merge.bed
Promoter_TSS1K=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_Promoter_TSS1K.bed_merge
UpDown5K=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_UpDown5K.sort.merge.bed
deglist_0=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/deglist_0fold_sorted.bed
deglist_2=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/deglist_2fold_sorted.bed
deglist_4=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/deglist_4fold_sorted.bed
GeneticRegion=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_genePlusUpDown5k.sort.merge.bed
FiveUTR=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_five_prime_UTR.sort.merge.bed
ThreeUTR=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_three_prime_UTR.sort.merge.bed
geneRegion=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_gene.sort.merge.bed
UpDown5K=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_UpDown5K.sort.merge.bed

#bedtools intersect -a $Poly_Freq_landrace -b ${GERP_all}_CDS > ${Poly_Freq_landrace}_CDS
#bedtools intersect -a $Poly_Freq_landrace -b $CDS > ${Poly_Freq_landrace}_CDS
#bedtools intersect -a $Poly_Freq_landrace -b $Intron > ${Poly_Freq_landrace}_Intron
#bedtools intersect -a $Poly_Freq_landrace -b $UTR > ${Poly_Freq_landrace}_UTR
#bedtools intersect -a $Poly_Freq_landrace -b $Promoter_TSS1K> ${Poly_Freq_landrace}_Promoter_TSS1K
#bedtools intersect -a $Poly_Freq_landrace -b $UpDown5K > ${Poly_Freq_landrace}_UpDown5K

#geneRegion=/home/wuyaoyao/03-Solanaceae/08_Deleterious/01_GenomicRegion/Sol_msa_4d_v6.1.repr_hc/Solanum_tuberosumDM_v6_repr_hc_gene.sort.merge.bed
#bedtools intersect -a $Poly_Freq_landrace -b $geneRegion > ${Poly_Freq_landrace}_GeneRegion
#bedtools intersect -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_Genetic
#bedtools subtract -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_InterGenetic
##genome sites stat
  Site_num=`cat  ${GERP_all} | wc -l`
  #awk '($4>=2){print}' ${GERP_all} > ${GERP_all}_Conserved2
  CDS_Site_num=`cat   ${GERP_all}_CDS | wc -l`
  UTR_Site_num=`cat   ${GERP_all}_UTR | wc -l`
  deglist_0_Site_num=`cat ${GERP_all}_deglist_0 | wc -l`
  deglist_4_Site_num=`cat   ${GERP_all}_deglist_4 | wc -l`
  Promoter_TSS1K_Site_num=`cat ${GERP_all}_Promoter_TSS1K | wc -l`
  Intron_Site_num=`cat  ${GERP_all}_Intron | wc -l`
  UpDown5K_Site_num=`cat  ${GERP_all}_UpDown5K | wc -l`
  InterGenetic_Site_num=`cat   ${GERP_all}_InterGenetic | wc -l`

#snp sites
  Poly_num=`cat $Poly_Freq_landrace | wc -l`
  CDS_Poly_num=`cat ${Poly_Freq_landrace}_CDS | wc -l`
  deglist_0_Poly_num=`cat ${Poly_Freq_landrace}_deglist_0 | wc -l`
 # deglist_2_Poly_num=`cat ${Poly_Freq_landrace}_deglist_2 | wc -l`
  deglist_4_Poly_num=`cat ${Poly_Freq_landrace}_deglist_4 | wc -l`
  Intron_Poly_num=`cat ${Poly_Freq_landrace}_Intron | wc -l`
  UTR_Poly_num=`cat ${Poly_Freq_landrace}_UTR | wc -l`
 # FiveUTR_Poly_num=`cat ${Poly_Freq_landrace}_FiveUTR | wc -l`
 # ThreeUTR_Poly_num=`cat ${Poly_Freq_landrace}_ThreeUTR | wc -l`
  Promoter_TSS1K_Poly_num=`cat ${Poly_Freq_landrace}_Promoter_TSS1K | wc -l`
  Genetic_Poly_num=`cat ${Poly_Freq_landrace}_Genetic | wc -l`
 # GeneRegion_Poly_num=`cat ${Poly_Freq_landrace}_GeneRegion | wc -l`
  UpDown5K_Poly_num=`cat ${Poly_Freq_landrace}_UpDown5K | wc -l`
  InterGenetic_Poly_num=`cat ${Poly_Freq_landrace}_InterGenetic | wc -l`

##snp different region
#bedtools intersect -a $Poly_Freq_landrace -b $CDS > ${Poly_Freq_landrace}_CDS
#bedtools intersect -a $Poly_Freq_landrace -b ${deglist_0}_just.bed > ${Poly_Freq_landrace}_deglist_0
#bedtools intersect -a $Poly_Freq_landrace -b ${deglist_2}_just.bed > ${Poly_Freq_landrace}_deglist_2
#bedtools intersect -a $Poly_Freq_landrace -b ${deglist_4}_just.bed > ${Poly_Freq_landrace}_deglist_4
#bedtools intersect -a $Poly_Freq_landrace -b $Intron > ${Poly_Freq_landrace}_Intron
#bedtools intersect -a $Poly_Freq_landrace -b $UTR > ${Poly_Freq_landrace}_UTR
#bedtools intersect -a $Poly_Freq_landrace -b  $FiveUTR> ${Poly_Freq_landrace}_FiveUTR
#bedtools intersect -a $Poly_Freq_landrace -b  $ThreeUTR> ${Poly_Freq_landrace}_ThreeUTR
#bedtools intersect -a $Poly_Freq_landrace -b $promoter > ${Poly_Freq_landrace}_Promoter_TSS1K
#bedtools intersect -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_Genetic
#bedtools subtract -a $Poly_Freq_landrace -b $GeneticRegion > ${Poly_Freq_landrace}_InterGenetic
#bedtools intersect -a $Poly_Freq_landrace -b $geneRegion > ${Poly_Freq_landrace}_GeneRegion
#bedtools intersect -a $Poly_Freq_landrace -b $UpDown5K > ${Poly_Freq_landrace}_UpDown5K

##constrainted and dele region
#bedtools intersect -a ${GERP_all}_Conserved2 -b $Poly_Freq_landrace -wb > ${GERP_all}_Conserved2_DeleLandrace
  #awk '($4>=2){print}' ${GERP_all}_CDS > ${GERP_all}_CDS_Conserved2
  # bedtools intersect -a ${GERP_all}_CDS_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_CDS_Conserved2
 #awk '($4>=2){print}' ${GERP_all}_deglist_0 > ${GERP_all}_deglist_0_Conserved2
 #bedtools intersect -a ${GERP_all}_deglist_0_Conserved2 -b $Poly_Freq_landrace -wb>${Poly_Freq_landrace}_deglist_0_Conserved2
 #deglist_2_Site_num=`cat ${GERP_all}_deglist_2 | wc -l`
 #deglist_2_Conserved_num=`awk '($4>=2){print}' ${GERP_all}_deglist_2 | wc -l`
 #awk '($4>=2){print}' ${GERP_all}_deglist_4 > ${GERP_all}_deglist_4_Conserved2
#bedtools intersect -a ${GERP_all}_deglist_4_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_deglist_4_Conserved2
 # awk '($4>=2){print}' ${GERP_all}_InterGenetic > ${GERP_all}_InterGenetic_Conserved2
 #bedtools intersect -a ${GERP_all}_InterGenetic_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_InterGenetic_Conserved2
 #awk '($4>=2){print}' ${GERP_all}_Promoter_TSS1K > ${GERP_all}_Promoter_TSS1K_Conserved2
 #bedtools intersect -a ${GERP_all}_Promoter_TSS1K_Conserved2 -b $Poly_Freq_landrace -wb >${Poly_Freq_landrace}_PromoterTSS1K_Conserved2
 #awk '($4>=2){print}' ${GERP_all}_Intron > ${GERP_all}_Intron_Conserved2
 #bedtools intersect -a ${GERP_all}_Intron_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_Intron_Conserved2
 #awk '($4>=2){print}' ${GERP_all}_UTR > ${GERP_all}_UTR_Conserved2
  #bedtools intersect -a ${GERP_all}_UTR_Conserved2 -b $Poly_Freq_landrace -wb> ${Poly_Freq_landrace}_UTR_Conserved2
  # awk '($4>=2){print}' ${GERP_all}_UpDown5K > ${GERP_all}_UpDown5K_Conserved2
#bedtools intersect -a ${GERP_all}_UpDown5K_Conserved2 -b $Poly_Freq_landrace -wb > ${Poly_Freq_landrace}_UpDown5K_Conserved

cat SelectedCutOff.txt | awk 'NR>1' | while read CutOff CutOff_ID
do
  Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_Conserved2 | wc -l`
  Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Poly_Freq_landrace}_Conserved2 | wc -l`
  Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | wc -l`

  CDS_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_CDS_Conserved2 | wc -l`
  CDS_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_CDS_Conserved2 | wc -l`
  CDS_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}'  | wc -l`

  deglist_0_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_deglist_0_Conserved2 | wc -l`
  deglist_0_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2 | wc -l`
  deglist_0_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

  deglist_4_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_deglist_4_Conserved2 | wc -l`
  deglist_4_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_deglist_4_Conserved2 | wc -l`
  deglist_4_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}'  |wc -l`

  Promoter_TSS1K_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_Promoter_TSS1K_Conserved2 | wc -l`
  Promoter_TSS1K_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2 | wc -l`
  Promoter_TSS1K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}'  | wc -l`

  Intron_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_Intron_Conserved2 | wc -l`
  Intron_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | wc -l`
  Intron_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_Intron_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

  UTR_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${GERP_all}_UTR_Conserved2 | wc -l`
  UTR_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_UTR_Conserved2 | wc -l`
  UTR_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' | wc -l`

  UpDown5K_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_UpDown5K_Conserved2 | wc -l`
  UpDown5K_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_UpDown5K_Conserved | wc -l`
  UpDown5K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_UpDown5K_Conserved| awk '($4>="'${CutOff}'"){print}'  | wc -l`

  InterGenetic_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${GERP_all}_InterGenetic_Conserved2 | wc -l`
  InterGenetic_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_InterGenetic_Conserved2 | wc -l`
  InterGenetic_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'  ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}'  | wc -l`

####
#Stat Burden Index
 Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_Conserved2 | awk '{sum+=$4};END{print sum}'`
 Dele_RareBurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Poly_Freq_landrace}_Conserved2 | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 CDS_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_CDS_Conserved2 | awk '{sum+=$4};END{print sum}'`
 CDS_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_CDS_Conserved2| awk '($4>="'${CutOff}'"){print}'  | awk '{sum+=$4};END{print sum}'`
 
 deglist_0_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_deglist_0_Conserved2 | awk '{sum+=$4};END{print sum}'`
 deglist_0_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Poly_Freq_landrace}_deglist_0_Conserved2| awk '($4>="'${CutOff}'"){print}'  | awk '{sum+=$4};END{print sum}'`

 deglist_4_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_deglist_4_Conserved2 | awk '{sum+=$4};END{print sum}'`
 deglist_4_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_deglist_4_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 Promoter_TSS1K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_PromoterTSS1K_Conserved2 | awk '{sum+=$4};END{print sum}'`
 Promoter_TSS1K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Poly_Freq_landrace}_PromoterTSS1K_Conserved2| awk '($4>="'${CutOff}'"){print}'  | awk '{sum+=$4};END{print sum}'`

 Intron_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Poly_Freq_landrace}_Intron_Conserved2 | awk '{sum+=$4};END{print sum}'`
 Intron_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_Intron_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 UTR_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_UTR_Conserved2 | awk '{sum+=$4};END{print sum}'`
 UTR_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_UTR_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 UpDown5K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_UpDown5K_Conserved | awk '{sum+=$4};END{print sum}'`
 UpDown5K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_UpDown5K_Conserved| awk '($4>="'${CutOff}'"){print}'  | awk '{sum+=$4};END{print sum}'`

 InterGenetic_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Poly_Freq_landrace}_InterGenetic_Conserved2 | awk '{sum+=$4};END{print sum}'`
 InterGenetic_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Poly_Freq_landrace}_InterGenetic_Conserved2| awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

echo WholeGenome WholeGenome All $Site_num $Conserved_num $Dele_num $Dele_Rare05_num $$Poly_num $Dele_BurdenIndex $Dele_RareBurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome CDS $CDS_Site_num $CDS_Conserved_num $CDS_Dele_num $CDS_Dele_Rare05_num $CDS_Poly_num $CDS_Dele_BurdenIndex $CDS_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome deglist_0 $deglist_0_Site_num $deglist_0_Conserved_num $deglist_0_Dele_num $deglist_0_Dele_Rare05_num $deglist_0_Poly_num $deglist_0_Dele_BurdenIndex $deglist_0_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome deglist_4 $deglist_4_Site_num $deglist_4_Conserved_num $deglist_4_Dele_num $deglist_4_Dele_Rare05_num $deglist_4_Poly_num $deglist_4_Dele_BurdenIndex $deglist_4_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome Intron $Intron_Site_num $Intron_Conserved_num  $Intron_Dele_num $Intron_Dele_Rare05_num $Intron_Poly_num $Intron_Dele_BurdenIndex $Intron_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome UTR $UTR_Site_num $UTR_Conserved_num $UTR_Dele_num $UTR_Dele_Rare05_num $UTR_Poly_num $UTR_Dele_BurdenIndex $UTR_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
#echo WholeGenome WholeGenome FivePrimerUTR $FiveUTR_Site_num $FiveUTR_Conserved_num $FiveUTR_Dele_num $FiveUTR_Dele_Rare05_num $FiveUTR_Poly_num $FiveUTR_Dele_BurdenIndex $FiveUTR_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
#echo WholeGenome WholeGenome ThreePrimerUTR $ThreeUTR_Site_num $ThreeUTR_Conserved_num  $ThreeUTR_Dele_num $ThreeUTR_Dele_Rare05_num $ThreeUTR_Poly_num $ThreeUTR_Dele_BurdenIndex $ThreeUTR_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome Promoter_TSS1K $Promoter_TSS1K_Site_num $Promoter_TSS1K_Conserved_num  $Promoter_TSS1K_Dele_num $Promoter_TSS1K_Dele_Rare05_num $Promoter_TSS1K_Poly_num $Promoter_TSS1K_Dele_BurdenIndex $Promoter_TSS1K_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome Genetic $Genetic_Site_num $Genetic_Conserved_num $Genetic_Dele_num $Genetic_Dele_Rare05_num $Genetic_Poly_num $Genetic_Dele_BurdenIndex $Genetic_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome UpDown5K $UpDown5K_Site_num $UpDown5K_Conserved_num  $UpDown5K_Dele_num $UpDown5K_Dele_Rare05_num $UpDown5K_Poly_num $UpDown5K_Dele_BurdenIndex $UpDown5K_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo WholeGenome WholeGenome InterGenetic $InterGenetic_Site_num $InterGenetic_Conserved_num $InterGenetic_Dele_num $InterGenetic_Dele_Rare05_num $InterGenetic_Poly_num $InterGenetic_Dele_BurdenIndex $InterGenetic_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
done


echo CutOff_ID CutOff Loci info Region Site_num Conserved_num Dele_num Dele_Rare05_num SNPNum Dele_BurdenIndex Dele_Rare0.05BurdenIndex  >  Loci_DM_V6_DeleStat_BurdenIndex.txt
#cat SD_DeleLoci_DM_V6.txt | awk 'NR>1' |  while read Loci	SD_population	V6_Chrom	V6_Start V6_End  others
cat Classified_Deleterious_Genes.txt | awk 'NR>1' |  while read Loci	Gene	V6_Chrom	V6_Start V6_End  Promoter_Start Promoter_End others
do
 GERP_chr_prefix=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/01_GERP/Sol_msa/SolMsa_100Species_${V6_Chrom}.fa.full.gerp.chrpos.rates.bed.withDepth.bed
 infor=${V6_Chrom}:${V6_Start}-${V6_End}_${Gene}

 awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix} >  ${Loci}_${infor}_Sol_msa.GERP.txt
 Site_num=`cat  ${Loci}_${infor}_Sol_msa.GERP.txt | wc -l`
  Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' $Poly_Freq_landrace | wc -l`

 awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_CDS > ${Loci}_${infor}_Sol_msa.GERP_CDS.txt
 CDS_Site_num=`cat ${Loci}_${infor}_Sol_msa.GERP_CDS.txt | wc -l`
  CDS_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_CDS | wc -l`

 awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_deglist_0 > ${Loci}_${infor}_Sol_msa.GERP_deglist_0.txt
 deglist_0_Site_num=`cat ${Loci}_${infor}_Sol_msa.GERP_deglist_0.txt | wc -l`
 deglist_0_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_deglist_0 | wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_deglist_4 > ${Loci}_${infor}_Sol_msa.GERP_deglist_4.txt
deglist_4_Site_num=`cat   ${Loci}_${infor}_Sol_msa.GERP_deglist_4.txt | wc -l`
deglist_4_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_deglist_4 | wc -l`

awk '($3>='${Promoter_Start}' && $3<='${Promoter_End}'){print}' ${GERP_chr_prefix} > ${Loci}_${infor}_Sol_msa.GERP_Promoter_TSS1KupStream
Promoter_TSS1K_Site_num=`cat  ${Loci}_${infor}_Sol_msa.GERP_Promoter_TSS1KupStream | wc -l`
Promoter_TSS1K_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${Promoter_Start}' && $3<='${Promoter_End}'){print}' ${Poly_Freq_landrace} | wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_Intron  > ${Loci}_${infor}_Sol_msa.GERP_Intron
Intron_Site_num=`cat ${Loci}_${infor}_Sol_msa.GERP_Intron  | wc -l`
Intron_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_Intron | wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_UTR > ${Loci}_${infor}_Sol_msa.GERP_UTR
UTR_Site_num=`cat   ${Loci}_${infor}_Sol_msa.GERP_UTR | wc -l`
UTR_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_UTR| wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_Genetic  > ${Loci}_${infor}_Sol_msa.GERP_Genetic
Genetic_Site_num=`cat ${Loci}_${infor}_Sol_msa.GERP_Genetic | wc -l`
Genetic_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_Genetic | wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_UpDown5K > ${Loci}_${infor}_Sol_msa.GERP_UpDown5K
UpDown5K_Site_num=`cat   ${Loci}_${infor}_Sol_msa.GERP_UpDown5K | wc -l`
UpDown5K_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_UpDown5K| wc -l`

awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_InterGenetic > ${Loci}_${infor}_Sol_msa.GERP_InterGenetic
InterGenetic_Site_num=`cat ${Loci}_${infor}_Sol_msa.GERP_InterGenetic | wc -l`
InterGenetic_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_InterGenetic | wc -l`

#gene dele
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_Conserved2 > ${Loci}_${infor}_Conserved2_DeleLandrace
 awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_CDS_Conserved2 > ${Loci}_${infor}_CDS_Conserved2_DeleLandrace
 awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}'${Poly_Freq_landrace}_deglist_0_Conserved2 > ${Loci}_${infor}_deglist_0_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_deglist_4_Conserved2 > ${Loci}_${infor}_deglist_4_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${Promoter_Start}' && $3<='${Promoter_End}'){print}' ${Poly_Freq_landrace}_Conserved2 > ${Loci}_${infor}_Promoter_TSS1kupStream_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_Intron_Conserved2 > ${Loci}_${infor}_Intron_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_UTR_Conserved2 > ${Loci}_${infor}_UTR_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_Conserved2 > ${Loci}_${infor}_Genetic_Conserved2_DeleLandrace
awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_UpDown5K_Conserved > ${Loci}_${infor}_UpDown5K_Conserved2_DeleLandrace
 awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_InterGenetic_Conserved2 > ${Loci}_${infor}_InterGenetic_Conserved2_DeleLandrace

cat SelectedCutOff.txt | awk 'NR>1' | while read CutOff CutOff_ID
do

  Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP.txt | wc -l`
   Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Conserved2_DeleLandrace | wc -l`
  Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Loci}_${infor}_Conserved2_DeleLandrace| awk '($4>="'${CutOff}'"){print}' | wc -l`

 CDS_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_CDS.txt | wc -l`
 CDS_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_CDS_Conserved2_DeleLandrace | wc -l`
 CDS_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_CDS_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}'  | wc -l`

  deglist_0_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_deglist_0.txt | wc -l`
 deglist_0_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_deglist_0_Conserved2_DeleLandrace | wc -l`
 deglist_0_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}' ${Loci}_${infor}_deglist_0_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}'  | wc -l`

deglist_4_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_deglist_4.txt | wc -l`
 deglist_4_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_deglist_4_Conserved2_DeleLandrace | wc -l`
 deglist_4_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_deglist_4_Conserved2_DeleLandrace |awk '($4>="'${CutOff}'"){print}'  | wc -l`

Promoter_TSS1K_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_Promoter_TSS1KupStream | wc -l`
 Promoter_TSS1K_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_Promoter_TSS1KupStream_Conserved2_DeleLandrace | wc -l`
 Promoter_TSS1K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Promoter_TSS1KupStream_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | wc -l`

Intron_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP_Intron  | wc -l`
 Intron_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Intron_Conserved2_DeleLandrace | wc -l`
 Intron_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Intron_Conserved2_DeleLandrace |awk '($4>="'${CutOff}'"){print}' | wc -l`

UTR_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP_UTR  | wc -l`
 UTR_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_UTR_Conserved2_DeleLandrace | wc -l`
 UTR_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_UTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | wc -l`

#awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_FivePrimeUTR > ${Loci}_${infor}_Sol_msa.GERP_FiveUTR
#FiveUTR_Site_num=`cat   ${Loci}_${infor}_Sol_msa.GERP_FiveUTR | wc -l`
#FiveUTR_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP_FiveUTR  | wc -l`
 #awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_all}_FiveUTR_Conserved2_DeleLandrace > ${Loci}_${infor}_FiveUTR_Conserved2_DeleLandrace
# FiveUTR_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_FiveUTR_Conserved2_DeleLandrace | wc -l`
# FiveUTR_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_FiveUTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | wc -l`
# FiveUTR_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_FiveUTR | wc -l`

#awk '($3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_chr_prefix}_ThreePrimeUTR > ${Loci}_${infor}_Sol_msa.GERP_ThreeUTR
#ThreeUTR_Site_num=`cat   ${Loci}_${infor}_Sol_msa.GERP_ThreeUTR | wc -l`
#ThreeUTR_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP_ThreeUTR  | wc -l`
#awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${GERP_all}_ThreeUTR_Conserved2_DeleLandrace > ${Loci}_${infor}_ThreeUTR_Conserved2_DeleLandrace
# ThreeUTR_Dele_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_ThreeUTR_Conserved2_DeleLandrace | wc -l`
# ThreeUTR_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_ThreeUTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' |wc -l`
#ThreeUTR_Poly_num=`awk '($1=="'${V6_Chrom}'" && $3>='${V6_Start}' && $3<='${V6_End}'){print}' ${Poly_Freq_landrace}_ThreeUTR | wc -l`

Genetic_Conserved_num=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Sol_msa.GERP_Genetic | wc -l`
 Genetic_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Genetic_Conserved2_DeleLandrace | wc -l`
 Genetic_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Genetic_Conserved2_DeleLandrace| awk '($4>="'${CutOff}'"){print}'  | wc -l`

UpDown5K_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_UpDown5K | wc -l`
 UpDown5K_Dele_num=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_UpDown5K_Conserved2_DeleLandrace | wc -l`
 UpDown5K_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_UpDown5K_Conserved2_DeleLandrace |awk '($4>="'${CutOff}'"){print}' | wc -l`

InterGenetic_Conserved_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_Sol_msa.GERP_InterGenetic | wc -l`
 InterGenetic_Dele_num=`awk '($4>="'${CutOff}'"){print}' ${Loci}_${infor}_InterGenetic_Conserved2_DeleLandrace | wc -l`
 InterGenetic_Dele_Rare05_num=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_InterGenetic_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | wc -l`

##Stat Burden Index
 Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}' ${Loci}_${infor}_Conserved2_DeleLandrace| awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

  CDS_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_CDS_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 CDS_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_CDS_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

deglist_0_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_deglist_0_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 deglist_0_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_deglist_0_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

 deglist_4_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_deglist_4_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 deglist_4_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_deglist_4_Conserved2_DeleLandrace |awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

Promoter_TSS1K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_Promoter_TSS1KupStream_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 Promoter_TSS1K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Promoter_TSS1KupStream_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

Intron_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_Intron_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 Intron_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Intron_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

UTR_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_UTR_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 UTR_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_UTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' | awk '{sum+=$4};END{print sum}'`

FiveUTR_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_FiveUTR_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 FiveUTR_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_FiveUTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

ThreeUTR_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_ThreeUTR_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 ThreeUTR_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_ThreeUTR_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

Genetic_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'  ${Loci}_${infor}_Genetic_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 Genetic_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_Genetic_Conserved2_DeleLandrace| awk '($4>="'${CutOff}'"){print}' |awk '{sum+=$4};END{print sum}'`

UpDown5K_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_UpDown5K_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 UpDown5K_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_UpDown5K_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}'  |awk '{sum+=$4};END{print sum}'`

InterGenetic_Dele_BurdenIndex=`awk '($4>="'${CutOff}'"){print}'   ${Loci}_${infor}_InterGenetic_Conserved2_DeleLandrace | awk '{sum+=$4};END{print sum}'`
 InterGenetic_Dele_Rare05_BurdenIndex=`awk '($12<0.05 || $14<0.05){print}'${Loci}_${infor}_InterGenetic_Conserved2_DeleLandrace | awk '($4>="'${CutOff}'"){print}'  |awk '{sum+=$4};END{print sum}'`

echo $CutOff $CutOff_ID ${Loci} ${infor}  All $Site_num $Conserved_num $Dele_num $Dele_Rare05_num $Poly_num $Dele_BurdenIndex $Dele_RareBurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  CDS $CDS_Site_num $CDS_Conserved_num $CDS_Dele_num $CDS_Dele_Rare05_num $CDS_Poly_num  $CDS_Dele_BurdenIndex $CDS_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  deglist_0 $deglist_0_Site_num $deglist_0_Conserved_num $deglist_0_Dele_num $deglist_0_Dele_Rare05_num $deglist_0_Poly_num $deglist_0_Dele_BurdenIndex $deglist_0_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  deglist_4 $deglist_4_Site_num $deglist_4_Conserved_num $deglist_4_Dele_num $deglist_4_Dele_Rare05_num $deglist_4_Poly_num $deglist_4_Dele_BurdenIndex $deglist_4_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  Intron $Intron_Site_num $Intron_Conserved_num  $Intron_Dele_num $Intron_Dele_Rare05_num $Intron_Poly_num $Intron_Dele_BurdenIndex $Intron_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  UTR $UTR_Site_num $UTR_Conserved_num $UTR_Dele_num $UTR_Dele_Rare05_num $UTR_Poly_num $UTR_Dele_BurdenIndex $UTR_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
#echo ${Loci} ${infor}  FivePrimerUTR $FiveUTR_Site_num $FiveUTR_Conserved_num $FiveUTR_Dele_num $FiveUTR_Dele_Rare05_num $FiveUTR_Poly_num $FiveUTR_Dele_BurdenIndex $FiveUTR_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
#echo ${Loci} ${infor}  ThreePrimerUTR $ThreeUTR_Site_num $ThreeUTR_Conserved_num  $ThreeUTR_Dele_num $ThreeUTR_Dele_Rare05_num $ThreeUTR_Poly_num $ThreeUTR_Dele_BurdenIndex $ThreeUTR_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  Promoter_TSS1K $Promoter_TSS1K_Site_num $Promoter_TSS1K_Conserved_num  $Promoter_TSS1K_Dele_num $Promoter_TSS1K_Dele_Rare05_num $Promoter_TSS1K_Poly_num $Promoter_TSS1K_Dele_BurdenIndex $Promoter_TSS1K_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  Genetic $Genetic_Site_num $Genetic_Conserved_num $Genetic_Dele_num $Genetic_Dele_Rare05_num $Genetic_Poly_num $Genetic_Dele_BurdenIndex $Genetic_Dele_Rare05_BurdenIndex >> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  UpDown5K $UpDown5K_Site_num $UpDown5K_Conserved_num $UpDown5K_Dele_num $UpDown5K_Dele_Rare05_num $UpDown5K_Poly_num $UpDown5K_Dele_BurdenIndex $UpDown5K_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
echo $CutOff $CutOff_ID ${Loci} ${infor}  InterGenetic $InterGenetic_Site_num $InterGenetic_Conserved_num $InterGenetic_Dele_num $InterGenetic_Dele_Rare05_num $InterGenetic_Poly_num $InterGenetic_Dele_BurdenIndex $InterGenetic_Dele_Rare05_BurdenIndex>> Loci_DM_V6_DeleStat_BurdenIndex.txt
done
 done
