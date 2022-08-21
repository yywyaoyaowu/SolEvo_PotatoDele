
#ln -s  /home/huyong/SolanaceaeGenomeAnalyze/Call_Variant/final_RH_E8669_C151_C10-20_NGS/20G_filter_merge/*20G_merge_DP4_0.3_0.7.g.vcf ./
accession_list=Landrace_187_Add6Accession.txt
accession367_list=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/02_RareAllele/01_VCF/Diploid.list_367.txt
for i in {01..12}
do
  chr=chr${i}
  vcf=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/02_RareAllele/01_VCF/filter1_no_cluster_missing_filterd.snp_${chr}.header.vcf
  vcf_367=DMV6_${chr}_GATKFilterd.snp_367DiploidJustDP4MR0.2

  cp slurm.sh src/vcf_combine_${chr}.sh
   echo "vcftools --vcf $vcf --minDP 4 --maxDP 100 --minGQ 10 --minQ 30 --keep $accession367_list --max-missing 0.2  --maf 0.001  --recode --recode-INFO-all --out $filter_prefix " >>src/vcf_filterFormat_${chr}_MR0.2.sh

   #add header, sort and index
  echo "cat header ${vcf_367}.recode.vcf  > Diploid367_chr${i}.h.vcf " >> src/vcf_combine_${chr}.sh
  echo "sleep 30s "  >> src/vcf_combine_${chr}.sh
  mkdir /home/wuyaoyao/tmp/bcftools.${i}_3 /home/wuyaoyao/tmp/bcftools.${i}_4
  echo "bcftools sort -T /home/wuyaoyao/tmp/bcftools.${i}_3 Diploid367_chr${i}.h.vcf -o Diploid367_chr${i}.sort.vcf"  >> src/vcf_combine_${chr}.sh
  echo "bcftools sort -T /home/wuyaoyao/tmp/bcftools.${i}_4 chr${i}_20G_merge_DP4_0.3_0.7.g.vcf -o chr${i}_20G_merge_DP4_0.3_0.7.sort.g.vcf" >> src/vcf_combine_${chr}.sh
  echo "bgzip -@ 20  Diploid367_chr${i}.sort.vcf " >> src/vcf_combine_${chr}.sh
  echo "bgzip -@ 20 chr${i}_20G_merge_DP4_0.3_0.7.sort.g.vcf" >> src/vcf_combine_${chr}.sh
  echo "tabix Diploid367_chr${i}.sort.vcf.gz" >> src/vcf_combine_${chr}.sh
  echo "tabix chr${i}_20G_merge_DP4_0.3_0.7.sort.g.vcf.gz" >> src/vcf_combine_${chr}.sh
  ###merge, missing is NA
  echo "bcftools merge Diploid367_chr${i}.sort.vcf.gz chr${i}_20G_merge_DP4_0.3_0.7.sort.g.vcf.gz > Diploid367_PlusNew_Combine_chr${i}_DP4_0.3_0.7_missing.vcf"  >> src/vcf_combine_${chr}.sh
 echo "sleep 30s "  >> src/vcf_combine_${chr}.sh
 ##Substr Landrance
  filter_prefix=Diploid367DP4Qaulity_PlusNewDP4_Combine_chr${i}_SubstrLandraceMR0.5
  echo "vcftools --vcf Diploid367_PlusNew_Combine_chr${i}_DP4_0.3_0.7_missing.vcf --keep $accession_list --max-missing 0.5  --maf 0.001  --recode --recode-INFO-all --out $filter_prefix " >> src/vcf_combine_${chr}.sh
  echo "vcftools --vcf ${filter_prefix}.recode.vcf --freq --out ${filter_prefix}" >> src/vcf_combine_${chr}.sh
  echo "sed 's/:/\t/g' ${filter_prefix}.frq >  ${filter_prefix}_Freq" >> src/vcf_combine_${chr}.sh
  echo "awk '{print \$1, \$2-1,\$2,\$4,\$5,\$6,\$7,\$8}' OFS='\t' ${filter_prefix}_Freq  | sed '1d' >  ${filter_prefix}_bed " >> src/vcf_combine_${chr}.sh
  echo "/home/wuyaoyao/miniconda3/bin/plink --vcf ${filter_prefix}.recode.vcf --make-bed --out $filter_prefix --allow-extra-chr " >> src/vcf_combine_${chr}.sh
  echo "cut -f 1-8 ${filter_prefix}.recode.vcf > ${filter_prefix}.recode.info.vcf " >> src/vcf_combine_${chr}.sh
  echo "java -jar /public/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar DM_v6.1 -ud 5000 ${filter_prefix}.recode.info.vcf > ${filter_prefix}.snpEff  " >> src/vcf_combine_${chr}.sh
  sed -i 's/job_name/JustDP4Combine'${chr}'/g'  src/vcf_combine_${chr}.sh
  sbatch src/vcf_combine_${chr}.sh
done














accession_list=Landrace_187_Add6Accession.txt
for i in {01..12}
do
  chr=chr${i}
vcf=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/05_AddNewAcession_MissRef/01_VCF/chr${i}_merge_missing.vcf
filter_prefix=Pop367_MR5_PlusNewAccession_SubstrLandrace_chr${i}_MR06DP4
cp slurm.sh src/vcf_filterFormat_${chr}_landrace.sh
echo "vcftools --vcf $vcf --minDP 4 --keep $accession_list --max-missing 0.6 --maf 0.001 --recode --recode-INFO-all --out $filter_prefix " >> src/vcf_filterFormat_${chr}_landrace.sh
echo "vcftools --vcf ${filter_prefix}.recode.vcf --freq --out ${filter_prefix}" >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "sed 's/:/\t/g' ${filter_prefix}.frq >  ${filter_prefix}_Freq" >> src/vcf_filterFormat_${chr}_landrace.sh
  echo " awk '{print \$1, \$2-1,\$2,\$4,\$5,\$6,\$7,\$8}' OFS='\t' ${filter_prefix}_Freq  | sed '1d' >  ${filter_prefix}_bed " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "/home/wuyaoyao/miniconda3/bin/plink --vcf ${filter_prefix}.recode.vcf --make-bed --out $filter_prefix --allow-extra-chr " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "cut -f 1-8 ${filter_prefix}.recode.vcf > ${filter_prefix}.recode.info.vcf " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "java -jar /public/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar DM_v6.1 -ud 5000 ${filter_prefix}.recode.info.vcf > ${filter_prefix}.vcf.snpEff  " >> src/vcf_filterFormat_${chr}_landrace.sh
  sed -i 's/job_name/SubstrlandrSNP'${chr}'/g'  src/vcf_filterFormat_${chr}_landrace.sh
  sbatch src/vcf_filterFormat_${chr}_landrace.sh
done


accession_list=Landrace_187_Add6Accession.txt
for i in {01..12}
do
  chr=chr${i}
  vcf=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/05_AddNewAcession_MissRef/01_VCF/chr${i}_merge_missing.vcf
 filter_prefix=Pop367_MR5_PlusNewAccession_SubstrLandrace_chr${i}_MR05
cp slurm.sh src/vcf_filterFormat_${chr}_landrace.sh
echo "vcftools --vcf $vcf --keep $accession_list --max-missing 0.5 --maf 0.001 --recode --recode-INFO-all --out $filter_prefix " >> src/vcf_filterFormat_${chr}_landrace.sh
echo "vcftools --vcf ${filter_prefix}.recode.vcf --freq --out ${filter_prefix}" >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "sed 's/:/\t/g' ${filter_prefix}.frq >  ${filter_prefix}_Freq" >> src/vcf_filterFormat_${chr}_landrace.sh
  echo " awk '{print \$1, \$2-1,\$2,\$4,\$5,\$6,\$7,\$8}' OFS='\t' ${filter_prefix}_Freq  | sed '1d' >  ${filter_prefix}_bed " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "/home/wuyaoyao/miniconda3/bin/plink --vcf ${filter_prefix}.recode.vcf --make-bed --out $filter_prefix --allow-extra-chr " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "cut -f 1-8 ${filter_prefix}.recode.vcf > ${filter_prefix}.recode.info.vcf " >> src/vcf_filterFormat_${chr}_landrace.sh
  echo "java -jar /public/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar DM_v6.1 -ud 5000 ${filter_prefix}.recode.info.vcf > ${filter_prefix}.vcf.snpEff  " >> src/vcf_filterFormat_${chr}_landrace.sh
  sed -i 's/job_name/SubstrlandrSNP'${chr}'/g'  src/vcf_filterFormat_${chr}_landrace.sh
  sbatch src/vcf_filterFormat_${chr}_landrace.sh
done



for i in {01..12}
do
  chr=chr${i}
cp slurm.sh src/vcfTogds_${chr}_landrace2.sh
echo "Rscript ChangeTogbs_Minor.R $i" >> src/vcfTogds_${chr}_landrace2.sh
sed -i 's/job_name/ChangeTogds'${chr}'/g'  src/vcfTogds_${chr}_landrace2.sh
  sbatch  src/vcfTogds_${chr}_landrace2.sh
done


nohup Rscript Burden_AllAccession_BothMinor.R > Burden_AllAccession_BothMinor.R.log &

