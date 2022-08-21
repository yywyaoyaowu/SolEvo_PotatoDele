#get the landrace genotype
sh 00_vcfFilter.sh

01_GetConservation_PercentileCutoff.R


sh 02_Summary_DiffCutOff_ConservedDele.sh

sh 03_CloneGene_DeleNum_BasedRareAllele.sh

sh 04_Loci_DeleNum_BasedRareAllele.sh

Rscript Burden_AllAccession_BothMinor.R

##plot the genome-wide burden
sh DeleBurden_GenomePlot.sh

##burden correlation
Rscript HeterBurden_GeneticBurden_AddRH.C1020.E8669.R

##repulsion phase dele
Rscript RHC151E8669_DeleStat.R

##plot figures
Rscript Figure2d_ConstraintedRatio_GenomeRegion_plot.R

Rscript Figure3c.d.S6_DeleProportion_DiffRegion.R
for i in {01..12}
do
  chr=chr${i}
   cp slurm.sh src/bed2gds_chr${i}.sh
  echo "Rscript ChangeTogbs_Minor.R $i " >> src/bed2gds_chr${i}.sh
  sed -i 's/job_name/JustDP4Combine'${chr}'/g' src/bed2gds_chr${i}.sh
  sbatch src/bed2gds_chr${i}.sh
done




for i in {01..12}
do
  chr=chr${i}
   cp slurm.sh  src/vcf2gds_NumRefchr${i}.sh
  echo "Rscript ChangeTogbs_RefNumber.R $i " >>  src/vcf2gds_NumRefchr${i}.sh
  sed -i 's/job_name/'${chr}'/g'  src/vcf2gds_NumRefchr${i}.sh
  sbatch  src/vcf2gds_NumRefchr${i}.sh
done

for i in {01..12}
do
    cp slurm.sh  src/vcfFreq_chr${i}.sh
  geno_prefix=/home/wuyaoyao/03-Solanaceae/SolEvo_Paper/08_Deleterious/04_WholeGenomeBurden/Pre_PopMinor_Dele/02_RHC151E8669_gvcf.367vcf/01_VCF/C151_RH_E8669_2hap_chr${i}_merge.g
  vcf=${geno_prefix}.vcf
  echo $vcf
  echo "vcftools --vcf $vcf --freq --out ${geno_prefix}" >> src/vcfFreq_chr${i}.sh
 sed -i 's/job_name/vcfFreqChr'${i}'/g'  src/vcfFreq_chr${i}.sh
  sbatch  src/vcfFreq_chr${i}.sh
 done


