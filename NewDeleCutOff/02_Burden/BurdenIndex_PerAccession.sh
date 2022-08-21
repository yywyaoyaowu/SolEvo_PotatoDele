
cat SelectedCutOff.txt | awk 'NR>2' | while read cutoff_id cutoff type
do
for i in {01..12}
do
cp slurm.sh src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
echo "Rscript BurdenIndex_OneCutoff_Sol100NonMajorDele.R $i $cutoff" >> src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
sed -i 's/job_name/Burden_chr'${i}'_Cutoff'${cutoff}'/g'  src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
sbatch src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
done
done


cat SelectedCutOff.txt | awk 'NR==2' | while read cutoff_id cutoff type
do
for i in {01..12}
do
cp slurm.sh src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
echo "Rscript BurdenIndex_OneCutoff_filteredSNP.R $i $cutoff" >> src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
sed -i 's/job_name/Burden_chr'${i}'_Cutoff'${cutoff}'/g' src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
sbatch src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
done
done


cat SelectedCutOff.txt | awk 'NR>2' | while read cutoff_id cutoff type
do
for i in {01..12}
do
cp slurm.sh src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
echo "Rscript BurdenIndex_OneCutoff_filteredSNP_TurnOverRareDele.R $i $cutoff" >> src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
sed -i 's/job_name/Burden_chr'${i}'_Cutoff'${cutoff}'/g' src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh
sbatch src/BurdenIndex_chr${i}_Cutoff${cutoff}_AddAcessions.sh

cp slurm.sh src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
echo "Rscript BurdenIndex_OneCutoff_Sol100NonMajorDele_TurnOverRareDele.R $i $cutoff" >> src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
sed -i 's/job_name/Burden_chr'${i}'_Cutoff'${cutoff}'/g'  src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
sbatch src/BurdenIndex_chr${i}_Cutoff${cutoff}.sh
done
done

Rscript catAllChrBurdenIndexTogether.R
for i in {01..12}
do
cp slurm.sh src/Landrance_chr${i}.sh
echo "Rscript AllLandrace_heterInfor_EachChr.R $i " >> src/Landrance_chr${i}.sh
sed -i 's/job_name/Landrace_chr'${i}'/g'  src/Landrance_chr${i}.sh
sbatch src/Landrance_chr${i}.sh
done
