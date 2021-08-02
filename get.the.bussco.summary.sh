dir1=/vol1/agis/huangsanwen_group/huyong/sol_assmebly/assmebly_result/fasta_result/busco_result
dir2=/vol3/agis/huangsanwen_group/wuyaoyao/03-Solanaceae/03-Assembly/busco
cat species_list | while read species
do
  cut -f 1-7 ${dir1}/${species}_hifiasm.p_ctg/${species}_hifiasm.p_ctg_busco/run_solanales_odb10/full_table.tsv > ${dir2}/${species}_hifiasm.p_ctg_busco_full_table.tsv
done
