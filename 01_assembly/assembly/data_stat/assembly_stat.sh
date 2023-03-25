#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 52
#SBATCH --qos=queue1
#SBATCH --error=err_%J_con.sh
#SBATCH --output=out_%J_con.sh
##################################################################
# @Author: huyong
# @Created Time : Sat Jul  9 22:41:48 2022

# @File Name: con.sh
# @Description:
##################################################################

#for i in $(cat hic_list)
#do
#       pigz -c -p 52 -k -d /home/wuyaoyao/03-Solanaceae/01_Assembly/02_HiC/01_Primary_Contig/Clean_contig/${i}_hifiasm.hic.p_ctg_clean.fasta.gz > ${i}.fa
#done
#
#for i in $(cat nohic_list)
#do
#       cp /home/wuyaoyao/03-Solanaceae/05_PartOne/01_Ref/${i}.fa .
#done
#
#
#ls *fa > list
#for i in $(cat list)
#do
#       seqkit stat -a -G "N" ${i} >> stat_results
#done
#
#for i in $(cat list)
#do
#       seq_n50.pl ${i} > ${i%.fa*}
#done
#for i in $(cat list); do echo -n "${i} "; grep "N90" ${i%.fa*} | cut -f2; done > N90_stat

#pigz -c -p 52 -k -d /home/wuyaoyao/03-Solanaceae/01_Assembly/02_HiC/01_Primary_Contig/Clean_contig/Solanum_ochranthum_hifiasm.l3.hic.hap1.p_ctg_clean.fasta.gz > Solanum_ochranthum.fa
#seqkit stat -a -G "N" Solanum_ochranthum.fa
#seq_n50.pl Solanum_ochranthum.fa > Solanum_ochranthum

#sed 's/list/xaa/g' busco.sh > busco_a.sh
#sed 's/list/xab/g' busco.sh > busco_b.sh
#sed 's/list/xac/g' busco.sh > busco_c.sh
#sed 's/list/xad/g' busco.sh > busco_d.sh

#grep "C:" *.fa_busco/short* | awk -F '10.' '{print $2}'
