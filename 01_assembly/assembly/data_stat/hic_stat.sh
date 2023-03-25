#!/bin/bash
#SBATCH --partition=queue1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J stat.sh
#SBATCH --qos=queue1
#SBATCH --error=err_%J_stat.sh
#SBATCH --output=out_%J_stat.sh
##################################################################
# @Author: huyong
# @Created Time : Thu May 12 14:07:04 2022

# @File Name: stat.sh
# @Description:
##################################################################

for i in $(cat list)
do
        echo """#!/bin/bash
#SBATCH --partition=smp01
#SBATCH -N 1
#SBATCH -c 1
seq_n50.pl <(fq2fa.pl <(zcat ../${i}/*_R*fq.gz)) > ${i}_stat
""" > ${i}.sh
sbatch ${i}.sh
done
