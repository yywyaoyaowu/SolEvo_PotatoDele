#! /bin/bash
#SBATCH -p low
#SBATCH -N 1
#SBATCH -c 10
#SBATCH -J seqkit_stats.sh

##################################################################
# @Author: huyong
# @Created Time : Sat 04 Sep 2021 08:51:18 AM CST

# @File Name: seqkit_stats.sh
# @Description:
##################################################################


for i in $(cat fastq_list)
do
        j=${i%.*}
        p=${j##*/}
        seqkit stats -a -j 10 -t dna -o ${p} ${i}
done
