#! /bin/bash
#SBATCH --partition=queue1
#SBATCH --qos=queue1
#SBATCH -N 1
#SBATCH -c 52
#SBATCH -J busco.sh

##################################################################
# @Author: huyong
# @Created Time : Sat Dec  4 21:08:37 2021

# @File Name: busco.sh
# @Description:
##################################################################

for i in $(cat list)
do
        singularity exec /home/huyong/contianer/busco_container \
                busco -m genome \
                        -i ${i} \
                        -o ${i}_busco \
                        -l /home/huyong/software/solanales_odb10 \
                        -c 52 \
                        --offline
done
