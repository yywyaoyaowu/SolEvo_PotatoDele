#!/bin/bash

awk 'BEGIN{OFS="\t"}{if($4>=2&&$6>=20){print "SNP0",$1,$2,$4}}' ../02_Pop367Diploid_Combined/Diploid367DP4Qaulity_PlusNewDP4_Combine_SubstrLandraceMR0.5_AllChrs.bed_Conserved2 > DMV6_AllChrs_GATKFilterd.snp_367MR0.5_GERP_Infor_4_plot.bed
#
rm window_count.sh
for i in {01..12}
do
    bedtools makewindows -g  DM_chr_len.xls -w 10000 | grep chr${i} > DM_chr_len_windows_chr${i}.xls
    grep chr${i} DMV6_AllChrs_GATKFilterd.snp_367MR0.5_GERP_Infor_4_plot.bed > DMV6_AllChrs_GATKFilterd.snp_367MR0.5_GERP_Infor_4_plot.bed_chr${i}.bed
    cat DM_chr_len_windows_chr${i}.xls | while read w
    do
        chr=`echo $w | awk '{print $1}'`
       start=`echo $w | awk '{print $2}'`
        end=`echo $w | awk '{print $3}'`
        echo -e """awk '\$3>="$start"&&\$3<="$end"' DMV6_AllChrs_GATKFilterd.snp_367MR0.5_GERP_Infor_4_plot.bed_chr${i}.bed | awk '{i+=\$4}END{if(NR!=0){print \$2,"$start","$end",i}else{print \""$chr"\","$start","$end",0}}' >> DMV6_AllChrs_GATKFilterd.snp_367MR0.5_GERP_Infor_4_plot.bed_chr${i}_window_count""" >> window_count.sh
    done
done
#
rm window_count.sh_*
split -a 3 -d -l 293 window_count.sh window_count.sh_
#
ls window_count.sh_* | while read i
do
    cat <(echo -e '#!/bin/bash') $i > tmp && mv tmp $i && chmod 755 $i
done
#
rm *_window_count
ls window_count.sh_* | while read i
do
    sbatch -J ${i} -p queue1 --qos=queue1 -N 1 --ntasks-per-node=1 -e %x.err -o %x.out "./${i}"
done


cat DMV6_AllChrs_LandraceMR0.5_GERP_Infor_4_plot.bed_chr*_window_count | sort -k1,1 -k2,2n | awk '{print "SNP0",$1,($3-$2)/2+$2,$4/10000}'  > DMV6_AllChrs_LandraceMR0.5_GERP_Infor_4_plot.bed_Allchr_window_count_4_plot
#
./plot_burden.R DMV6_AllChrs_LandraceMR0.5_GERP_Infor_4_plot.bed_Allchr_window_count_4_plot

rm functional_gene_burden.xls
cat functional_gene.xls | while read i
do
    chr=`echo $i | awk '{print $1}'`
    pos=`echo $i | awk '{print $2}'`
    gene=`echo $i | awk '{print $3}'`
    cat DMV6_AllChrs_LandraceMR0.5_GERP_Infor_4_plot.bed_Allchr_window_count_4_plot | awk '{if($2=="'"$chr"'"&&$3>=int("'"$pos"'")-10000&&$3<=int("'"$pos"'")+10000){print $0"\t""'"$gene"'"}}' >> functional_gene_burden.xls
done


