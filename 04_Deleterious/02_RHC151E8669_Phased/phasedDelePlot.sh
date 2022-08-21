cat RH_RH1015_heter_stat_50k_0.02_250k_merge.bed | grep 74385000 | awk '{print $1"_"$2"_"$3,$1,$1,$2,$3}' | while read l
do
chr=`echo $l | awk '{print $3}'`
locus=`echo $l | awk '{print $1}'`
start=`echo $l | awk '{print $4}'`
end=`echo $l | awk '{print $5}'`
length=`expr $end - $start`

#    for g in 0.749 2 3.48 4.13 4.9
        for g in  0.01 # 2
        do
    len=$length
    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_RH_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-2)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_RH_${chr}_GERP${g}_H1.txt
    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_RH_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-1)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_RH_${chr}_GERP${g}_H2.txt

    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_E8669_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-2)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_E8669_${chr}_GERP${g}_H1.txt
    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_E8669_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-1)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_E8669_${chr}_GERP${g}_H2.txt

    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_C151_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-2)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_C151_${chr}_GERP${g}_H1.txt
    cat <(echo -e "chr\tpos") <(awk '$3>=int("'"$start"'")&&$3<=int("'"$end"'")' JustPhasedDele/PhasedHap_JustPhasedDele_C151_${chr}_GERPTop${g}_GERP*txt | awk '$(NF-1)==2' | awk '{print $2,$3-int("'"$start"'")+1}') > PhasedHap_JustPhasedDele_C151_${chr}_GERP${g}_H2.txt

    ./plot.R PhasedHap_JustPhasedDele_RH_${chr}_GERP${g}_H1.txt PhasedHap_JustPhasedDele_RH_${chr}_GERP${g}_H2.txt PhasedHap_JustPhasedDele_C151_${chr}_GERP${g}_H1.txt PhasedHap_JustPhasedDele_C151_${chr}_GERP${g}_H2.txt PhasedHap_JustPhasedDele_E8669_${chr}_GERP${g}_H1.txt PhasedHap_JustPhasedDele_E8669_${chr}_GERP${g}_H2.txt $len E8669_C151_RH_onlyPhased_${locus}_GERP${g}.pdf
    done
done