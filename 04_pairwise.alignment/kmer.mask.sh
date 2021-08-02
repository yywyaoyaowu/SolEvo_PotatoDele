#need kat
#source /home/yw2326/software/miniconda3/etc/profile.d/conda.sh
#conda activate kat

Species_name=$1
QUAIR_NAME=$2
QUERY_PREFIX=$Species_name
gunzip ${QUAIR_NAME}.gz
kat hist -t 20 ${QUAIR_NAME} -m 20 -o ${QUERY_PREFIX}.kat.m20.hist
#plot the k-mer frequency and select the
Rscript /home/yw2326/scripts/kat.break.point_double.slope_baoxing.r ${QUERY_PREFIX}.kat.m20.hist

# kat jellyfish quary genome sequence
kat_jellyfish count -m 20 -s 100M -t 20 -C ${QUAIR_NAME} -o ${QUERY_PREFIX}_k20_count.js
kat_jellyfish dump ${QUERY_PREFIX}_k20_count.js > ${QUERY_PREFIX}_count_dumps.fa
echo $Species_name $f >> genome_kmer_cutoff
f=`cat ${QUERY_PREFIX}.kat.m20.hist_double.slope.break.point.txt`
echo $f
###hard mask
/workdir/yw2326/andropogoneae-conservation/SSW/and_CNS maskGenome -i ${QUAIR_NAME} -o Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa -k  ${QUERY_PREFIX}_count_dumps.fa -f $f  > ${Species_name}_Kmerrm.log &
#soft mask
 /workdir/yw2326/andropogoneae-conservation/SSW/and_CNS maskGenome -i ${QUAIR_NAME} -d -o Kmersm/masked_${QUERY_PREFIX}_k20_f${f}_sm.fa -k  ${QUERY_PREFIX}_count_dumps.fa -f $f > ${Species_name}_Kmersm.log &

perl /home/yw2326/usual.code/fa.length.pl ${QUAIR_NAME} ${QUAIR_NAME}_length.txt
Original_len=`cat ${QUAIR_NAME}_length.txt | awk 'NR>1' | awk '{sum+=$2};END{print sum}'`
perl /home/yw2326/usual.code/N50Stat.pl -i ${QUAIR_NAME} -o N50/${QUAIR_NAME}_N50
N50=`grep 'N50 length' N50/${QUAIR_NAME}_N50`
N90=`grep 'N90 length' N50/${QUAIR_NAME}_N50`
echo $Species_name $N50 $N90

perl /home/yw2326/usual.code/fa.length.pl Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa_length.txt

length=`cat Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa_length.txt | awk 'NR>1' | awk '{sum+=$2};END{print sum}'`
N_bp=`cat Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa_length.txt | awk 'NR>1' | awk '{sum+=$3};END{print sum}'`
n_bp=`cat Kmerrm/masked_${QUERY_PREFIX}_k20_f${f}_rm.fa_length.txt | awk 'NR>1' | awk '{sum+=$4};END{print sum}'`
echo Species_name Genome_size N50 N90 Kmer_freq mask_geenome_length N_bp Kmer.mask.bp > ${Species_name}.txt
echo $Species_name Original_len $N50 $N90 $f $length $N_bp $n_bp >> ${Species_name}.txt
echo $Species_name Original_len $N50 $N90  $f $length $N_bp $n_bp >> Kmerrm.mask.stasstic.txt
