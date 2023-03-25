################## step1 Generate circular consensus sequences (ccs) from subreads #####################

## pbccs (https://github.com/PacificBiosciences/pbbioconda)
prefix=$1
#conda install -c bioconda pbccs
ccs -j 88 ${prefix}.subreads.bam > ${prefix}.ccs.bam;
samtools fasta -@ 88 ${prefix}.ccs.bam > ${prefix}.ccs.fasta
samtools fastq -@ 88 ${prefix}.ccs.bam > ${prefix}.ccs.fastq

################## step2-A hifiasm assembly #####################
## hifiasm assembly (https://github.com/chhylp123/hifiasm)
hifiasm -o ${prefix}_hifiasm -t 88 ${prefix}.ccs.fasta.gz -l3 2> ${prefix}_hifiasm.log

#get the contig seq
awk '/^S/{print ">"$2;print $3}' ${prefix}_hifiasm.a_ctg.gfa > ${prefix}_hifiasm.a_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${prefix}_hifiasm.r_utg.gfa > ${prefix}_hifiasm.r_utg.fasta
awk '/^S/{print ">"$2;print $3}' ${prefix}_hifiasm.p_utg.gfa > ${prefix}_hifiasm.p_utg.fasta
awk '/^S/{print ">"$2;print $3}' ${prefix}_hifiasm.p_ctg.gfa > ${prefix}_hifiasm.p_ctg.fasta

seq_n50.pl ${prefix}_hifiasm.a_ctg.fasta > ${prefix}_hifiasm.a_ctg.fasta.n50
seq_n50.pl ${prefix}_hifiasm.r_utg.fasta > ${prefix}_hifiasm.r_utg.fasta.n50
seq_n50.pl ${prefix}_hifiasm.p_utg.fasta > ${prefix}_hifiasm.p_utg.fasta.n50
seq_n50.pl ${prefix}_hifiasm.p_ctg.fasta > ${prefix}_hifiasm.p_ctg.fasta.n50

#step3- KAT check (https://github.com/TGAC/KAT)
kat comp -o p_ctg -t 30 -H 10000000000 -I 10000000000 -m 31 -h ${prefix}.ccs.fasta.gz ${prefix}_hifiasm.p_ctg.fasta


