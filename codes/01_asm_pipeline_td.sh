################## step1 Generate circular consensus sequences (ccs) from subreads #####################

## pbccs (https://github.com/PacificBiosciences/pbbioconda)
#conda install -c bioconda pbccs
ccs -j 88 C151.subreads.bam > C151.ccs.bam

samtools fasta -@ 88 C151.ccs.bam > C151.ccs.fasta
samtools fastq -@ 88 C151.ccs.bam > C151.ccs.fastq

################## step2-A hifiasm assembly #####################

## hifiasm assembly (https://github.com/chhylp123/hifiasm)

hifiasm -o C151_hifiasm -t 88 C151.ccs.fasta.gz 2> C151_hifiasm.log

awk '/^S/{print ">"$2;print $3}' C151_hifiasm.a_ctg.gfa > C151_hifiasm.a_ctg.fasta 
awk '/^S/{print ">"$2;print $3}' C151_hifiasm.r_utg.gfa > C151_hifiasm.r_utg.fasta 
awk '/^S/{print ">"$2;print $3}' C151_hifiasm.p_utg.gfa > C151_hifiasm.p_utg.fasta 
awk '/^S/{print ">"$2;print $3}' C151_hifiasm.p_ctg.gfa > C151_hifiasm.p_ctg.fasta 

seq_n50.pl C151_hifiasm.a_ctg.fasta > C151_hifiasm.a_ctg.fasta.n50
seq_n50.pl C151_hifiasm.r_utg.fasta > C151_hifiasm.r_utg.fasta.n50
seq_n50.pl C151_hifiasm.p_utg.fasta > C151_hifiasm.p_utg.fasta.n50
seq_n50.pl C151_hifiasm.p_ctg.fasta > C151_hifiasm.p_ctg.fasta.n50

################## step2-B hicanu assembly #####################

## hicanu assembly (https://canu.readthedocs.io/en/latest/)
canu -d C151_hicanu_asm -p C151_hicanu genomeSize=1.6g maxThreads=88 -pacbio-hifi C151.ccs.fasta.gz 1>&C151_hicanu.log 2>&1


################## step3 purge_dups #####################

minimap2 -t 40 -x map-pb C151_hifiasm.p_ctg.fasta C151.ccs.fasta.gz > C151_aln.paf
/lustre/home/tangdie/software/purge_dups/bin/pbcstat C151_aln.paf
/lustre/home/tangdie/software/purge_dups/bin/calcuts PB.stat > cutoffs
/lustre/home/tangdie/software/purge_dups/bin/split_fa C151_hifiasm.p_ctg.fasta > C151_split.fa
minimap2 -xasm5 -DP -t 40 C151_split.fa C151_split.fa > C151_split_self_aln.paf
/lustre/home/tangdie/software/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov C151_split_self_aln.paf > dups.bed
/lustre/home/tangdie/software/purge_dups/bin/get_seqs -p C151_hifiasm dups.bed C151_hifiasm.p_ctg.fasta

# rename contigs and merge results
cat purged.fa hap.fa a_ctg.fa > final.fa
#purged.fa + hap.fa a_ctg.fa + homolog.fa; 周姚分型，但会10倍降低N50。
################## step4 KAT check for every assemblies  #####################

# KAT check (https://github.com/TGAC/KAT)

kat comp -o p_ctg -t 30 -H 10000000000 -I 10000000000 -m 31 -h C151.ccs.fasta.gz C151_hifiasm.p_ctg.fasta
kat comp -o purged -t 30 -H 10000000000 -I 10000000000 -m 31 -h C151.ccs.fasta.gz purged.fa
kat comp -o final -t 30 -H 10000000000 -I 10000000000 -m 31 -h C151.ccs.fasta.gz final.fa

################## step5 Repeat Annotation and Gene Structure Annotation #####################

## EDTA pipeline (https://github.com/oushujun/EDTA#quick-installation-using-conda)

## one round braker pipeline (https://github.com/Gaius-Augustus/BRAKER), two rounds maker pipeline (http://www.yandell-lab.org/software/maker.html)

# /public/agis/huangsanwen_group/baozhigui/Pipeline/Annotation/Snakefile
















