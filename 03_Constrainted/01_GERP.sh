
  msa_fasta=$1
  neutral_tre=$2
  REF=$3
  export PATH=/home/wuyaoyao/software/gerp++KRT/:$PATH
  gerpcol -t $neutral_tre -f ${msa_fasta} -a -e ${REF} -j -z -x .gerp.rates
  source ~/miniconda3/bin/activate python2

python src/Add_N_Sites_GERP.py ${msa_fasta} ${msa_fasta}.gerp.rates $Ref ${msa_fasta}.full.gerp.rates
lines=`cat  ${msa_fasta}.full.gerp.rates | wc -l`
paste -d "\t" <(yes "chr01" | head -n ${lines}) <(seq 1 ${lines}) <(cat ${msa_fasta}.full.gerp.rates) | tr -s ' '| tr ' ' '\t'  > ${msa_fasta}.full.gerp.chrpos.rates
gerpelem -f ${msa_fasta}.full.gerp.chrpos.rates