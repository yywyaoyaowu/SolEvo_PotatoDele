cp -r hifiasm_l3_hic_workflow ./${1}_l3

cd ${1}_l3

ln -s ../../../../../ccs/${1}.ccs.fasta.gz rawdata/02_ccs

cd rawdata/03_hic
ln -s ../../../../../hic_raw/${1}/${1}_* .
cd ../..

sed -i "s/sequence/${1}/g" config.yaml

mv qsub.sh ${1:0:1}_${1##*_}.sh

ln -s ../GenomeAssemblyContainer .

cd ..
