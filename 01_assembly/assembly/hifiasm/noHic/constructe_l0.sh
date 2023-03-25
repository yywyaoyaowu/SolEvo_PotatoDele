cp -r hifi_assembly_l0_workflow ./${1}_l0

cd ${1}_l0

ln -s ../../../../../ccs/${1}.ccs.fasta.gz rawdata/02_ccs

sed -i "s/sequence/${1}/g" config.yaml

mv qsub.sh ${1:0:1}_${1##*_}.sh

ln -s ../GenomeAssemblyContainer_v0.2 .

cd ..
