dir=/workdir/yw2326/
cp -r /programs/busco-3.1.0/sample_data/example ./
cp  /programs/busco-3.1.0/sample_data/target.fa ./
cp /programs/busco-3.1.0/config/config.ini ./
cp -r /programs/Augustus-3.3.2/config ./

export AUGUSTUS_CONFIG_PATH=${dir}/config
export BUSCO_CONFIG_FILE="${dir}/config.ini"
export PYTHONPATH=/programs/busco-3.1.0/lib/python3.6/site-packages
export PATH=/programs/busco-3.1.0/scripts:/programs/Augustus-3.3.2/bin:/programs/Augustus-3.3.2/scripts:$PATH

QUARY_DIR=$1
query_genome=$2
Species_name=$3
output=Busco_${Species_name}

echo ${QUARY_DIR}/${Species_name}
#gunzip ${QUARY_DIR}/${QUAIR_NAME}.gz
perl src/N50Stat.pl -i ${QUARY_DIR}/${query_genome} -o ${QUARY_DIR}/${query_genome}_N50
N50=`grep 'N50 length' ${QUARY_DIR}/${query_genome}_N50`
N90=`grep 'N90 length' ${QUARY_DIR}/${query_genome}_N50`
echo $Species_name $N50 $N90

echo "run_BUSCO.py --in ${QUARY_DIR}/${QUAIR_NAME}  --lineage_path ./solanaceae_odb10 --mode genome --out ${output} --cpu 10 --species tomato -f -c 1 > ${Species_name}.log"
#gunzip ${QUARY_DIR}/${QUAIR_NAME}.gz
nohup run_BUSCO.py --in ${QUARY_DIR}/${QUAIR_NAME}  --lineage_path ./solanaceae_odb10 --mode genome --out ${output} --cpu 5 --species tomato -f -c 1 > ${Species_name}.log &
