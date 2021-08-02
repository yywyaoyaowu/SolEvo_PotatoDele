mkdir N50
mkdir Kmerrm Kmersm
cat Sol.hifi.species.info.txt | while read species genome
do
sh kmer.mask.sh $species $genome
done
