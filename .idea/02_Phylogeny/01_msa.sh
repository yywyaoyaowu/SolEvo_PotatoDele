source /home/wuyaoyao/software/cactus-bin-v2.0.3/venv/bin/activate
export PATH=/home/wuyaoyao/software/cactus-bin-v2.0.3/bin/:$PATH
export PYTHONPATH=/home/wuyaoyao/software/cactus-bin-v2.0.3/lib:$PYTHONPATH

cactus jobStore_Cactus Cactus_Species.txt Sol_Cactus.hal  --binariesMode local  --stats
