export PATH=/public/agis/huangsanwen_group/huyong/software/anaconda3/envs/busco-5.0.0/bin:$PATH
export BUSCO_CONFIG_FILE="/public/agis/huangsanwen_group/huyong/software/anaconda3/envs/busco-5.0.0/config/config.ini"
export AUGUSTUS_CONFIG_PATH="/public/agis/huangsanwen_group/huyong/software/anaconda3/envs/busco-5.0.0/config"

python3 /public/agis/huangsanwen_group/huyong/software/anaconda3/envs/busco-5.0.0/bin/generate_plot.py -wd /vol3/agis/huangsanwen_group/wuyaoyao/03-Solanaceae/03-Assembly/Busco_summary
