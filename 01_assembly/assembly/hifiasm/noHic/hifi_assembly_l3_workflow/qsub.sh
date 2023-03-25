#! /bin/bash
#SBATCH --exclusive
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

export PATH=~/software/Anaconda/envs/snakemake/bin/:$PATH

snakemake --stats snake.stats --latency-wait 120 -k -j 64

