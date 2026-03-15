#!/bin/bash
#Python script scheduler
#Change cpus as needed (still must pass same number of threads to python script)

#SBATCH --job-name=pyscheduler
#SBATCH --output=/home/ffalk/links/scratch/scheduler.out
#SBATCH --error=/home/ffalk/links/scratch/scheduler.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --time=12:00:00
#SBATCH --account=(fill)

module load python/3.10
module load trimmomatic
module load fastqc/0.12.1
module load bwa-mem2
module load samtools
module load StdEnv/2023
module load gcc/12.3
module load bcftools
module load bbmap

python3 "$@"
