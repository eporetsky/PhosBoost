#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="embed"
#SBATCH -N1
#SBATCH -n2
#SBATCH --mem=8GB
#SBATCH -t 1:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date

# $1 - Name of fasta file data/fasta/ folder, excluding ".fasta"
# $2 - Number of split batches to run as separate slurm jobs 

date

for i in $(seq 1 $2)
do
    echo embed.sh $1 $2 $i
	sbatch embed.sh $1 $2 $i
done

date