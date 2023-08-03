#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="embed" 
#SBATCH -N1
#SBATCH -n48
#SBATCH --mem=350GB
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date
 
# $1 - Name of the output embedding files saved to the data/tmp/ folder

python embed.merge.py $1

date
