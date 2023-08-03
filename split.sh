#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="split" 
#SBATCH -N1
#SBATCH -n4
#SBATCH --mem=350GB
#SBATCH -t 48:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date

# $1 - Name of the embedding file located in the data/embeddings/ folder
# $2 - Value of train ratio (0.6)
# $3 - Value of validation ratio (0.2)
# $4 - Value of test ratio (0.2)
# Note: Ratio values must add to 1

python split.py $1 $2 $3 $4

date
