#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="infer"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=360GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 3-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# $1 - Name of the fasta file to be queried against the database
# $2 - Name of the database queired, requires both fasta and site files
# $3 - Window size (recommended 15 residues on each side of the phosphosite)

date

echo python infer.py $1 $2 $3 
python infer.py $1 $2 $3

date
