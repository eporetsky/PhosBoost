#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="gff3"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=360GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 3-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# $1 - Name of GFF3 file
# $2 - Name of prediction file
# $3 - Name of inferred csv file

date

echo python gff3.py $1 $2 $3 
python gff3.py $1 $2 $3

date