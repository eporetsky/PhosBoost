#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="label" 
#SBATCH -N1
#SBATCH -n4
#SBATCH --mem=350GB
#SBATCH -t 48:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date
 
# $1 - Name of the output embedding files saved to the data/tmp/ folder
# $2 - Name of the positive label site list file in the data/sites/ folder

python label.py $1 $2

date
