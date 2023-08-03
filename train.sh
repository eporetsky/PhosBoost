#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="train"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=300GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 3-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# $1 - Name of the embedding file located in the data/embeddings/ folder
# $2 - Residues to train on ("ST" or "Y")
# $3 - Name of the param files to use located in data/params/
# $4 - Name of the output model pkl file that will be saved to the data/models/ folder

date

echo python train.py $1 $2 $3 $4
python train.py $1 $2 $3 $4

date