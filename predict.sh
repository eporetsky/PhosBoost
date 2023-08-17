#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="predict"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=360GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# $1 - Name of the embedding file located in the data/embeddings/ folder
# $2 - Residues to train on ("ST" or "Y")
# $3 - Type "test" to predict only the test sites or "all" to predict all sites 
# $4 - Name of the model to use located in the data/models
# $5 - Name of the output prediction csv files saved to data/preds/ folder
# Note: Results will not include true labels if there is no "y" column in emb file 

python predict.py $1 $2 $3 $4 $5
