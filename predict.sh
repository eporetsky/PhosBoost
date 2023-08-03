#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="predict"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=360GB           	#number of memory
#SBATCH --ntasks=1              #number of nodes
#SBATCH --cpus-per-task=48	#number of cores
#SBATCH -t 2-00:00:00         		#maximum runtime
#SBATCH -o "./log/stdout.%j.%N" 	# standard output
#SBATCH -e "./log/stderr.%j.%N" 	#standard error

# $1 - Name of the embedding file located in the data/embeddings/ folder
# $2 - Residues to train on ("ST" or "Y")
# $3 - Type "test" to predict only the test sites or "all" to predict all sites 
# $4 - Name of the output prediction csv files saved to data/preds/ folder
# Note: Results will not include true labels if there is no "y" column in emb file 

python predict.py $1 $2 $3 $4
