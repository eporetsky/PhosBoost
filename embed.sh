#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="embed"   	#name of this job
#SBATCH -N1                  		#number of nodes
#SBATCH -n2                  		#number of cores
#SBATCH --mem=360GB           	#number of memory
#SBATCH -t 5-00:00:00         		#maximum runtime
#SBATCH -o "./log/stdout.%j.%N" 	# standard output
#SBATCH -e "./log/stderr.%j.%N" 	#standard error

date

# $1 - Name of the fasta file to embed in the data/fasta/ folder (.fasta)
# $2 - Number of split batches to run as separate slurm jobs (use "1" to not split) 
# $3 - The specific number of the batch 

date

echo embed.py $1 $2 $3
python embed.py $1 $2 $3

date