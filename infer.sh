#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="infer"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=360GB           	#number of memory
#SBATCH --ntasks=1              #number of nodes
#SBATCH --cpus-per-task=48	#number of cores
#SBATCH -t 3-00:00:00         		#maximum runtime
#SBATCH -o "./log/stdout.%j.%N" 	# standard output
#SBATCH -e "./log/stderr.%j.%N" 	#standard error

# $1 - Name of the fasta file to be queried against the database
# $2 - Name of the database queired, requires both fasta and site files
# $3 - Window size (recommended 15 residues on each side of the phosphosite)

date

echo python infer.py $1 $2 $3 
python infer.py $1 $2 $3

date