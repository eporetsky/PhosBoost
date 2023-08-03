#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="tune"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=300GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 5-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# $1 - Name of the experiment
# $2 - Which residues to include ("ST" or "Y")
# $3 - Number of iterations for the BayesianOptimization
# $4 - Set to optimize using "balanced" or "equal" weights
# $5 - Name of the output parameter file



# Run the python hyperparameter tuning function
# Output files containing the best parameters is saved to:
# ./data/hyperparameters/stdout.%j.%N
# If errors occured, they are saved to:
# ./log/stderr.%j.%N

date

echo python tune.py $1 $2 $3 $4 $5
python tune.py $1 $2 $3 $4 $5

date