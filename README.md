# PhosBoost

## Introduction

PhosBoost is a machine learning approach that leverages protein language models and gradient boosting trees to predict protein phosphorylation from experimentally derived data. PhosBoost offers improved performance when recall is prioritized while consistently providing more confident probability scores.

The following files are provided for reproducing results:

Phosphoprotein fasta files and site list:
* Arabidopsis.qPTMplants.fasta
* qPTMplants.fasta
* qPTMplants.sites

Trained models:
* Benchmarking: Arabidopsi.qPTMplants.ST.joblib.pkl
* Benchmarking: Arabidopsi.qPTMplants.Y.joblib.pkl
* Genome-wide predictions: qPTMplants.ST.joblib.pkl
* Genome-wide predictions: qPTMplants.Y.joblib.pkl

# Conda environment

Creating the PhosBoost conda environment
```
conda create --prefix /project_folder/phosboost
conda activate /project_folder/phosboost
cd /project_folder/phosboost
conda env update --file environments.yml

# Optional: Specify a project folder to save large pre-trained models
export TORCH_HOME=/project_folder/models/
export TRANSFORMERS_CACHE=/project_folder/models/
```

# Table of Contents

* [1. Data preparation](#preparation)
* [2. Model training](#training)
* [3. Model prediction](#prediction)
* [4. DIAMOND inference](#inference)
* [5. Genome browser](#browser)
* [6. Citation](#citation)

## Methods

<a name="preparation"></a>
## 1. Data preparation

### 1.1. Generate emedding vector data
PhosBoost was trained using the combined S/T/Y residue embedding vector data and average protein embedding vector data, generated using the pretrained ProtT5-XL-U50 pLM. A modified version of the [embed_ProtT5.ipynb](https://colab.research.google.com/drive/1TUj-ayG3WO52n5N50S7KH9vtt6zRkdmj?usp=sharing) from the [ProtTrans](https://github.com/agemagician/ProtTrans/tree/master) GitHub repository was used to obtain the residue and protein embedding data for all S/T/Y residues in an input FASTA file. The main benefit of this modified version is that it doesn't require the extraction of the residue and protein embedding data from the h5 files, and the generated compressed python pickle file is ready for training and prediction by PhosBoost. While it is possible to use GPUs to generate the embedding data, we used multiple CPU nodes (48 threads, 350GB RAM) to process the input FASTA file in multiple batches. Use the large amount of available RAM meant that we were not limited to protein sequences of certain length, while generating the embedding data in reasonable time.

Below is an example for running the modified embedding code using SLURM

```
# Activate the phosboost conda environment
conda activate phosboost

# Recommended: Specifiy the path for the large pre-trained pLM models to be saved at
export TORCH_HOME=/project_folder/models/
export TRANSFORMERS_CACHE=/project_folder/models

########### Generate embeddings in batch mode ###########
# $1 - Name of the fasta file to embed in the data/fasta/ folder
# $2 - Name of the output embedding file(s) saved to the data/tmp/ folder
# $3 - Number of split batches to run as separate slurm jobs 
sbatch embed.batch.sh name.fasta name 30

# Merge the batch output pickle files
# $1 - Name of the output embedding files located in the data/tmp/ folder
sbatch embed.merge.sh name
rm data/tmp/*

########### Generate embeddings using a single slurm node mode ###########
# $1 - Name of the fasta file to embed in the data/fasta/ folder
# $2 - Name of the output embedding file(s) saved to the data/embeddings/ folder
sbatch embed.sh name name 1
```

### 1.2. Label embedding data with positive and negative phosphosites

The embedding vector file generated in the previous step can be used directly for PhosphoSite prediction but it lacks the label (y) column needed for hyperparameter tuning, training, and testing. Using a text file containing all positive phosphosites, use the python script to add a label column ("y"), labeling positive phosphosites as 1s and negative as 0s. 

```
# $1 - Name of embedding pickle file located in data/embeddings/
# $2 - Name of sites file to use as positive labels located in data/sites/
# Note: This overwrites the existing embedding file with one containing a "y" column
sbatch label.sh name qPTMplants.sites
```

### 1.3. Splitting the embedding data into training, validation, and test sets
The labeled embedding data file contains the data for the Ser/Thr/Tyr residues. We use the train_test_split (shuffle=True, stratify=y) function from scikit-learn to split the data of each residue individually while maintaing the same label frequency in each split. The Ser/Thr embeddings are combined. This should generate six sites files saved to the data/sites/ folder, training, validation and test splits files for both the ST and Y models.

```
# $1 - Name of the embedding and positive sites files located in the data/embeddings/ and data/sites/ folders
# $2 - Value of train ratio split (0.6)
# $3 - Value of validation split (0.2)
# $4 - Value of test split (0.2)
# Note: Ratio values must add to 1
```

<a name="training"></a>
## 2. Model training

### 2.1. Hyper-parameter tuning
We used the BayesianOptimization python package to perform Bayesian hyperparameter tuning. Because PhosBoost uses a stacking classifier to predict an S/T and Y models, we performed hyperparameter tuning for each stacked classifier separately. Three hyperparameters for each of the CatBoost classifiers, one trained with balanced weights and one trained with equal weights, were tuned: n_estimators, depth and learning_rate, and the Bayesian hyperparameter tuning ran over 100 iterations, optimizing for the F1 score at each iteration.

```
sbatch tune.sh name ST 100 balanced
sbatch tune.sh name ST 100 equal
sbatch tune.sh name Y 100 balanced
sbatch tune.sh name Y 100 equal
```

The optimal hyperparameters can be obtained from the log/ subfolder, where nnnnn stands for the batch ID, as to not over-write previous results.

### 2.2. Training the PhosBoost model
We provide to options for training a PhosBoost model: 1. Simultaneously trains a PhosBoost model using the combined training and validation sets and then predicts the result using a  Users will have to copy optimal hyperparameters generated in the previous step and manually place them in the PhosBoost.train.predict.py or PhosBoost.train.py  python scripts. 


For training a model, there should be an embedding file and site list files 
```
# $1 - Name of the embedding file located in the data/embeddings/ folder
# $2 - Residues to train on ("ST" or "Y")
# $3 - Name of the param files to use located in data/params/
# $4 - Name of the output model pkl file that will be saved to the data/models/ folder
sbatch train.sh name ST name name
sbatch train.sh name Y name name
```

<a name="prediction"></a>
## 3. Model prediction
If you already have a mode and just want to predict 

Alternatively, if you have a trained model and you just want to predict protein phosphorylation on a given dataset, you can use the following method:
```
# $1 - Name of the embedding file located in /data/embeddings/ (if $2=test provide a sites file in /data/sites)
# $2 - "ST" to predict Ser/Thr and "Y" to predict Tyr
# $3 - "all" to predict all sites in embedding file or "test" if to use the test set located in /data/sites 
# $4 - Name of the model to use located in the data/models
sbatch predict.sh name ST test name
sbatch predict.sh name Y test name
```

<a name="inference"></a>
## 4. DIAMOND inference

To conduct the DIAMOND pairwise sequence alignment analysis step, we need the subject (phosphosite database) fasta file and positive site list and the query (predicted phosphosites) fasta file.

```
# $1 - Name of the fasta file to be queried against the database
# $2 - Name of the database queired, requires both fasta and site files
# $3 - Window size (recommended 15 residues on each side of the phosphosite producing a window size of 31)
sbatch infer.sh name qPTMplants 15
``` 

<a name="browser"></a>
## 5. Genome browser

This part generates a new gff3 file from the genomic gff3 and annotate with the PhosBoost prediction results and the DIAMOND pairwise alignment analysis results. The output gff3 can be uploaded directly to a JBrowse genome browser to visualize the phosphosites on specific genes. The coordinates of each protein were mapped to the genomic coordiantes of the encoding three-basepair codon.

```
# $1 - Name of GFF3 file
# $2 - Name of prediction file
# $3 - Name of inferred csv file
sbatch gff3.sh name name name name_qPTMplants_15 
```

We generated the phosphosite predictions for the maize, wheat, oat and barley that can be uploaded and visualized directly in the genome browser. As a working example, the PhosBoost predictions for the [oat Sang genome](https://graingenes.org/jb/?data=%2Fggds%2Foat-sang&loc=chr1A%3A279100401..279355200&tracks=phospho&highlight=) have been added as a track in the GrainGenes database genome browser.

<a name="citation"></a>
## 6. Citation
The PhosBoost manuscript by Poretsky <i>et al.</i> is currently under review.