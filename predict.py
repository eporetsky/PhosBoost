import sys
import pickle
import joblib
import pandas as pd
import numpy as np

# Name of the embedding file(s) located in the data/embeddings/ folder for model training
# Optional: To combine multiple embedding files comma-separate each input file name
# (For example: combining training and validation data after hyperparameter tuning)
name = sys.argv[1]
# Load the embedding data file(s) for the model to be trained on
X = pd.read_pickle("data/embeddings/{}.pkl".format(name))
y = True if "y" in X.columns else False

# Load all the site names for the residue in all three splits
residues = sys.argv[2]

test = sys.argv[3]
if test == "test":
    with open("data/sites/{}.{}.test.pkl".format(name, residues), 'rb') as f:
        sites = pickle.load(f)
    # Prepare the training data
    X = X.loc[sites]
    if y:
        y = X["y"].tolist()
        X = X.drop("y", axis=1)
elif test == "all":
    X["res"] = X.index.str.split("_").str[-1].str[0]
    X = X[X["res"].isin(list(residues))]
    X = X.drop("res", axis=1)
    if y:
        y = X["y"].tolist()
        X = X.drop("y", axis=1)
else:
    print("The test variable should be either 'test' or 'all'")

# Load the model pickle file
model_name = sys.argv[4]
model = joblib.load("data/models/{}.{}.model".format(model_name, residues))

# Name of the output results saved to data/preds/
pred_name = sys.argv[5]

preds_proba = model.predict_proba(X)

preds = pd.DataFrame()
preds["site"] = X.index.tolist()
if y:
    preds["label"] = y

preds["prob"]  = preds_proba[:,1]
preds["pred"]  = preds["prob"].apply(lambda x: 1 if x >= 0.5 else 0)
preds.to_csv("data/preds/{}.{}.{}.csv".format(pred_name, residues, test), index=False)

print("Prediction results for model file:", model_name)
print("Prediction results for input file:", name)
print("Name of output file:", pred_name)
print("Residues used:", residues)
