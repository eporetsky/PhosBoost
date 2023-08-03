import sys
import pickle
import joblib
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression

# Name of the embedding file(s) located in the data/embeddings/ folder for model training
# Optional: To combine multiple embedding files comma-separate each input file name
# (For example: combining training and validation data after hyperparameter tuning)
name = sys.argv[1]
# Load the embedding data file(s) for the model to be trained on
X_train = pd.read_pickle("data/embeddings/{}.pkl".format(name))

# Load all the site names for the residue in all three splits
residues = sys.argv[2]
# Name of the training site list located in the data/sites/ folder
with open("data/sites/{}.{}.train.pkl".format(name, residues), 'rb') as f:
    train_sites = pickle.load(f)
# Name of the validation site list located in the data/sites/ folder
with open("data/sites/{}.{}.val.pkl".format(name, residues), 'rb') as f:
    val_sites = pickle.load(f)
# Prepare the training data
X_train = X_train.loc[train_sites+val_sites]
y_train = X_train["y"].tolist()
X_train = X_train.drop("y", axis=1)

# Load the hyperparameter dictionaries for the balanced to equal weights
param_name = sys.argv[3]
with open('data/params/{}.{}.balanced.params'.format(param_name, residues), 'rb') as f:
    balanced_param = pickle.load(f)
with open('data/params/{}.{}.equal.params'.format(param_name, residues), 'rb') as f:
    equal_param = pickle.load(f)
# Convert depth and n_estimators parameters to integers
make_int = ["depth", "n_estimators"]
for m in make_int:
    balanced_param[m] = int(balanced_param[m])
    equal_param[m] = int(equal_param[m])

# Name of the output model pkl file that will be saved to the data/models/ folder
model_name = sys.argv[4]

# Create the balanced and equal CatBoost classifier objects
catboost_balanced = CatBoostClassifier(auto_class_weights="Balanced", random_state=42, verbose = 0)
catboost_balanced.set_params(**balanced_param)
catboost_equal = CatBoostClassifier(random_state=42, verbose = 0)
catboost_equal.set_params(**equal_param)

# Make a list of the stacked classifiers
estimators = [
    ('catboost_balanced', catboost_balanced),
    ('catboost_equal', catboost_equal),
]

# Create the stacking classifier object
clf = StackingClassifier(
    estimators=estimators,
    final_estimator=LogisticRegression(class_weight="balanced", max_iter=1000),
    cv=5
)

print("Starting to train the mode:")
print("Name of input embedding file:", ",".join(name))
print("Residues used:", residues)
print("Name of output file:", output_name)
print("Name of params:", param_name)
print("Balanced params:", balanced_param)
print("Equal hyperparameters:", equal_param)

# Train the model using training data
clf.fit(np.array(X_train), np.array(y_train))

# Save the model as a python pickle file to the data/models/ folder
filename = "data/models/{}.{}.model".format(model_name, residues)
_ = joblib.dump(clf, filename, compress=9)