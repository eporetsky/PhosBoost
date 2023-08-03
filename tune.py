#!pip install pandas numpy catboost pip install scikit-learn bayesian-optimization

import sys
import pickle
import pandas as pd
import numpy as np
from catboost import CatBoostClassifier
from bayes_opt import BayesianOptimization
from sklearn.metrics import f1_score

name = sys.argv[1]


residues = sys.argv[2]
assert residues=="ST" or residues=="Y", "Need to specifiy either 'ST' or 'Y' residues"

# Name of the training site list located in the data/sites/ folder
with open("data/sites/{}.{}.train.pkl".format(name, residues), 'rb') as f:
    train_sites = pickle.load(f)
# Name of the validation site list located in the data/sites/ folder
with open("data/sites/{}.{}.val.pkl".format(name, residues), 'rb') as f:
    val_sites = pickle.load(f)

# Number of iterations for the BayesianOptimization
n_iter = int(sys.argv[3]) # number of optimization iterations

balanced_input = sys.argv[4]
# Set to optimize using "balanced" or "equal" weights
if balanced_input == "balanced":
    balanced = "Balanced"
elif balanced_input == "equal":
    balanced = "None"
else:
    print("Need to specifiy either 'balanced' or 'equal'")

# For running different tunings
output_name = sys.argv[5]


# Prepare the training data
X_train = pd.read_pickle("data/embeddings/{}.pkl".format(name))
X_train = X_train.loc[train_sites]
y_train = X_train["y"].tolist()
X_train = X_train.drop("y", axis=1)

# Prepare the validationdata
X_val = pd.read_pickle("data/embeddings/{}.pkl".format(name))
X_val = X_val.loc[train_sites]
y_val = X_val["y"].tolist()
X_val = X_val.drop("y", axis=1)

# Based on the following example: https://ai.plainenglish.io/catboost-cross-validated-bayesian-hyperparameter-tuning-91f1804b71dd
def bayesian_opt(n_estimators, depth, learning_rate): 
    cls = CatBoostClassifier(verbose = 0,
                            n_estimators = int(n_estimators),
                            max_depth = int(depth),
                            learning_rate=learning_rate,
                            random_state = 42,
                            grow_policy = "SymmetricTree",
                            use_best_model = True, 
                            auto_class_weights=balanced
                            )

    # Bayesian optimization over the F1 score
    cls.fit(np.array(X_train), np.array(y_train), eval_set = (np.array(X_val), np.array(y_val)))
    return f1_score(np.array(y_val), cls.predict(np.array(X_val)))

# Set the hyperparameters to tune
pbounds = {'n_estimators': (50, 2000),
           'depth': (2, 10),
           'learning_rate': (0.05, 0.5)}

# Set the optimization function
optimizer = BayesianOptimization(
    f = bayesian_opt,
    pbounds = pbounds,
    verbose = 2,
    random_state = 42)

print("Starting hyperparameter tuning:")
print("Embedding name:", name)
print("Residue used:", residues)
print("Number of iterations:", n_iter)
print("Mode:", balanced)
print("Parameters and ranges:", pbounds)

optimizer.maximize(init_points = 2, n_iter = n_iter)

with open("data/params/{}.{}.{}.params".format(output_name, residues, balanced_input), 'wb') as f:
    pickle.dump(optimizer.max["params"], f, protocol=pickle.HIGHEST_PROTOCOL)

# Print the optimization results
print(optimizer.max)