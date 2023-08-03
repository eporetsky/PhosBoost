import sys
import pickle
import pandas as pd
from sklearn.model_selection import train_test_split

name = sys.argv[1]

X = pd.read_pickle("data/embeddings/{}.pkl".format(name))

# Make a column for making the stratified split on each residue individually 
X["aa"] = X.index.str.split("_").str[-1].str[0].astype(str) + X["y"].astype(str)

# The specific ratio for each split 
train_ratio = float(sys.argv[2])
val_ratio = float(sys.argv[3])
test_ratio = float(sys.argv[4])
assert train_ratio+val_ratio+test_ratio == 1, "Sum of split ratios should be 1"

sites_dict = {}

for res in ["S", "T", "Y"]: 
    X_tmp = X[X["aa"]==res]
    y_tmp = X_tmp["y"].tolist()
    X_tmp = X_tmp.drop("aa", axis=1)
    X_tmp = X_tmp.drop("y", axis=1)

    # Split the complete dataset into training and test (temporary) sets 
    X_train, X_test, y_train, y_test = train_test_split(X_tmp, y_tmp, test_size=1-train_ratio,
                                                        random_state=42, shuffle=True, stratify=y_tmp)

    # Leave the training data and split the test data into test and validation sets
    X_val, X_test, y_val, y_test = train_test_split(X_test, y_test, test_size=test_ratio/(test_ratio+val_ratio), 
                                                    random_state=42, shuffle=True, stratify=y_test)
    
    sites_dict[res] = {"train": X_train.index.tolist(), 
                       "val": X_val.index.tolist(), 
                       "test": X_test.index.tolist()}

# Save the stratified site lists as python pickle files
for split in ["train", "val", "test"]:
    with open("data/sites/{}.ST.{}.pkl".format(name, split), "wb") as f:
        pickle.dump(sites_dict["S"][split]+sites_dict["T"][split], f)

    with open("data/sites/{}.Y.{}.pkl".format(name, split), "wb") as f:
        pickle.dump(sites_dict["Y"][split], f)