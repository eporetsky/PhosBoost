import sys
import pandas as pd

name = sys.argv[1]

emb_file = "data/embeddings/{}.pkl".format(name)

pos_sites= "data/sites/"+sys.argv[2]
pos_sites = pd.read_csv(pos_sites, header=None)[0].tolist()

df = pd.read_pickle(emb_file)    
    
df["y"] = df.apply(lambda row: 1 if row.name in pos_sites else 0, axis=1)

# Overwrites the original file without the y column
df.to_pickle(emb_file)