import pandas as pd
import glob
import sys

# $1 - Name of the output embedding files saved to the data/tmp/ folder

name = sys.argv[1]

df = pd.DataFrame()
for fl in glob.glob("data/tmp/*"):
    if name in fl:
        df = pd.concat([df, pd.read_pickle(fl)], axis=0)

df.to_pickle("data/embeddings/{}.pkl".format(name))