import sys
import multiprocessing
import pandas as pd
import numpy as np

gff_name = sys.argv[1]
gff_file = "data/gff3/{}.gff3".format(gff_name)
gff_output = "data/gff3/{}.phosboost.gff3".format(gff_name)
# Load the GFF file and get all CDSs
gff = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
gff = gff[gff[2]=="CDS"]
### Can manually add chromosome names for the next steps in the analysis
### gff = gff[gff[0].isin(['1H', '2H', '3H', '4H', '5H', '6H', '7H'])]
# Splits the attribute column and keep the gene ID part
gff[8] = gff[8].str.split(";").str[1].str.replace("Parent=","").str.replace("transcript:","").str.replace("Name=","")

# This is wheat GFF3 specific but can be used to replace other parts of string if added to gene ID
gff[8] = gff[8].str.replace(".v2.1", "")

# Load all prediction results
prediction_name = sys.argv[2]

preds = pd.concat([pd.read_csv("data/preds/{}.ST.all.csv".format(prediction_name)),
                   pd.read_csv("data/preds/{}.Y.all.csv".format(prediction_name))],
                   axis=0)
preds["aa"] = preds["site"].str.split("_").str[-1].str[0]
preds["loc"] = preds["site"].str.split("_").str[-1].str[1:].astype(int)
preds["gene_id"] = preds["site"].str.split("_").str[:-1].str.join("_")

# These columns from the GFF3 file are needed for the final GFF3
cds = gff[[0,1,6,8]].drop_duplicates(8)
cds.columns = ["seqid", "source","strand", "gene_id"]
preds = pd.merge(preds, cds, on="gene_id", how="left")
preds = preds.dropna()

inferred_name = sys.argv[3]
inferred = pd.read_csv("data/inferred/{}".format(inferred_name), sep="\t")
inferred_dict = inferred.groupby('q_site').apply(lambda x: '|'.join(x['h_site'])).to_dict()

# For all predictions (pos and neg) check if conserved and keep either predicted or conserved sites
preds["conserved"] = preds.apply(lambda x: inferred_dict[x["site"]] if x["site"] in inferred_dict.keys() else "None", axis=1)
preds = preds[(preds["conserved"]!="None") | (preds["pred"]==1)]

def aa_to_genomic_coord(row):
    gene_id = row[0]
    loc = row[1]
    # If gene is on negative strand then we want to iterate over CDSs in reverse
    # and reverse the column order to get the appropriate ranges
    cds = gff[gff[8]==gene_id][[0,3,4,6]]
    #print(cds)
    if len(cds)==0:
        return(-2)
    
    # Reoder the column in reverse if the strand is reversed
    strand = cds.iloc[-1,-1]
    if strand == "-":
        cds = cds.iloc[::-1][[0,3,4,6]]
    
    prot_len = 0
    for _, row in cds.iterrows():
        row = list(row)
        chrom = row[0]
        start = row[1]-1
        end = row[2]
        
        cds_len = (end-start)/3 # length in aminoa-acids of the CDS currently iterating over
        prot_len += cds_len     # sum total length of CDSs as they are iterated over
        
        # if the protein falls within the CDS then we just want to extract the genomic coordinates
        # Simplest way is to subtract the distance from the end of the CDS because we have that info
        if loc < prot_len:
            if strand == "+":
                return(int(1+end-(prot_len-loc+1)*3))
            else:
                return(int(start+(prot_len-loc+1)*3))
    return(-1)

pool = multiprocessing.Pool(processes=48)
genomic_coords = pool.map(aa_to_genomic_coord, zip(preds["gene_id"].tolist(), preds["loc"].tolist()))

preds["start"] = genomic_coords
preds = preds.dropna()

# For all predictions (pos and neg) check if conserved and keep either predicted or conserved sites
preds["conserved"] = preds.apply(lambda x: inferred_dict[x["site"]] if x["site"] in inferred_dict.keys() else "None", axis=1)
preds = preds[(preds["conserved"]!="None") | (preds["pred"]==1)]

preds["short"] = preds["site"].str.split("_").str[-1]
preds["type"] = "phosphosite"
preds["end"] = preds.apply(lambda x: x["start"]+2 if x["strand"]=="+" else x["start"], axis=1)
preds["start"] = preds.apply(lambda x: x["start"]-2 if x["strand"]=="-" else x["start"], axis=1)
preds["phase"] = "."
preds["prob"] = preds["prob"].round(3)
# Might want to turn this into a cons dictionary for faster access

# Use one of three annotations do describe each site: Conserved, Predicted or Both
preds["status"] = preds.apply(lambda x: "Conserved" if x["pred"]==0 else ("Predicted" if x["conserved"]=="None" else "Both"), axis=1)
preds["attributes"] = preds.apply(lambda x: f"ID={x['site']};Phosphoprotein={x['gene_id']};probability={x['prob']};label={x['short']};DIAMOND={x['conserved']};status={x['status']}", axis=1)

# This is the score column, which is "." in the original genomic GFF file
preds["prob"] = "."

preds = preds[["seqid", "source", "type", "start", "end", "prob", "strand", "phase", "attributes"]]

# All start positions should be larger than end positions
assert len(preds[preds["start"]>preds["end"]])==0, "Some of the coordinates are flipped"

# Not necessary but copy the original comment lines from the origianl gff3 file
gff_chr_list = []
with open(gff_file, 'r') as f:
    for line in f:
        # Can stop at the first contig if present, manually edit
        if line.startswith('#'):
            gff_chr_list.append(line)
        else:
            break

with open(gff_output, 'w') as f:
    for gff_chr in gff_chr_list:
        f.write(gff_chr)
    preds["seqid"] = preds["seqid"].astype("str")
    preds.to_csv(f,  sep="\t", header=None, index=False)
