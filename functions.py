import gzip
import glob
import os
import re

from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import numpy as np
import h5py

import matplotlib.pyplot as plt
import seaborn as sns

def fasta_to_dict(fasta_path):
    """
    Load a fasta file and parse it using biopython
    
    - Input
    fasta_path : str
        relative path of the fasta file

    - Output
    sec_dict : pandas dataframe
        dictionary with gene IDs as keys and sequences as string values
    """
    
    seq_dict = OrderedDict()
    #print(fasta_path)
    record_iterator = SeqIO.parse(fasta_path, "fasta")
    for record in record_iterator:
        seq_dict[str(record.id)] = record
        #seq_dict[str(record.id)] = str(record.seq)
    return(seq_dict)

def fasta_to_dataframe(fasta_path):
    """
    Load a fasta file and parse it using biopython
    
    - Input
    fasta_path : str
        relative path of the fasta file

    - Output
    sec_dict : pandas dataframe
        dictionary with gene IDs as keys and sequences as string values
    """
    
    seq_dict = fasta_to_dict(fasta_path)

    # Convert fasta dict to PDF for a simpler merge  
    seq_df  = pd.Series(seq_dict, index=seq_dict.keys())
    seq_df = pd.DataFrame(seq_df)
    seq_df = seq_df.reset_index()
    seq_df.columns = ["gene_id", "seq"]
    return(seq_df)


def extract_fasta_sites(fasta_path, sites):
    record_iterator = SeqIO.parse(fasta_path, "fasta")
    site_list = []
    for record in record_iterator:
        for pattern in sites:
            indices = re.finditer(pattern=pattern, string=str(record.seq))
            for ix in indices:
                ix = ix.start()
                site_list.append(str(record.id)+"_"+pattern+str(ix+1))
    return(site_list)

def make_phospholingo_prediction_fasta(fasta_path, sites, output_path):
    seq_dict = fasta_to_dict(fasta_path)
    for key, val in seq_dict.items():
        for site in sites:
            seq_dict[key].seq = val.seq.replace(site, site+"#")
    
    with open(output_path, 'w') as handle:
        SeqIO.write(seq_dict.values(), handle, 'fasta')

def sites_table(input_path, input_type, output_path=None):
    """
    Load different types of input data and extract the sites table
    Assumes that the index of the table is of the form geneid_Xnum
    Where "X" is the aa at the "num" position in the protein

    TODO: 1) Replace "_" separator since often used in geneids
          2) Test that the indicated aa matches the protein position
    
    - Input
    input_path : str
        relative path to data table containing site-informative 
        index (geneid_Xnum) and assigned labels in the y-column
    input_type : str
        the type of format the data table is saved in (csv, tsv, pickle)
    output_path : str
        relative path to save the output sites-table

    - Output
    None : written csv file
        dictionary with gene IDs as keys and sequences as string values
    """

    if input_type == "pickle":
        # Load the pre-split pickles
        df = pd.read_pickle(input_path)
        
        # Generate the column containing just the gene IDs
        df["gene_id"] = df.index.str.split("_").str[0]
        
        # \D strips all strings from the AA position and \d strips all digits
        df["pos"] = df.index.str.split("_").str[1].str.replace("\D", "", regex=True)
        df["aa"] = df.index.str.split("_").str[1].str.replace("\d", "", regex=True)

        # Get the columns that will be used to generate the final file
        df = df[["gene_id", "pos", "aa", "y"]]

        if output_path:
            # Save the sites-table as a CSV file
            df.to_csv(output_path, index=False)
        else:
            # Or return the sites dataframe
            return(df)

def label_fasta(sites, fasta):
    return(None)


def label_table(sites, fasta, output_path):
    # seq_df = sites_table()
    # Add the protein sequence and PTM info and save as csv files
    # df = df.merge(pd.DataFrame(seq_df), how="inner", on="gene_id")
    #df.to_csv(output_path, index=False, header=False)
    return(None)