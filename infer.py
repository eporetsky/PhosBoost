import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from ptmtools.functions import extract_fasta_sites
import re
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

##### Load the files that you want to query against a database
name = sys.argv[1]
seq_dict = {}
record_iterator = SeqIO.parse("data/fasta/{}.fasta".format(name), "fasta")
for record in record_iterator:
    record_id = str(record.id)
    record_seq = str(record.seq)
    seq_dict[record_id] = SeqRecord(seq=Seq(record_seq), id=record_id, name="", description="")

# Get all phosphosites in the fasta file of interest
site_list = extract_fasta_sites("data/fasta/{}.fasta".format(name), "STY")
#####


##### Load phosphorylation database files
# Get the path for the database fasta file 
database_name = sys.argv[2]
database_fasta = "data/fasta/{}.fasta".format(database_name)
tmp_diamond_db = "data/tmp/tmp.window.dmnd.db"
# Prepare the blast output data to parse
os.system("diamond makedb --in {} --db {}".format(database_fasta, tmp_diamond_db))

# Load the list of positive sites in the queried database
pos_sites = pd.read_csv("data/sites/{}.sites".format(database_name), header=None)[0].tolist()
#####

# Set size of window from each side of phosphosite
window = int(sys.argv[3])

# Set name of the output csv file with inferred sites
output_name = 'data/inferred/{}_{}_W{}.tsv'.format(name, database_name, str(window))

##### Create the temporary fasta files with the window around phosphosites
tmp_windows_fasta = open("data/tmp/tmp.window.fasta", 'a')
for site in site_list:
    prot_id = site.split("_")
    aa = prot_id[-1][0]
    pos = int(prot_id[-1][1:])
    prot_id = "_".join(prot_id[:-1])

    pos = pos - 1 # shifting to python indexing
    
    aa_ix = window # the index of the modified site
    left = pos - window
    right = pos + window + 1
    seq = str(seq_dict[prot_id].seq)
    # deal with edge cases and adjust windows and index
    if len(seq) < window * 2 + 1:
        left = 0
        right = len(seq)
        aa_ix = pos
    if len(seq) < right:
        left -= right - len(seq) 
        aa_ix += right - len(seq)
        right = len(seq)
    if left < 0:
        right += abs(left)
        aa_ix -= abs(left)
        left = 0
    # Site name is followed by the coordinate of the phosphosite in the seq
    tmp_windows_fasta.write(">"+site+"|"+str(aa_ix+1)+"\n")
    tmp_windows_fasta.write(seq[left:right]+"\n")
tmp_windows_fasta.close()

tmp_windows_fasta = "data/tmp/tmp.window.fasta"
output_xml = "data/tmp/tmp.dmnd.xml"
os.system("diamond blastp --masking none --ultra-sensitive --max-target-seqs 100 --db {} --out {} --query {} --outfmt 5".format(tmp_diamond_db, output_xml, tmp_windows_fasta))

# Parse the resulting XML file
from Bio.Blast import NCBIXML

blast_records = NCBIXML.parse(open(output_xml))
count = 0
result_df = []
for blast_record in blast_records:
    query = blast_record.query
    q_site = query.split("|")[0] # assuming no other "| in ID"
    q_id = "_".join(q_site.split("_")[:-1])
    ix = int(query.split("|")[-1])
    for alignment in blast_record.alignments:
        h_id = alignment.hit_id
        
        # As long as the peptide sequence is short it shouldn't have multiple hits
        assert len(alignment.hsps) == 1, "Not expecting multiple hits per query"
        
        # Load the first (and only) HSP which is the pairwise alignment result
        hsp = alignment.hsps[0]

        # Sometime the beginning of the query is clipped and not aligned
        hsp_ix = ix - hsp.query_start

        # Find the position of interest in the alignment, accounting for gaps in the query
        # Basically: Move along the pairwise alignment and check for "-" on both the
        #    query and hit alignments. Add to the counter of each as needed. 
        matching_coord = None
        q_gap_count = 0
        h_gap_count = 0
        for i, (q_char, h_char) in enumerate(zip(hsp.query, hsp.sbjct)):
            if q_char == "-":
                q_gap_count += 1
                continue
            if h_char == "-":
                h_gap_count += 1
            # Step when the phosphosite in the query sequence is reached
            if i - q_gap_count == hsp_ix:
                h_site = h_id + "_" + h_char + str(hsp.sbjct_start + i - h_gap_count)
                # Keep only hits that art S/T/Y
                if h_char == "S" or h_char == "T" or h_char == "Y":
                    result_df.append([q_id, h_id, q_site, h_site, q_char, h_char])
                
result_df = pd.DataFrame(result_df)
result_df.columns = ["q_id", "h_id", "q_site", "h_site","q_aa", "h_aa"]
result_df_pos = result_df[result_df["h_site"].isin(pos_sites)]
result_df_pos = result_df_pos[result_df_pos["q_id"]!=result_df_pos["h_id"]]
result_df_pos.to_csv(output_name, sep="\t", index=False)
