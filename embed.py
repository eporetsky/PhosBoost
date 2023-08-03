#@title Import dependencies and check whether GPU is available. { display-mode: "form" }
from transformers import T5EncoderModel, T5Tokenizer
import torch
import time
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

#device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
device = 'cpu'
print("Using {}".format(device))

name = sys.argv[1]

fasta_name = "data/fasta/{}.fasta".format(name)

# These batches are for processing the FASTA file in multiple SLURM jobs
# and have nothing to do with the batches in the get_embeddings function
batch_total = int(sys.argv[2])
if batch_total > 1:
    batch_num = int(sys.argv[3])
    emb_name = "data/tmp/{}.batch{}".format(name, str(batch_num))
else:
    batch_num = 1
    emb_name = "data/embeddings/{}.emb.pkl".format(name)
    
print("Running ProtT5 embedding analysis on", fasta_name)

# Based on Rost lab notebook. For more information you can view my notebook at:
# http://localhost:8888/notebooks/GrainGenes/Embeddings/embed_ProtT5.ipynb

#@title Network architecture for secondary structure prediction. { display-mode: "form" }
# Convolutional neural network (two convolutional layers) to predict secondary structure
class ConvNet( torch.nn.Module ):
    def __init__( self ):
        super(ConvNet, self).__init__()
        # This is only called "elmo_feature_extractor" for historic reason
        # CNN weights are trained on ProtT5 embeddings
        self.elmo_feature_extractor = torch.nn.Sequential(
                        torch.nn.Conv2d( 1024, 32, kernel_size=(7,1), padding=(3,0) ), # 7x32
                        torch.nn.ReLU(),
                        torch.nn.Dropout( 0.25 ),
                        )
        n_final_in = 32
        self.dssp3_classifier = torch.nn.Sequential(
                        torch.nn.Conv2d( n_final_in, 3, kernel_size=(7,1), padding=(3,0)) # 7
                        )
        
        self.dssp8_classifier = torch.nn.Sequential(
                        torch.nn.Conv2d( n_final_in, 8, kernel_size=(7,1), padding=(3,0))
                        )
        self.diso_classifier = torch.nn.Sequential(
                        torch.nn.Conv2d( n_final_in, 2, kernel_size=(7,1), padding=(3,0))
                        )

    def forward( self, x):
        # IN: X = (B x L x F); OUT: (B x F x L, 1)
        x = x.permute(0,2,1).unsqueeze(dim=-1) 
        x         = self.elmo_feature_extractor(x) # OUT: (B x 32 x L x 1)
        d3_Yhat   = self.dssp3_classifier( x ).squeeze(dim=-1).permute(0,2,1) # OUT: (B x L x 3)
        d8_Yhat   = self.dssp8_classifier( x ).squeeze(dim=-1).permute(0,2,1) # OUT: (B x L x 8)
        diso_Yhat = self.diso_classifier(  x ).squeeze(dim=-1).permute(0,2,1) # OUT: (B x L x 2)
        return d3_Yhat, d8_Yhat, diso_Yhat
            
#@title Load encoder-part of ProtT5 in half-precision. { display-mode: "form" }
# Load ProtT5 in half-precision (more specifically: the encoder-part of ProtT5-XL-U50) 
def get_T5_model():
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    #model.full()
    model = model.to(device) # move model to GPU
    model = model.eval() # set model to evaluation model
    
    # https://github.com/agemagician/ProtTrans/blob/cc432db8ce0b28754a4383f53b0a387267574f1c/Embedding/PyTorch/Use_Case:NetSurfp_Dataset_Secondary_Structure_Prediction_(3_States)_Feature_Extraction.ipynb
    #if torch.cuda.is_available():
    #    model = model.half()
    
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)

    return model, tokenizer
    
import re
def get_ixs_sites(name, seq, sites):
    """
    Returns 
    """
    ix_list = []
    site_list = []
    for pattern in sites:
        indices = re.finditer(pattern=pattern, string=seq)
        for ix in indices:
            ix = ix.start()
            ix_list.append(ix)
            site_list.append(name+"_"+pattern+str(ix+1))
            
    return(site_list, ix_list)
    
    
#@title Generate embeddings. { display-mode: "form" }
# Generate embeddings via batch-processing
# per_residue indicates that embeddings for each residue in a protein should be returned.
# per_protein indicates that embeddings for a whole protein should be returned (average-pooling)
# max_residues gives the upper limit of residues within one batch
# max_seq_len gives the upper sequences length for applying batch-processing
# max_batch gives the upper number of sequences per batch
def get_embeddings( model, tokenizer, seqs, 
                   max_residues=50000, max_seq_len=20000, max_batch=1000 ):
# original: max_residues=4000, max_seq_len=1000, max_batch=100
    
    site_list = []
    results = np.empty((0, 2048))
    
    # sort sequences according to length (reduces unnecessary padding --> speeds up embedding)
    seq_dict   = sorted( seqs.items(), key=lambda kv: len( seqs[kv[0]] ), reverse=True )
    start = time.time()
    batch = list()
    
    for seq_idx, (pdb_id, seq) in enumerate(seq_dict,1):
                
        seq = seq
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id,seq,seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed 
        n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len 
        if len(batch) >= max_batch or n_res_batch>=max_residues or seq_idx==len(seq_dict) or seq_len>max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()

            # add_special_tokens adds extra token at the end of each sequence
            token_encoding = tokenizer.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)
            
            try:
                with torch.no_grad():
                    # returns: ( batch-size x max_seq_len_in_minibatch x embedding_dim )
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                print("RuntimeError during embedding for {} (L={})".format(pdb_id, seq_len))
                continue

            print("START")
            for batch_idx, identifier in enumerate(pdb_ids): # for each protein in the current mini-batch
                s_len = seq_lens[batch_idx]
                
                # Get the site names and the ixs in the numpy array
                sites, ixs = get_ixs_sites(identifier, seqs[batch_idx].replace(" ", ""), "STY")
                
                # slice off padding --> batch-size x seq_len x embedding_dim
                emb = embedding_repr.last_hidden_state[batch_idx,:s_len]
                                      
                tmp_emb = np.hstack([
                    emb.detach().cpu().numpy().squeeze()[ixs,:],
                    np.tile(emb.mean(dim=0).detach().cpu().numpy().squeeze(), (len(ixs), 1))
                ])

                # Add the generated ProtT5 batch embedding results 
                results = np.vstack([results, tmp_emb])
                
                site_list += sites
                    
    results = pd.DataFrame(results)
    results.index = site_list
    return results

################################################# PhosBoost read_fasta #################################################
# Custom fasta reading function that returns a dictionary of sequences for the specified batch 
def split_fasta(fasta, batch_total, batch_num):
    batch_size = len(fasta) // batch_total
    last = len(fasta) % batch_total
    start = 0
    batch = []

    for i in range(batch_total):
        end = start + batch_size
        
        if i == batch_total-1 and len(fasta)%batch_total:
            batch.append(fasta[start:])
        else:
            batch.append(fasta[start:end])
        
        start = end
    
    batch = batch[batch_num-1]
    seq_dict = {}
    for seq in batch:
        seq_dict[str(seq.id)] = str(seq.seq).replace('U','X').replace('Z','X').replace('O','X').replace('*', '')
        
    return(seq_dict)

################################################# Beginning #################################################

# We replaced the original read_fasta function with a BioPython parser
# If batches > 1 then read_fasta will process the input fasta in batches

fasta_list = list(SeqIO.parse(open(fasta_name), "fasta"))
seqs = split_fasta(fasta_list, batch_total, batch_num)

# Load the encoder part of ProtT5-XL-U50 in half-precision (recommended)
model, tokenizer = get_T5_model()

# Compute embeddings and/or secondary structure predictions
results = get_embeddings(model, tokenizer, seqs)

results.to_pickle(emb_name)