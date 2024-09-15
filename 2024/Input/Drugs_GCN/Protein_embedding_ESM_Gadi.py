import torch
import pandas as pd
import re
import requests
import esm
import os
import gc
from tqdm import tqdm


#in Gadi I created a virtual environment
#module load intel-mkl/2020.3.304
#module load python3/3.9.2
#source /g/data/yr31/pm5363/2024/Env/esm/bin/activate
#python3 ESM.py

# Set environment variable for PyTorch memory allocation
#os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


#This function only convert the Fasta format to only sequence of aminoacids 
def fasta_to_sequence(fasta_string):
    lines = fasta_string.strip().split('\n')
    sequence_lines = [line.strip() for line in lines if not line.startswith('>')]
    sequence = ''.join(sequence_lines)
    return sequence

#for computational purpose I set the batch_size to 1(It crashed for long sequences)
def esm_process(protein, file_path, batch_size=1):
    dictt = {}
    
    # Process in batches
    for i in tqdm(range(0, len(protein), batch_size)):
        batch = protein[i:i+batch_size]
        
        data = [(pid, seq) for pid, seq in batch]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        
        # Move batch_tokens to GPU
        batch_tokens = batch_tokens.to(device)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
        
        # Extract per-residue representations (on GPU)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33]
        
        # Generate per-sequence representations via averaging
        for j, (protein_id, _) in enumerate(batch):
            tokens_len = batch_lens[j]
            sequence_representation = token_representations[j, 1 : tokens_len - 1].mean(0)
            dictt[protein_id] = sequence_representation.cpu().numpy().tolist()
        
        # Clear GPU memory
        del results, token_representations
        torch.cuda.empty_cache()
        gc.collect()
    
    dataf = pd.DataFrame.from_dict(dictt, orient='index')
    dataf.reset_index(inplace=True)
    dataf.columns = ['ID'] + [f'column_{i}' for i in range(1, len(dataf.columns))]
    
    path = os.path.join(file_path, 'esm_protein_650M.csv')
    dataf.to_csv(path, index=False)

# Load the model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
print("Model downloaded")

# Move model to GPU and set to eval mode
model = model.to(device)
model.eval()

batch_converter = alphabet.get_batch_converter()
print("batch_converter Done")

# Read data and prepare protein list, This data is in the format of .csv file, where first column is the uniprot ID("Uniprot ID") and the second column("Fasta") is the fasta sequence of the protein
p = pd.read_csv("/g/data/yr31/pm5363/2024/Codes/Drug_Target_Fasta.csv")

#The path that I want the embedding to be stored
file_path = "/g/data/yr31/pm5363/2024/Codes/"

protein = []
for _, row in p.iterrows():
    current_uniprotid = row["UniProt ID"]    
    fasta_input = row["Fasta"]

    sequence = fasta_to_sequence(fasta_input)
    protein.append((current_uniprotid, sequence))

# Process the data
r = esm_process(protein, file_path)

# Clear GPU memory after processing
del model
torch.cuda.empty_cache()
gc.collect()