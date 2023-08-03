## libraries 
import anndata
import scanpy as sc
import pandas as pd
import os
import scEntropy.scEntropy as scEntropy
import csv
import numpy as np
import anndata
from matplotlib import pyplot as plt 



## functions 

# calculating E, S, C
def calculating_ESC(adata, 
                    Ige_label="I_ge", 
                    I_go_label="I_go"):
    
    E = adata.obs[I_go_label]/adata.obs[Ige_label]
    S = adata.obs[Ige_label] - adata.obs[I_go_label]
    C = E * S
    
    adata.obs['E'] = E 
    adata.obs['S'] = S
    adata.obs['C'] = C 
    
    
    return adata

## start of script 

## finding Ige 

## data imports
df_data = pd.read_csv("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/gene_expression.csv", delimiter=",")

# assign row names
df_data = df_data.set_index("Unnamed: 0") 
df_data = df_data.T

# assigning df_data
gene_expression = df_data 

#running scEntropy
I_ge = scEntropy.scEntropy(gene_expression, option='RCSA')
I_ge_df = pd.DataFrame(I_ge, columns=["I_ge"])

# Saving the array 
I_ge_df.to_csv("Pyige.csv")
#firstarray = np.genfromtxt("PyIge.csv", delimiter=",")

# Sanity Check 
print(I_ge_df) 


## finding Igo



## computing metrics

# loading data
adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")
I_ge = pd.read_csv("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/Pyige.csv")
I_go = pd.read_csv("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/Rigo.csv")

I_ge.shape 
I_go.shape
adata.shape

# inscribing I_go & I_ge into adata.obs 
pd.concat([adata.obs, I_go ,I_ge], axis = 1)



adata = calculating_ESC(adata, Ige_label="I_ge", 
                    I_go_label="I_go")

# saving adata.with new variables 
adata.write_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/Paul15ESC.h5ad")
adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/Paul15ESC.h5ad")

# sanity check 
adata.obs

