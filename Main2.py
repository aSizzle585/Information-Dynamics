## libraries 
# from yaml import scanpy
import yaml 

with open('EntropyV2.yaml', 'r') as file:
    data = yaml.safe_load(file)
    
import subprocess
import anndata
import scanpy as sc
import pandas as pd
import os
import scEntropy.scEntropy as scEntropy
import csv
import numpy as np
from matplotlib import pyplot as plt 

adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")

#####################
## functions 
def adataX2pandas(adata):
    gene_expression = pd.DataFrame(adata.X)
    gene_expression.columns = adata.var.index 
    gene_expression.index = adata.obs.index
    gene_expression = gene_expression.T
    return(gene_expression)

# calculating E, S, C
def calculating_ESC(adata, I_ge_label="I_ge", I_go_label="I_go"):
    
    E = adata.obs[I_go_label]/adata.obs[I_ge_label]
    S = adata.obs[I_ge_label] - adata.obs[I_go_label]
    C = E * S
    
    adata.obs['E'] = E 
    adata.obs['S'] = S
    adata.obs['C'] = C 
    
    
    return adata

#####################
## start of script 

## finding Ige 

## data imports
# df_data = pd.read_csv("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/gene_expression.csv", delimiter=",")

# assign row names
# df_data = df_data.set_index("Unnamed: 0") 
# df_data = df_data.T

# assigning df_data
# gene_expression = adata.X 
# gene_expression
# print(gene_expression.shape)

#running scEntropy
I_ge = scEntropy.scEntropy(adataX2pandas(adata), option='RCSA')
# I_ge_df = pd.DataFrame(I_ge, columns=["I_ge"])

# adata.write("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")
# Saving the array 
# I_ge_df.to_csv("Pyige.csv")
#firstarray = np.genfromtxt("PyIge.csv", delimiter=",")

# Sanity Check 
# print(I_ge_df) 


## finding Igo
subprocess.call("R CMD BATCH REntropyIgo.r", shell=True)


## computing metrics

# loading data
# adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")
# I_ge = pd.read_csv("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/Pyige.csv")
I_go = pd.read_csv("data/Rigo.csv")

# inscribing I_go & I_ge into adata.obs 
adata.obs["I_ge"] = I_ge
adata.obs["I_go"] = I_go.I_go
# adata.obs = pd.concat([adata.obs, I_ge], axis=1)

# adata.write("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")
 

adata = calculating_ESC(adata, I_ge_label="I_ge", I_go_label="I_go")
print(adata.obs)

# saving adata.with new variables 
adata.write_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/Paul15ESC.h5ad")
# adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/Paul15ESC.h5ad")

# sanity check 
# adata.obs

