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



#####################
## functions 
def adataX2pandas(adata):
    gene_expression = pd.DataFrame(adata.X)
    gene_expression.columns = adata.var.index 
    gene_expression.index = adata.obs.index
    gene_expression = gene_expression.T
    return(gene_expression)

def Generating_data_for_I_go(adata):
    gene_expression_df = pd.DataFrame()
    gene_expression_df = pd.concat([gene_expression_df, pd.DataFrame(adata.X)], axis=1)
    gene_expression_df.columns = [adata.var.index]
    gene_expression_df.to_csv('data/gene_expression.csv')
    adata.obs.to_csv("data/adataobs.csv")
    adata.var.to_csv("data/adatavar.csv")

# calculating E, S, C
def calculating_ESC(adata, I_ge_label="I_ge", I_go_label="I_go"):
    
    E = adata.obs[I_go_label]/adata.obs[I_ge_label]
    S = adata.obs[I_ge_label] - adata.obs[I_go_label]
    C = E * S
    
    adata.obs['E'] = E 
    adata.obs['S'] = S
    adata.obs['C'] = C 
    
    
    return adata

def Returning_I_go(adata):

    # 1. Generating relavent files
    Generating_data_for_I_go(adata)

    # 2. Running R script
    subprocess.call("R CMD BATCH REntropyIgo.r", shell=True)

    # 3. loading data
    I_go = pd.read_csv("data/Rigo.csv")

    return(I_go)

#####################
## start of script 

## finding Ige 

## loading adata object
adata = sc.read_h5ad("data/TrajectInfer.h5ad")
print("loading adata object")

# Getting I_ge
I_ge = scEntropy.scEntropy(adataX2pandas(adata), option='RCSA')
print("Getting I_ge")

# Getting I_go 
I_go = Returning_I_go(adata)
print("Getting I_go")

# inscribing I_go & I_ge into adata.obs 
adata.obs["I_ge"] = I_ge
adata.obs["I_go"] = I_go.I_go
# adata.obs = pd.concat([adata.obs, I_ge], axis=1)

# adata.write("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/TrajectInfer.h5ad")
 

adata = calculating_ESC(adata, I_ge_label="I_ge", I_go_label="I_go")
print("Calculating E,S,C")

# saving adata.with new variables 
adata.write_h5ad("data/Paul15ESC.h5ad")
# adata = sc.read_h5ad("C:/Users/Dr. Hamad/Downloads/SRSI/SRSI projects/data/Paul15ESC.h5ad")

# sanity check 
# adata.obs

