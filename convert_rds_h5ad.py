#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 12:39:49 2023

@author: balqees
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as pl
import scanpy as sc 
import igraph
import scvelo as scv 
import loompy as lmp 
import anndata
from scipy import io 
from  scipy.sparse import coo_matrix, csr_matrix
import os 

X = io.mmread("/home/balqees/Documents/convert seurat to anndata/matrix.mtx")

adata = anndata.AnnData(X=X.transpose().tocsr())
adata

metadata = pd.read_csv("/home/balqees/Documents/convert seurat to anndata/metadata.csv")

with open("/home/balqees/Documents/convert seurat to anndata/gene_names.csv","r") as f:
   gene_names = f.read().splitlines()
   
adata.obs = metadata 
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

adata

pca = pd.read_csv("/home/balqees/Documents/convert seurat to anndata/pca.csv")
pca.index = adata.obs.index

adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs["UMAP_1"].to_numpy(),adata.obs["UMAP_2"].to_numpy())).T

sc.pl.umap(adata , color= ['seurat_clusters'], frameon =False , save = True)

adata.write("/home/balqees/Documents/convert seurat to anndata/query_integrated.h5ad")

adata
   

   
   

