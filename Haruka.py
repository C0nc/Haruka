import os


import sys

import scanpy as sc 
import anndata  as ad 
import squidpy as sq
import pandas as pd 
import numpy as  np


import scenvi
import anndata as ad 
import numpy as np

import torch
import os 

import scvi


import cellcharter as cc
import anndata as ad 

import numpy as np

import scipy.sparse

from pytorch_lightning import seed_everything

import pytorch_lightning as pl

class Haruka():
    def __init__(self, adata, cov_gene_num=32, seed=42):
        self.adata = adata
        self.input_adata = None
        self.cov_gene_num = cov_gene_num
        self.background_indices = None 
        self.target_indices = None
        self.model = None
        self.output_adata = None

    def construct_cov(self):

        g = self.cov_gene_num


        adata = self.adata.copy()

        if scipy.sparse.issparse(adata.X):

            adata.X  = adata.X.toarray()

        covet, covet_sqrt, gene = scenvi.utils.compute_covet(adata, g=g, batch_key='slice_id', batch_size=100000)
        
        #print(adata.X.max())

        if g == -1 :
            g = adata.X.shape[1]

        new_adata = ad.AnnData(X=np.concatenate([adata.X, covet_sqrt.reshape(-1, g*g)], axis=1), obs=adata.obs)

        new_adata.obsm = adata.obsm

        self.input_adata = new_adata
    
    def setup_model(self, latent_dim=256, gene_likelihood='poisson', me_weight=0.3, batch_key='slice_id', wasserstein_penalty=0.2, condition_key='condition', condition_label='disease', seed=42):
        

        
        scvi.external.ContrastiveVI.setup_anndata(self.input_adata, batch_key=batch_key)

        self.model = scvi.external.ContrastiveVI(
            self.input_adata, n_me=self.cov_gene_num * self.cov_gene_num,  n_salient_latent=latent_dim, n_background_latent=latent_dim, use_observed_lib_size=False, wasserstein_penalty=wasserstein_penalty, 
            gene_likelihood=gene_likelihood, me_weight=me_weight
        )

        self.background_indices = np.where(self.input_adata.obs[condition_key] != condition_label)[0]
        self.target_indices = np.where(self.input_adata.obs[condition_key] == condition_label)[0]

    def train(self, max_epochs=200,  batch_size=512, **trainer_kwargs):

        #current_seed = torch.random.initial_seed()
        #print(current_seed) 

        self.model.train(
            background_indices=self.background_indices,
            target_indices=self.target_indices,
            early_stopping=True,
            max_epochs=max_epochs,
            batch_size=batch_size,
            **trainer_kwargs,
        )

        self.output_adata = self.input_adata[self.target_indices].copy()

    def extract_rep(self, n_layer=3):

        self.output_adata.obsm["salient_rep"] = self.model.get_latent_representation(
            self.output_adata, representation_kind="salient"
        ).astype(np.float32)


        self.output_adata.obsm["background_rep"] = self.model.get_latent_representation(
            self.output_adata, representation_kind="background"
        ).astype(np.float32)


        sq.gr.spatial_neighbors(self.output_adata, library_key='slice_id', spatial_key='spatial')\

        cc.gr.remove_long_links(self.output_adata)

        cc.gr.aggregate_neighbors(self.output_adata, n_layers=n_layer, use_rep='salient_rep', out_key='X_cellcharter_salient', sample_key='slice_id')

        cc.gr.aggregate_neighbors(self.output_adata, n_layers=n_layer, use_rep='background_rep', out_key='X_cellcharter_background', sample_key='slice_id')
    
    def cluster(self, salient_cluster_num=4, background_cluster_number=4, seed=42):
        


        model = cc.tl.Cluster(n_clusters=salient_cluster_num, trainer_params=dict(accelerator='gpu', devices=1), random_state=seed)

        pl.seed_everything(seed)

        model.fit(self.output_adata, use_rep='X_cellcharter_salient')

        self.output_adata.obs['cluster_haruka_salient'] = model.predict(self.output_adata, use_rep='X_cellcharter_salient')

        model = cc.tl.Cluster(n_clusters= background_cluster_number, trainer_params=dict(accelerator='gpu', devices=1), random_state=seed)
        
        pl.seed_everything(seed)

        model.fit(self.output_adata, use_rep='X_cellcharter_background')

        self.output_adata.obs['cluster_haruka_background'] = model.predict(self.output_adata, use_rep='X_cellcharter_background')


