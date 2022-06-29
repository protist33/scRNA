#!/usr/bin/env python
# coding: utf-8

# In[94]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[26]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


# In[27]:


dat = sc.read_10x_mtx(r"C:\Users\Виталик\Documents\RNAproject\10x genomics", var_names="gene_symbols")


# In[28]:


dat.var_names_make_unique()  


# In[29]:


sc.pl.highest_expr_genes(dat, n_top=20)


# In[30]:


sc.pp.filter_cells(dat, min_genes=200)


# In[31]:


sc.pp.filter_genes(dat, min_cells=3)


# In[32]:


dat.var['mt'] = dat.var_names.str.startswith('mt-') #annotate the group of mitochondrial genes as 'mt'


# In[33]:


dat


# In[34]:


sc.pp.calculate_qc_metrics(dat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[35]:


sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[36]:


dat.obs


# In[37]:


#Remove cells that have too many mitochondrial genes expressed or too many total counts:
sc.pl.scatter(dat, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(dat, x='total_counts', y='n_genes_by_counts')


# In[38]:


#filter out genes with counts < 4500 and percent of mitochondrial genes < 25
adata = dat[dat.obs.n_genes_by_counts < 4500, :]
adata = dat[dat.obs.pct_counts_mt < 25, :]


# In[39]:


#Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)


# In[40]:


sc.pp.normalize_total


# In[41]:


#Identify highly-variable gen
sc.pp.log1p(adata)
#Identify highly-variable genes with some conditions
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[42]:


#make the plot with highly variable genes
sc.pl.highly_variable_genes(adata)


# In[43]:


#select only highly variable genes
adata = adata[:, adata.var.highly_variable]


# In[18]:


#look at the plot
sc.pl.highly_variable_genes(adata)


# In[19]:


#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


# In[44]:


#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)


# In[45]:


sc.pl.highest_expr_genes(adata, n_top=20)


# In[46]:


#pca and observe some genes(B2m,Ly86,Kif5c)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=['B2m', 'Ly86', 'Kif5c'])


# In[47]:


#take a look at adata object names
adata.var_names
adata


# In[48]:


#elbow plot with PCAs
sc.pl.pca_variance_ratio(adata, log=True)


# In[49]:


#neighbors 
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)


# In[50]:


#umap method for dimencion reduction
sc.tl.umap(adata)


# In[51]:


sc.pl.umap(adata, color= ['B2m', 'Ly86', 'Kif5c'] )


# In[52]:


#using leiden algoritm for clusterization
sc.tl.leiden(adata)


# In[53]:


#observe clusters with umap
sc.pl.umap(adata, color='leiden')


# In[31]:


#Let us compute a ranking for the highly differential genes in each cluster.
#and select the firs 25
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[56]:


sc.settings.verbosity = 2  # reduce the verbosity


# In[57]:


#As an alternative, let us rank genes using logistic regression
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[58]:


#Show the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)


# In[66]:


#make dataframe with ranked genes
df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)


# In[67]:


df


# In[61]:


df.info()


# In[70]:


#compare cluster1 genes, only stores top 100 by default
tt = sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')


# In[71]:


#use wilcoxon method
wil = sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')


# In[72]:


ttov = sc.tl.rank_genes_groups(adata, 'leiden', method='t-test_overestim_var')


# In[79]:


#observe the genes with heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True)


# In[80]:


#another type of plot which can show as the cell fractions by size of circle
#and gene expression by color
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)


# In[81]:


sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5)


# In[83]:


#select top 20 genes from 1 cluster
sc.pl.rank_genes_groups(adata, groups=['1'], n_genes=20)


# In[87]:


#genes with same expression
# convert numpy.recarray to list
mynames = [x[0] for x in adata.uns['rank_genes_groups']['names'][:20]]
sc.pl.stacked_violin(adata, mynames, groupby="leiden")


# In[99]:


import gseapy as gp


# In[101]:



gene_set_names = gp.get_library_name(organism='Mouse')
print(gene_set_names)


# In[ ]:





# In[35]:


#Compare to a single cluster:
sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)


# In[ ]:


get_ipython().system('jt -r')


# In[36]:


sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)


# In[64]:


#If you want to compare a certain gene across groups, use the following.
sc.pl.violin(adata, ['Palm', 'B2m'], groupby='leiden')


# In[65]:





# In[39]:


sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')


# In[40]:


sc.pl.dotplot(adata,['Lyz2', 'B2m'], groupby='leiden')


# In[42]:


sc.pl.stacked_violin(adata,['Lyz2', 'B2m', 'Atf3'] , groupby='leiden', rotation=90)


# In[ ]:




