##scRNA
I am studying some bioinformatics methods, so I've been making the project  with single cell RNA sequencing data, which I found here 
https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSM4914027
I take only one part of the experiment, adrenal gland's tissue and make the Seurat pipeline, GSEA, enrichment analysis with R and Python.
First time I have made the directory with 10x genomic files(matrix.mtx, barcodes.tsv and genes.tsv)
Let's observe the files in directory. With R 
```
dir = "C:\\Users\\Виталик\\Documents\\RNAproject\\10x"
```
look at the files in dir
```
list.files(dir)
```
With python: 
```
import os
os.listdir(r"C:\Users\Виталик\Documents\RNAproject\10x")
```
i try to make some tables for the single cell pipeline

| R | Python |
| --- | --- |
|libraries:                         |packages:
|library(dplyr)                     |import numpy as np
|library(Seurat)                    |import pandas as pd
|library(ggplot2)                   |import matplotlib.pyplot as plt
|library(patchwork)                 |import seaborn as sns
|                                   |import scanpy as sc

After  I need to make a Seurat object in R  and AnnData object in Python(scanpy), from the files in directory
using R
read data from the dir(barcodes, features/genes, matrix)
```
data <- Read10X(data.dir = dir)
```
make a seurat object from the file.
```
sob <- CreateSeuratObject(counts = data, project = "scdat", min.cells = 3,
                          min.features = 200)
```
With python  
```
dat = sc.read_10x_mtx(r"C:\Users\Виталик\Documents\RNAproject\10x", var_names="gene_symbols")
```
In this place the pipelines are different. Here I try to the structure of a Seurat object, and AnnData
Seurat
![seur](https://user-images.githubusercontent.com/90727271/178681897-5cdb1904-0472-4ac1-9785-f4662d7688b5.jpg)
AnnData
<img width="429" alt="AnnData" src="https://user-images.githubusercontent.com/90727271/178682481-147f011f-e804-4f0b-adc9-d0e9cd23f7a3.png">
In the scanpy I can observe some differentially expression genes
I should make variant names unique, for no count one gene as two or more name. 
```
dat.var_names_make_unique()
```
After, the pipeline tell me to look at 20 highest expression genes
The function estimates percent of counts exists for each gene and show me a plot
```
sc.pl.highest_expr_genes(dat, n_top=20)
```
![exp1](https://user-images.githubusercontent.com/90727271/178683345-8756d50c-2009-43bc-b02d-8d176649dc3d.png)
I could see a lot of genes with prefix “mt-” are  mitochondrial genes which usually went from disrupted cells, and the lider is malat1 it is a long noncoding RNA.
During the both pipelines I need to filter out the cells with a lot of mitochondrial genes because of its tells as about cell damage. 
with R
add another column with persentage of mitochondrial genes
the genes have prefix "mt-"
```
sob[["percent.mit"]] <- PercentageFeatureSet(sob, pattern = '^mt-')
```
with python
Lets filter the set of genes with at least 200 genes
```
sc.pp.filter_cells(dat, min_genes=200)
```
And at least 3 cells
```
sc.pp.filter_genes(dat, min_cells=3)
```
And make a set of mitochondrial genes
```
dat.var['mt'] = dat.var_names.str.startswith('mt-')
```
And calculate quality control metrics 
```
sc.pp.calculate_qc_metrics(dat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```
After it, I make some plots for observe distributions, remember what I used the same data
with R
```
VlnPlot(sob, features = c("nFeature_RNA", "nCount_RNA", "percent.mit"),
        ncol = 3)
```
![vlnpl](https://user-images.githubusercontent.com/90727271/178684636-ebce1417-c268-4e6c-a7cd-d71e8da21a89.jpg)
Python
```
sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```
![pythonviol](https://user-images.githubusercontent.com/90727271/178686434-1b721533-b433-48a8-9e34-bf090b3797f7.png)
So before the process of cutting some cells i should to observe counts and mitochondrial genes by plots
R decision
```
pl1 <- FeatureScatter(sob, feature1 = "nCount_RNA", 
                      feature2 = "nFeature_RNA")

pl2 <- FeatureScatter(sob, feature1 = "nCount_RNA",
                      feature2 = "percent.mit")

pl1 + pl2
```
![scatterplot](https://user-images.githubusercontent.com/90727271/178687268-b219fa8e-b3a9-4d11-b4a6-c5c0ace6319c.jpg)

Python decision
```
sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```
|---|---|
|![pythonmit](https://user-images.githubusercontent.com/90727271/178688255-3c3a8047-285a-4779-a144-afb782109500.png) |![pythonfeat](https://user-images.githubusercontent.com/90727271/178688343-0a9d38d9-ff06-43d9-9ab3-09170d0e648d.png)
From these plots, I'm trying to understand a threshold and filter out some genes with enormous mitochondrial genes expression and small gene expression (nFeature_RNA is the number of genes detected in each cell) , 
nCount_RNA total number of molecules (reads) detected within a cell shouldn't be too large. 
R decision
after i can subset some genes
```
sob <- subset(sob, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mit < 15)
```
Python
filter out genes with counts < 4500 and percent of mitochondrial genes < 25
```
adata = dat[dat.obs.n_genes_by_counts < 4500, :]
adata = dat[dat.obs.pct_counts_mt < 25, :]
```
For avoid some replication problems i need use log normalisation. One good video about log normalization 
https://www.youtube.com/watch?v=VSi0Z04fWj0&t=81s.
R decision with lognormalization
```
sob1 <- NormalizeData(sob, normalization.method = "LogNormalize", scale.factor = 10000)
```
and scanpy total-count normalization
Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
```
sc.pp.normalize_total(adata, target_sum=1e4)
```
Next step is highly variable gene observing and outlying, because of these are our points of  interest 
R function and plot, for finding more variable genes, select top 10 and make a plot
```
sob2 <- FindVariableFeatures(sob1, selection.method = "vst",
                             nfeatures = 2000)
top10 <- head(VariableFeatures(sob2),10)
plot1 <- VariableFeaturePlot(sob2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2
```
![higvar](https://user-images.githubusercontent.com/90727271/178699518-d6919bc9-e9d9-40e5-abdb-f7f48d63581e.jpg)
and scanpy decision and plot
#Identify highly-variable gen
```
sc.pp.log1p(adata)
```
#Identify highly-variable genes with some conditions
```
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
```
make the plot with highly variable genes
```
sc.pl.highly_variable_genes(adata)
```
![pythonhigh](https://user-images.githubusercontent.com/90727271/178700735-50feeac3-cc7b-4449-a8b1-d885f94dd9e5.png)
select only highly variable genes
adata = adata[:, adata.var.highly_variable]
Step the data scaling might be make with some different ways, in the Seurat pipeline by the way
making the list of genes names
after scaling data
```
all.genes <- rownames(sob2)
sob3 <- ScaleData(sob2, features = all.genes)
```
For understand some about scaling process I made one picture( it is only one way for me, because of I am bad in math) , 
mean for every cell is 0, and sd = 1 after the process of normalization.
![scale](https://user-images.githubusercontent.com/90727271/178701731-a338dbdc-0b38-467d-ba74-6ec993021b8a.jpg)

And now it is time for PCA analysis, one perfect video about PCA https://www.youtube.com/watch?v=FgakZw6K1QQ
#R
```
sob4 <- RunPCA(sob3, features = VariableFeatures(object = sob3))
```
Examine and visualize PCA results a few different ways
```
print(sob4[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sob4, dims = 1:2, reduction = 'pca')
DimPlot(sob4, reduction = "pca")
DimHeatmap(sob4, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sob4, dims = 1:5, cells = 500, balanced = TRUE)
```
![rpca2](https://user-images.githubusercontent.com/90727271/178702989-70c41f2c-4c0c-440b-8dc2-49166d7d3f1a.jpg)
#scanpy
pca and observe use two genes as color(B2m,Ly86)
```
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=['B2m','Ly86'])
```
![pythonpca](https://user-images.githubusercontent.com/90727271/178703418-8591c169-af16-4644-b4b8-00403e122142.png)
With R I make JackStraw method this is significance test, in detail described, in scanpy it doesnt util, I think because can not find it
http://cran.nexr.com/web/packages/jackstraw/vignettes/jackstraw.pdf
making JackStraw method and the plot
sob5 <- JackStraw(sob4, num.replicate = 100)
sob5 <- ScoreJackStraw(sob5, dims = 1:20)
JackStrawPlot(sob5, dims = 1:15)
![jack](https://user-images.githubusercontent.com/90727271/178703953-f1e0d99a-789a-4291-84ea-0f7e8b489162.jpg)
Observe elbow plot for PCA selection, with R
ElbowPlot(sob5)
![elbow](https://user-images.githubusercontent.com/90727271/178704453-04bee1b8-a940-40ae-9c22-d20c1a11cf8d.jpg)
and scanpy
sc.pl.pca_variance_ratio(adata, log=True)
![elbowpyt](https://user-images.githubusercontent.com/90727271/178704718-a6a631e6-f98b-4a1b-9ba4-b983e9662c4f.png)

My final aim is to find some clusters(which show us some cells ), so, in both cases I use KNN method, and use the results for UMAP dimentional reduction.
So, with R I can make it by
```
sob6 <- FindNeighbors(sob5, dims = 1:10)
sob6 <- FindClusters(sob6, resolution = 0.5)
head(Idents(sob6), 5)
um <- RunUMAP(sob6, dims = 1:10)
DimPlot(um, reduction = 'umap', label = TRUE)
```
![umapr](https://user-images.githubusercontent.com/90727271/178706655-d32687a5-47e8-4aeb-8f64-73ccc31898ea.jpg)

#scanpy
neighbors 
```
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
```
umap method for dimencion reduction
```
sc.tl.umap(adata)
```
using leiden algoritm for clusterization
```
sc.tl.leiden(adata)
```
using leiden algoritm for clusterization
```
sc.pl.umap(adata, color='leiden')
```
![pytumap](https://user-images.githubusercontent.com/90727271/178707263-e827f016-3ed6-403f-ade6-85acdaee2025.png)
I am able to observe highly expressed genes in every cluster and used the top for manual annotation
with R
find markers for every cluster compared to all remaining cells, report only the positive ones
```
clmarkall <- FindAllMarkers(um, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25)
 ```
so i group top 5 genes for each cluster(ordered by avg_log2FC) and read the table
```
marks <- clmarkall %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC)
View(marks)
```
![clust](https://user-images.githubusercontent.com/90727271/178708702-3e1680e5-31f9-407a-bfbe-bec28a8e4735.jpg)
With scanpy I can make a plot like this, let us compute a ranking for the highly differential genes in each cluster and select the firs 25
```
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```
![geneselbp](https://user-images.githubusercontent.com/90727271/178709401-5eacdea9-44a8-42c6-95ef-9ffd94ea46c8.png)
I might using different statistics methods, like t-test or logistic regression, with scanpy it is possible to take a look at the table.
```
sc.settings.verbosity = 2 
df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
```
![pclustrun](https://user-images.githubusercontent.com/90727271/178710256-0ff70ac2-79e4-4dd5-9a52-19d1c122fe7e.jpg)
For manual annotation with R make list of top genes for some cluster, for example 3 and 5.
from cluster 3
```
features3 <- c("Lyz2","C1qb","","Ctss","Apoe","C1qa")
```
from cluster 5
```
features5 <- c("Il1b","Ccl4","Ccl3","Cxcl2","Cxcl1")
```
Feature and Violine plots make it possible to look at the expression
```
FeaturePlot(um, features = features3)
```
![exprplot](https://user-images.githubusercontent.com/90727271/178710664-bff106b8-7c70-4f01-a699-c4f7fbdd3670.jpg)
```
VlnPlot(um, features = features3)      
```

After, it is possible to use many databases and papers to find information about cluster genes and their association with cell types. I have found what my clusters 3 and 5 can be associated with macrophages. I stopped on it and didn't find which type of macrophages. 
With Scanpy I used, for example, some genes from my third cluster in R, and look at the expression
```
sc.pl.umap(adata, color= ['Lyz2',"C1qb","Ctss","C1qa"])
```
![pyt3genes](https://user-images.githubusercontent.com/90727271/178711919-d2784018-6805-46d4-bed6-5d2c91fdf44e.png)
See what they are highly expressed on the “island”, so i might think about this island how about macrophages

Now I'm trying to research and understand automated single cell annotation in Python and don't see many packages for it. From another way, I have found some packages for annotation with R. Here I used SingleR package.
```
library(SingleR)
BiocManager::install("celldex")
```
celldex is the package what have some dataset for using as a reference, my data are about mice, so, I used the dataset
```
library(celldex)
ref <- celldex::MouseRNAseqData()
```
Before I need to convert my seurat object to a single cell experiment, i need make the seurat object less using the ‘diet’ function
```
library(SingleCellExperiment)
```
DietSeurat can save only sertain seurat aspects
```
diet <- DietSeurat(um)
```
after i should convert it to a single cell object
```
sce <- as.SingleCellExperiment(diet, assay = "RNA")
```
After I use two way for annotation with singleR function, ‘main’ option does it with less detail, and ‘fine’ method for more detailed clustering and annotation
data and the MouseSingleCellData(ref) like a reference data, lets look at the heatmap plot
```
main <- SingleR(test = sce, assay.type.test = 1, ref = ref,
                labels = ref$label.main)
fain <- SingleR(test = sce, assay.type.test = 1, ref = ref,
                labels = ref$label.fine)
SingleR::plotScoreHeatmap(main)
```
![Rplot02](https://user-images.githubusercontent.com/90727271/178713044-295d35a2-8c5f-44ff-b7a9-c8ee3184e397.png)

After I can add the results of annotation to my seurat object and use for another plots

add the labels to metadata in my seurat object.
```
um@meta.data$main <- main$pruned.labels
um@meta.data$fine <- fain$pruned.labels
```
take the labels of cells
```
um@meta.data$mainlab <- main$labels
um@meta.data$fainlab <- fain$labels
```
in the end making the labeled umap plot
```
maindim <- DimPlot(um, reduction = "umap", group.by = "mainlab", label = TRUE, label.size = 1.5,
                   label.box = TRUE)
```
![dimplotlm](https://user-images.githubusercontent.com/90727271/178713707-8035a1e9-0360-4d46-8f51-755cb74268a8.png)


































































































































































             







        






























































