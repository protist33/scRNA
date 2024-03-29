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
##With python: 
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
##using R
read data from the dir(barcodes, features/genes, matrix)
```
data <- Read10X(data.dir = dir)
```
make a seurat object from the file.
```
sob <- CreateSeuratObject(counts = data, project = "scdat", min.cells = 3,
                          min.features = 200)
```
##With python  
```
dat = sc.read_10x_mtx(r"C:\Users\Виталик\Documents\RNAproject\10x", var_names="gene_symbols")
```
In this place the pipelines are different. Here I try to the structure of a Seurat object, and AnnData
#Seurat
![seur](https://user-images.githubusercontent.com/90727271/178681897-5cdb1904-0472-4ac1-9785-f4662d7688b5.jpg)
#AnnData
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
##with R
add another column with persentage of mitochondrial genes
the genes have prefix "mt-"
```
sob[["percent.mit"]] <- PercentageFeatureSet(sob, pattern = '^mt-')
```
##with python
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
##with R
```
VlnPlot(sob, features = c("nFeature_RNA", "nCount_RNA", "percent.mit"),
        ncol = 3)
```
![vlnpl](https://user-images.githubusercontent.com/90727271/178684636-ebce1417-c268-4e6c-a7cd-d71e8da21a89.jpg)
##Python
```
sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```
![pythonviol](https://user-images.githubusercontent.com/90727271/178686434-1b721533-b433-48a8-9e34-bf090b3797f7.png)
So before the process of cutting some cells i should to observe counts and mitochondrial genes by plots
##R decision
```
pl1 <- FeatureScatter(sob, feature1 = "nCount_RNA", 
                      feature2 = "nFeature_RNA")

pl2 <- FeatureScatter(sob, feature1 = "nCount_RNA",
                      feature2 = "percent.mit")

pl1 + pl2
```
![scatterplot](https://user-images.githubusercontent.com/90727271/178687268-b219fa8e-b3a9-4d11-b4a6-c5c0ace6319c.jpg)

##Python decision
```
sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```
|---|---|
|![pythonmit](https://user-images.githubusercontent.com/90727271/178688255-3c3a8047-285a-4779-a144-afb782109500.png) |![pythonfeat](https://user-images.githubusercontent.com/90727271/178688343-0a9d38d9-ff06-43d9-9ab3-09170d0e648d.png)
From these plots, I'm trying to understand a threshold and filter out some genes with enormous mitochondrial genes expression and small gene expression (nFeature_RNA is the number of genes detected in each cell) , 
nCount_RNA total number of molecules (reads) detected within a cell shouldn't be too large. 
##R decision
after i can subset some genes
```
sob <- subset(sob, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mit < 15)
```
##Python
filter out genes with counts < 4500 and percent of mitochondrial genes < 25
```
adata = dat[dat.obs.n_genes_by_counts < 4500, :]
adata = dat[dat.obs.pct_counts_mt < 25, :]
```
For avoid some replication problems i need use log normalisation. One good video about log normalization 
https://www.youtube.com/watch?v=VSi0Z04fWj0&t=81s.
##R decision with lognormalization
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
##and scanpy decision and plot
###Identify highly-variable gen
```
sc.pp.log1p(adata)
```
###Identify highly-variable genes with some conditions
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
##R
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
##scanpy
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
##and scanpy
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

##scanpy
nearest neighbors method 
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
##with R
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

Exist a lot of packages for automated cluster annotation, which work with some different approaches
![clustering](https://user-images.githubusercontent.com/90727271/185595697-9cbc6639-00a2-4801-9cb4-63193170a5a4.jpg)
I am going to try some methods and packages with the same data for better understanding. Above I made analysis with singleR,
and now trying scCATCH package from 
https://github.com/ZJUFanLab/scCATCH
```
#scCATCH method
install.packages("scCATCH")
library('scCATCH')
#clusters should be character
clust <- as.character(um@meta.data$seurat_clusters)
cat1 <- createscCATCH(um@assays$RNA@data, cluster = clust)
cat2 <- findmarkergene(object = cat1, species = 'Mouse', 
                       marker = cellmatch, tissue = 'Kidney',
)
cat3 <- findcelltype(cat2)
```
I made the picture with single cell object structure
<img width="618" alt="catch object" src="https://user-images.githubusercontent.com/90727271/185597900-3338562f-93e7-43b3-bb60-64deafe0d634.png">
I couldnt find the adrenal gland tissue, so I used kidney tissue references. I like what the package find some pubmed articles that are used
for cluster annotation, so I could  make my own reference data with papeles manually. My 'macrophage clusters' are identifined as macrophage cells,
as in my manually annotation and singleR
reference 
##So I wanna add here some trajectory analysis, I hadn't made it before because have found one good tutorial and 
finally understand how to do it, but need work for understand it deeper, ok, going to make the trajectory of my macrophages.
#Another way I try to understand the package "Clustifyr", have found some good tutorials and examples,
and do it by two ways: one with separate matrix and metadata and another with the whole Seurat object
```
BiocManager::install("clustifyr")
library(clustifyr)
library(ggplot2)
library(cowplot)
# this package can it one whole Seurat object,
#but, firstly I make it by another way, separate matrix
#and metadata from my Seurate object
matr <- um@assays$RNA@data
meta <- um@meta.data
#my metadata hasn't umap, so I need extract it
#and add manually 
umap1 <- um@reductions$umap@cell.embeddings
umap11 <- umap1[,1]
umap22 <- umap1[,2]
meta$UMAP_1 <- umap11
meta$UMAP_2 <- umap22
#I use my genes, ordered by log_FC
genes <- marks$gene[1:500]
#as a reference matrix I use the datahub
BiocManager::install("clustifyrdatahub")
#look at the useful datasets
knitr::kable(dplyr::select(
  read.csv(system.file("extdata", "metadata.csv", package = "clustifyrdatahub")),
  c(1, 9, 2:7)))
library(ExperimentHub)
eh <- ExperimentHub()
#wanna try two reference datasets
refs <- query(eh, 'clustifyrdatahub')
#mouse cell atlas
refs_mous <- refs[['EH3444']]
#ref mouse data from singleR Package 
refs_mous2 <- refs[['EH3447']]
#lets make the clustify function
res <- clustify(input = matr, 
                metadata = meta,
                cluster_col = 'seurat_clusters',
                ref_mat = refs_mous,
                query_genes = genes)
                
#take it best for each cluster
rescor <- cor_to_call(cor_mat = res,
                    cluster_col = 'seurat_clusters')
#add it to metadata
meta2 <- call_to_metadata(res = rescor,
                         metadata = meta,
                         cluster_col = 'seurat_clusters')
                         
#plot_cor give me some errors I gonna fixed it later 
plot_cor_heatmap(cor_mat = res)
#plot similarity measures with umap
corr_um <- plot_cor(cor_mat = res,
                    metadata = meta,
                    data_to_plot = colnames(res)[6:7],
                    cluster_col = 'seurat_clusters')
```
![makrplot](https://user-images.githubusercontent.com/90727271/185801888-bb63209c-fdfd-4ea9-aa87-87de58bcb849.png)
```
#best types
clust_types <- plot_best_call(cor_mat = res,
                              metadata = meta,
                              do_label = TRUE,
                              do_legend = FALSE,
                              cluster_col = 'seurat_clusters')   
   #dimention plot(umap), colored by type
known_types <- plot_dims(data = meta2,
                         feature = 'type',
                         do_label = TRUE,
                         do_legend = FALSE
                         ) + ggtitle("cell type")
```
#clust_type plot
![clusttype](https://user-images.githubusercontent.com/90727271/185802270-61401773-75f0-42c1-91dd-776ee31a48b2.jpg)
```
plot_grid(known_types, clust_types)
#and the simple way for using Seurat object 
resseu <- clustify(input = um,
                   ref_mat = refs_mous,
                   cluster_col = 'seurat_clusters',
                   obj_out = TRUE)
                 
resseu2 <- clustify(input = um,
                    ref_mat = refs_mous2,
                    cluster_col = 'seurat_clusters',
                    obj_out = TRUE
                    )   
plotres1 <- DimPlot(resseu, reduction = 'umap', group.by = 'type', 
                  label = TRUE, label.box = TRUE, label.size = 2)
plotres2 <- DimPlot(resseu2, reduction = 'umap', group.by = 'type',
                    label = TRUE, label.box = TRUE)
```
#With cell atlas annotatyon, I think what something wrong a lot of cells that could not be 
![plotrefseq1](https://user-images.githubusercontent.com/90727271/185802352-950d32d1-3b9e-4ad0-af2a-d167be568f5d.png)
and with another reference, it could see my lovely macrophages but don't understant the huge part of cells
![plotrefseq2](https://user-images.githubusercontent.com/90727271/185802652-e109c38a-40b8-4195-8c34-a81da05f54c9.png)

#monocle
```
library(SingleCellExperiment)
library(monocle3)
```
I use the monocle3 package and pipeline in the end of my Seurat process, 
firstly subset only macrophages from the Seurat.
```
macro <- subset(um, main == 'Macrophages')
```
I make the object for monocle3 pipeline with as.cell_data_set
```
unique(macro@meta.data$seurat_clusters)
cds1 <- as.cell_data_set(macro)
```
Start preparation for the graph 
```
fData(cds1)
rownames(fData(cds1))[1:10]
fData(cds1)$short_gene_name <- rownames(fData(cds1))
counts(cds1)
```
Counts give me a matrix with counts, cells and genename,
after I make partition as factor. 
```
recreate.partition <- c(rep(1, length(cds1@colData@rownames)))
names(recreate.partition) <- cds1@colData@rownames
recreate.partition <- as.factor(recreate.partition)
red <- reducedDim(cds1, 'UMAP')
cds1@clusters$UMAP <- red
cds1@clusters$UMAP$partitions <- recreate.partition
list_cl <- macro@active.ident
cds1@clusters$UMAP$clusters <- list_cl
cds1@int_colData@listData$reducedDims$UMAP <- macro@reductions$umap@cell.embeddings
cl_name <- macro@meta.data$fine

plot1 <- plot_cells(cds1,
                    color_cells_by = 'cluster',
                    label_groups_by_cluster = FALSE,
                    group_label_size = 5) +
  theme(legend.position = 'right')

cds3 <- learn_graph(cds1, use_partition = FALSE)
plot_cells(cds3,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           group_label_size = 5)

```
![traj2](https://user-images.githubusercontent.com/90727271/185802973-7b96f812-774f-460a-a6e4-98e8673d3798.png)

Lets look at the enrichment analysis process
Install packages for enrichment analysis

```
library(biomaRt)
library(escape)
```
Also, I need one gene set and one gene list 
```
gen_set <- idsub$entrezgene_id
gen_set <- as.character(gen_set)
```
Durind the tutorial reading i had some problem with gen list sense understanding
it should be the list of genes from my data, ordered by decreasing avg_log2FC
Fist step I use bioMart for adding some others gene names from Ensembl, with getBM function
and merge the results with my dataset of clusters, 
```
library(biomaRt)
library(escape)
GS.hallmark <- getGeneSets(library = "H")
gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                "entrezgene_id"), 
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")))
id_marks <- merge(clmarkall, gene_id[,c(1,2,3)], by.x="gene", by.y="external_gene_name")
```
For the gen list I use not all my genes, only because of my laptop isn't very powerful, 
so make a subset with log2FC some more than median value
```
id_marksrg <- id_marks %>%
  arrange(desc(id_marks$avg_log2FC))
summary(id_marksrg)
idsub <- subset(id_marksrg, subset = avg_log2FC > 0.75)
dim(id_marksrg)
dim(idsub)
```
Ok, it is time to make gene set
```
gen_set <- idsub$entrezgene_id
gen_set <- as.character(gen_set)
```
And the gene list 
```
genlog <- idsub$avg_log2FC
names(genlog) <- as.character(idsub$entrezgene_id)
genlog <- na.omit(genlog)
genlist <- sort(genlog, decreasing = TRUE)
```
Need install two databases and make enrichment analysis with mice database
```
library(clusterProfiler)
hsdb <- BiocManager::install("org.Hs.eg.db", force = TRUE)
musdb <- BiocManager::install("org.Mm.eg.db", force = TRUE)
engo <- enrichGO(gene = gen_set, OrgDb = musdb, pvalueCutoff = 0.05,
                 ont = "all", readable = T)
```
I describe some about the GO object structure
![go](https://user-images.githubusercontent.com/90727271/179626367-062b496d-bb04-40ba-b67c-b4ce138df05a.jpg)
#gene enrichment analysis can be visualyse with the dotplot
```
barplot(engo, showCategory = 20)
```
![barplotenrich](https://user-images.githubusercontent.com/90727271/179952392-9c453040-5c17-480b-b2b6-80176d18f127.png)
#select go terms and other
```
engo2 <- engo@result[2,]
```
groupGO function
```
ggo1 <- groupGO(gene = gen_set, OrgDb = musdb, ont = 'CC',
                level = 3, readable = TRUE)
ggo1@result[1:3,]
```
After, I can use the site http://geneontology.org/
for observing enrichGO results, for example I take the second object
<img width="515" alt="enrichscr" src="https://user-images.githubusercontent.com/90727271/179960011-18fad418-9d19-4f20-8c35-510d32d7ff04.png">
And paste it on the site, and observe result, it tell me about ubiquitin ligase complex, it is a cellular component, I had made groupGO analysis about cellular components(ont = 'CC'). 
<img width="940" alt="gosite" src="https://user-images.githubusercontent.com/90727271/179960929-198e9f85-2838-4215-841f-5848685a938e.png">
```
gsego <- gseGO(geneList = genlist, OrgDb = get('org.Mm.eg.db'), ont = "CC", minGSSize = 100,
              maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)

gsego1 <- gseGO(geneList = genlist, OrgDb = get('org.Mm.eg.db'), ont = "BP", minGSSize = 100,
               maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)

dotplot(gsego, showCategory = 20)
```
![gsegobp](https://user-images.githubusercontent.com/90727271/180742661-d1f28bf9-7538-49c1-a00b-d5b63abc8d93.png)
```
gseaplot(gsego, geneSetID = 1, by = 'runningScore')
```
![gseaplot11](https://user-images.githubusercontent.com/90727271/180743963-a1b3f63e-d13b-485e-83bc-3a36a5d3a6cf.png)

Now, after I show and understand some tools and steps lets back to my data, from research paper
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
they used two group of mice, control and stressed(n=15)
something about the stress test
![stresscondition](https://user-images.githubusercontent.com/90727271/182120451-2b33ee77-f97b-4a46-8319-82cae5aadfd3.jpg)
Their mice were stressed durig 21 day and had 48 days for recovery, after they measure level of cortisol it was significantly higher for the
stressed mice. Fur score(it means grouming conditions) is different too. 
They were using Scanpy and defined 16 clusters in adrenal data, I had found twenty with Scanpy and eleven with the Seurat package
from adrenal glands they outlined Abcb1b like one new molecular marker. It is significantly unique into the zona fasciculata cells.
Three top genes were found as which what define this cells: Abcb1b, Sbsn, Srd5a2. These genes previously associated with GC transport, cell poliferation
and glucose metabolism. Cyp11b1 expressed in all zona fasciculata cells(from zFasc1 to zFasc2). 

I look at the three genes from the original paper
```
features6 <- c('Abcb1a','Abcb1b', 'Sbsn', 'Srd5a2')
FeaturePlot(um, features = features6)
VlnPlot(um, features = features6)
````
I see some expression in the 7 cluster
![abcb1a](https://user-images.githubusercontent.com/90727271/182140515-02fcaa5e-1f66-42cb-bd48-691dc98df1b9.png)
and violine expression plot
![abcbexp](https://user-images.githubusercontent.com/90727271/182140583-fe7959bd-1ded-44d3-a22a-9683f52ff99f.jpg)











               

























































































































































































































































             







        






























































