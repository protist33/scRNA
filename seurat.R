library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
#I can observe the methods what could be used with seurat object 
utils::methods(class = 'Seurat')
#read the directory with 10X genomics files
dir = "C:\\Users\\Виталик\\Documents\\RNAproject\\10x genomics"
#look at the files in dir
list.files(dir)
#read data from the dir(barcodes, features/genes, matrix)
data <- Read10X(data.dir = dir)
#make a seurat object from the file.
sob <- CreateSeuratObject(counts = data, project = "scdat", min.cells = 3,
                          min.features = 200)

#add another column with persentage of mitochondrial genes
# the genes have prefix "mt-"
sob[["percent.mit"]] <- PercentageFeatureSet(sob, pattern = '^mt-')
#Violine plot help me look at the percentage of mitochondrial genes
#features, count distributions
VlnPlot(sob, features = c("nFeature_RNA", "nCount_RNA", "percent.mit"),
        ncol = 3)
#and scatter plot
pl1 <- FeatureScatter(sob, feature1 = "nCount_RNA", 
                      feature2 = "nFeature_RNA")

pl2 <- FeatureScatter(sob, feature1 = "nCount_RNA",
                      feature2 = "percent.mit")

pl1 + pl2
#nFeature_RNA is the number of genes detected in each cell
#nCount_RNA total number of molecules(reads) detected within a cell.
#after i can subset some genes
sob <- subset(sob, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mit < 15)
#for avoid some replication problems i need use log normalisation
#sob1 <- NormalizeData(sob, normalization.method = "LogNormalize", 
#about log,  https://www.youtube.com/watch?v=VSi0Z04fWj0&t=81s
sob1 <- NormalizeData(sob, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#find more variable genes, select top 10 and make a plot
#
sob2 <- FindVariableFeatures(sob1, selection.method = "vst",
                             nfeatures = 2000)

top10 <- head(VariableFeatures(sob2),10)
plot1 <- VariableFeaturePlot(sob2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

#making the list of genes names
all.genes <- rownames(sob2)
#after scaling data
sob3 <- ScaleData(sob2, features = all.genes)
head(sob3[["RNA"]]@scale.data)
#make PCA
#one perfect video about PCA https://www.youtube.com/watch?v=FgakZw6K1QQ
sob4 <- RunPCA(sob3, features = VariableFeatures(object = sob3))
# Examine and visualize PCA results a few different ways
print(sob4[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sob4, dims = 1:2, reduction = 'pca')
DimPlot(sob4, reduction = "pca")
DimHeatmap(sob4, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sob4, dims = 1:5, cells = 500, balanced = TRUE)
sob4@reductions$pca[1:10, 1:3]
#making JackStraw method and the plot
sob5 <- JackStraw(sob4, num.replicate = 100)
sob5 <- ScoreJackStraw(sob5, dims = 1:20)
JackStrawPlot(sob5, dims = 1:15)
#and elbow plot for select PC
ElbowPlot(sob5)
#after i used KNN method an find clusters
#use UMAP method for dimentional reduction
#one video about UMAP function
#https://www.youtube.com/watch?v=eN0wFzBA4Sc&t=974s
sob6 <- FindNeighbors(sob5, dims = 1:10)
sob6 <- FindClusters(sob6, resolution = 0.5)
head(Idents(sob6), 5)
um <- RunUMAP(sob6, dims = 1:10)
DimPlot(um, reduction = 'umap', label = TRUE)
#i can  find all markers(differencially expressed genes) of cluster 2
clmark <- FindMarkers(um, ident.1 = 3, min.pct = 0.25)
head(clmark)
features <- c('Vsnl1', '1500015O10Rik','Rgs4','Wnt4','Smpx')
#is a lot of methods to observe the distribution
#for example the ridge plot of PC_1 by clusters
#by contrast of PC_2
r1 <- RidgePlot(um, features = 'PC_1')
r2 <- RidgePlot(um, features = 'PC_2')
r1+r2


#also it is possible to observe some genes expression by cluster
RidgePlot(um, features = features)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
clmarkall <- FindAllMarkers(um, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25)
#so i group top 5 genes for each cluster(ordered by avg_log2FC)
marks <- clmarkall %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC)
#and read the table
View(marks)
#from the marks table i manually take five genes from three clusters
#and made the vectors and 
#cell annotation is especial problem(subjectivity and other)
#top 5 genes from cluster 0
features0 <- c("1500015O10Rik", "Wnt4", "Hmgb3", "Vsnl1", "Steap1")
#from cluster 3
features3 <- c("Lyz2","C1qb","","Ctss","Apoe","C1qa")
#from cluster 5
features5 <- c("Il1b","Ccl4","Ccl3","Cxcl2","Cxcl1")
#i take a look at with some different plots
#feature plots
FeaturePlot(um, features = features0)
FeaturePlot(um, features = features3)
FeaturePlot(um, features = features5)
#violine plots
VlnPlot(um, features = features0)
VlnPlot(um, features = features3)
VlnPlot(um, features = features5)

#Existen multiple tools for the automated annotation
#i trying to use SingleR package
library(SingleR)
BiocManager::install("celldex")
#celldex library could make connection for some datasets
#it should be used for reference
library(celldex)
ref <- celldex::MouseRNAseqData()


library(SingleCellExperiment)
#DietSeurat can save only sertain seurat aspects
diet <- DietSeurat(um)
#after i should convert it to a single cell object
sce <- as.SingleCellExperiment(diet, assay = "RNA")
#on the next step i try the SingleR function using sce like a test 
#data and the MouseSingleCellData(ref) like a reference data
main <- SingleR(test = sce, assay.type.test = 1, ref = ref,
                labels = ref$label.main)
fain <- SingleR(test = sce, assay.type.test = 1, ref = ref,
                labels = ref$label.fine)
#have fogot to notice what it is possible to observe some statistics
#about single cell object
colMeans(counts(sce))
colData(sce)

cell_info <- as.data.frame(colData(sce))
SingleR::plotScoreHeatmap(main)
SingleR::plotScoreHeatmap(fain)
SingleR::plotDeltaDistribution(main)


table(main$pruned.labels)
table(fain$pruned.labels)

um@meta.data$main <- main$pruned.labels
um@meta.data$fine <- fain$pruned.labels
um@meta.data$mainlab <- main$labels
um@meta.data$fainlab <- fain$labels

head(um@meta.data)
meta <- um@meta.data
meta <- meta[order(meta$seurat_clusters),]
View(meta)
cl3 <- meta[meta$seurat_clusters == 3,]
cl5 <- meta[meta$seurat_clusters == 5,]

maindim <- DimPlot(um, reduction = "umap", group.by = "mainlab", label = TRUE)
faindim <- DimPlot(um, reduction = "umap", group.by = "fainlab", label = TRUE)
maindim + faindim
faindim

BiocManager::install("TabulaMurisData")
library(TabulaMurisData)
BiocManager::install("scRNAseq")
library("scRNAseq")
out <- listDatasets()

BiocManager::install("ExperimentHub", force = TRUE)

