#分析计划
#1.读入矩阵 ----
rm(list=ls(all=TRUE))
.libPaths()
library (Seurat)
library(dplyr)
library(ggplot2)
packageVersion('Seurat')
data_dir1 <- ".\\expressions\\Flox"
data_dir2 <- ".\\expressions\\KO"
###step2
bat.data <- Read10X(data.dir = data_dir1)
#2.更改矩阵列名----
#3.创建seurat对象，简单过滤取子集，到标准化之前----
# 改细胞名
head(colnames(bat.data))
colnames(bat.data) <- paste0('Flox_', colnames(bat.data))
#基因至少在1个细胞表达，每个细胞至少表达10个基因
bat1 <- CreateSeuratObject(counts = bat.data, project = "Flox",min.cells = 1, min.features =10)
#bat1 <- RenameCells(bat1,levels(bat1@active.ident))
bat1@meta.data$orig.ident<-as.factor(bat1@meta.data$orig.ident)

###step3 KO 读入
bat2.data <- Read10X(data.dir = data_dir2)
# 改细胞名
colnames(bat2.data) <- paste0('KO_', colnames(bat2.data))
bat2 <- CreateSeuratObject(counts = bat2.data, project = "KO",min.cells = 1, min.features =10)
#bat2 <- RenameCells(bat2,levels(bat2@active.ident))
bat2@meta.data$orig.ident<-as.factor(bat2@meta.data$orig.ident)

### step5 分别取子集，简单过滤，Flox
bat1[["percent.mt"]] <- PercentageFeatureSet(bat1, pattern = "^mt-")
#VlnPlot(bat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =-1)
VlnPlot(bat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0)
VlnPlot(bat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0.1)
dim(bat1)
### QC
bat1.qc <- subset(bat1, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)
dim(bat1.qc)

#bat2 <- subset(bat1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 15000 & nCount_RNA > 10 & percent.mt<2)
VlnPlot(bat1.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)###pt.size =-1
VlnPlot(bat1.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)###pt.size =-1

plot1 <- FeatureScatter(bat1.qc, feature1 = "nCount_RNA", feature2 = "percent.mt")
NoLegend()
list(plot1)
plot2 <- FeatureScatter(bat1.qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
NoLegend()
list(plot2)

### step5 分别取子集，简单过滤，KO
bat2[["percent.mt"]] <- PercentageFeatureSet(bat2, pattern = "^mt-")
#VlnPlot(bat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =-1)
VlnPlot(bat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0)
VlnPlot(bat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0.1)
dim(bat2)
### QC
bat2.qc <- subset(bat2, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)
dim(bat2.qc)
#bat2 <- subset(bat1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 15000 & nCount_RNA > 10 & percent.mt<2)
VlnPlot(bat2.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.1)###pt.size =-1
VlnPlot(bat2.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)###pt.size =-1

plot3 <- FeatureScatter(bat2.qc, feature1 = "nCount_RNA", feature2 = "percent.mt")
NoLegend()
list(plot3)
plot4 <- FeatureScatter(bat2.qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
NoLegend()
list(plot4)

##保存数据 seurat_qc

#4.创建样本列表，利用张强老师的代码，过滤双胞----
bat_Flox <- bat1.qc
bat_KO <- bat2.qc
samples.list <- list(bat_Flox, bat_KO)
save(bat_Flox, bat_KO, file = 'after_seurat_QC.rdata')

##5.Doublet detection: which needs to run through the Seurat routine up to PCA/TSNE/UMAP without merging samples-----------------------

samples.list <- list(bat_Flox, bat_KO)
# Normalize and scale datasets independently--------------------------
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
})

# Run dimensionality reduction and identify nuclei clusters ---------------------------
# A clustering resolution of 0.2 was selected based on preliminary evaluations which produced similar number of clusters as identified hepatic cell types from previous studies.
samples.list <- lapply(X = samples.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'  
  x <- RunPCA(x, features = VariableFeatures(object = x), verbose = FALSE)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20, verbose = FALSE)
  x <- FindClusters(x, resolution = 0.8, verbose = FALSE)
  #x <- RunTSNE(x, dims = 1:30, max_iter = 2000, verbose = FALSE)
  #x <- RunUMAP(x, dims = 1:30, verbose = FALSE, umap.method = "umap-learn",metric = "correlation")
  # x <- RunUMAP(x, dims = 1:30, verbose = FALSE) #From Nault originally
})

#DimPlot(object = x, reduction = 'umap')

#Show some results for each sample separately

# Examine and visualize PCA results a few different ways
print(samples.list[[1]][["pca"]], dims = 1:5, nfeatures = 10)
VizDimLoadings(samples.list[[1]], dims = 1:4, reduction = "pca")
DimPlot(samples.list[[1]], reduction = "pca")

print(samples.list[[2]][["pca"]], dims = 1:5, nfeatures = 10)
VizDimLoadings(samples.list[[2]], dims = 1:4, reduction = "pca")
DimPlot(samples.list[[2]], reduction = "pca")

#Displaying expression levels (Qiang's intepretation) of top genes contributing to a PC. Each column is a cell. Setting "cells" parameter to a number plots the ‘extreme’ cells on both ends of the spectrum
DimHeatmap(samples.list[[1]], dims = 1:2, cells = 500, balanced = TRUE) 

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
#sample1_PCA <- JackStraw(samples.list[[1]], num.replicate = 100)
#sample1_PCA <- ScoreJackStraw(sample1_PCA, dims = 1:20)
#The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.
#JackStrawPlot(sample1_PCA, dims = 1:15)

#An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(samples.list[[1]])

PCs_num <- 20 #Number of PCs used for clustering
#Cluster the cells
sample1_cluster <- FindNeighbors(samples.list[[1]], dims = 1:PCs_num)
sample1_cluster <- FindClusters(sample1_cluster, resolution = 0.8)

#Run non-linear dimensional reduction (UMAP/tSNE). 
#As input to the UMAP and tSNE, use the same PCs as input to the clustering analysis.
sample1_cluster <- RunUMAP(sample1_cluster, dims = 1:20)
DimPlot(sample1_cluster, label=TRUE, reduction = "umap")

sample1_cluster <- RunTSNE(sample1_cluster, dims = 1:20, max_iter = 2000)
DimPlot(sample1_cluster, label=TRUE, reduction = "tsne")



#Doublet detection. Doublet detection is performed using DoubletFinder assuming doublet rate to be 0.7666% (~0.8%) per 1000 cells:
#(https://bioinformatics.stackexchange.com/questions/3165/what-are-doublets-in-single-cell-rna-seq-data) (https://www.biotech.wisc.edu/services/gec/services/10X-Genomics-single-cell-services)
#Need to run each sample separately because samples are not mixed during the library generation process forming gel beads.
#https://rdrr.io/github/chris-mcginnis-ucsf/DoubletFinder/
#https://github.com/chris-mcginnis-ucsf/DoubletFinder

# pK Identification (with no ground-truth)
library(DoubletFinder)
sweep.res.list_Flox <- paramSweep_v3(samples.list[[1]], PCs = 1:20, sct = FALSE) #Pcs should be equal to the number of PCs used above to obtain the number of clusters in Seurat
sweep.stats_Flox <- summarizeSweep(sweep.res.list_Flox, GT = FALSE) #GT: ground truth
bcmvn_Flox <- find.pK(sweep.stats_Flox)

#optimal pK of the PC neighborhood size for KNN classification
pk_optimal = as.numeric(levels(bcmvn_Flox$pK))[which.max(bcmvn_Flox$BCmetric)]  

plot(bcmvn_Flox$pK, bcmvn_Flox$BCmetric, type="b", xlab="pK", ylab="BCmvn") #Plot BCmvn (mean-variance-normalized bimodality coefficient) against pK. 

#Plotting histograms of distributions of pANN (proportion of artificial nearest neighbors, which is = number of ANN / total artificial doublets added) of each cell for the optimal pK case
pN <- seq(0.05, 0.3, by=0.05) #This is the range and values of pN swept by DoubletsFinder
for( i in pN )
{
  pNpK_string <- paste0("pN_",i,"_pK_",pk_optimal) #paste0() will not leave a white space between each element as "paste()" does
  hist(sweep.res.list_Flox[[pNpK_string]][["pANN"]], xlab="pANN", main=pNpK_string, breaks=seq(0,1,by=0.02)) #
}


## Heterotypic Doublet Proportion Adjustment ------------------------
homotypic.prop <- modelHomotypic(samples.list[[1]]@meta.data$seurat_clusters) # Homotypic Doublet Proportion Estimate. The number of identified cell clusters each containing cells of same similar type is important here.   #May use "RNA_snn_res.0.2" to replace "seurat_clusters" as these two metadata are the same at this stage

num_of_cells_Flox <-length(rownames(samples.list[[1]]@meta.data))#Number of cells in the sample

nExp_poi <- round(0.007666*num_of_cells_Flox/1000 *num_of_cells_Flox) # EXP_poi:expected number of doublets in the sample. Assuming 0.7666% doublet formation rate per 1000 cells - tailor for your dataset.
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop)) #Expected number of heterotypic doublets after adjustment for homotypic doublets

## Run DoubletFinder with parameters defined above to identify doublet cells or nuclei -----------------------------

pN_value <- 0.25 #recommended to just use 0.25 by the DoubletFinder authors 

#Use un-adjusted expected number of doublets
#Sample_01 <- doubletFinder_v3(samples.list[[1]], PCs = 1:PCs_num, pN = pN_value, pK = pk_optimal, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#Use adjusted expected number of doublets, i.e., heterotypic doublets. (recommended)
Sample_01 <- doubletFinder_v3(samples.list[[1]], PCs = 1:PCs_num, pN = pN_value, pK = pk_optimal, nExp = nExp_poi.adj,  sct = FALSE)

#Removing doublets

DF_subset_string <- paste0("DF.classifications_",pN_value,"_",pk_optimal,"_",nExp_poi.adj)

# Flox
table(Sample_01$DF.classifications_0.25_0.005_1724)
Sample_01_Double_removed <- subset(Sample_01, subset = DF.classifications_0.25_0.005_1724 == "Singlet")


## KO从这里开始
#An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(samples.list[[2]])

PCs_num <- 20 #Number of PCs used for clustering
#Cluster the cells
sample2_cluster <- FindNeighbors(samples.list[[2]], dims = 1:PCs_num)
sample2_cluster <- FindClusters(sample2_cluster, resolution = 0.8)

#Run non-linear dimensional reduction (UMAP/tSNE). 
#As input to the UMAP and tSNE, use the same PCs as input to the clustering analysis.
sample2_cluster <- RunUMAP(sample2_cluster, dims = 1:20)
DimPlot(sample2_cluster, label=TRUE, reduction = "umap")

sample2_cluster <- RunTSNE(sample2_cluster, dims = 1:20, max_iter = 2000)
DimPlot(sample2_cluster, label=TRUE, reduction = "tsne")



#Doublet detection. Doublet detection is performed using DoubletFinder assuming doublet rate to be 0.7666% (~0.8%) per 1000 cells: (https://bioinformatics.stackexchange.com/questions/3165/what-are-doublets-in-single-cell-rna-seq-data) (https://www.biotech.wisc.edu/services/gec/services/10X-Genomics-single-cell-services)
#Need to run each sample separately because samples are not mixed during the library generation process forming gel beads.
#https://rdrr.io/github/chris-mcginnis-ucsf/DoubletFinder/
#https://github.com/chris-mcginnis-ucsf/DoubletFinder

# pK Identification (with no ground-truth)
sweep.res.list_KO <- paramSweep_v3(samples.list[[2]], PCs = 1:20, sct = FALSE) #Pcs should be equal to the number of PCs used above to obtain the number of clusters in Seurat
sweep.stats_KO <- summarizeSweep(sweep.res.list_KO, GT = FALSE) #GT: ground truth
bcmvn_KO <- find.pK(sweep.stats_KO)

#optimal pK of the PC neighborhood size for KNN classification
pk_optimal = as.numeric(levels(bcmvn_KO$pK))[which.max(bcmvn_KO$BCmetric)]  

plot(bcmvn_KO$pK, bcmvn_KO$BCmetric, type="b", xlab="pK", ylab="BCmvn") #Plot BCmvn (mean-variance-normalized bimodality coefficient) against pK. 

#Plotting histograms of distributions of pANN (proportion of artificial nearest neighbors, which is = number of ANN / total artificial doublets added) of each cell for the optimal pK case
pN <- seq(0.05, 0.3, by=0.05) #This is the range and values of pN swept by DoubletsFinder
for( i in pN )
{
  pNpK_string <- paste0("pN_",i,"_pK_",pk_optimal) #paste0() will not leave a white space between each element as "paste()" does
  hist(sweep.res.list_KO[[pNpK_string]][["pANN"]], xlab="pANN", main=pNpK_string, breaks=seq(0,1,by=0.02)) #
}


## Heterotypic Doublet Proportion Adjustment ------------------------
homotypic.prop <- modelHomotypic(samples.list[[2]]@meta.data$seurat_clusters) # Homotypic Doublet Proportion Estimate. The number of identified cell clusters each containing cells of same similar type is important here.   #May use "RNA_snn_res.0.2" to replace "seurat_clusters" as these two metadata are the same at this stage

num_of_cells_KO <-length(rownames(samples.list[[2]]@meta.data))#Number of cells in the sample

nExp_poi <- round(0.007666*num_of_cells_KO/1000 *num_of_cells_KO) # EXP_poi:expected number of doublets in the sample. Assuming 0.7666% doublet formation rate per 1000 cells - tailor for your dataset.
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop)) #Expected number of heterotypic doublets after adjustment for homotypic doublets

## Run DoubletFinder with parameters defined above to identify doublet cells or nuclei -----------------------------

pN_value <- 0.25 #recommended to just use 0.25 by the DoubletFinder authors 

#Use un-adjusted expected number of doublets
#Sample_02 <- doubletFinder_v3(samples.list[[2]], PCs = 1:PCs_num, pN = pN_value, pK = pk_optimal, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#Use adjusted expected number of doublets, i.e., heterotypic doublets. (recommended)
Sample_02 <- doubletFinder_v3(samples.list[[2]], PCs = 1:PCs_num, pN = pN_value, pK = pk_optimal, nExp = nExp_poi.adj,  sct = FALSE)




#Removing doublets

DF_subset_string <- paste0("DF.classifications_",pN_value,"_",pk_optimal,"_",nExp_poi.adj)
# KO
table(Sample_02$DF.classifications_0.25_0.005_1082)
Sample_02_Double_removed <- subset(Sample_02, subset = DF.classifications_0.25_0.005_1082 == "Singlet")

#5.取单包子集，进行后续的分析,到亚群鉴定。-----
##2021.6.20
save(Sample_01_Double_removed, Sample_02_Double_removed,file = '11-10_doublet_filtered.rdata')


##合并样本,降维聚类
DefaultAssay(Sample_01_Double_removed) <- 'RNA'
DefaultAssay(Sample_02_Double_removed) <- 'RNA'
Flox_matrix <- GetAssayData(Sample_01_Double_removed, slot = 'counts')
as.matrix(Flox_matrix)[1:6, 1:6]
KO_matrix <- GetAssayData(Sample_02_Double_removed, slot = 'counts')
as.matrix(KO_matrix)[1:6, 1:6]
bat_Flox_filtered <- CreateSeuratObject(counts = Flox_matrix, project = "Flox")
bat_KO_filtered <- CreateSeuratObject(counts = KO_matrix, project = "KO")
bat_total<-merge(bat_Flox_filtered, bat_KO_filtered)
dim(bat_total)
dim(bat_Flox_filtered)
dim(bat_KO_filtered)

#### 6 合并降维聚类
### step5
bat_total[["percent.mt"]] <- PercentageFeatureSet(bat_total, pattern = "^mt-")
#VlnPlot(bat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =-1)
VlnPlot(bat_total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0)
VlnPlot(bat_total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0.1)


plot1 <- FeatureScatter(bat_total, feature1 = "nCount_RNA", feature2 = "percent.mt")
NoLegend()
list(plot1)
plot2 <- FeatureScatter(bat_total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
NoLegend()
list(plot2)

### 
CombinePlots(plots = list(plot1, plot2))
#表达量数据标准化： LogNormalize 的算法： A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
bat_total_3<-NormalizeData(bat_total,normalization.method = "LogNormalize",scale.factor = 10000)

#鉴定表达高变基因(2000 个），用于下游分析，如 PCA；
bat_total_4 <- FindVariableFeatures(bat_total_3, selection.method = "vst", nfeatures = 2000)
# 提取表达量变变化最高的 10 个基因；
top10 <- head(VariableFeatures(bat_total_4), 10)
top10
plot3 <- VariableFeaturePlot(bat_total_4)
NoLegend()
list(plot3)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE,xnudge=0,ynudge=0)
plot4
#PCA 分析数据准备，使用 ScaleData()进行数据归一化；默认只是标准化高变基因（2000 个），速度更快，不影响 PCA 和分群，但影响热图的绘制。
#all_genes <- rownames(bat_total_4)
bat_total_5 <- ScaleData(bat_total_4, vars.to.regress = "percent.mt")
#而对所有基因进行标准化的方法如下：
#all.genes <- rownames(pbmc4)
#pbmc5 <- ScaleData(pbmc4,features = all.genes,vars.to.regress = "percent.mt")
#线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定；
bat_total_6 <- RunPCA(bat_total_5, features = VariableFeatures(object = bat_total_5))
# 检查 PCA 分群结果， 这里只展示前 12 个 PC,每个 PC 只显示 3 个基因；
print(bat_total_6[["pca"]], dims = 1:12, nfeatures = 3)
#
VizDimLoadings(bat_total_6, dims = 1:2, reduction = "pca")
#绘制 pca 散点图；
DimPlot(bat_total_6, reduction = "pca")
#画前 2 个主成分的热图；
DimHeatmap(bat_total_6, dims = 1:2, cells = 500, balanced = TRUE)
#方法 1：肘部图（碎石图），基于每个主成分对方差解释率的排名；
#识别数据集的真实维度 - 对用户来说可能具有挑战性/不确定性。因此，我们建议考虑这三种方法。第一个是更多的监督，探索PC以确定异质性的相关来源，并且可以与GSEA一起使用。第二个实现基于随机空模型的统计测试，但对于大型数据集来说是耗时的，并且可能不会返回明确的PC截止。第三种是常用的启发式算法，可以立即计算。在这个例子中，所有三种方法产生了类似的结果，但我们可能有理由在PC 7-12之间选择任何作为截止值的东西。
#我们在这里选择了10，但鼓励用户考虑以下内容：
#1.树突细胞和NK研究者可以认识到与PC12和13强烈相关的基因定义了罕见的免疫亚群（即MZB1是浆细胞样DC的标记）。然而，这些组是如此罕见，如果没有先验知识，很难将这种大小的数据集与背景噪声区分开来。
#2.我们鼓励用户使用不同数量的PC（10,15甚至50！）重复下游分析。正如您将观察到的，结果通常没有显着差异。
#3.在选择此参数时，我们建议用户在较高的一侧犯错。例如，仅使用5台PC进行下游分析确实会对结果造成严重影响。
ElbowPlot(bat_total_6)
#方法 2：Jackstraw 置换检验算法；重复取样（原数据的 1%），重跑 PCA,鉴定 p-value 较小的 PC；计算‘null distribution’(即零假设成立时)时的基因 scores;
#该JackStrawPlot功能提供了一种可视化工具，用于比较每个PC的p值分布和均匀分布（虚线）。“重要的”PC将显示具有低p值的特征的强烈丰富（虚线上方的实线）。在这种情况下，看起来在前10-12个PC之后，重要性急剧下降
#bat_total_7 <- JackStraw(bat_total_6, num.replicate = 100)
#bat_total_8 <- ScoreJackStraw(bat_total_7, dims = 1:20)
#JackStrawPlot(bat_total_8, dims = 1:20)
#基于 PCA 空间中的欧氏距离计算 nearest neighbor graph，优化任意两个细胞间的距离权重（输入上一步得到的 PC 维数）；

bat_total_9 <- FindNeighbors(bat_total_6, dims = 1:20)
save(bat_total_9, file = 'bat_total_9.rdata')
#接着优化模型，resolution 参数决定下游聚类分析得到的分群数，对于 3K 左右的细胞，设为 0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大，该参数也应该适当增大；
bat_total_10 <- FindClusters(bat_total_9, resolution = 1)
#使用 Idents（）函数可查看不同细胞的分群；
head(Idents(bat_total_10), 8)
#Seurat 提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起），比如 UMAP 和 t-SNE，运行 UMAP 需要先安装'umap-learn'包,这里不做介绍。
bat_total_11 <- RunTSNE(bat_total_10, dims = 1:20)
tsneplot<-TSNEPlot(bat_total_11,label = TRUE, pt.size = 1)+ NoLegend()
tsneplot
#用 DimPlot()函数绘制散点图， reduction = "tsne"，指定绘制类型；如果不指定，默认先从搜索 umap， 然后 tsne, 再然后 pca； 也可以直接使用这 3 个函数 PCAPlot()、 TSNEPlot()、UMAPPlot()； cols， pt.size 分别调整分组颜色和点的大小小
#library(umap)
bat_total_12 <- RunUMAP(bat_total_11, dims = 1:20) 
DimPlot(bat_total_12, reduction = "umap", label = TRUE, pt.size = 1)+ NoLegend()
DimPlot(bat_total_12, reduction = "tsne", label = TRUE, pt.size = 1)+ NoLegend()
#tsneplot<-TSNEPlot(bat_total_12,label = TRUE, pt.size = 1)+ NoLegend()

###分组绘图
Flox <- bat_total_12[,bat_total_12@meta.data[["orig.ident"]]=='Flox']
dim(Flox)
tsneplot_Flox<-TSNEPlot(Flox,label = TRUE, pt.size = 1,label.size = 6)+ NoLegend()
tsneplot_Flox
umapplot_Flox<-UMAPPlot(Flox,label = TRUE, pt.size = 1,label.size = 6)+ NoLegend()
umapplot_Flox


KO <- bat_total_12[,bat_total_12@meta.data[["orig.ident"]]=='KO']
dim(KO)
tsneplot_KO<-TSNEPlot(KO,label = TRUE, pt.size = 1,label.size = 6)+ NoLegend()
tsneplot_KO
umapplot_KO<-UMAPPlot(KO,label = TRUE, pt.size = 1)+ NoLegend()
umapplot_KO

saveRDS(bat_total_12,'clustered_res_1.Rds')
bat_total_12 <- clustered_res_1
###caculate diff cluster cell number and proportion -----
##all cell
table(bat_total_12@meta.data[["seurat_clusters"]])
a <- table(bat_total_12@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'cluster_number.csv')

###Flox cell number
table(Flox@meta.data[["seurat_clusters"]])
a <- table(Flox@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'Flox_cluster_number.csv')


###KO cell number
table(KO@meta.data[["seurat_clusters"]])
a <- table(KO@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'KO_cluster_number.csv')




### FindAllMarkers -----
##找到各亚群前十marker，进行鉴定
library(dplyr)
All_cluster.markers <- FindAllMarkers(object = bat_total_12, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(All_cluster.markers, file = 'all-cluster-diff.csv')
saveRDS(All_cluster.markers,'All_cluster_res1.markers.Rds')
Gene <- All_cluster.markers$gene
dim(All_cluster.markers)
top10 <- All_cluster.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top10,file = 'cluster_res1_diff_top10.csv')
DoHeatmap(object = bat_total_12, features = top10$gene, label = FALSE)
##### 从文献找marker 亚群鉴定-----
## cluster_20
cluster_20_top10 <- top10$gene[top10$cluster=='20']
FeaturePlot(object =  bat_total_12,reduction = 'tsne', features =cluster_20_top10 )
DotPlot(object =  bat_total_12, features =cluster_20_top10)
VlnPlot(object = bat_total_12,features = cluster_20_top10,pt.size = 0,adjust = 3)

# marker not adipo from paper RBB
library(readxl)
marker.from.paper. <- read_xlsx(path = 'marker-not-adipo.xlsx')
marker.from.paper. <- marker.from.paper.[2:21,]
marker.fb <- marker.from.paper.$FB
marker <- marker.fb
## macro certain
marker <- marker.from.paper.$MAC[1:10]
FeaturePlot(object =  bat_total_12,reduction = 'tsne', features = marker )+NoAxes()+
  theme(plot.title = element_text(size = 20, face = 'italic'))

#### marker_from article-adipo
marker <- 'Ppara'
marker <- 'Pparg'
marker <- 'Ucp1'
marker <- 'Prdm16'
marker <- 'Cidea'
## brown adipo
marker <- 'Pdk4'

## marker wat
marker <- 'Retn'
marker <- 'Tle3'
marker <- 'Tcf21'
marker <- 'Zfp423'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## progenitor 
marker <- 'En1'
marker <- 'Pax7'
marker <- 'Myf5'

## brown preadipo
marker <- 'Ebf2'
marker <- 'Pdgfra'
marker <- c('Pdgfra', 'Sca1', 'Cd34', 'Pref1', 'Cd29')


## macrophage
marker <- 'Adgre1'
marker <- 'Cd68'
marker <- 'Tnf'

marker <- marker.from.paper.$MAC[1:10]
marker <- marker.from.paper.$MAC[11:20]

## NKT
marker <- marker.from.paper.$NKT[1:10]
marker <- marker.from.paper.$NKT[11:20]

## VEC
marker <- marker.from.paper.$VEC[1:10]

## FB
marker <- marker.from.paper.$FB[1:10]
marker <- marker.from.paper.$FB[11:20]

## beige adipo
marker <- 'Pat2'
marker <- 'P2rx5'

## interested
marker <- 'Me1'

## marker from rensuping
rensu.marker <- read.csv('article/关注的基因.txt',header = F)

## 蛋白酶体系列 验证敲除
grep(pattern = '^Psm', rownames(bat_total_12))
rownames(bat_total_12)[grep(pattern = '^Psm', rownames(bat_total_12))]

psm.series <- rownames(bat_total_12)[grep(pattern = '^Psm', rownames(bat_total_12))]
pse.diff <- FindMarkers(bat_total_12, features = psm.series, ident.1 = 'KO', 
                        ident.2 = 'Flox', group.by = 'orig.ident',subset.ident = c(0,1,2,3,4,11,18))

marker <- rensu.marker$V1[1:8]
marker <- rensu.marker$V1[9:16]
marker <- rensu.marker$V1[17:24]
FeaturePlot(object =  bat_total_12,reduction = 'tsne', features = marker, label = TRUE, label.size = 5 )+NoAxes()+
  theme(plot.title = element_text(size = 20, face = 'italic'))

VlnPlot(object = bat_total_12,features = marker,pt.size = 0,adjust = 1)+NoLegend()+
  theme(plot.title = element_text(size = 20, face = 'italic'))
ggsave(filename = 'vln-mac-marker-from-paper-2.pdf', width = 20, height = 15)



##interested gene qpcr_test -----
# 炎症
marker <- c('Adgre1', 'Cd68')
# 坏死性凋亡 升高
marker <- c('Ripk1', 'Ripk3', 'Mlkl')
# 升高
marker <- c('Tnf', 'Casp1', 'Nlrp3', 'Pycard')
# 新生与蛋白分解 20周升高，4无
marker <- c('Myog', 'Desmin', 'Myh6'  )
# 4 20 60 全部下降
marker <- c('Fbxo32')
# 线粒体损伤 4无（第一个上升），60周上升
marker <- c('Cpt1a', 'Cpt1b')
# 代谢物改变 4无，20周上升
marker <- c( 'B4galt5', 'B4galt6')



#####给群加id-------
new.cluster.ids <-c()
names(new.cluster.ids) <- levels(bat_total_12)
bat_total_13 <- RenameIdents(bat_total_12, new.cluster.ids)
tsneplot<-TSNEPlot(bat_total_13,label = TRUE, pt.size = 1.5)+ NoLegend()
tsneplot
umapplot<-UMAPPlot(bat_total_13,label = TRUE, pt.size = 1.5)+ NoLegend()
umapplot
##
save(bat_total_13, file = 'gas_13_identified.rdata')
tsneplot<-TSNEPlot(bat_total_13,label = TRUE, pt.size = 1.5,
                   label.size=4.9,repel=TRUE)+ NoLegend()
tsneplot

### 给群加id方法二-----
head(bat_total_12@meta.data)
bat_total_12$cell.type <- NA
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(0,1,2,3,4,18,11)] <- 'Adipocytes'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(8,17)] <- 'MAC'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(14)] <- 'NKT'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(6, 10, 16)] <- 'VEC'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(5,7,9,19)] <- 'Preadipocytes'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(12)] <- 'Pericytes'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(13)] <- 'Myocytes'
bat_total_12$cell.type[bat_total_12$seurat_clusters %in% c(15)] <- 'NC'

head(bat_total_12@meta.data)


Idents(bat_total_12) <- 'cell.type'
DimPlot(bat_total_12, reduction = 'tsne', label = TRUE, pt.size = 1.5, label.size = 5)+NoLegend()+NoAxes()

###分组绘图
Flox <- bat_total_13[,bat_total_13@meta.data[["orig.ident"]]=='Flox']
dim(Flox)
tsneplot_Flox<-TSNEPlot(Flox,label = TRUE, pt.size = 1)+ NoLegend()
tsneplot_Flox
KO <- bat_total_13[,bat_total_13@meta.data[["orig.ident"]]=='KO']
dim(KO)
tsneplot_KO<-TSNEPlot(KO,label = TRUE, pt.size = 1)+ NoLegend()
tsneplot_KO

###caculate diff cluster cell number and proportion -----
##all cell
table(bat_total_13@meta.data[["seurat_clusters"]])
a <- table(bat_total_13@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'cluster_number.csv')

table(bat_total_13@active.ident)
a <- table(bat_total_13@active.ident)
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'cell_type_cluster_number.csv')
###Flox cell number
table(Flox@meta.data[["seurat_clusters"]])
a <- table(Flox@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'Flox_cluster_number.csv')

table(Flox@active.ident)
a <- table(Flox@active.ident)
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'Flox_cell_type_cluster_number.csv')
###KO cell number
table(KO@meta.data[["seurat_clusters"]])
a <- table(KO@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'KO_cluster_number.csv')

table(KO@active.ident)
a <- table(KO@active.ident)
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'KO_cell_type_cluster_number.csv')

save(Flox,KO,file = 'Flox_KO.Rdata')

### chiquare test

###save


#### 提取barcode 信息----------
library(tidyverse)
library(tidyr)
library(stringr)
x <- adipo.scrna@meta.data
head(x)
count(x,adipo.type)
colnames(x)
x <- x[,c('orig.ident', 'adipo.type')]
x$barcode <- rownames(x)
x <- separate(x, barcode, into = c('group','barcode'), sep = '_')
x <- x[,-1]
write.csv(x,file = 'subadipo_cluster_barcode.csv')
## 验证  barcode是否一致
y <- row.names(x[x$cell.type == 'Adipocytes',])
y2 <- rownames(bat_total_12@meta.data[bat_total_12@meta.data$cell.type == 'Adipocytes',])
identical(y,y2)

### 计算某基因原始reads的比例 --------
x <- bat_total_13@assays$RNA@counts
x <- as.matrix(x)
x[1:6,1:6]
x_Flox <- x[,bat_total_13@meta.data$orig.ident =='Flox']
x_KO <- x[,bat_total_13@meta.data$orig.ident =='KO']
yki <- rowSums(x_Flox)
yki['Nfe2l1']
yko <- rowSums(x_KO)
yko['Nfe2l1']
sum(yki)
sum(yko)
9294 / 41967944 #矫正
1830 / 34886178 #矫正
0.0002214547 / 5.245631e-05 #4.221698 ko 是ki的这么多倍，跟bam的reads对应上~。
# 我的bam的reads没有进行标准化（两组不可比），只进行了归一化（同组内，不同reads的缩放）

## 
DimPlot(bat_total_12, reduction = 'tsne')







