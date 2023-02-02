DimPlot(bat_total_12, reduction = 'tsne')
DimPlot(bat_total_12)
DimPlot(adipo.Flox)
DimPlot(adipo.Flox, reduction = 'tsne')
DimPlot(adipo.KO)
## 
Cells(bat_total_12)
Cells(adipo.Flox)
head(Cells(adipo.Flox))
head(bat_total_12@meta.data)
head(adipo.Flox@meta.data)
flox <- bat_total_12[, bat_total_12@meta.data$orig.ident == 'Flox']
dim(flox)
table(flox@meta.data$orig.ident)

## use adipo.flox and adipo.type
## cell ID
write.csv(Cells(adipo.Flox), file = "cellID_obs.csv", row.names = FALSE)

## UMAP
write.csv(Embeddings(adipo.Flox, reduction = "tsne"), file = "cell_embeddings.csv")

## seurat clustering
write.csv(Idents(adipo.Flox), "clusters.csv")

## 
# Cells(bat_total_12)
head(bat_total_12@meta.data)
ko <- bat_total_12[, bat_total_12@meta.data$orig.ident == 'KO']
dim(ko)
table(ko@meta.data$orig.ident)

## 
## cell ID
write.csv(Cells(adipo.KO), file = "ko_cellID_obs.csv", row.names = FALSE)
## UMAP
write.csv(Embeddings(adipo.KO, reduction = "tsne"), file = "ko_cell_embeddings.csv")
## seurat clustering
write.csv(Idents(adipo.KO), "ko_clusters.csv")

## 调取tsne图颜色，用于速率图
adipo.tsne <- TSNEPlot(adipo.scrna)
adipo.tsne

data1 <- ggplot_build(adipo.tsne)$data[[1]]
data1$celltype <- adipo.scrna$adipo.type

cell_color = data1[,c('colour', 'celltype')]

cell_color = cell_color[!duplicated(cell_color$celltype),]

cell_color
cell_color$celltype <- factor(cell_color$celltype, 
                              levels = c("A0" , "A1" ,"A2" , "A3" , "A4","A11" ,"A18"  ))

cell_color <- cell_color[order(cell_color$celltype, decreasing = FALSE),]
cell_color$colour

