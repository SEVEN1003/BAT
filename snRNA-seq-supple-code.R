## package load ----
library(Seurat)
library(ggplot2)
library(tidyverse)
## fs1 cluster identified----
head(bat_total_12)
DimPlot(bat_total_12, reduction = 'tsne', group.by = 'seurat_clusters', label = TRUE)+
  NoLegend()
#### marker_from article-adipo
marker <- 'Ppara'
marker <- 'Pparg'
marker <- 'Ucp1'
marker <- 'Prdm16'
marker <- 'Cidea'
marker <- 'Adipoq'

FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5, label = FALSE)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 18),
                 legend.text = element_text(size = 15))
ggsave(filename = paste0('feature_',marker, '_.pdf'), width = 3.5, height = 2.8)


FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))


## brown preadipo
marker <- 'Ebf2'
marker <- 'Pdgfra'
marker <- 'Cd34'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## macrophage
marker <- 'Adgre1'
marker <- 'Cd68'
marker <- 'Tnf'
marker <- 'Ctss'
marker <- 'Ly86'
marker <- marker.from.paper.$MAC[1:10]
marker <- marker.from.paper.$MAC[11:20]
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## NKT
marker <- marker.from.paper.$NKT[1:10]
marker <- marker.from.paper.$NKT[11:20]
marker <- 'Skap1'
marker <- 'Gimap3'
marker <- "Il7r"
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## VEC
marker <- marker.from.paper.$VEC[1:10]
marker <-"Flt1"
marker <-"Ptprb"
marker <- "Cyyr1"
marker <- 'Nos3'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## PERICYTES/Smooth muscle cell SMC
marker <- 'Notch3'
marker <- 'Abcc9'
marker <- 'Kcnq5'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

### MYPCYTES
marker <- 'Ttn'
marker <- 'Atp2a1'
marker <- 'Tnnt3'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## nc nerve cell
marker <- 'Cadm2'
marker <- 'Csmd1'
marker <- 'Nkain2'
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## fs2 group validation----
## protesam violin
marker <- psm.diff.gene
marker <- 'Psmd1'
marker <- "Psmd11"
marker <- "Psmd14"
marker <- "Psmc6"
marker <- "Psmd4"
VlnPlot(bat_total_12, marker, split.by = 'orig.ident', pt.size = 0, cols = c('#078992', '#f8766d'))+
  labs(x = '', y= 'Expression level')+
  theme(plot.title = element_text(face = 'italic', size = 20), 
        axis.ticks.x = element_blank())
FeaturePlot(Flox, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))
FeaturePlot(KO, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## psm deg analysis
marker <- psm.series
psm.deg <- FindMarkers(bat_total_12, features = marker, ident.1 = 'KO', ident.2 = 'Flox', 
                       subset.ident = 'Adipocytes', group.by = 'orig.ident' )
write.csv(psm.deg, file = 'psm-adipo-diff.csv')

## Nfe2l1 
marker <- 'Nfe2l1'
VlnPlot(adipo.scrna, marker, split.by = 'orig.ident', pt.size = 0, cols = c('#078992', '#f8766d'))+
  labs(x = '', y= 'Expression level')+
  theme(plot.title = element_text(face = 'italic', size = 20), 
        axis.ticks.x = element_blank())

FeaturePlot(Flox, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))
FeaturePlot(KO, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

### t1 chi squre of all cell type propor-----
library(readxl)
cell.type.chi <- read_xlsx(path = 'cell_type_cluster_number.xlsx')
table(bat_total_12@active.ident, bat_total_12$orig.ident)
###写个循环做卡方
mytable <- cell.type.chi
for (x in 1:nrow(mytable)) {
  a=mytable$KO.Freq[x]
  b=sum(mytable$KO.Freq)-a
  c=mytable$Flox.Freq[x]
  d=sum(mytable$Flox.Freq)-c
  table_for_chi=matrix(c(a,b,c,d),nrow = 2,byrow = TRUE)
  e=chisq.test(table_for_chi)
  mytable$p_value[x]=e[["p.value"]]
  mytable$X_squared[x]=e[["statistic"]][["X-squared"]]
  
}
write.csv(mytable,'all-cell-type-pvalue.csv')
## tsne
DimPlot(bat_total_12, reduction = 'tsne', split.by = 'orig.ident', 
        pt.size = 1.5, label = TRUE, label.size = 5 )+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))+NoLegend()

## ADIPO SUbSet chi test----
adipo.type.chi <- read_xlsx(path = 'cluster_number.xlsx', sheet = 2)
mytable <- adipo.type.chi
mytable <- mytable[-nrow(mytable),]
mytable$Var1 <- paste0('A', mytable$Var1)
for (x in 1:nrow(mytable)) {
  a=mytable$KO.Freq[x]
  b=sum(mytable$KO.Freq)-a
  c=mytable$Flox.Freq[x]
  d=sum(mytable$Flox.Freq)-c
  table_for_chi=matrix(c(a,b,c,d),nrow = 2,byrow = TRUE)
  e=chisq.test(table_for_chi)
  mytable$p_value[x]=e[["p.value"]]
  mytable$X_squared[x]=e[["statistic"]][["X-squared"]]
  
}
write.csv(mytable,'ADIPO-type-pvalue.csv')

## a4 explore high concered TF----
marker <- 'Cebpa'
marker <- 'Cebpb'
marker <- 'Esrrg'
marker <- 'Pparg'
adipo.Flox <- subset(adipo.scrna, subset = orig.ident == 'Flox')
adipo.KO <- subset(adipo.scrna, subset = orig.ident == 'KO')
FeaturePlot(adipo.Flox, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))
FeaturePlot(adipo.KO, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))


VlnPlot(adipo.scrna, marker, split.by = 'orig.ident', pt.size = 0, cols = c('#078992', '#f8766d'),
        adjust = 2)+
  labs(x = '', y= 'Expression level')+
  theme(plot.title = element_text(face = 'italic', size = 20), 
        axis.ticks.x = element_blank())
DotPlot(adipo.scrna, features = marker)
dot.df <- DotPlot(adipo.scrna, features = marker)
dot.df$data


DotPlot(adipo.scrna, features = marker, split.by = 'orig.ident')
dot.df <- DotPlot(adipo.scrna, features = marker, split.by = 'orig.ident')
dot.df <- DotPlot(preadipo.scrna, features = marker, split.by = 'orig.ident')
dot.df$data
str(adipo.scrna$orig.ident)

marker <- 'Ucp1'
marker <- 'Ppargc1a'
# 
library(tidyverse)
library(tidyr)
library(stringr)
df <- DotPlot(adipo.scrna, features = marker, split.by = 'orig.ident')
df$data
df <- df$data
df <- df %>% select(1, 2, 3, 4) %>% separate(id, into = c('Adipo.type', 'Group'), sep = '_') 
###长列表变宽列表
df_KO <- df[df$Group == 'KO', ]
row.names(df_KO) <- df_KO$Adipo.type
df_Flox <- df[df$Group == 'Flox', ]
row.names(df_Flox) <- df_Flox$Adipo.type
df_adj <- cbind(df_Flox, df_KO)
head(df_adj)
df_adj <- df_adj[,c(1,2,6, 7)]
write.csv(df_adj, file = 'ppargc1a_expr.csv')

## pseud analysis data prepare
levels(adipo.scrna)
save(adipo.scrna, file = 'adipo-scrna.rdata')
## go to pseud file
### fs3 group deg a 0, 1, 2-----
deg.a0 <- FindMarkers(adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident',
                      subset.ident = 'A0')
write.csv(deg.a0, file = 'deg-a0-group.csv')
## remove nfe2l1 and plot heatmap
head(deg.a0)
str(deg.a0)
nrow(deg.a0)
deg.a0 <- deg.a0[rownames(deg.a0) != 'Nfe2l1', ]
deg.a0.sig <- deg.a0[deg.a0$p_val_adj < 0.05,]
deg.a0.sig <- deg.a0.sig[order(deg.a0.sig$avg_log2FC, decreasing = TRUE),]
a0.up.20 <- rownames(deg.a0.sig)[1:20]
nrow(deg.a0.sig)
a0.down.20 <- rownames(deg.a0.sig)[320:339]
a0.top20 <- c(a0.down.20, a0.up.20)
a0.scrna <- adipo.scrna[,adipo.scrna@active.ident =='A0']
table(a0.scrna@active.ident)
a0.scrna <- ScaleData(a0.scrna, features = a0.top20)
DoHeatmap(a0.scrna, features = a0.top20, group.by = 'orig.ident')+
  theme(axis.text.y = element_text(face = 'italic'),
        text = element_text(size = 15))+guides(color = FALSE )
# ## 尝试修改字体 yyplot, 先不安装了
# install.packages("qrcode")
# install.packages("ggimage")
# install.packages("yulab.utils")
# devtools::install_github("GuangchuangYu/yyplot")
# ## https://stackoverflow.com/questions/68794484/loading-magick-package-after-tidyverse
deg.a1 <- FindMarkers(adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident',
                      subset.ident = 'A1')
deg.a1 <- deg.a1[rownames(deg.a1) != 'Nfe2l1', ]
write.csv(deg.a1, file = 'deg-a1-group.csv')
## remove nfe2l1 and plot heatmap
head(deg.a1)
str(deg.a1)
nrow(deg.a1)
deg.a1.sig <- deg.a1[deg.a1$p_val_adj < 0.05,]
deg.a1.sig <- deg.a1.sig[order(deg.a1.sig$avg_log2FC, decreasing = TRUE),]
a1.up.20 <- rownames(deg.a1.sig)[1:20]
nrow(deg.a1.sig)
a1.down.20 <- rownames(deg.a1.sig)[727:746]
a1.top20 <- c(a1.down.20, a1.up.20)
a1.scrna <- adipo.scrna[,adipo.scrna@active.ident =='A1']
table(a1.scrna@active.ident)
a1.scrna <- ScaleData(a1.scrna, features = a1.top20)
DoHeatmap(a1.scrna, features = a1.top20, group.by = 'orig.ident')+
  theme(axis.text.y = element_text(face = 'italic'),
        text = element_text(size = 15))+guides(color = FALSE )
#rownames(adipo.scrna)[grep(pattern = 'TUSC3', rownames(adipo.scrna))]

deg.a2 <- FindMarkers(adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident',
                      subset.ident = 'A2')
deg.a2 <- deg.a2[rownames(deg.a2) != 'Nfe2l1', ]
write.csv(deg.a2, file = 'deg-a2-group.csv')
## remove nfe2l1 and plot heatmap
head(deg.a2)
str(deg.a2)
nrow(deg.a2)
deg.a2.sig <- deg.a2[deg.a2$p_val_adj < 0.05,]
deg.a2.sig <- deg.a2.sig[order(deg.a2.sig$avg_log2FC, decreasing = TRUE),]
a2.up.20 <- rownames(deg.a2.sig)[1:20]
nrow(deg.a2.sig)
a2.down.20 <- rownames(deg.a2.sig)[273:292]
a2.top20 <- c(a2.down.20, a2.up.20)
a2.scrna <- adipo.scrna[,adipo.scrna@active.ident =='A2']
table(a2.scrna@active.ident)
a2.scrna <- ScaleData(a2.scrna, features = a2.top20)
DoHeatmap(a2.scrna, features = a2.top20, group.by = 'orig.ident')+
  theme(axis.text.y = element_text(face = 'italic'),
        text = element_text(size = 15))+guides(color = FALSE )
#rownames(adipo.scrna)[grep(pattern = 'TUSC3', rownames(adipo.scrna))]


deg.a3 <- FindMarkers(adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident',
                      subset.ident = 'A3')
deg.a3 <- deg.a3[rownames(deg.a3) != 'Nfe2l1', ]
write.csv(deg.a3, file = 'deg-a3-group.csv')
## remove nfe2l1 and plot heatmap
head(deg.a3)
str(deg.a3)
nrow(deg.a3)
deg.a3.sig <- deg.a3[deg.a3$p_val_adj < 0.05,]
deg.a3.sig <- deg.a3.sig[order(deg.a3.sig$avg_log2FC, decreasing = TRUE),]
a3.up.20 <- rownames(deg.a3.sig)[1:20]
nrow(deg.a3.sig)
a3.down.20 <- rownames(deg.a3.sig)[134:153]
a3.top20 <- c(a3.down.20, a3.up.20)
a3.scrna <- adipo.scrna[,adipo.scrna@active.ident =='A3']
table(a3.scrna@active.ident)
a3.scrna <- ScaleData(a3.scrna, features = a3.top20)
DoHeatmap(a3.scrna, features = a3.top20, group.by = 'orig.ident')+
  theme(axis.text.y = element_text(face = 'italic'),
        text = element_text(size = 15))+guides(color = FALSE )
#rownames(adipo.scrna)[grep(pattern = 'TUSC3', rownames(adipo.scrna))]


deg.a4 <- FindMarkers(adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident',
                      subset.ident = 'A4')
deg.a4 <- deg.a4[rownames(deg.a4) != 'Nfe2l1', ]
write.csv(deg.a4, file = 'deg-a4-group.csv')
## remove nfe2l1 and plot heatmap
head(deg.a4)
str(deg.a4)
nrow(deg.a4)
deg.a4.sig <- deg.a4[deg.a4$p_val_adj < 0.05,]
deg.a4.sig <- deg.a4.sig[order(deg.a4.sig$avg_log2FC, decreasing = TRUE),]
a4.up.20 <- rownames(deg.a4.sig)[1:20]
nrow(deg.a4.sig)
a4.down.20 <- rownames(deg.a4.sig)[132:151]
a4.top20 <- c(a4.down.20, a4.up.20)
a4.scrna <- adipo.scrna[,adipo.scrna@active.ident =='A4']
table(a4.scrna@active.ident)
a4.scrna <- ScaleData(a4.scrna, features = a4.top20)
DoHeatmap(a4.scrna, features = a4.top20, group.by = 'orig.ident')+
  theme(axis.text.y = element_text(face = 'italic'),
        text = element_text(size = 15))+guides(color = FALSE )
#rownames(adipo.scrna)[grep(pattern = 'TUSC3', rownames(adipo.scrna))]


## figs data mining figure-----
## nature
m <- 'Cyp2e1'
m <- 'Slc7a10'
m <- 'Auts2'
m <- 'Retn'

FeaturePlot(adipo.scrna, features = m, reduction = 'tsne', pt.size = 1.5, label = TRUE)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))
## gse
m <- 'Iigp1'
m <- 'Xdh'
m <- 'Gbp7'
m <- 'Gbp2'
m <- 'Irf1'
m <- 'Ucp1'
m <- 'Cidea'
m <- 'Cox7a1'
m <- 'Lipe'
m <- 'Ctsb'

m <- NOD_A3
AverageExpression(seu.gse.qc.adip, features = m)
FeaturePlot(seu.gse.qc.adip, features = m, reduction = 'tsne', 
            pt.size = 1.5, label = TRUE)+NoAxes()+theme(plot.title  = element_text(face = 'italic'))

FeaturePlot(adipo.Flox, features = m, reduction = 'tsne', 
            pt.size = 1.5, label = TRUE)+NoAxes()+
  theme(plot.title  = element_text(face = 'italic')
  )

## 更改配色范围
FeaturePlot(adipo.Flox, features = m, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))+
  scale_color_gradient(limits = c(0, 2), # 数据上下限
                       breaks = c(0, 0.5, 1,1.5, 2), # 分段点
                       low = "grey", # 下限颜色
                       high = "blue") # 上限颜色



