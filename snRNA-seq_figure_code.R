## package loading----
library(Seurat)
library(ggplot2)
library(tidyverse)

## figure 1a tsne----
DimPlot(bat_total_12, reduction = 'tsne',label = TRUE, label.size = 5,pt.size = 1.5)+NoLegend()+
  NoAxes()

marker <- 'Adipoq'
## feature plot
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))
FeaturePlot(bat_total_12, features = marker, reduction = 'tsne', pt.size = 1.5,cols = c('lightgrey', 'red'))+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))

## figure1b proportin barplot excel data-----
## all cell type
## caculate diff cluster cell number and proportion
##all cell
table(bat_total_12@meta.data[["seurat_clusters"]])
a <- table(bat_total_12@meta.data[["seurat_clusters"]])
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'cluster_number.csv')

table(bat_total_12@active.ident)
a <- table(bat_total_12@active.ident)
b <- data.frame(a)
b$proportion <- b$Freq/sum(b$Freq)
write.csv(b,file = 'cell_type_cluster_number.csv')

Flox <- bat_total_12[,bat_total_12$orig.ident == 'Flox']
KO <- bat_total_12[,bat_total_12$orig.ident == 'KO']
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

## 读入数据
library(readxl)
cell.type.barplot <- read_xlsx(path = 'cell_type_cluster_number.xlsx')
cell.type.barplot <- cell.type.barplot[,c(1,3,5,7)]
cell.type.barplot <- cell.type.barplot[,-2]
cell.type.barplot.flox <- cell.type.barplot[,c(1,2)]
colnames(cell.type.barplot.flox)[2] <- 'Freq'
cell.type.barplot.flox$Group <- 'Flox'

cell.type.barplot.ko <- cell.type.barplot[,c(1,3)]
colnames(cell.type.barplot.ko)[2] <- 'Freq'
cell.type.barplot.ko$Group <- 'KO'

cell.type.barplot.rbind <- rbind(cell.type.barplot.flox, cell.type.barplot.ko)
colnames(cell.type.barplot.rbind)[2] <- 'Proportion'
colnames(cell.type.barplot.rbind)[1] <- 'Cluster.number'
### 
cell.type.barplot.rbind$Cluster.number <- factor(cell.type.barplot.rbind$Cluster.number)
str(cell.type.barplot.rbind)
cell.type.barplot.rbind$Group <- factor(cell.type.barplot.rbind$Group, levels = c('KO', 'Flox'))
ggplot(data= cell.type.barplot.rbind, 
       aes(x = Group, y = Proportion, fill = Cluster.number)) +
  geom_bar(stat="identity")+theme_classic(base_size = 15)+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_fill_manual(values=c('#f8766d','#00a9ff','#00bfc4','#ff61cc', '#c77cff','#7cae00','#cd9600','#00be67'))+
  NoLegend()+coord_flip()
ggsave(filename = 'barplot_bat_横.pdf', width = 4, height = 1.7)

## f1b2 get adipo and preadipo----
colnames(bat_total_12@meta.data)
adipo.scrna <- subset(bat_total_12, subset = cell.type == 'Adipocytes')
preadipo.scrna <- subset(bat_total_12, subset = cell.type == 'Preadipocytes')

## adipo barplot
adipo.scrna$adipo.type <- NA
head(adipo.scrna@meta.data)
adipo.scrna$adipo.type <- paste0('A',adipo.scrna$seurat_clusters)
table(adipo.scrna$orig.ident, adipo.scrna$seurat_clusters)
Idents(adipo.scrna) <- 'adipo.type'
adipo.scrna$adipo.type <- factor(adipo.scrna$adipo.type, 
                                 levels = c('A0','A1','A11','A18', 'A2', 'A3', 'A4'))
adipo.scrna@active.ident<- factor(adipo.scrna@active.ident, 
                                 levels = c('A0','A1','A11','A18', 'A2', 'A3', 'A4'))

DimPlot(adipo.scrna, reduction = 'tsne',label = TRUE, label.size = 5,pt.size = 1.5, split.by = 'orig.ident')+NoLegend()+
  NoAxes()

adipo.bar <- table(adipo.scrna$adipo.type, adipo.scrna$orig.ident)
adipo.bar <- data.frame(adipo.bar)
colnames(adipo.bar)[1:2] <- c('Adipo.type', 'Group')
table(adipo.scrna$orig.ident)
adipo.bar$Propor[1:7] <- adipo.bar$Freq[1:7]/8931
adipo.bar$Propor[8:14] <- adipo.bar$Freq[8:14]/8268
colnames(adipo.bar)
##
adipo.bar$Group <- factor(adipo.bar$Group,levels = c('KO', 'Flox'))

ggplot(data = adipo.bar, aes(x = Group, y = Propor ,fill = Adipo.type))+
  geom_bar(stat = 'identity')+theme_classic(base_size = 15)+
  labs( y = 'Proportion', fill = 'Adipocyte\nsubclusters')+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(size = 15),
        text = element_text(size = 15),
        legend.position = 'bottom')+
  coord_flip()
ggsave(filename = 'barplot_adipo_横.pdf', width = 4.7,height = 2.8)
## preadipo bar
preadipo.scrna$preadipo.type <- NA
head(preadipo.scrna@meta.data)
preadipo.scrna$preadipo.type <- paste0('PA',preadipo.scrna$seurat_clusters)
table(preadipo.scrna$orig.ident, preadipo.scrna$seurat_clusters)
Idents(preadipo.scrna) <- 'preadipo.type'
DimPlot(preadipo.scrna, reduction = 'tsne',label = TRUE, label.size = 5,pt.size = 1.5)+NoLegend()+
  NoAxes()
preadipo.bar <- table(preadipo.scrna$preadipo.type, preadipo.scrna$orig.ident)
preadipo.bar <- data.frame(preadipo.bar)
colnames(preadipo.bar)[1:2] <- c('preadipo.type', 'Group')
table(preadipo.scrna$orig.ident)
preadipo.bar$Propor[1:4] <- preadipo.bar$Freq[1:4]/2140
preadipo.bar$Propor[5:8] <- preadipo.bar$Freq[5:8]/1496
colnames(preadipo.bar)
ggplot(data = preadipo.bar, aes(x = Group, y = Propor ,fill = preadipo.type))+
  geom_bar(stat = 'identity')+theme_minimal()+labs( y = 'Proportion', fill = 'Preadipo type ')+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(size = 15),
        text = element_text(size = 15))


## f1c heatmap of  cell marker-------
## 计算亚群鉴定后的marker
all.bat.cell.marker <- FindAllMarkers(bat_total_12, logfc.threshold = 0.25,
                                      only.pos = TRUE)
## 使用前5个marker 每种细胞抽样500个，绘制热图
head(all.bat.cell.marker)

top5 <- all.bat.cell.marker %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)
bat_total_12 <- ScaleData(bat_total_12, top5$gene)
bat.subset <- subset(bat_total_12, downsample = 500)

DoHeatmap(bat.subset, top5$gene, group.bar = FALSE)+
  theme(axis.text.y.left = element_text(face = 'italic', 
                                        size = 10,
                                        colour = 'black'))

ggsave(filename = 'heatmap-bat.pdf',width = 8,height = 20, units = 'cm')
levels(bat.subset)
## 4.48 5.55


## f1f a3 a4 top 10 heatmap----
## adipo subcluster top10
adipo.deg <- FindAllMarkers(object = adipo.scrna, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(adipo.deg, file = 'adipo-deg.csv')
colnames(adipo.deg)
adipo.top20 <- adipo.deg %>% group_by(cluster) %>% top_n(20, avg_log2FC)
adipo.top10 <- adipo.deg %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(adipo.top20, file = 'adipo-top20-deg.csv')

adipo.top10.gene <- adipo.top10$gene
adipo.scrna <- ScaleData(adipo.scrna, features = adipo.top10.gene)
## gene order
adipo.top10$cluster <- factor(adipo.top10$cluster, levels = c('A0','A1','A11','A18', 'A2', 'A3', 'A4'))
adipo.top10 <- adipo.top10[order(adipo.top10$cluster, decreasing = FALSE),]
adipo.top10.gene <- adipo.top10$gene

DoHeatmap(adipo.scrna, features = adipo.top10.gene)

## f1f去掉a18的热图-------
## 排序基因
adipo.top10 <- adipo.deg %>% group_by(cluster) %>% top_n(10, avg_log2FC)
adipo.top10$cluster <- factor(adipo.top10$cluster, levels = c('A0','A1','A2', 'A3', 'A4','A11','A18'))
adipo.top10 <- adipo.top10[order(adipo.top10$cluster, decreasing = FALSE),]
adipo.top10.gene <- adipo.top10$gene[adipo.top10$cluster != 'A18']

## 去除18号亚群
adipo.scrna <- ScaleData(adipo.scrna, features = adipo.top10.gene)
adipo.scrna.no.a18 <- subset(adipo.scrna, subset = adipo.type != 'A18')

levels(adipo.scrna.no.a18) <- c("A0" , "A1"  , "A2" , "A3"  ,"A4","A11")
## 每个亚群只取五百个展示
adipo.scrna.no.a18.downsampel <- subset(adipo.scrna.no.a18, downsample = 500 )

DoHeatmap(adipo.scrna.no.a18.downsampel, features = adipo.top10.gene, group.bar = FALSE)+
  theme(axis.text.y.left = element_text(face = 'italic', 
                                        size = 10,
                                        colour = 'black'))

ggsave(filename = 'heatmap-adipo.pdf',width = 4.4,height = 7.1)


## top 20 heatmap a3 and a4
adipo.a34.deg.gene <- adipo.top20$gene[adipo.top20$cluster %in% c('A3', 'A4')]
adipo.scrna <- ScaleData(adipo.scrna, features = adipo.a34.deg.gene)
adipo.a34 <- rownames(adipo.scrna@meta.data[adipo.scrna$adipo.type %in% c('A3', 'A4'),])
DoHeatmap(adipo.scrna, features = adipo.a34.deg.gene, cells = adipo.a34 )+
 

## a3 a4 go and kegg ----
a3.all.marker <- adipo.deg$gene[adipo.deg$cluster == 'A3']
a4.all.marker <- adipo.deg$gene[adipo.deg$cluster =='A4']
## a0 deg 富集分析
a0.group.down <- a0.adipo.diff$gene[which(a0.adipo.diff$p_val_adj < 0.05 & a0.adipo.diff$avg_log2FC < 0)]
a0.group.up <- a0.adipo.diff$gene[which(a0.adipo.diff$p_val_adj < 0.05 & a0.adipo.diff$avg_log2FC > 0)]
## a1 deg 富集分析
a1.group.down <- a1.adipo.diff$gene[which(a1.adipo.diff$p_val_adj < 0.05 & a1.adipo.diff$avg_log2FC < 0)]
a1.group.up <- a1.adipo.diff$gene[which(a1.adipo.diff$p_val_adj < 0.05 & a1.adipo.diff$avg_log2FC > 0)]
## a2 deg 富集分析
a2.group.down <- a2.adipo.diff$gene[which(a2.adipo.diff$p_val_adj < 0.05 & a2.adipo.diff$avg_log2FC < 0)]
a2.group.up <- a2.adipo.diff$gene[which(a2.adipo.diff$p_val_adj < 0.05 & a2.adipo.diff$avg_log2FC > 0)]
## all deg 
all.group.down <- all.adipo.diff$gene[which(all.adipo.diff$p_val_adj < 0.05 & all.adipo.diff$avg_log2FC < 0)]
all.group.up <- all.adipo.diff$gene[which(all.adipo.diff$p_val_adj < 0.05 & all.adipo.diff$avg_log2FC > 0)]
## deg 012
head(deg.a012)
dim(deg.a012)
deg.a012 <- rename(deg.a012, 'gene' = 'X')

down.gene <- deg.a012$gene[deg.a012$p_val_adj < 0.05 & deg.a012$avg_log2FC < 0]
up.gene <- deg.a012$gene[deg.a012$p_val_adj < 0.05 & deg.a012$avg_log2FC >0]

library(clusterProfiler)
library(org.Mm.eg.db)

# dir.create('enrich-scrna-deg-between-group')
gene.for.go <- all.group.down
gene.for.go <- all.group.up
gene.for.go <- a0.group.down
gene.for.go <- a0.group.up
gene.for.go <- a2.group.down
gene.for.go <- a2.group.up

gene.for.go <- down.gene
gene.for.go <- up.gene[up.gene !='Nfe2l1']
# ego_ALL <- enrichGO(gene          = gene.for.go,
#                     #universe     = ,
#                     OrgDb         = 'org.Mm.eg.db',
#                     keyType       = 'SYMBOL',
#                     ont           = "ALL",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.05,
#                     qvalueCutoff  = 0.2)
# ego_all <- data.frame(ego_ALL)
# write.csv(ego_all,'enrich-scrna-deg/enrichGO-all-adip-a4.csv')           
ego_BP <- enrichGO(gene          = gene.for.go,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,## 矫正p
                   qvalueCutoff  = 0.2)
ego_BP.df <- data.frame(ego_BP)
## up
write.csv(ego_BP.df, file = 'enrich-scrna-deg-between-group/a012-up-bp.csv')
## down
write.csv(ego_BP.df, file = 'enrich-scrna-deg-between-group/a012-down-bp.csv')


# ## 这里可以对term的长度进行修改，使图片更好看，遍历没一个term,取其1到70的字符（不超过35个字母，包括空格）
# ego_BP@result$Description <- substring(ego_BP@result$Description,1,40)
# ## 这里绘图，注意参数,另外他支持ggplot2的函数
# p_BP <- barplot(ego_BP,showCategory = 20) + 
#   ggtitle("Biological process of down-regulated genes")+
#   labs(x = 'Gene numbers')
# p_BP

## marker基因kegg富集分析
## 这一部,是必须的，要转化成entrezid。
library(tidyverse)
genelist <- bitr(gene.for.go, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')
#write.csv(genelist, file = 'enrich-scrna-deg/a4-gene-entreid.csv')
genelist.ENTREZID <- pull(genelist,ENTREZID)  
## 注意默认参数，最后只剩下两条, 默认参数p 0.05 ,q 0.2
ekegg <- enrichKEGG(gene = genelist.ENTREZID, organism = 'mmu')
ekegg.loose <- enrichKEGG(gene = genelist.ENTREZID, organism = 'mmu', pvalueCutoff = 1, qvalueCutoff = 1)
ekegg.loose.df <- data.frame(ekegg.loose)
## 尝试把基因ID变成gene symbol 之后再做吧
ekegg.loose.df$gene.symbol <- NA
for (i in 1:nrow(ekegg.loose.df)) {
  gene.entreid <- str_split(ekegg.loose.df$geneID[i], pattern = '/')[[1]]
  gene.entyid.symble <- bitr(gene.entreid, fromType="ENTREZID",toType="SYMBOL", OrgDb='org.Mm.eg.db')
  gene.symble <- pull(gene.entyid.symble, SYMBOL)
  gene.symble <-  paste(gene.symble, collapse = '/')
  ## 循环跑不通时，就要思考了，一步一步运行，这里由于gene symbol是个多元素的向量，转到
  ## 数据框中的一个元素时，就会出错，将其转化成单元素即可
  ekegg.loose.df$gene.symbol[i] <- gene.symble
}
## up
write.csv(ekegg.loose.df, file = 'enrich-scrna-deg-between-group/a012-up-kegg-add-symbol.csv')
## down
write.csv(ekegg.loose.df, file = 'enrich-scrna-deg-between-group/a012-down-kegg-add-symbol.csv')

## f2画图KEGG 
## dotplot 
library(ggplot2)
library(tidyverse)
rt <- read.csv(file = 'enrich-scrna-deg/down-bp.csv')
rt <- ekegg.loose.df
rt <- ego_BP.df
#
keggSig = rt[rt$p.adjust < 0.05,]
rownames(keggSig) <- 1:nrow(keggSig)
# # 前20条
keggSig <- keggSig[1:20,]

# # 前10条
# keggSig <- keggSig[1:10,]
# # 关注的通路选择
# keggSig <- keggSig[1,]
# str(keggSig)
## 写个循环 得到fold.enrich
for (i in 1:nrow(keggSig)) {
  gene.number <- str_split(keggSig$GeneRatio[i], pattern  = '/', simplify = TRUE)[1]
  gene.all <- str_split(keggSig$GeneRatio[i], pattern  = '/', simplify = TRUE)[2]
  
  gene.ratio.value <- as.numeric(gene.number)/as.numeric(gene.all)
  
  bg.number <- str_split(keggSig$BgRatio[i], pattern  = '/', simplify = TRUE)[1]
  bg.all <- str_split(keggSig$BgRatio[i], pattern  = '/', simplify = TRUE)[2]
  
  bg.ratio.value <- as.numeric(bg.number)/as.numeric(bg.all)
  keggSig$fold.enrich[i] <- gene.ratio.value/bg.ratio.value
}
###  改fold enrich的范围, 最大fe为10
# for (i in 1:nrow(keggSig)) {
#   ifelse(keggSig$fold.enrich[i] < 10, 
#          keggSig$fold.enrich[i] <- keggSig$fold.enrich[i], keggSig$fold.enrich[i] <- 10)
# }

##构建顺序因子
##先按fc排序
keggSig <- keggSig[order(keggSig$fold.enrich, decreasing = TRUE), ]
##改变rowname,
row.names(keggSig) <- seq(1:length(row.names(keggSig)))
## 截取长度
keggSig$Description_short <- substr(keggSig$Description, 1, 40)
## bp首字母没大写
keggSig$Description_short <- str_c(str_to_upper(str_sub(keggSig$Description_short, 1,1)), 
                                   str_sub(keggSig$Description_short, 2, str_length(keggSig$Description_short)))


porder = factor(as.integer(rownames(keggSig)),labels=rev(keggSig$Description_short))

pathbar = ggplot(keggSig,aes(x = Description_short,y = fold.enrich))

## 4.75x3.59
pathbar+ geom_point(stat="identity",
                    aes(x=rev(porder),
                        color=-log10(p.adjust), size=Count))+
  scale_colour_gradient(low="#0000ff",high="#f8766d")+
  labs(
    color=expression(-log[10](FDR)),
    size="Gene numbers",
    y="Fold enrichment"
    # x="Pathway name",
    #title="KEGG enrichment"
  )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3), colour = 'black'),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = 'black', size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )+coord_flip()

ggsave(filename = 'down-kegg-a012.pdf', width = 15, height = 10, units = 'cm')
## 另一种排序
### 按基因数目排序，点的颜色为p.adj,点的大小为fc
##构建顺序因子
##先按count排序
keggSig <- keggSig[order(keggSig$Count, decreasing = TRUE), ]
##改变rowname,
row.names(keggSig) <- seq(1:length(row.names(keggSig)))
porder = factor(as.integer(rownames(keggSig)),labels=rev(keggSig$Description))
colnames(keggSig)
pathbar = ggplot(keggSig,aes(x=Description,y=Count))



##绘制气泡图, 成功啦，果然以后用别人的代码修改，必须新建一个脚本，一句一句读
pathbar+ 
  geom_point(stat="identity",aes(x=rev(porder),color=-log10(p.adjust), size=fold.enrich))+
  scale_colour_gradient(low="#0000ff",high="#ff0000")+
  labs(
    color=expression(-log[10](p.adjust)),
    size="Fold enrichment",
    y="Gene numbers"
    # x="Pathway name",
    #title="KEGG enrichment"
    )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )+coord_flip()
### 画到下调的bp了。

### 验证

psm.diff <- FindMarkers(adipo.scrna, features = psm.series, ident.1 = 'KO', ident.2 = 'Flox', group.by = 'orig.ident')
psm.diff.gene <- rownames(psm.diff)
marker <- psm.diff.gene
DotPlot(adipo.scrna, features = psm.diff.gene, group.by = 'orig.ident')
DimPlot(adipo.scrna, reduction = 'tsne', split.by = 'orig.ident')
marker <- 'Nfe2l1'
marker <- 'Adipoq'
marker <- ''
marker <- 'Ppara'
marker <- 'Pparg'
marker <- 'Ucp1'
marker <- 'Prdm16'
marker <- 'Cidea'
DotPlot(adipo.scrna, features = marker , group.by = 'orig.ident')
FeaturePlot(adipo.scrna, features = marker, split.by = 'orig.ident', reduction = 'tsne')
VlnPlot(adipo.scrna, marker, split.by = 'orig.ident')
DotPlot(preadipo.scrna, features = marker , group.by = 'orig.ident')
FeaturePlot(preadipo.scrna, features = marker, split.by = 'orig.ident', reduction = 'tsne')
VlnPlot(preadipo.scrna, marker, split.by = 'orig.ident')

dot.df <- DotPlot(adipo.scrna, features = marker , group.by = 'orig.ident')
dot.df$data

dot.df <- DotPlot(preadipo.scrna, features = marker , group.by = 'orig.ident')
dot.df$data
dot.df <- DotPlot(bat_total_12, features = marker , group.by = 'orig.ident')

FeaturePlot(bat_total_12, features = marker, split.by = 'orig.ident', reduction = 'tsne')

## f2 差异分析A0 1 2 KO/FLOX-----
## 
head(bat_total_12@meta.data)
levels(bat_total_12)
levels(adipo.scrna)
deg.a012 <- FindMarkers(object = adipo.scrna, ident.1 = 'KO', ident.2 = 'Flox',
                        subset.ident = c('A0', 'A1', 'A2'), group.by = 'orig.ident')

write.csv(deg.a012, 'deg.012.csv')
deg.a012 <- read.csv('deg.012.csv')
head(deg.a012)

## 进行差异基因富集分析

## 蛋白酶体评分------
## 添加到list
genelist <- list()

genelist$Proteasome <- deg.a012$gene[grepl('^Psm', deg.a012$gene) & deg.a012$p_val_adj < 0.05]
genelist
## nod
NOD <- c('Stat2/Nampt/Gbp7/Sugt1/Ctsb')
NOD <- str_split(NOD, pattern =   '/') %>% unlist()

genelist <- list()
genelist$NOD <- NOD

## Thermogenesis 
## Acsl1/Gnas/Mgll/Cox8b/Slc25a20 是否可以用更多基因？
Thermogenesis <- c('Acsl1/Gnas/Mgll/Cox8b/Slc25a20')
Thermogenesis <- str_split(Thermogenesis, pattern =   '/') %>% unlist()

genelist <- list()
genelist$Thermogenesis <- Thermogenesis
## Thermogenesis all gene 
## kegg rest for gene-------
library(KEGGREST)
listDatabases()
## 获得某一条
## human hsa mouse mmu 
## KEGG_Thermogenesis bing搜索后点击Pathway menu 查看编号 | Organism menu查看物种
gs <- keggGet('mmu04714')## 产热
gs <- keggGet('mmu04621')## nod 

#获取通路中gene信息 
# head(gs[[1]]$GENE)
# # hsa 查找所有基因 
# # genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
# # genelist <- genes[1:length(genes)%%3 ==2] 

## mmu

genes <- gs[[1]]$GENE 
genes <- genes[1:length(genes) %% 2 == 0]
genes <- str_split(genes, pattern = ';',simplify = TRUE)[,1]
genes <- sort(genes)


Thermogenesis_all <- genes
NOD_all <- genes
genelist <- list(Thermogenesis_all = Thermogenesis_all,
                 NOD_all = NOD_all)
head(rownames(adipo.scrna))
intersect(rownames(adipo.scrna), Thermogenesis_all)
length(Thermogenesis_all)
## a012 up
fatty_meta <- c('Acsl1/Acaca/Hacd2/Elovl6/Ehhadh/Acox1')
fatty_meta <- str_split(fatty_meta, '/') %>% unlist()
fatty_meta
## a3 up 
NOD_A3 <- c('Gbp7/Gbp2/Gbp3/Stat1/Gbp5/Stat2/Nampt/Nod1/Txn1/Ctsb/Irf9/Sugt1')
NOD_A3 <- str_split(NOD_A3, '/') %>% unlist()
NOD_A3

genelist <- list(fatty_meta = fatty_meta,
                 NOD_A3 = NOD_A3)
## more score-------
## A3 necro and antigen 评分
## A4 ppar 通路评分
## A012 AMPK 评分
##  基因集准备见绘图代码汇总
head(genelist)

## 基因集的准备见绘图代码汇总
names(genelist)

adipo.scrna <-  AddModuleScore(adipo.scrna, 
                               features = genelist, ctrl = 500, seed = 1)
head(adipo.scrna@meta.data)


## 修改名称
names(genelist)
adipo.scrna@meta.data <- rename(adipo.scrna@meta.data, 
                                Propanoate_meta = Cluster1,
                                Propanoate_meta.all = Cluster2
                                )
head(adipo.scrna@meta.data)

# FeaturePlot(object =  adipo.scrna, split.by = 'orig.ident',reduction = 'tsne',
#             features = 'Proteasome',cols = c('lightgrey', 'red'),
#             pt.size =1.5 )+NoAxes()+
#   theme(plot.title = element_text(face = 'italic'))
# f2b  去掉离群点画箱式图--------
metadata <- adipo.scrna@meta.data
str(metadata$orig.ident)
head(metadata)
metadata$orig.ident <- factor(metadata$orig.ident, levels = c('KO', 'Flox') )
metadata <- metadata %>% filter(adipo.type != 'A18')

## Proteasome Fatty_acid_metabolism NOD_A3 Thermogenesis_all
## APP Necroptosis PPAR AMPK PPAR_all AMPK_all Propanoate metabolism
## base size 指的是坐标标题的size
colnames(metadata)
ggplot(metadata, aes(y= adipo.type, x=Propanoate_meta.all))+
  geom_boxplot(aes(color = orig.ident), outlier.colour = NA)+theme_classic(base_size = 12)+
  labs(title = 'Propanoate metabolism (all) score', x = 'Gene set score', y= 'Adipocyte subcluster', color = 'Group')+
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.ticks.y = element_blank(),
        axis.text =  element_text(color = 'black', size = 10),
        axis.title = element_text(size = 12))+
  scale_color_manual(values =  c("#f8766d","#0000ff"))

#  
ggsave(filename = 'f3-boxplot-Propanoate_meta_all.pdf',width = 3.0, height = 3.31)

## 控制坐标轴的范围~
ggplot(metadata, aes(y= adipo.type, x=Thermogenesis_all))+
  geom_boxplot(aes(color = orig.ident), outlier.colour = NA)+theme_classic(base_size = 12)+
  labs(title = 'Thermogenesis score', x = 'Gene set score', y= 'Adipocyte subcluster', color = 'Group')+
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.ticks.y = element_blank(),
        axis.text =  element_text(color = 'black', size = 10),
        axis.title = element_text(size = 12))+
  scale_color_manual(values =  c("#f8766d","#0000ff"))+
  xlim(c(0, 0.5))
ggsave(filename = 'f2b-boxplot-therm.pdf',width = 3.0, height = 3.31)

## f2e 散点图------
## 
head(metadata)

avg.score <- metadata %>% group_by(group.cluster) %>% 
  summarise(nod.avg = median(NOD),
            pro.avg = median(Proteasome),
            them.avg = median(Thermogenesis),
            ada_them.avg = median(Adaptive_thermogenesis),
            them.all.avg = median(Thermogenesis_all),
            nod.all.avg = median(NOD_all),
            nod.a3.avg = median(NOD_A3),
            fatty_meta.avg = median(Fatty_acid_metabolism),
            app.avg = median(APP),
            Necroptosis.avg = median(Necroptosis),
            PPAR.avg = median(PPAR),
            PPAR_all.avg = median(PPAR_all),
            AMPK.avg = median(AMPK),
            AMPK_all.avg = median(AMPK_all),
            Proteasome.all = median(Proteasome.all),
            APP.all = median(APP.all), 
            Necroptosis.all = median(Necroptosis.all),
            Propanoate_meta = median(Propanoate_meta),
            Propanoate_meta.all = median(Propanoate_meta.all))


avg.score <- avg.score %>% separate(col = group.cluster, into = c('Group', 'Subcluster'), sep = '_')
head(avg.score)

library(ggrepel)
## y = protea, x = noda3; them.all.avg, fatty_meta.avg, xintercept
## 
colnames(avg.score)

ggplot(avg.score, 
       aes(y = pro.avg, x = Propanoate_meta, 
           color = Group, label = Subcluster))+theme_bw(base_size = 12)+
  geom_point()+
  geom_text_repel()+
  geom_vline(xintercept = 0.4, linetype = "dashed",alpha = 0.5)+
  geom_hline(yintercept = 0,linetype = "dashed", alpha = 0.5)+
  labs(y = 'Proteasome score (median)', x = 'Propanoate metabolism score (median)')+
  scale_color_manual(values =  c("#0000ff","#f8766d"))+
  theme(axis.text = element_text(color = 'black', size = 10))
## 
ggsave(filename = 'f3g-scatter-Propanoate_meta.pdf', width = 3.91, height = 3.26)
## 使用ampk Propanoate metabolism Thermogenesis

head(avg.score)
colnames(avg.score)

ggplot(avg.score, 
       aes(y = AMPK.avg, x = Propanoate_meta.all, 
           color = Group, label = Subcluster))+theme_bw(base_size = 12)+
  geom_point()+
  geom_text_repel()+
  geom_vline(xintercept = 0.45, linetype = "dashed",alpha = 0.5)+
  geom_hline(yintercept = 0.8,linetype = "dashed", alpha = 0.5)+
  labs(y = 'AMPK score (median)', x = 'Propanoate metabolism (all) score (median)')+
  scale_color_manual(values =  c("#0000ff","#f8766d"))+
  theme(axis.text = element_text(color = 'black', size = 10))

#
ggsave(filename = 'f4-scatter--AMPK-prop-meta-all.pdf', width = 3.91, height = 3.26)

## 
## 
# Ctsz/Ctsb/Npc1
m <- 'Ctsb'
m <- 'Ctsz'
m <- 'Npc1'
m <- 'Gadd45g'
VlnPlot(adipo.scrna, m, split.by = 'orig.ident', pt.size = 0)
DotPlot(adipo.scrna,features = m)
FeaturePlot(adipo.scrna, m, reduction = 'tsne', split.by = 'orig.ident')

dot.df <- DotPlot(adipo.scrna,features = m, split.by = 'orig.ident')
dot.df$data

## f2 fam 气泡图--------
fam <- c('Acsl1/Acaca/Hacd2/Elovl6/Ehhadh/Acox1')
fam <- str_split(fam, '/') %>% unlist()
gene.dot <- c('Cyp2e1','Slc7a10','Auts2','Retn', fam)
## ampk
ampk <- c('Pparg/Igf1r/Ppp2r3a/Eef2k/Camkk2/Prkag2')
ampk <- str_split(ampk,'/') %>% unlist()
gene.dot <- ampk
AverageExpression(adipo.scrna, features = gene.dot)
DotPlot(adipo.scrna,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks.x =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

DotPlot(seu.rt.qc,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks.x =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggsave(filename = 'f2h-data-mine-adipo-dot-ampk.pdf', width = 4.3, height = 3.3)
ggsave(filename = 'f2h-data-mine-nature-dot-ampk.pdf', width = 4.3, height = 3.3)


df <- DotPlot(adipo.scrna,features = gene.dot,split.by = 'orig.ident')+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
df
df$data








## 使用自己画的小提琴图，它的太蠢了，展示各亚群的阳性比例和阳性的表达提琴

## 阳性比例条图------
df <- dot.df$data
df <- df %>% 
  dplyr::select(1:4) %>% 
  separate(col = 'id', into = c('subcluster', 'group')) %>% 
  mutate(pct.exp = pct.exp/100,
         pct.noexp = 1-pct.exp)
df <- df %>% filter(subcluster != 'A18')
df




ggplot(df, aes(x = subcluster, y = pct.exp))+
  geom_bar(stat = 'identity', aes(fill = group),alpha = .85,position = position_dodge())+
  theme_classic(base_size = 12)+
  scale_fill_manual(values = c("#0000ff","#f8766d"))+
  labs(title = 'Percentage of Ctsz-positive cells',
       x = '', y = 'Percentage', fill = 'Group')+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks.x = element_blank(),
        title = element_text(size = 12))

# table(adipo.scrna@meta.data$orig.ident, adipo.scrna$adipo.type)
## 继续画 阳性表达的箱式图
head(adipo.scrna@meta.data)
## 
ctsb <- FetchData(adipo.scrna, vars = c('orig.ident','adipo.type','Ctsb'))
head(ctsb)
ctsb <- ctsb %>% filter(Ctsb > 0,adipo.type !='A18' )
str(ctsb)
ggplot(ctsb, aes(x = adipo.type, y = Ctsb, fill = orig.ident))+
  theme_classic(base_size = 12)+
  geom_boxplot(outlier.shape =  NA, alpha = 0.85, color = 'black')+
  labs(x = '', y = 'Normalized  expression', fill = 'Group', title = 'Ctsb')+
  theme(plot.title = element_text(face = 'italic',hjust = 0.5),
        axis.text = element_text(colour = 'black',size = 10),
        axis.ticks.x =  element_blank())+
  scale_fill_manual(values = c("#0000ff","#f8766d"))
## ctsz
ctsb <- FetchData(adipo.scrna, vars = c('orig.ident','adipo.type','Ctsz'))
head(ctsb)
ctsb <- ctsb %>% filter(Ctsz > 0,adipo.type !='A18' )
str(ctsb)
ggplot(ctsb, aes(x = adipo.type, y = Ctsz, fill = orig.ident))+
  theme_classic(base_size = 12)+
  geom_boxplot(outlier.shape =  NA, alpha = 0.85, color = 'black')+
  labs(x = '', y = 'Normalized  expression', fill = 'Group', title = 'Ctsz')+
  theme(plot.title = element_text(face = 'italic',hjust = 0.5),
        axis.text = element_text(colour = 'black',size = 10),
        axis.ticks.x =  element_blank())+
  scale_fill_manual(values = c("#0000ff","#f8766d"))



## flox A3/012-----
levels(Flox)
levels(adipo.Flox)
head(adipo.Flox@meta.data)
flox.a3_v_012 <- FindMarkers(adipo.Flox, ident.1 = 'A3', 
                             ident.2 = c('A0', 'A1', 'A2'), 
                             group.by = 'adipo.type', logfc.threshold = 0.25
                            )
write.csv(flox.a3_v_012, file = 'flox.3v012.csv')
dim(flox.a3_v_012)
## 下调前20的气泡
down.gene
DotPlot(adipo.Flox, features = )

## 富集分析
## flox合并ko a3/a012
a3_v_012 <- FindMarkers(adipo.scrna, ident.1 = 'A3', 
                             ident.2 = c('A0', 'A1', 'A2'), 
                             group.by = 'adipo.type', logfc.threshold = 0.25)

write.csv(a3_v_012, file = 'a3_v_012.csv')
a3_v_012$gene <- rownames(a3_v_012)





## koA4/012
ko.a4_v_012 <- FindMarkers(adipo.KO, ident.1 = 'A4', 
                             ident.2 = c('A0', 'A1', 'A2'), 
                             group.by = 'adipo.type', logfc.threshold = 0.25)

ko.a4_v_012$gene <- rownames(ko.a4_v_012)

m <- 'Grina'

df <- DotPlot(adipo.scrna, features = m, split.by = 'orig.ident')
df$data

## koa4/ flox 012
head(adipo.scrna@meta.data)
adipo.scrna$group.cluster <- paste0(adipo.scrna$orig.ident,'_', adipo.scrna$adipo.type)
table(adipo.scrna$group.cluster)

ko.a4_v_flox_012 <- FindMarkers(adipo.scrna, ident.1 = 'KO_A4', 
                           ident.2 = c('Flox_A0', 'Flox_A1', 'Flox_A2'), 
                           group.by = 'group.cluster', logfc.threshold = 0.25)

ko.a4_v_flox_012$gene <- rownames(ko.a4_v_flox_012)

## a4/a012
a4_v_012 <- FindMarkers(adipo.scrna, ident.1 = 'A4', 
                           ident.2 = c('A0', 'A1', 'A2'), 
                           group.by = 'adipo.type', logfc.threshold = 0.25)

a4_v_012$gene <- rownames(a4_v_012)

## preadipo analys------
## pparg
preadipo.scrna
head(preadipo.scrna@meta.data)
levels(preadipo.scrna)
m <- 'Pparg'
m <- 'Psmd1'
df <- DotPlot(preadipo.scrna, features = m, split.by = 'orig.ident')
df$data

str(preadipo.scrna@meta.data)
preadipo.scrna$orig.ident <- factor(preadipo.scrna$orig.ident, levels = c('Flox', 'KO'))

VlnPlot(preadipo.scrna, features = m, group.by = 'orig.ident', pt.size = 0)+
  scale_fill_manual(values = c('#00bfc4','#f8766d'))
VlnPlot(preadipo.scrna, features = m, split.by ='orig.ident' , pt.size = 0)+
  scale_fill_manual(values = c('#00bfc4','#f8766d'))





## f3 a3 a4 a11 分析------
## 速率分析加组内差异基因富集分析+数据挖掘

## paper figure caption---
dim(bat_total_12)
table(bat_total_12$orig.ident)
median(bat_total_12$nFeature_RNA)

## f 数据挖掘气泡图
gene.dot <- c('Iigp1','Xdh','Irf1','Gbp7','Gbp2','Fabp4', 'Ucp1', 'Cidea','Cox7a1','Ndufa12' )
AverageExpression(adipo.scrna, features = gene.dot)
DotPlot(adipo.scrna,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks.x =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

DotPlot(adipo.Flox,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks.x =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggsave(filename = 'f3-data-mine-flox-A3-dot.pdf', width = 4.3, height = 3.3)

DotPlot(seu.gse.qc.adip,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks.x =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggsave(filename = 'f3-data-mine-BAL-dot.pdf', width = 4.3, height = 3.3)


df <- DotPlot(adipo.scrna,features = gene.dot)+coord_flip()+
  labs(x = '', y = '')+
  theme(axis.text.y = element_text(face = 'italic'),
        axis.ticks =  element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
df
df$data

## 二代测序 四周分析------
## 蛋白酶体和ampk gsea分析
## 读入差异分析结果
library(tidyverse)
library(readxl)
bulk_4_week <- read_xlsx(path = 'bulk_4week/Gene_expr_with_pval.xlsx')
head(bulk_4_week)
bulk_4_week

## 读入两个样本的组间分析
bulk_4_week_2v2 <- read.csv('bulk_4week/bulk_deg_2v2_deseq2.csv')
head(bulk_4_week_2v2)

## venn 图
## 上，下调与蛋白酶体
Proteasome.all
bulk_down <- bulk_4_week$Geneid[bulk_4_week$PValue < 0.05 & bulk_4_week$logFC < -0.25]
intersect(Proteasome.all, bulk_down)
## 还是使用gsea吧

## bulk 14 周分析报告------
bulk_14 <- read.xlsx('bulk_14week/deg_14wks.xlsx', rowNames = TRUE)
head(bulk_14)

##去往绘图代码做gsea分析

## 脂肪酸代谢基因查看-------
# a012 上调 脂肪酸代谢
library(tidyverse)
genes <- c('Acsl1/Acaca/Hacd2/Elovl6/Ehhadh/Acox1')
genes <- str_split(genes,'/') %>% unlist()
str_c(genes, collapse = ', ')
## 4
bulk_4_week_2v2[bulk_4_week_2v2$gene %in% genes,]
deg.a012[deg.a012$gene %in% genes ,]

# bulk 4 week fam 下调
genes_4_week <- c('Eci2/Aldh7a1/Echs1/Acsl1/Acat1/Eci1/Acads/Hadha/Acadm/Acadvl/Acadl/Acaa2/Hadhb/Hadh/Cpt2')
genes_4_week <- str_split(genes_4_week,'/') %>% unlist()
genes_4_week_df <- bulk_4_week_2v2[bulk_4_week_2v2$gene %in% genes_4_week,]


intersect(genes, genes_4_week_down)

genes_4_week_down <- genes_4_week_df$gene[genes_4_week_df$sig == 'Down']
str_c(genes_4_week_down, collapse = ', ')
## 
deg.fam_bulk4_down <- FindMarkers(object = adipo.scrna, 
                                  features = genes_4_week_down,
                                  ident.1 = 'KO', ident.2 = 'Flox',
                                  logfc.threshold = 0,
                        subset.ident = c('A0', 'A1', 'A2'), group.by = 'orig.ident')

deg.fam_bulk4_down
openxlsx::write.xlsx(deg.fam_bulk4_down, 
                     file = 'fam_4_week_down_in_sc.xlsx',
                     rowNames  = TRUE)
## 倍数变化 不设阈值
deg.a012.all <- FindMarkers(object = adipo.scrna, 
                            ident.1 = 'KO', ident.2 = 'Flox',
                            logfc.threshold = 0,
                            subset.ident = c('A0', 'A1', 'A2'), group.by = 'orig.ident')


## 同样 ko a4/a012
## koA4/012
ko.a4_v_012_fam_4down <- FindMarkers(adipo.KO, ident.1 = 'A4', 
                           ident.2 = c('A0', 'A1', 'A2'), 
                           features = genes_4_week_down,
                           
                           group.by = 'adipo.type', logfc.threshold = 0)

## 取交集
y = bulk_4_week_2v2$gene[bulk_4_week_2v2$log2FoldChange < 0 & 
                           bulk_4_week_2v2$pvalue < 0.05]
x = deg.a012$gene[deg.a012$avg_log2FC < 0]
intersect(x,y)





