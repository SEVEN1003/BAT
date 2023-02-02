## get all gene of a kegg pathway------
## kegg rest
library(KEGGREST)
listDatabases()
## 获得某一条
## human hsa mouse mmu 
## KEGG_Thermogenesis bing搜索后点击Pathway menu 查看编号 | Organism menu查看物种
gs <- keggGet('mmu04714') ## 产热
gs<-keggGet('mmu04621') ## nod 
gs <- keggGet('mmu03320') ## pparg
gs <- keggGet('mmu04152') ## AMPK

gs <- keggGet('mmu03050') ## protea
gs<-keggGet('mmu04621') ## nod 
gs <- keggGet('mmu04612') ##app
gs <- keggGet('mmu04217') ## necro
gs <- keggGet('mmu00640') ## propan
gs <- keggGet('mmu01212')
#获取通路中gene信息 
# head(gs[[1]]$GENE)
# hsa 查找所有基因 
# genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
# genelist <- genes[1:length(genes)%%3 ==2] 

## mmu
genes <- gs[[1]]$GENE 
genes <- genes[1:length(genes) %% 2 == 0]
genes <- str_split(genes, pattern = ';',simplify = TRUE)[,1]
genes
genes <- sort(genes)

## 
Propanoate_meta.all <- genes
Proteasome.all <- genes
APP.all <- genes
NOD_all
Necroptosis.all <- genes

head(adipo.scrna)

genelist <- list(Proteasome.all = Proteasome.all,
                 APP.all = APP.all,
                 Necroptosis.all = Necroptosis.all)

## propan丙酸 meta
Propanoate_meta <- c('Acaca/Ehhadh/Acox1/Acss1')
Propanoate_meta <- str_split(Propanoate_meta, pattern = '/') %>% unlist()
genelist <- list(Propanoate_meta = Propanoate_meta,
                 Propanoate_meta.all = Propanoate_meta.all)


PPAR_all <- genes
AMPK_all <- genes
## feature plot 映射图-----
m <- 'Iigp1'
m <- 'Xdh'
m <- 'Gbp7'
m <- 'Gbp2'
m <- 'Irf1'
m <- 'Ucp1'

FeaturePlot(seu.gse.qc.adip, features = m, reduction = 'tsne', 
            pt.size = 1.5, label = TRUE)+NoAxes()+
  theme(plot.title  = element_text(face = 'italic')
  )

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



                       
FeaturePlot(seu.gse.qc.adip, features = marker, reduction = 'tsne', pt.size = 1.5)+
  NoAxes()+theme(plot.title = element_text(face = 'italic', size = 20))


m <- 'Cyp2e1'
FeaturePlot(bat_total_12, features = m, reduction = 'tsne')

## more score-------
## A3 necro and antigen 评分
## A4 ppar 通路评分
## A012 AMPK 评分

## Antigen processing and presentation
APP <- c('Tap1/Tap2/H2-T22/Tapbp/H2-T23/Psme1/Ctsb/B2m/H2-K1')
Necroptosis <- c('Stat1/Stat2/Mlkl/Stat3/Glul/Eif2ak2/Bid/Irf9/Ifngr1/Camk2d/Fas')

## PPAR 做个全部的
PPAR <- c('Angptl4/Gk/Fabp3/Fabp4/Ehhadh/Slc27a2/Lpl')
PPAR_all
## AMPK 做个全部的
AMPK <- c('Pparg/Igf1r/Ppp2r3a/Eef2k/Camkk2/Prkag2')
AMPK_all

genelist <- list()
geneset <- c(APP, Necroptosis, PPAR, AMPK)
for (i in seq_along(geneset)) {
  genelist[[i]] <- str_split(geneset[[i]], '/') %>% unlist()
}
names(genelist) <- c('APP', 'Necroptosis', 'PPAR', 'AMPK')
## 评分
genelist$PPAR_all <- PPAR_all
genelist$AMPK_all <- AMPK_all

length(PPAR_all)
intersect(rownames(adipo.scrna), PPAR_all) %>% length()

intersect(rownames(adipo.scrna), AMPK_all) %>% length()

## all protea nod app necro

## 箱式图秩和检验---------
head(metadata)
Proteasome_wilcox <- wilcox.test(Proteasome ~ orig.ident, data = metadata,
                      exact = FALSE,correct = FALSE, subset = adipo.type == 'A0')

## 写个循环做做全部
colnames(metadata)
table(metadata$adipo.type)
metadata$adipo.type <- droplevels(metadata$adipo.type)
cell_type_list <- levels(metadata$adipo.type)
# 提取统计量和p值，保存在数据库框
wilcox_test <- data.frame(adipo_type = levels(metadata$adipo.type) )
## 
colnames(metadata)
## Proteasome Thermogenesis_all; Fatty_acid_metabolism; NOD_A3 ; APP; Necroptosis
## PPAR; AMPK

## protea
for (x in 1:length(cell_type_list)) {
  wilcox_results <- wilcox.test(PPAR ~ orig.ident, data = metadata,
                        exact = FALSE, correct = FALSE, subset = adipo.type == cell_type_list[x])
  wilcox_test$W[x] <- wilcox_results$statistic
  wilcox_test$p_value[x] <- wilcox_results$p.value
}

Proteasome_wilcox <- wilcox_test
AMPK_wilcox <- wilcox_test
fam_wilcox <- wilcox_test

NOD_A3_wilcox <- wilcox_test
app_wilcox <- wilcox_test
necro_wilcox <- wilcox_test

therm_wilcox <- wilcox_test
ppar_wilcox <- wilcox_test

## gsea 4week -------
## 准备gmt文件
## pack load-----
library(GSEABase) 
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(tidyverse)
library(enrichplot)

## 准备差异基因列表----
## 使用2v2的分析
genelist <- bulk_4_week_2v2$log2FoldChange
#head(genelist)

names(genelist) <- bulk_4_week_2v2$gene

#head(names(genelist))
genelist <- sort(genelist, decreasing = TRUE)
head(genelist)

# 已制作
kegg.gmt <- read.gmt('bulk_4week/mouse-kegg.gmt')
kegg.gmt$term <- str_remove(kegg.gmt$term, 'KEGG_')

## 添加AMPK NCRO THERM

## ampk
AMPK_all
Necroptosis.all
Thermogenesis_all
APP.all
df <- data.frame(term = c(rep('AMPK signaling pathway',length(AMPK_all)),
                          rep('Necroptosis', length(Necroptosis.all)),
                          rep('Thermogenesis', length(Thermogenesis_all))),
                 gene = c(AMPK_all, Necroptosis.all, Thermogenesis_all))
## rbind
kegg.gmt <- rbind(kegg.gmt, df)

## 分析
egmt.4wks <- GSEA(geneList = genelist, TERM2GENE = kegg.gmt, verbose = FALSE,
             pvalueCutoff = 0.5, eps = 0)

egmt.df.4wks <- data.frame(egmt.4wks)

# write.csv(egmt.df, file = 'gsea.csv')
write.xlsx(egmt.df.4wks, file = 'bulk_4week/gsea_2v2.xlsx')
# 读入 4week 2v2
# egmt.df.4wks <- read.xlsx('bulk_4week/gsea_2v2.xlsx', rowNames = TRUE)

# par(pin=c(3,3))
## 下调
setid_down <- c('OXIDATIVE_PHOSPHORYLATION', 'CITRATE_CYCLE_TCA_CYCLE',
                'FATTY_ACID_METABOLISM','Thermogenesis',
                'PROTEASOME', 'NITROGEN_METABOLISM')

## 上调 
setid_up <- c('GRAFT_VERSUS_HOST_DISEASE','ASCORBATE_AND_ALDARATE_METABOLISM',
              'CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
              'NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',
              'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY')
## 绘图 
# 下调
for (setid in setid_down) {
  gseaplot2(x = egmt.4wks, geneSetID = setid, title = gsub('_',' ',str_to_title(setid)), base_size = 15, 
            subplots = 1:2, pvalue_table = FALSE, color = 'green')
  ggsave(filename = paste0('down_',setid, '_gsea.pdf'), width = 4.0, height = 2.8, path = 'bulk_4week/')
  
}

# 上
for (setid in setid_up) {
  gseaplot2(x = egmt.4wks, geneSetID = setid, title = gsub('_',' ',str_to_title(setid)), base_size = 15, 
            subplots = 1:2, pvalue_table = FALSE, color = 'red')
  ggsave(filename = paste0('up_',setid, '_gsea.pdf'), width = 4.0, height = 2.8, path = 'bulk_4week/')
  
}

## gsea 14 week -------
## 准备gmt文件
## pack load-----
library(GSEABase) 
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(tidyverse)
library(enrichplot)

## 准备差异基因列表----
head(bulk_14)

genelist <- bulk_14$log2FoldChange
#head(genelist)

names(genelist) <- bulk_14$gene

#head(names(genelist))
genelist <- sort(genelist, decreasing = TRUE)
head(genelist)

## kegggmt已经准备好了
head(kegg.gmt)
kegg.gmt <- read.gmt('bulk_4week/mouse-kegg.gmt')
kegg.gmt$term <- str_remove(kegg.gmt$term, 'KEGG_')

## 添加AMPK NCRO THERM
## ampk
AMPK_all
Necroptosis.all
Thermogenesis_all
APP.all
df <- data.frame(term = c(rep('AMPK signaling pathway',length(AMPK_all)),
                          rep('Necroptosis', length(Necroptosis.all)),
                          rep('Thermogenesis', length(Thermogenesis_all))),
                 gene = c(AMPK_all, Necroptosis.all, Thermogenesis_all))
## rbind
kegg.gmt <- rbind(kegg.gmt, df)

## gsea分析~
egmt <- GSEA(geneList = genelist, TERM2GENE = kegg.gmt, verbose = FALSE,
             pvalueCutoff = 0.5, eps = 0)

egmt.df <- data.frame(egmt)
write.xlsx(egmt.df, file = 'bulk_14week/gsea_results.xlsx')

# par(pin=c(3,3))
## 下调前五条
# 1
setid <- 'OXIDATIVE_PHOSPHORYLATION'
setid <- 'Thermogenesis'
setid <- 'NITROGEN_METABOLISM'
setid <- 'PROTEASOME'
setid <- 'CITRATE_CYCLE_TCA_CYCLE' # 显著性不高
setid_down <- c('OXIDATIVE_PHOSPHORYLATION','Thermogenesis','NITROGEN_METABOLISM',
                'PROTEASOME','CITRATE_CYCLE_TCA_CYCLE')

# 上调前五条
setid <- 'ECM_RECEPTOR_INTERACTION'
setid <- 'LYSOSOME'
setid <- 'ANTIGEN_PROCESSING_AND_PRESENTATION'
setid <- 'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'
setid <- 'B_CELL_RECEPTOR_SIGNALING_PATHWAY'

setid_up <- c('ECM_RECEPTOR_INTERACTION','LYSOSOME','ANTIGEN_PROCESSING_AND_PRESENTATION',
                'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY','B_CELL_RECEPTOR_SIGNALING_PATHWAY')

## 绘图 
# 下调
for (setid in setid_down) {
  gseaplot2(x = egmt, geneSetID = setid, title = gsub('_',' ',str_to_title(setid)), base_size = 15, 
            subplots = 1:2, pvalue_table = FALSE, color = 'green')
  ggsave(filename = paste0('down_',setid, '_gsea.pdf'), width = 4.0, height = 2.8, path = 'bulk_14week/')
  
}

# 上
for (setid in setid_up) {
  gseaplot2(x = egmt, geneSetID = setid, title = gsub('_',' ',str_to_title(setid)), base_size = 15, 
            subplots = 1:2, pvalue_table = FALSE, color = 'red')
  ggsave(filename = paste0('up_',setid, '_gsea.pdf'), width = 4.0, height = 2.8, path = 'bulk_14week/')
  
}

## 4周差异基因 饼图--------
bulk_4_week_2v2 %>% head()
bulk_4_week_2v2$sig <- ifelse(bulk_4_week_2v2$pvalue < 0.05,
                              ifelse(bulk_4_week_2v2$log2FoldChange > 1,'Up',
                                     ifelse(bulk_4_week_2v2$log2FoldChange < -1,'Down','No sig.')),'No sig.')
table(bulk_4_week_2v2$sig)

## 14周差异基因 饼图
bulk_14$sig <- ifelse(bulk_14$pvalue < 0.05,
                      ifelse(bulk_14$log2FoldChange > 1,'Up',
                             ifelse(bulk_14$log2FoldChange < -1,'Down','No sig.')),'No sig.')
table(bulk_14$sig)

## 四周饼图
library(scales)
pie.df <- as.data.frame(table(bulk_14$sig))


colnames(pie.df)[1] <- 'Type'

percentage <- scales::percent(pie.df$Freq/sum(pie.df$Freq))
percentage

pie.df$labs <- paste0(pie.df$Type,' (', pie.df$Freq,', ',percentage,')')
pie.df

library(ggpubr)
p <- ggpie(pie.df, 'Freq',  #绘图，只用写频数就行，切记不用再写分组
            fill = 'Type', palette  = c('#00ff00', '#868686','#ff2525'), #按照Cylinders填充，颜色板为jco.
            label = pie.df$labs, lab.pos = 'out', lab.font = c(4, 'white'))+NoLegend() #设置标签，标签的位置在图的内部，标签的大小为4， 颜色为白色.
p

ggsave(filename = 'pie_bulk_14_wks.pdf', width = 5, height = 3)

## 饼图 方法二------
library(ggstatsplot)

bulk_4_week_2v2$sig <- factor(bulk_4_week_2v2$sig, levels = c('Down', 'Up', 'No sig.'))

p2 <- ggpiestats(bulk_4_week_2v2, 'sig', 
                 results.subtitle = F, #标题中不显示统计结果
                 label = 'both', #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
                 perc.k = 2, #百分比的小数位数为2
                 direction = 1, #1为顺时针方向，-1为逆时针方向
                  #设置调色板
                 title = '',
                 legend.title = 'Type',
                 label.repel = TRUE)+
  scale_fill_manual(values = c('#868686', '#ff2525','#00ff00'))
p2

ggsave(filename = 'pie_bulk_4_wks.pdf', width = 5, height = 3)

## gsea 气泡图------
bulk_14 %>% head()
egmt.df %>% head()
## 14 Week dot
setid_down <- c('OXIDATIVE_PHOSPHORYLATION','Thermogenesis','NITROGEN_METABOLISM',
                'PROTEASOME','CITRATE_CYCLE_TCA_CYCLE')
setid_up <- c('ECM_RECEPTOR_INTERACTION','LYSOSOME','ANTIGEN_PROCESSING_AND_PRESENTATION',
              'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY','B_CELL_RECEPTOR_SIGNALING_PATHWAY')

## 4 WKS DOT
## 下调
setid_down <- c('OXIDATIVE_PHOSPHORYLATION', 'CITRATE_CYCLE_TCA_CYCLE',
                'FATTY_ACID_METABOLISM','Thermogenesis',
                'PROTEASOME', 'NITROGEN_METABOLISM')

## 上调 
setid_up <- c('GRAFT_VERSUS_HOST_DISEASE','ASCORBATE_AND_ALDARATE_METABOLISM',
              'CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
              'NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',
              'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY')


df <- egmt.df.4wks[c(setid_down, setid_up),]
## 按nes排序
df <- df %>% arrange(desc(NES))
rownames(df) <- 1:nrow(df)
# colnames(df)
df$ID <- gsub('_',' ',str_to_title(df$ID))
df$change <- ifelse(df$NES > 0, 'Up', 'Down')
table(df$change)
## 构建索引
porder <- factor(as.integer(rownames(df)), labels = rev(df$ID))
# porder

pathbar = ggplot(df,aes(y = NES, x = ID))


pathbar+ 
  geom_point(stat="identity",
             aes(x=rev(porder),
             size =-log10(p.adjust),color = change))+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_color_manual(values = c('#00ff00', '#ff2525'))+
  labs(
    color= 'Type',
    y="NES",
    # x="Pathway name",
    #title="KEGG enrichment"
  )+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = rel(1.3), colour = 'black'),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )+coord_flip()

ggsave(filename = 'bulk_4_dot_gsea.pdf', width = 5.5, height = 3.2)

## venn 图用于通路的交集------
# 4
term <- egmt.df.4wks$ID
term <- gsub('_', ' ',str_to_title(term))
# down
term_4wks_down <- term[egmt.df.4wks$NES < 0]
term_4wks_up <- term[egmt.df.4wks$NES > 0]

intersect(term_4wks_down, term_14wks_down)
setdiff(term_4wks_down, term_14wks_down)
setdiff(term_14wks_down,term_4wks_down)

# 14
term <- egmt.df$ID
term <- gsub('_', ' ',str_to_title(term))

term_14wks_down <- term[egmt.df$NES < 0]
term_14wks_up <- term[egmt.df$NES > 0]

intersect(term_4wks_up, term_14wks_up)
setdiff(term_4wks_up, term_14wks_up)
setdiff(term_14wks_up,term_4wks_up)


library(ggvenn)
library(ggtext)
## down
venn.down <- list(week_4 = term_4wks_down,
     week_14 = term_14wks_down) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "white",
         fill_color = c("#1E90FF","#E41A1C"),
         set_name_color = c("#1E90FF","#E41A1C"))
venn.down
ggsave(filename = 'bulk.venn.down.pdf', width = 5.4, height = 3.3)
## up
venn.up <- list(week_4 = term_4wks_up,
                  week_14 = term_14wks_up) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "white",
         fill_color = c("#1E90FF","#E41A1C"),
         set_name_color = c("#1E90FF","#E41A1C"))

venn.up

ggsave(filename = 'bulk.venn.up.pdf', width = 5.4, height = 3.3)


## 绘制交集热图-----
# 下调
intersect(term_4wks_down, term_14wks_down)

setdiff(term_4wks_down, term_14wks_down)
setdiff(term_14wks_down,term_4wks_down)

## 要得到两个列表
## 4周的term 及对应的基因
common_down <- c('OXIDATIVE_PHOSPHORYLATION','Thermogenesis','NITROGEN_METABOLISM',
                'PROTEASOME','CITRATE_CYCLE_TCA_CYCLE')

term_gene_list <- str_split(egmt.df$core_enrichment, pattern = '/')
names(term_gene_list) <- egmt.df$ID
## 绘制14周这五条term的热图,取这些通路中p < 0.05的
gene <- term_gene_list[common_down] %>% unlist()
gene

gene <- gene[gene %in% bulk_14$gene[ bulk_14$pvalue< 0.05]]

## nitro.gen
gene <- gene[grepl('NITROGEN_METABOLISM', names(gene))]
## 读入rpkm文件
library(openxlsx)

rpkm_14wks <- read.xlsx('bulk_14week/All_samples_rpkm.xlsx', rowNames = TRUE)
colnames(rpkm_14wks)
rpkm_14wks <- rpkm_14wks %>% select( Flox_1:Flox_3,KO_1:KO_3)

gene
pheatmap::pheatmap(rpkm_14wks[gene, ], scale = 'row', 
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_rownames = TRUE,gaps_col = 3,cellwidth = 20)



term_gene_list_4_week <- str_split(egmt.df.4wks$core_enrichment, pattern = '/')
names(term_gene_list_4_week) <- egmt.df.4wks$ID










