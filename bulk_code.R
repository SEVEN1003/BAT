## pack load ------
library(readxl)
library(tidyverse)
if(!require(ggplotify))install.packages("ggplotify")
if(!require(patchwork))install.packages("patchwork")
if(!require(cowplot))install.packages("cowplot")
if(!require(DESeq2))BiocManager::install('DESeq2')
if(!require(edgeR))BiocManager::install('edgeR')
if(!require(limma))BiocManager::install('limma')
rm(list = ls())

## data prep----
proj <- 'bat_14wks'
exp <- read_xlsx(path = 'bulk_14wks/All_reads_counts.xlsx', sheet = 1)
colnames(exp)
## data clean
exp <- exp %>% column_to_rownames(var = "Geneid") %>% 
  select(Flox_1:KO_3)
head(exp)

## qc PCA cor heatmap--------
# if(!require(tinyarray))devtools::install_local("tinyarray-master.zip",upgrade = F) 
# BiocManager::install('clusterProfiler')
# BiocManager::install('org.Hs.eg.db')
library(ggplot2)
library(tinyarray)
# cpm counts per million 
dat = log2(cpm(exp)+1) #dat for plot
cor(dat)
Group <- factor(c(rep('Flox',3), rep('KO',3)),levels = c('Flox','KO'))
# install.packages('FactoMineR')
# install.packages('factoextra')
pca.plot = draw_pca(dat,Group);pca.plot
## save
ggsave(filename = paste0(proj, '_pca.pdf'), width = 4.0, height = 3.2)
# cor heatmap
library(pheatmap)
cor_heat <- pheatmap(cor(dat))
save_pheatmap_pdf <- function(x, filename, width=4.5, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(cor_heat, filename = 'cor_heat.pdf',4.5,4)
# 
colData <- data.frame(row.names =colnames(exp), 
                      condition = Group) 

colData
## deseq2 -----
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = exp,  
  colData = colData, 
  design = ~ condition) 

dds <- DESeq(dds)  

save(dds,file = paste0(proj,"_dd.Rdata"))

# get results
res <- results(dds, contrast = c("condition",rev(levels(Group)))) 
resOrdered <- res[order(res$pvalue),] 
DEG <- as.data.frame(resOrdered)
DEG = na.omit(DEG) 
DEG$gene <- rownames(DEG)
library(openxlsx)
write.xlsx(DEG, file = paste0(proj,'_deseq.xlsx'), rowNames = TRUE)

## 
exp_14wks <- exp
dat_14wks <- dat
deg_14wks <- DEG

deg <- deg_14wks
deg$sig <- ifelse(deg$pvalue < 0.05,
                  ifelse(deg$log2FoldChange > 1,'Up',
                         ifelse(deg$log2FoldChange < -1, 'Down', 'No sig.')),'No sig.')
table(deg$sig)
deg$sig <- factor(deg$sig, levels = c('No sig.', 'Up','Down'))

## GSEA------
## pack load
library(GSEABase) 
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(tidyverse)
library(enrichplot)

head(deg)
genelist <- deg$log2FoldChange
names(genelist) <- deg$gene
genelist <- sort(genelist, decreasing = TRUE)
head(genelist)

## kegg gmt 
## go gmt
kegg.gmt <- read.gmt('public-data/mouse-kegg.gmt')
kegg.gmt$term <- str_remove(kegg.gmt$term, 'KEGG_')


## gsea
egmt <- GSEA(geneList = genelist, TERM2GENE = kegg.gmt, verbose = FALSE,
             pvalueCutoff = 0.5, eps = 0)
egmt.df <- data.frame(egmt)
write.xlsx(egmt.df, file = 'enrich_results/gsea_results.xlsx')
save.image(file = paste0(proj,'.Rdata'))

# plot ------
# pie
library(ggstatsplot)
pie <- ggpiestats(
  deg,
  'sig',
  results.subtitle = F,
  #标题中不显示统计结果
  label = 'both',
  #设置饼图中标签的类型（默认percentage：百分比, 还有counts：频数, both : 频数和百分比都显示）
  perc.k = 2,
  #百分比的小数位数为2
  direction = 1,
  #1为顺时针方向，-1为逆时针方向
  #设置调色板
  title = '',
  legend.title = 'Type',
  label.repel = TRUE
) +
  scale_fill_manual(values = c('#00ff00','#ff2525','#868686' ))# 绿红灰: 下上无
pie

ggsave(filename = paste0(proj,'_pie.pdf'),
       width = 5,
       height = 3)

# GSEA 气泡图------
setid_down <- c('OXIDATIVE_PHOSPHORYLATION','NITROGEN_METABOLISM',
                'PROTEASOME','CITRATE_CYCLE_TCA_CYCLE')
setid_up <- c('ECM_RECEPTOR_INTERACTION','LYSOSOME','ANTIGEN_PROCESSING_AND_PRESENTATION',
              'NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY','B_CELL_RECEPTOR_SIGNALING_PATHWAY')
df <- egmt.df
df <- df[c(setid_down, setid_up),]
## 按nes排序
df <- df %>% arrange(desc(NES))
rownames(df) <- 1:nrow(df)
# colnames(df)
df$ID <- gsub('_',' ',str_to_title(df$ID))
df$change <- ifelse(df$NES > 0, 'Up', 'Down')
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
    y="NES")+
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

ggsave(filename = 'gsea.pdf', width = 5.5, height = 3.2)


# venn图------
library(ggvenn)
library(ggtext)
deg2 <- deg_4wks
deg2$sig <- ifelse(deg2$pvalue < 0.05,
                   ifelse(deg2$log2FoldChange > 1,'Up',
                          ifelse(deg2$log2FoldChange < -1, 'Down', 'No sig.')),'No sig.')
# 上调交集
up_1 <- deg$gene[deg$sig == 'Up']
up_2 <- deg2$gene[deg2$sig == 'Up']
venn <- list(up_1  = up_1,up_2 = up_2) %>%
  ggvenn(
    show_percentage = T,
    show_elements = F,
    label_sep = ",",
    digits = 1,
    stroke_color = "white",
    fill_color = c("#1E90FF", "#E41A1C"),
    set_name_color = c("#1E90FF", "#E41A1C")
  )
venn
ggsave(filename = 'venn.pdf',
       width = 5.4,
       height = 3.3)

























