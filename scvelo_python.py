# import
import matplotlib
matplotlib.use('TkAgg')
import anndata
import scvelo as scv
import pandas as pd
import numpy as np

## setting
# scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# check the figure
# flox final figure
flox_filtered_data = scv.read('flox_veloed_data.h5ad')
scv.pl.velocity_embedding_stream(flox_filtered_data, basis="tsne", color="seurat_clusters")
# KO final figure
ko_final_data = scv.read('ko_veloed_data.h5ad')
scv.pl.velocity_embedding_stream(ko_final_data, basis="tsne", color="seurat_clusters")

## flox analysis
ldata = scv.read('Flox_KEZW5.loom', cache=True)
ldata.var_names_make_unique()

ldata.obs_names

## 替换x
obs_name_rep = ldata.obs_names.str.replace('x', '')  ## 这有点蠢，万一需要x呢
obs_name_rep = ldata.obs_names.str[11:27]  ## 要头不要尾
obs_name_rep  ## 成功，有没有更好的方法呢，显然flox要加上才对，毕竟要分组，不然barcode容易重复
obs_name_rep = 'Flox_' + obs_name_rep
obs_name_rep
ldata.obs_names = obs_name_rep
ldata.obs_names

## 读取barcode、tsne以及seurat的分群结果
cellID_obs = pd.read_csv("./excel/cellID_obs.csv")
tsne_cord = pd.read_csv("./excel/cell_embeddings.csv")
cell_clusters = pd.read_csv("./excel/clusters.csv")
## 把seurat导出文件里barcode中的'-1'去掉
cellID_obs['x'] = cellID_obs['x'].str.replace('-1', '')
## 根据barcode选出我们要的细胞
filtered_ldata = ldata[cellID_obs['x']].copy()
## 添加tsne信息
ldata_index = pd.DataFrame(filtered_ldata.obs.index)
tsne_cord = tsne_cord.rename(columns={'Unnamed: 0': 'CellID'})
tsne_cord['CellID'] = tsne_cord['CellID'].str.replace('-1', '')
tsne_ordered = ldata_index.merge(tsne_cord, on="CellID")
tsne_ordered = tsne_ordered.iloc[:, 1:]
filtered_ldata.obsm['X_tsne'] = tsne_ordered.values
## 添加分群信息
cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'CellID'})
cell_clusters['CellID'] = cell_clusters['CellID'].str.replace('-1', '')
cell_clusters = ldata_index.merge(cell_clusters, on="CellID")
cell_clusters = cell_clusters.iloc[:, 1:]
filtered_ldata.obs['seurat_clusters'] = cell_clusters.values
## 保存过滤后还未速率分析的flox
# filtered_ldata.write('flox_filtered_unvelo_data.h5ad')

## 进行速率分析~
scv.pp.filter_and_normalize(filtered_ldata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(filtered_ldata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(filtered_ldata)
scv.tl.velocity_graph(filtered_ldata)
scv.pl.velocity_embedding_stream(filtered_ldata, basis="tsne", color="seurat_clusters")

## 保存flox filtered data
flox_filtered_data = filtered_ldata
# flox_filtered_data.write('flox_veloed_data.h5ad')

flox_filtered_data = scv.read('flox_veloed_data.h5ad')
scv.pl.velocity_embedding_stream(flox_filtered_data, basis="tsne", color="seurat_clusters")

#########################
## 分析KO

ldata = scv.read('KO_LX6HU.loom', cache=True)
ldata.var_names_make_unique()
ldata.obs_names

## 替换x
# obs_name_rep = ldata.obs_names.str.replace('x','') ## 这有点蠢，万一需要x呢
obs_name_rep = ldata.obs_names.str[9:25]
obs_name_rep  ## 成功，有没有更好的方法呢，显然flox要加上才对，毕竟要分组，不然barcode容易重复
obs_name_rep = 'KO_' + obs_name_rep
obs_name_rep
ldata.obs_names = obs_name_rep

## 读取barcode、tsne以及seurat的分群结果
cellID_obs = pd.read_csv("./excel/ko_cellID_obs.csv")
tsne_cord = pd.read_csv("./excel/ko_cell_embeddings.csv")
cell_clusters = pd.read_csv("./excel/ko_clusters.csv")
## 把seurat导出文件里barcode中的'-1'去掉
cellID_obs['x'] = cellID_obs['x'].str.replace('-1', '')
## 根据barcode选出我们要的细胞
filtered_ldata = ldata[cellID_obs['x']].copy()
## 添加tsne信息
ldata_index = pd.DataFrame(filtered_ldata.obs.index)
tsne_cord = tsne_cord.rename(columns={'Unnamed: 0': 'CellID'})
tsne_cord['CellID'] = tsne_cord['CellID'].str.replace('-1', '')
tsne_ordered = ldata_index.merge(tsne_cord, on="CellID")
tsne_ordered = tsne_ordered.iloc[:, 1:]
filtered_ldata.obsm['X_tsne'] = tsne_ordered.values
## 添加分群信息
cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'CellID'})
cell_clusters['CellID'] = cell_clusters['CellID'].str.replace('-1', '')
cell_clusters = ldata_index.merge(cell_clusters, on="CellID")
cell_clusters = cell_clusters.iloc[:, 1:]
filtered_ldata.obs['seurat_clusters'] = cell_clusters.values
## 保存未velo的 ko
ko_filtered_unvelo = filtered_ldata
ko_filtered_unvelo.write('ko_filtered_unvelo.h5ad')

## 进行速率分析~
scv.pp.filter_and_normalize(filtered_ldata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(filtered_ldata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(filtered_ldata)
scv.tl.velocity_graph(filtered_ldata)
scv.pl.velocity_embedding_stream(filtered_ldata, basis="tsne", color="seurat_clusters")

## 保存ko
ko_final_data = filtered_ldata
ko_final_data = scv.read('ko_veloed_data.h5ad')
scv.pl.velocity_embedding_stream(ko_final_data, basis="tsne", color="seurat_clusters")
ko_final_data.write('ko_veloed_data.h5ad')