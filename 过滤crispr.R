library(Seurat)
library(infercnv)
library(glue)
library(ggsci)
library(ComplexHeatmap)
col_map_func <- function(dat){
  res = list()
  for(n in colnames(dat)){
    tmp_name = sort(unique(dat[, n]))
    tmp_col = pal_igv()(length(tmp_name))
    names(tmp_col) = tmp_name
    res[[n]] = tmp_col
  }
  return(res)
}
wx_integer_CNV <- function(cnv, method='dist'){
  if(method=='dist'){
    #min cn is 0 ploid, cn=1 is 2 ploid
    min_cn = min(cnv)
    cnv = round((cnv-min_cn) / (1-min_cn) * 2, 0)
  }
  else{
    cnv[cnv > 1.1] = 2
    cnv[cnv < 0.9] = 0
    cnv[cnv!=2&cnv!=0] = 1
  }
  return(cnv)
}
process_seurat <-function(obj, dim=30, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = row.names(obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,n.neighbors = 10L,dims = 1:dim,reduction = "pca", min.dist = 0.3, n.components=n.components)
  return(obj)
}
setwd('/Volumes/WX_extend/细胞轨迹推断/RNAseq/')
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_tree/'


############ Huh7 ############
prefix = 'Huh7'
srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
srt = process_seurat(srt, resolution = 1)
a=DimPlot(srt, pt.size = 4) + scale_color_igv()
b=DimPlot(srt, pt.size = 4, group.by = 'feature_call') + scale_color_igv()
a+b
# 每个转录组的亚群，选择top的gRNA，保证不出现在其他亚群
res = srt@meta.data %>%
  group_by(seurat_clusters, feature_call) %>%
  summarise(num = n()) %>%
  arrange(-num) %>%
  as.data.frame()
res$seurat_clusters = as.vector(res$seurat_clusters)
unique_sc = as.vector(unique(res$seurat_clusters))
keep_res = c()
while (length(unique_sc)!=0) {
  sc = res[1, 'seurat_clusters']
  fc = res[1, 'feature_call']
  res = filter(res, seurat_clusters!=sc&feature_call!=fc)
  unique_sc = setdiff(unique_sc, sc)
  keep_res = rbind(keep_res, c(sc, fc))
}
keep_srt = list()
for(i in 1:nrow(keep_res)){
  tmp_srt = subset(srt, seurat_clusters==keep_res[i,1]&feature_call==keep_res[i,2])
  keep_srt[[i]] = tmp_srt
}
keep_srt = merge(keep_srt[[1]], keep_srt[2:length(keep_srt)],merge.dr='umap' )
#keep_srt = process_seurat(keep_srt, resolution = 0.1)

a=DimPlot(keep_srt, pt.size = 4) + scale_color_igv()
b=DimPlot(keep_srt, pt.size = 4, group.by = 'feature_call') + scale_color_igv()
a+b

# run infercnv_02_重新过滤.R
cnv_res = read.table(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_重新过滤/{prefix}/expr.infercnv.21_denoised.dat'),
                     check.names = F)
colnames(cnv_res) = gsub('_1', '',colnames(cnv_res))
# cnv_res = cnv_res-1 # state to int cnv, state6 as cnv 5
int_cnv_res = wx_integer_CNV(t(cnv_res))
write.table(int_cnv_res[colnames(keep_srt),], paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_infercnv_int_gene_filter.txt'), quote = F, sep='\t')
saveRDS(keep_srt, paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_crispr_srt_filter.txt'))



############ A549 ############
prefix = 'A549'
srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
srt = process_seurat(srt, resolution = 1)
a=DimPlot(srt, pt.size = 2) + scale_color_igv()
b=DimPlot(srt, pt.size = 2, group.by = 'feature_call') + scale_color_igv() + NoLegend()
a+b
# 每个转录组的亚群，选择top的gRNA，保证不出现在其他亚群
res = srt@meta.data %>%
  group_by(seurat_clusters, feature_call) %>%
  summarise(num = n()) %>%
  arrange(-num) %>%
  as.data.frame()
res$seurat_clusters = as.vector(res$seurat_clusters)
unique_sc = as.vector(unique(res$seurat_clusters))
keep_res = c()
while (length(unique_sc)!=0) {
  sc = res[1, 'seurat_clusters']
  fc = res[1, 'feature_call']
  res = filter(res, seurat_clusters!=sc&feature_call!=fc)
  unique_sc = setdiff(unique_sc, sc)
  keep_res = rbind(keep_res, c(sc, fc))
}
keep_srt = list()
for(i in 1:nrow(keep_res)){
  tmp_srt = subset(srt, seurat_clusters==keep_res[i,1]&feature_call==keep_res[i,2])
  keep_srt[[i]] = tmp_srt
}
keep_srt = merge(keep_srt[[1]], keep_srt[2:length(keep_srt)],merge.dr='umap' )
#keep_srt = process_seurat(keep_srt, resolution = 0.1)

a=DimPlot(keep_srt, pt.size = 4) + scale_color_igv()
b=DimPlot(keep_srt, pt.size = 4, group.by = 'feature_call') + scale_color_igv()
a+b

# run infercnv_02_重新过滤.R
cnv_res = read.table(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_重新过滤/{prefix}/expr.infercnv.21_denoised.dat'),
                     check.names = F)
colnames(cnv_res) = gsub('_1', '',colnames(cnv_res))
# cnv_res = cnv_res-1 # state to int cnv, state6 as cnv 5
int_cnv_res = wx_integer_CNV(t(cnv_res))
write.table(int_cnv_res[colnames(keep_srt),], paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_infercnv_int_gene_filter.txt'), quote = F, sep='\t')
saveRDS(keep_srt, paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_crispr_srt_filter.txt'))


############ U251 ############
prefix = 'U251'
srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
srt = process_seurat(srt, resolution = 1)
a=DimPlot(srt, pt.size = 2) + scale_color_igv()
b=DimPlot(srt, pt.size = 2, group.by = 'feature_call') + scale_color_igv() + NoLegend()
a+b
# 每个转录组的亚群，选择top的gRNA，保证不出现在其他亚群
res = srt@meta.data %>%
  group_by(seurat_clusters, feature_call) %>%
  summarise(num = n()) %>%
  arrange(-num) %>%
  as.data.frame()
res$seurat_clusters = as.vector(res$seurat_clusters)
unique_sc = as.vector(unique(res$seurat_clusters))
keep_res = c()
while (length(unique_sc)!=0) {
  sc = res[1, 'seurat_clusters']
  fc = res[1, 'feature_call']
  res = filter(res, seurat_clusters!=sc&feature_call!=fc)
  unique_sc = setdiff(unique_sc, sc)
  keep_res = rbind(keep_res, c(sc, fc))
}
keep_srt = list()
for(i in 1:nrow(keep_res)){
  tmp_srt = subset(srt, seurat_clusters==keep_res[i,1]&feature_call==keep_res[i,2])
  keep_srt[[i]] = tmp_srt
}
keep_srt = merge(keep_srt[[1]], keep_srt[2:length(keep_srt)],merge.dr='umap' )
#keep_srt = process_seurat(keep_srt, resolution = 0.1)

a=DimPlot(keep_srt, pt.size = 4) + scale_color_igv()
b=DimPlot(keep_srt, pt.size = 4, group.by = 'feature_call') + scale_color_igv()
a+b

# run infercnv_02_重新过滤.R
cnv_res = read.table(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_重新过滤/{prefix}/expr.infercnv.21_denoised.dat'),
                     check.names = F)
colnames(cnv_res) = gsub('_1', '',colnames(cnv_res))
# cnv_res = cnv_res-1 # state to int cnv, state6 as cnv 5
int_cnv_res = wx_integer_CNV(t(cnv_res))
write.table(int_cnv_res[colnames(keep_srt),], paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_infercnv_int_gene_filter.txt'), quote = F, sep='\t')
saveRDS(keep_srt, paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_crispr_srt_filter.txt'))



############ U251 ############
prefix = '786O'
srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
srt = process_seurat(srt, resolution = 1)
a=DimPlot(srt, pt.size = 2) + scale_color_igv()
b=DimPlot(srt, pt.size = 2, group.by = 'feature_call') + scale_color_igv() + NoLegend()
a+b
# 每个转录组的亚群，选择top的gRNA，保证不出现在其他亚群
res = srt@meta.data %>%
  group_by(seurat_clusters, feature_call) %>%
  summarise(num = n()) %>%
  arrange(-num) %>%
  as.data.frame()
res$seurat_clusters = as.vector(res$seurat_clusters)
unique_sc = as.vector(unique(res$seurat_clusters))
keep_res = c()
while (length(unique_sc)!=0) {
  sc = res[1, 'seurat_clusters']
  fc = res[1, 'feature_call']
  res = filter(res, seurat_clusters!=sc&feature_call!=fc)
  unique_sc = setdiff(unique_sc, sc)
  keep_res = rbind(keep_res, c(sc, fc))
}
keep_srt = list()
for(i in 1:nrow(keep_res)){
  tmp_srt = subset(srt, seurat_clusters==keep_res[i,1]&feature_call==keep_res[i,2])
  keep_srt[[i]] = tmp_srt
}
keep_srt = merge(keep_srt[[1]], keep_srt[2:length(keep_srt)],merge.dr='umap' )
#keep_srt = process_seurat(keep_srt, resolution = 0.1)

a=DimPlot(keep_srt, pt.size = 4) + scale_color_igv()
b=DimPlot(keep_srt, pt.size = 4, group.by = 'feature_call') + scale_color_igv()
a+b

# run infercnv_02_重新过滤.R
cnv_res = read.table(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_重新过滤/{prefix}/expr.infercnv.21_denoised.dat'),
                     check.names = F)
colnames(cnv_res) = gsub('_1', '',colnames(cnv_res))
# cnv_res = cnv_res-1 # state to int cnv, state6 as cnv 5
int_cnv_res = wx_integer_CNV(t(cnv_res))
write.table(int_cnv_res[colnames(keep_srt),], paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_infercnv_int_gene_filter.txt'), quote = F, sep='\t')
saveRDS(keep_srt, paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/', prefix,'_crispr_srt_filter.txt'))



