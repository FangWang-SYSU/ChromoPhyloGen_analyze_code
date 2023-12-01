library(Seurat)
seurat_process <-function(obj, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:20)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:20,reduction = "pca", n.components=n.components)
  return(obj)
}
srt = readRDS('/Volumes/WX_extend/CRC项目/tumors_filter_by_cnv_man.rds')


scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC/'

file_name = c('Pt01', 'Pt02', 'Pt03', 'Pt04', 'Pt05', 'Pt06', 'Pt07', 'Pt08', 'Pt09', 'Pt10',
              'Pt11', 'Pt13', 'Pt15', 'Pt16', 'Pt17', 'Pt18')

prefix = 'Pt16'
sub_srt = subset(srt, patient==prefix)
sub_srt = seurat_process(sub_srt)
DimPlot(sub_srt, group.by = 'lesions') + ggsci::scale_color_igv()


count_data = t(as.matrix(sub_srt@assays$RNA@data)[VariableFeatures(sub_srt), ])


# 过滤低表达基因
count_data = count_data[, colSums(count_data!=0)>=0.1, drop=F]

count_data = matrix(c(2,2,2,2,3,3,3,
                      2,1,1,2,5,5,5,
                      3,3,3,3,3,2,2,
                      3,3,3,3,3,2,2,
                      2,1,1,2,5,5,5), nrow=5, byrow = T)
colnames(count_data) = c('1', '2', '3', '4','5','6','7')
geneid = 'S100A6'
library(mcga)

fit_fun <- function(x){
  x = round(x)
  to_x = matrix(rep(x, length(x)), nrow=length(x), byrow=T)
  px=to_x / (x+1e-6)
  # px[is.na(px)] = 6
  #px[is.infinite(px)] = 6
  loss = sqrt(sum(real_ratio-px)^2)
  print(loss)
  return(loss)
}
all_cna = c()
for(geneid in colnames(count_data)[1:6]){
  print(geneid)
  to_data = matrix(rep(count_data[, geneid], nrow(count_data)), nrow=nrow(count_data), byrow=T)
  real_ratio=to_data / (count_data[, geneid]+1e-6)
  m <- mcga(  popsize=200,
              chsize=as.numeric(nrow(count_data)),
              minval=0.0,
              maxval=6,
              maxiter=2500,
              crossprob=1,
              mutateprob=0.01,
              evalFunc=fit_fun)
  predict_cna = round(m$population[1,])
  all_cna = cbind(all_cna, predict_cna)
}
all_cna[all_cna>6] = 6
pheatmap(all_cna, cluster_rows = F, cluster_cols = F)
pheatmap(count_data[,1:6], cluster_rows = F, cluster_cols = F)

x=all_cna[,1]
to_x = matrix(rep(x, length(x)), nrow=length(x), byrow=T)
px=to_x / (x+1e-6)


