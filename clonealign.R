setwd('/Volumes/WX_extend/细胞轨迹推断/RNAseq/')
library(reticulate)
use_python("/Users/lab/.virtualenvs/clonealign_env/bin/python")
use_virtualenv("clonealign_env")
library(tensorflow)
#install_tensorflow()
#library(tensorflow)
# tf$constant("Hello Tensorflow!")
tf$print(
  tf$constant('Hello, TensorFlow!')
)
# devtools::install_github("kieranrcampbell/clonealign")
library(clonealign)
################
library(Seurat)
srt = readRDS('./output/cell_line_srt_with_SC3_cytotrace_new.RDS')

expdata = FetchData(srt, vars = VariableFeatures(srt), slot = 'counts')
# cna
cna_dir = '/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/'
file_names = c(
  'huh7.txt',
  'plc.txt',
  #'hep3b.txt',
  'smc7721.txt',
  'hepg2.txt',
  'mhcc97h.txt',
  'mhcc97l.txt'
  #'lm3.txt'
)
DNA_sctc = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')

cna_data = list()
for(file_name in file_names){
  print(file_name)
  tmp_dir = paste0(cna_dir, '/', file_name, '_gene_data.csv')
  tmp_cna = read.csv(tmp_dir, sep=' ', header = T, check.names = F)
  tmp_cna[is.na(tmp_cna)] = 2
  rownames(tmp_cna) = gsub('\\.', '-', rownames(tmp_cna))
  tmp_obj = DNA_sctc[[strsplit(file_name, '\\.')[[1]][1]]]
  tmp_clone = tmp_obj$map_obj$cell_map_pos#tmp_obj$clone_graph_obj$cell_pos_reduction
  tmp_cna$clone = tmp_clone[rownames(tmp_cna), 'clone']
  tmp_cna_clone = tmp_cna %>%
    group_by(clone) %>%
    summarise_all(list(mean)) %>% as.data.frame()
  rownames(tmp_cna_clone) = paste0(file_name,'_',tmp_cna_clone[,1])
  tmp_cna_clone = tmp_cna_clone[,-1]
  tmp_cna_clone = as.data.frame(tmp_cna_clone)
  same_gene = intersect(colnames(expdata), colnames(tmp_cna))

  cna_data[[file_name]] = tmp_cna_clone[, same_gene]#round(colMeans(tmp_cna[, same_gene]),0)
}
all_same_gene = c()
for(i in cna_data){
  if(length(all_same_gene)==0){
    all_same_gene = colnames(i)
  }else{
    all_same_gene = intersect(all_same_gene, colnames(i))
  }
}

srt[['low_CL']] = tolower(srt$cell_line)
clone_res = c()
for(i in cna_data){
  tmp_cnadata = i[,all_same_gene]
  CL = strsplit(rownames(tmp_cnadata), '\\.')[[1]][1]
  tmp_cnadata = t(as.data.frame(tmp_cnadata))
  cnadata_filter = tmp_cnadata[rowSums(tmp_cnadata==0)==0, ]

  #tmp_cells = rownames(srt@meta.data[srt$low_CL==CL,])
  #tmp_exp = as.matrix(expdata[tmp_cells, rownames(cnadata_filter)])
  tmp_srt = subset(srt, low_CL==CL)
  #expdata = AverageExpression(tmp_srt, group.by = 'sc3_cluster')$RNA
  #
  expdata = FetchData(tmp_srt, vars = VariableFeatures(tmp_srt), slot = 'counts')
  expdata = as.matrix(expdata[,rownames(cnadata_filter)])

  cal <- run_clonealign(expdata, cnadata_filter)
  clones <- cal$clone
  names(clones) = rownames(expdata)
  clone_res = c(clone_res, clones)
}



srt[['clonealign']] = clone_res[colnames(srt)]
##
DimPlot(srt, group.by = 'clonealign')#+scale_color_igv()

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

file_names = c(
  'huh7.txt',
  'plc.txt',
  #'hep3b.txt',
  'smc7721.txt',
  'hepg2.txt',
  'mhcc97h.txt',
  'mhcc97l.txt'
  #'lm3.txt'
)
cl_list = list()
cl_plot_list = list()
for(i in file_names){
  CL = strsplit(i, '\\.')[[1]][1]
  tmp_srt = subset(srt, low_CL==CL&clonealign!='unassigned')
  tmp_srt  = seurat_process(tmp_srt)
  cl_list[[i]] = tmp_srt
  cl_plot_list[[i]] = DimPlot(tmp_srt, group.by = 'clonealign')+scale_color_igv()
}
cowplot::plot_grid(plotlist = cl_plot_list, nrow=2)

saveRDS(srt,'./output/cell_line_srt_with_SC3_cytotrace_clonealign.RDS')




