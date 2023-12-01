source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
library(Seurat)
process_srt = function(obj, dim=30, n.components=2,resolution=0.5, npcs=50){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "pca", n.components=n.components)
  return(obj)
}

# load pathway gene list
###### CRC #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv_infer/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc_all.rds')


## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv_infer/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc_all.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_score_sctc_all2.rds')

###### ccRCC_res #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv_infer/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_sctc_all.rds')


## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/tumor.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv_infer/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_sctc_all.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_score_sctc_all2.rds')

###### Lung_res #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv_infer/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_sctc_all.rds')

## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv_infer/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_sctc_all.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_score_sctc_all2.rds')




############ 第二次 #######
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
library(Seurat)
process_srt = function(obj, dim=30, n.components=2,resolution=0.5, npcs=50){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "pca", n.components=n.components)
  return(obj)
}

# load pathway gene list
###### CRC #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv2/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv_infer2/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc_all2.rds')


## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor_filter.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv_infer2/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc_all2.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_score_sctc_all2.rds')

###### ccRCC_res #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv2/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv_infer2/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_sctc_all2.rds')


## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/tumor_filter.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv_infer2/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_sctc_all2.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ccRCC_score_sctc_all2.rds')

###### Lung_res #####
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv2/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv_infer2/'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix%in%names(BC_list)){

  }else{
    sctc = tryCatch({
      scTraceClass(scTrace_dir,
                   prefix,
                   layout = 'slanted',
                   rate = 0.98,
                   min_cell_num=NULL,
                   min_cell_num2=NULL)
    },error = function(cond){
      return(NA)})
    BC_list[[prefix]] = sctc
  }
}
saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_sctc_all2.rds')

## score
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor_filter.rds')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv_infer2/'
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_sctc_all2.rds')
all_mechnism = c()
BC_list_with_score = list()

for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  if(is.na(sctc)){

  }else{
    CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
    all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
    all_limit_prop[all_limit_prop<0] = 0
    pvalue_thr = 0.001
    CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
    CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
    CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
    CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
    CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
    CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
    #
    all_rearrange_score2 = all_rearrange_score
    all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
    CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

    CNA_mechnism[is.na(CNA_mechnism)] = 0
    pde = predict_PDE(sctc)
    CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
    CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

    sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

    BC_list_with_score[[prefix]] = sctc

    CNA_mechnism$file = prefix
    all_mechnism = rbind(all_mechnism, CNA_mechnism)
  }
}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor_with_score2.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Lung_score_sctc_all2.rds')
















