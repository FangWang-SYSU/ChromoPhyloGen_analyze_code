library(Seurat)
library(infercnv)
library(ggplot2)
library(ggsci)
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

##### LUNG cancer ####
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung.rds')
DefaultAssay(srt) = 'RNA'
srt = process_srt(srt)
DimPlot(srt, group.by = c('seurat_clusters', 'patients', 'Cell_type.refined'), cols = pal_igv()(50))


# 加载CNV
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg38.txt', header = T)
cnv_b1 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_Lung/infercnv_outs/GSE123902/run.final.infercnv_obj')
cnv_b1_obs = cnv_b1@expr.data[, unlist(cnv_b1@observation_grouped_cell_indices)]
cnv_b2 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_Lung/infercnv_outs/GSE127465/run.final.infercnv_obj')
cnv_b2_obs = cnv_b2@expr.data[, unlist(cnv_b2@observation_grouped_cell_indices)]
cnv_b3 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_Lung/infercnv_outs/GSE131907/run.final.infercnv_obj')
cnv_b3_obs = cnv_b3@expr.data[, unlist(cnv_b3@observation_grouped_cell_indices)]
cnv_b4 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_Lung/infercnv_outs/PRJNA591860/run.final.infercnv_obj')
cnv_b4_obs = cnv_b4@expr.data[, unlist(cnv_b4@observation_grouped_cell_indices)]

cnv_b1_data = cnv_b1_obs[, intersect(colnames(srt), colnames(cnv_b1_obs))]
cnv_b2_data = cnv_b2_obs[, intersect(colnames(srt), colnames(cnv_b2_obs))]
cnv_b3_data = cnv_b3_obs[, intersect(colnames(srt), colnames(cnv_b3_obs))]
cnv_b4_data = cnv_b4_obs[, intersect(colnames(srt), colnames(cnv_b4_obs))]

same_gene = intersect(rownames(cnv_b1_data), rownames(cnv_b2_data))
same_gene = intersect(same_gene, rownames(cnv_b3_data))
same_gene = intersect(same_gene, rownames(cnv_b4_data))
same_gene = intersect(same_gene, gene_info$gene_name)

all_cnv = cbind(cnv_b1_data[same_gene,], cnv_b2_data[same_gene,], cnv_b3_data[same_gene,], cnv_b4_data[same_gene,])

epi_bc = colnames(all_cnv)
srt = subset(srt, cells = epi_bc)
srt = process_srt(srt)
DimPlot(srt, group.by = c('seurat_clusters', 'patients', 'Cell_type.refined'), cols = pal_igv()(50))
FeaturePlot(srt, features = c('PTPRC', 'CD3D', 'CD3E', 'KRT19', 'EPCAM'))
srt = FindClusters(srt, resolution = 0.1)
DimPlot(srt, group.by = c('seurat_clusters', 'patients', 'Cell_type.refined'), cols = pal_igv()(50))
srt = subset(srt, seurat_clusters%in%c(0,5,12), invert=T)
srt = process_srt(srt)
DimPlot(srt, group.by = c('seurat_clusters', 'patients', 'Cell_type.refined'), cols = pal_igv()(50))

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor.rds')
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/tumor.rds')
# 新的上皮barcode
epi_bc = colnames(srt)
all_cnv = all_cnv[, epi_bc]


gene_info = gene_info[gene_info$gene_name%in%same_gene, ]
gene_info = gene_info[!duplicated(gene_info$gene_name), ]
gene_info$seqnames = gsub('chr', '', gene_info$seqnames)
gene_info$id = paste0(gene_info$seqnames, '_', gene_info$start, '_',gene_info$end)
gene_info$len = gene_info$end-gene_info$start
rownames(gene_info) = gene_info$gene_name


for(i in unique(srt$patients)){
  print(i)
  tmp_bc = rownames(srt@meta.data[srt$patients==i, ])
  tmp_cnv = all_cnv[,tmp_bc]
  knum=10
  x = kmeans(t(tmp_cnv), centers = knum)
  xc = x$cluster

  new_cnv = c()
  for(xx in 1:knum){
    tmp_k_bc = xc[xc==xx]
    if(length(tmp_k_bc)>=5){
      tmp_k_cnv = t(tmp_cnv[,names(tmp_k_bc)])
      tmp_k_cnv2 = apply(tmp_k_cnv, 2, tmp_fun)
      # tmp_k_cnv2 = apply(tmp_k_cnv, 2, my_segment_sort, num=1)

      rownames(tmp_k_cnv2) = rownames(tmp_k_cnv)
      new_cnv = rbind(new_cnv, tmp_k_cnv2)
    }
  }
  new_cnv = wx_integer_CNV(t(new_cnv))
  # tmp_cnv = new_cnv
  tmp_cnv = tmp_cnv[, colnames(new_cnv)]
  tmp_cnv[new_cnv==2] = 1
  tmp_cnv = wx_integer_CNV(tmp_cnv)
  tmp_var = apply(tmp_cnv, 2, var)
  tmp_cnv = tmp_cnv[, tmp_var>quantile(tmp_var,0.25)]

  tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)

  #print(i)
  #tmp_bc = rownames(srt@meta.data[srt$patients==i, ])
  #tmp_cnv = all_cnv[,tmp_bc]
  #tmp_cnv = wx_integer_CNV(tmp_cnv)
  #tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)
  write.table(tmp_cnv, paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/Lung_res/patient_cnv/',i, '.txt'), quote=F, row.names = F)
}


#### ccRCC ######
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC.rds')
DefaultAssay(srt) = 'RNA'
# srt = process_srt(srt)

aa = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC.rds')
DefaultAssay(aa) = 'RNA'
aa = process_srt(aa)
aa2 = aa %>% subset(nFeature_RNA > 200 & percent.mt < 25 & nFeature_RNA<6000)

DimPlot(aa2)
DimPlot(aa, group.by = 'infercnv_type')


# 加载CNV

cnv_b1 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNV_ccRCC/infercnv_outs/GSE159115//run.final.infercnv_obj')
cnv_b1_obs = cnv_b1@expr.data[, unlist(cnv_b1@observation_grouped_cell_indices)]
cnv_b2 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNV_ccRCC/infercnv_outs/SCP1288//run.final.infercnv_obj')
cnv_b2_obs = cnv_b2@expr.data[, unlist(cnv_b2@observation_grouped_cell_indices)]
cnv_b3 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNV_ccRCC/infercnv_outs/SRZ190804//run.final.infercnv_obj')
cnv_b3_obs = cnv_b3@expr.data[, unlist(cnv_b3@observation_grouped_cell_indices)]

cnv_b1_data = cnv_b1_obs[, intersect(colnames(srt), colnames(cnv_b1_obs))]
cnv_b2_data = cnv_b2_obs[, intersect(colnames(srt), colnames(cnv_b2_obs))]
cnv_b3_data = cnv_b3_obs[, intersect(colnames(srt), colnames(cnv_b3_obs))]

same_gene = intersect(rownames(cnv_b1_data), rownames(cnv_b2_data))
same_gene = intersect(same_gene, rownames(cnv_b3_data))
same_gene = intersect(same_gene, gene_info$gene_name)

all_cnv = cbind(cnv_b1_data[same_gene,], cnv_b2_data[same_gene,], cnv_b3_data[same_gene,])

epi_bc = colnames(all_cnv)
srt = subset(srt, cells = epi_bc)
srt = process_srt(srt)
FeaturePlot(srt, features = c('PTPRC', 'CD3D', 'CD3E', 'KRT19', 'EPCAM'))
DimPlot(srt, group.by = c('seurat_clusters', 'patients', 'celltype_final'), cols = pal_igv()(50))

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/tumor.rds')

gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg38.txt', header = T)
gene_info = gene_info[gene_info$gene_name%in%same_gene, ]
gene_info = gene_info[!duplicated(gene_info$gene_name), ]
gene_info$seqnames = gsub('chr', '', gene_info$seqnames)
gene_info$id = paste0(gene_info$seqnames, '_', gene_info$start, '_',gene_info$end)
gene_info$len = gene_info$end-gene_info$start
rownames(gene_info) = gene_info$gene_name


for(i in unique(srt$Patient)){
  print(i)
  tmp_bc = rownames(srt@meta.data[srt$Patient==i, ])
  tmp_cnv = all_cnv[,tmp_bc]
  knum=10
  x = kmeans(t(tmp_cnv), centers = knum)
  xc = x$cluster

  new_cnv = c()
  for(xx in 1:knum){
    tmp_k_bc = xc[xc==xx]
    if(length(tmp_k_bc)>=5){
      tmp_k_cnv = t(tmp_cnv[,names(tmp_k_bc)])
      tmp_k_cnv2 = apply(tmp_k_cnv, 2, tmp_fun)
      # tmp_k_cnv2 = apply(tmp_k_cnv, 2, my_segment_sort, num=1)

      rownames(tmp_k_cnv2) = rownames(tmp_k_cnv)
      new_cnv = rbind(new_cnv, tmp_k_cnv2)
    }
  }
  new_cnv = wx_integer_CNV(t(new_cnv))
  # tmp_cnv = new_cnv
  tmp_cnv = tmp_cnv[, colnames(new_cnv)]
  tmp_cnv[new_cnv==2] = 1
  tmp_cnv = wx_integer_CNV(tmp_cnv)
  tmp_var = apply(tmp_cnv, 2, var)
  tmp_cnv = tmp_cnv[, tmp_var>quantile(tmp_var,0.25)]

  tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)
  #tmp_cnv = wx_integer_CNV(tmp_cnv)
  #tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)
  write.table(tmp_cnv, paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/ccRCC_res/patient_cnv/',i, '.txt'), quote=F, row.names = F)
}



#### CRC ######
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC.rds')
DefaultAssay(srt) = 'RNA'

all_cnv_list = list()
for(i in list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_CRC/cnv_PID/')){
  print(i)
  tmp_path = paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_CRC/cnv_PID/', i,'/run.final.infercnv_obj')
  tmp_cnv = readRDS(tmp_path)
  tmp_cnv = tmp_cnv@expr.data[, unlist(tmp_cnv@observation_grouped_cell_indices)]
  all_cnv_list[[i]] = tmp_cnv[, intersect(colnames(srt), colnames(tmp_cnv))]
}

same_gene = rownames(all_cnv_list$C105)
for(i in all_cnv_list){
  same_gene = intersect(same_gene, rownames(i))
}

all_cnv = c()
for(i in all_cnv_list){
  all_cnv = cbind(all_cnv, i[same_gene, ])
}

srt = subset(srt, cells=colnames(all_cnv))
srt = process_srt(srt)
FeaturePlot(srt, features = c('PTPRC', 'CD3D', 'CD3E', 'KRT19', 'EPCAM'))
DimPlot(srt, group.by = c('seurat_clusters', 'PID', 'celltype'), cols = pal_igv()(50))

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor.rds')

srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/tumor.rds')

gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg38.txt', header = T)
gene_info = gene_info[gene_info$gene_name%in%same_gene, ]
gene_info = gene_info[!duplicated(gene_info$gene_name), ]
gene_info$seqnames = gsub('chr', '', gene_info$seqnames)
gene_info$id = paste0(gene_info$seqnames, '_', gene_info$start, '_',gene_info$end)
gene_info$len = gene_info$end-gene_info$start
rownames(gene_info) = gene_info$gene_name

all_cnv = all_cnv[gene_info$gene_name, ]
for(i in unique(srt$PID)){
  print(i)
  tmp_bc = rownames(srt@meta.data[srt$PID==i, ])
  tmp_cnv = all_cnv[,tmp_bc]
  knum=10
  x = kmeans(t(tmp_cnv), centers = knum)
  xc = x$cluster

  new_cnv = c()
  for(xx in 1:knum){
    tmp_k_bc = xc[xc==xx]
    if(length(tmp_k_bc)>=5){
      tmp_k_cnv = t(tmp_cnv[,names(tmp_k_bc)])
      tmp_k_cnv2 = apply(tmp_k_cnv, 2, tmp_fun)
      # tmp_k_cnv2 = apply(tmp_k_cnv, 2, my_segment_sort, num=1)

      rownames(tmp_k_cnv2) = rownames(tmp_k_cnv)
      new_cnv = rbind(new_cnv, tmp_k_cnv2)
    }
  }
  new_cnv = wx_integer_CNV(t(new_cnv))
  # tmp_cnv = new_cnv
  tmp_cnv = tmp_cnv[, colnames(new_cnv)]
  tmp_cnv[new_cnv==2] = 1
  tmp_cnv = wx_integer_CNV(tmp_cnv)
  tmp_var = apply(tmp_cnv, 2, var)
  tmp_cnv = tmp_cnv[, tmp_var>quantile(tmp_var,0.25)]
  tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)

  #tmp_cnv = wx_integer_CNV(tmp_cnv)
  #tmp_cnv = cbind(gene_info[rownames(all_cnv), c('seqnames', 'start', 'end')], tmp_cnv)
  write.table(tmp_cnv, paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CRC_res/patient_cnv/',i, '.txt'), quote=F, row.names = F)
}


#### 更新拷贝数谱 #####
all_srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/pan_cancer.rds')

for(ct in c('Lung', 'CRC', 'ccRCC')){
  print(ct)
  tmp_meta = all_srt@meta.data[all_srt$CancerType==ct, ]
  for(pt in unique(tmp_meta$patient)){
    tmp_bc = rownames(tmp_meta[tmp_meta$patient==pt, ])
    tmp_cnv = read.table(paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/',ct,'_res/patient_cnv/',pt, '.txt'), header = T, check.names = F)
    tmp_cnv = tmp_cnv[, c('seqnames', 'start', 'end', tmp_bc)]
    write.table(tmp_cnv,
                paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/',ct,'_res/patient_cnv2/',pt, '.txt'),
                quote=F, row.names = F)
  }
}

