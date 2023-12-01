library(Seurat)
library(infercnv)
library(ggplot2)
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
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC12.final_wy-230420.rds')

srt = subset(srt, cell_annotion =='Epithelial')
srt = process_srt(srt)
DimPlot(srt)
DimPlot(srt, group.by = c('patient'))
# cytotrace
library(CytoTRACE)
srt[['cytotrace_score2']] = 0#results$CytoTRACE[colnames(srt)]

for(i in unique(srt$patient)){
  print(i)
  tmp_srt = subset(srt, patient==i)
  if(sum(!is.na(tmp_srt$clone))>=2){
    tmp_srt = subset(tmp_srt, cells = rownames(tmp_srt@meta.data[!is.na(tmp_srt$clone), ]))
    count_mat = as.matrix(GetAssayData(tmp_srt, 'counts'))
    count_mat = AverageExpression(tmp_srt, group.by = 'clone')$RNA
    results <- CytoTRACE(count_mat[same_gene, ])

    srt@meta.data[colnames(tmp_srt), 'cytotrace_score2'] = results$CytoTRACE[tmp_srt$clone]
  }

  #srt@meta.data[colnames(tmp_srt), 'cytotrace_score2'] = results$CytoTRACE[colnames(tmp_srt)]
}
srt[['tmp_group']] = paste0(srt$file, '_',srt$clone)
count_mat = AverageExpression(srt, group.by = 'tmp_group')$RNA
results <- CytoTRACE(count_mat)
srt@meta.data[, 'cytotrace_score2'] = results$CytoTRACE[srt$tmp_group]

FeaturePlot(srt, features = 'cytotrace_score2') +
  scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))

##
count_mat = as.matrix(GetAssayData(srt, 'counts'))
results <- CytoTRACE(count_mat)
srt[['cytotrace_score']] = results$CytoTRACE[colnames(srt)]
FeaturePlot(srt, features = 'cytotrace_score') +
  scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))

srt = AddModuleScore(srt, features=list(c(cc.genes$s.genes, cc.genes$g2m.genes)), name='CellCycle')

FeaturePlot(srt, features = 'CellCycle1')
###
saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor.rds')


srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor.rds')


# 加载CNV
library(changepoint)
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg38.txt', header = T)
cnv_b1 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_BreastCancer/inferCNV/b1/run.final.infercnv_obj')
cnv_b2 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_BreastCancer/inferCNV/b2/run.final.infercnv_obj')

cnv_b1_data = cnv_b1@expr.data[, intersect(colnames(srt), colnames(cnv_b1@expr.data))]
cnv_b2_data = cnv_b2@expr.data[, intersect(colnames(srt), colnames(cnv_b2@expr.data))]

# cnv_b1_data = wx_integer_CNV(cnv_b1_data)
# cnv_b1_data = t(cnv_b1_data)
#
# pheatmap::pheatmap(cnv_b1_data[1:100, ], cluster_rows = F, cluster_cols = F)

my_segment <- function(x,num=50){
  ansvar=cpt.meanvar(x, method='BinSeg', Q=num)
  seg_data = c()#gene_fc_data[ansvar@cpts, ]
  cpts_res = c(0,ansvar@cpts)
  for(i in 1:(length(cpts_res)-1)){
    s = cpts_res[i]
    e = cpts_res[i+1]
    xx = table(x[s:e])
    # seg_data = c(seg_data, rep(as.numeric(names(xx)[which.max(xx)]), e-s))
    seg_data = c(seg_data, rep(ansvar@param.est[["mean"]][i], e-s))
  }
  return(seg_data)
}
my_segment_sort <- function(x,num=5){
  # old_x = x
  x_order = names(x)
  x = sort(x)
  new_names = names(x)
  ansvar=cpt.meanvar(x, method='BinSeg', Q=num)
  seg_data = c()#gene_fc_data[ansvar@cpts, ]
  cpts_res = c(0,ansvar@cpts)
  for(i in 1:(length(cpts_res)-1)){
    s = cpts_res[i]
    e = cpts_res[i+1]
    xx = table(x[s:e])
    # seg_data = c(seg_data, rep(as.numeric(names(xx)[which.max(xx)]), e-s))
    seg_data = c(seg_data, rep(ansvar@param.est[["mean"]][i], e-s))
  }
  names(seg_data) = new_names
  seg_data = seg_data[x_order]
  return(seg_data)
}
# seg_res = t(apply(cnv_b1_data[1:100, ], 1, my_segment))
#
# pheatmap::pheatmap(seg_res, cluster_rows = F, cluster_cols = F)
#
# x = cnv_b1_data[200, ]
# ansvar=cpt.mean(x, method='BinSeg',penalty='MBIC',pen.value=1, Q=500)
# plot(ansvar)
#

same_gene = intersect(rownames(cnv_b1_data), rownames(cnv_b2_data))
same_gene = intersect(same_gene, gene_info$gene_name)




all_cnv = cbind(cnv_b1_data[same_gene,], cnv_b2_data[same_gene,])
gene_info = gene_info[gene_info$gene_name%in%same_gene, ]
gene_info = gene_info[!duplicated(gene_info$gene_name), ]
gene_info$seqnames = gsub('chr', '', gene_info$seqnames)
gene_info$id = paste0(gene_info$seqnames, '_', gene_info$start, '_',gene_info$end)
gene_info$len = gene_info$end-gene_info$start
rownames(gene_info) = gene_info$gene_name

library(DNAcopy)
library(dplyr)

tmp_fun <-function(x){
  rep(mean(x), length(x))
}

for(i in unique(srt$patient)){
  print(i)
  tmp_bc = rownames(srt@meta.data[srt$patient==i, ])
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

  write.table(tmp_cnv, paste0('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/patient_cnv_seg/',i, '.txt'),
              quote=F, row.names = F)
}

FeaturePlot(srt, features = 'CNV_SD')


x = tmp_cnv[3,]
x=sort(as.vector(as.matrix(tmp_cnv)))

plot(density(x))

aa = sample(x,1000)
ansvar=cpt.meanvar(sort(aa), method='BinSeg', Q=5)
plot(ansvar)
#x2=my_segment_sort(x)
#plot(x2)
#cpm.res = processStream(x, cpmType = "Cramer-von-Mises")
# 可视化变点
#plot(x, type = "l", col = "steelblue", lwd = 2)
#abline(v = cpm.res$changePoints, lwd = .5, col = "red")



CNA.object <- CNA(cbind(log(tmp_cnv)),
                  gene_info[rownames(all_cnv), 'seqnames'],
                  gene_info[rownames(all_cnv), 'start'],
                  data.type="logratio",
                  sampleid=colnames(tmp_cnv))
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=0)
seg_cnv_res = cbind(segment.smoothed.CNA.object$output, segment.smoothed.CNA.object$segRows)

new_cnv_ratio = matrix(0, nrow=nrow(tmp_cnv), ncol=5)
colnames(new_cnv_ratio) = c('chr', 'start',colnames(tmp_cnv))
for(i in 1:nrow(seg_cnv_res)){
  s = seg_cnv_res$startRow[i]
  e = seg_cnv_res$endRow[i]
  bc = seg_cnv_res$ID[i]
  new_cnv_ratio[s:e, bc] = seg_cnv_res$seg.mean[i]
  new_cnv_ratio[s:e, 'chr'] = as.numeric(seg_cnv_res[i, c('chrom')])
  new_cnv_ratio[s:e, 'start'] = as.numeric(seg_cnv_res[i, c('loc.start')])
}
new_cnv_ratio = new_cnv_ratio[order(as.numeric(new_cnv_ratio[,1]), as.numeric(new_cnv_ratio[,2])), ]
new_cnv_ratio = new_cnv_ratio[,-c(1,2)]
# new_cnv_ratio = apply(new_cnv_ratio, 2, my_segment)
# order_new = order(as.numeric(seg_cnv_res$chrom), as.numeric(seg_cnv_res$loc.start))
# tt = rbind(t(new_cnv_ratio),
#            t(log(tmp_cnv[,1:3]))
#            )
tt = wx_integer_CNV(exp(tt))

pheatmap::pheatmap(t(tmp_cnv)[1:100,],cluster_cols = F, cluster_rows = F,show_colnames = F)

min_seg = mean(gene_info$len)
new_gene_info = c()
for(i in 1:nrow(gene_info)){
  print(i)
  tmp_len = gene_info[i, 'len']
  if(tmp_len<=min_seg){
    tmp_num = 1
  }else{
    tmp_num = round(tmp_len/min_seg)
  }
  tmp = rep(unlist(gene_info[i, c('gene_name', 'id')]), tmp_num)
  tmp = matrix(tmp, ncol = 2, byrow = T)
  new_gene_info = rbind(new_gene_info, tmp)
}



