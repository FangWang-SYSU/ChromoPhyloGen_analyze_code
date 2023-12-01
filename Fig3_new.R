source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')

dna_file = c(
  'uber.t12.seg.txt',
  'uber.t11.seg.txt',
  'uber.t10.seg.txt',
  'uber.t9.seg.txt',
  'uber.t8.seg.txt',
  'uber.t7.seg.txt',
  'uber.t6.seg.txt',
  'uber.t5.seg.txt',
  'uber.t4.seg.txt',
  'uber.t3.seg.txt',
  'uber.t2.seg.txt',
  'uber.t1.seg.txt',
  'KTN615_cleaned_uber.seg.txt',
  'KTN302_cleaned_uber.seg.txt',
  'KTN206_cleaned_uber.seg.txt',
  'KTN152_cleaned_uber.seg.txt',
  'KTN132_cleaned_uber.seg.txt',
  'KTN129_cleaned_uber.seg.txt',
  'KTN126_cleaned_uber.seg.txt',
  'KTN102_cleaned_uber.seg.txt')

#
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sub/'
TNBC_list = list()
for(prefix in dna_file){
  print(prefix)
  if(prefix=="uber.t3.seg.txt"){
    min_cell_num = 50
  }else if(prefix=='KTN615_cleaned_uber.seg.txt'){
    min_cell_num=2
  }
  else{
    min_cell_num=30
  }
  #aa = load_scTrace(scTrace_dir, prefix)
  sctc = scTraceClass(scTrace_dir,
                      prefix,
                      layout = 'slanted',
                      rate = 0.95,
                      min_cell_num=min_cell_num,
                      min_cell_num2=min_cell_num)
  TNBC_list[[prefix]] = sctc
}

saveRDS(TNBC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all.rds')
TNBC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all.rds')
library(ggnewscale)
library(ggtreeExtra)
library(ggpubr)

#### 计算score ######
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sub2/'

TNBC_list_with_score = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = TNBC_list[[prefix]]
  CNA_mechnism = read.table(glue('{scTrace_dir}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop[all_limit_prop<0] = 0
  pvalue_thr = 0.001
  CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
  CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr& all_rearrange_score>0)
  CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5& all_rearrange_score>0)
  CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5& all_rearrange_score>0)
  CNA_mechnism$BFB = log1p(CNA_mechnism$BFB)# / ncol(sctc$orig.data$all_node_data)
  #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5& all_rearrange_score>0)] = NA
  CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
  #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5& all_rearrange_score>0)] = NA
  CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

  CNA_mechnism[is.na(CNA_mechnism)] = 0
  pde = predict_CEE(sctc)
  CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
  CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

  #CNA_mechnism = define_CNA_mechnism(sctc,
  #                                   cancer_type='BRCA',
  #                                   chromothripsis_chr_thr=4,
  #                                   chromoplexy_chr_thr=7,
  #                                   quantile_thr=0.75,
  #                                   min_length=5*1000*1000,
  #                                   random = FALSE,
  #                                   random_num = 1000)
  # pheatmap::pheatmap(-log(CNA_mechnism$orig.rearrange_score$all_rearrange_score_pvalue))
  # pheatmap::pheatmap(-log(CNA_mechnism$orig.rearrange_score$all_limit_prop_pvalue))

  sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]
  ##### 单细胞进化树 rearrange着色 ####
  colorby = 'rearrange_score'
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
  rearrange_p=gtree_res +
    aes(color=group, size=0.2)+
    scale_color_gradientn(colours = c('yellow','red'))+
    guides(color=guide_colourbar(title='Rearrange_score'))+
    theme_tree(legend.position='right')+
    #coord_flip()+
    #scale_x_reverse()+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(2, 'mm'),
          legend.title = element_text(size=6)
    )
  ##### 单细胞进化树 BFB着色 ####
  colorby = 'BFB'
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
  BFB_p=gtree_res +
    aes(color=group, size=isTip)+
    scale_color_gradientn(colours = c('yellow', 'red'))+
    theme_tree(legend.position='right')+
    guides(color=guide_colourbar(title='BFB_prop'))+
    #coord_flip()+
    #scale_x_reverse()+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(2, 'mm'),
          legend.title = element_text(size=6)
    )
  ##### 单细胞进化树 limit_num ####
  colorby = 'limit_num'
  #cell_info = sctc$map_obj$cell_map_pos[, colorby, drop=F]
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
  limit_p=gtree_res +
    aes(color=group, size=isTip)+
    scale_color_gradientn(colours = c('yellow', 'red'))+
    guides(color=guide_colourbar(title='limit_num'))+
    theme_tree(legend.position='right')+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(1, 'mm'),
          legend.title = element_text(size=6)
    )
  ##### 单细胞进化树 seismic_num ####
  colorby = 'seismic_num'
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
  seismic_p=gtree_res +
    aes(color=group, size=isTip)+
    scale_color_gradientn(colours = c('yellow', 'red'))+
    guides(color=guide_colourbar(title='seismic_num'))+
    theme_tree(legend.position='right')+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(1, 'mm'),
          legend.title = element_text(size=6)
    )

  ## pde ##
  colorby = 'pde'
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]

  pde_p=gtree_res +
    aes(color=group, size=isTip)+
    scale_color_gradientn(colours = c('yellow', 'red'))+
    theme_tree(legend.position='right')+
    guides(color=guide_colourbar(title='CEE'))+
    #coord_flip()+
    #scale_x_reverse()+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(2, 'mm'),
          legend.title = element_text(size=6)
    )

  # clone
  colorby = 'clone'
  cell_info = sctc$map_obj$cell_map_pos[!grepl('root|virtual', rownames(sctc$map_obj$cell_map_pos)), colorby, drop=F]
  gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
  gtree_res$data$group = cell_info[gtree_res$data$label, colorby]

  clone_col = setNames(pal_igv()(length(unique(cell_info$clone))), unique(cell_info$clone))
  clone_col = c(clone_col)
  clone_p=gtree_res +
    aes(color=group, size=0.2)+
    scale_color_manual(values=clone_col)+
    # scale_color_gradientn(colours = c('blue','yellow', 'red'))+
    theme_tree(legend.position='right')+
    guides(color=guide_legend(title='Clone'))+
    #coord_flip()+
    #scale_x_reverse()+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(1, 'mm'),
          legend.title = element_text(size=6)
    )


  p = cowplot::plot_grid(clone_p, rearrange_p, BFB_p, limit_p, seismic_p, pde_p, nrow = 1, align = 'v')

  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/进化树不同指标着色/tree_', prefix,'.pdf'),
         p, width=220, height = 60, units='mm', dpi = 600)

  # p = clone_p

  gtree_res = ggtree(sctc$orig.data$tree)
  gtree_res$data$group = cell_info[gtree_res$data$label, 'clone']
  gtree_res = ggtree(gtree_res$data, aes(color=group), size=0.2)

  clone_col = setNames(pal_igv()(length(unique(cell_info$clone))), unique(cell_info$clone))
  clone_col = c(clone_col)
  p=gtree_res +
    scale_color_manual(values=clone_col)+
    # scale_color_gradientn(colours = c('blue','yellow', 'red'))+
    theme_tree(legend.position='right')+
    guides(color=guide_legend(title='Clone'))+
    #coord_flip()+
    #scale_x_reverse()+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(1, 'mm'),
          legend.title = element_text(size=6)
    )




  score_data = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'wgd', 'CPS_num','BFB',
                                             'limit_num', 'seismic_num', 'limit_score','seismic_score', 'pde')]
  score_data$nameid = rownames(score_data)
  score_data$wgd = factor(score_data$wgd, c('WGD0', 'WGD1', 'WGD2'))
  #tile_data = reshape2::melt(score_data[, c(2,5,6,7,8)], id='nameid')
  score_data = score_data[p$data$label[p$data$isTip], ]
  score_data$x1=0
  #score_data$x2=score_data$x1+1
  bar_data = reshape2::melt(score_data[, c('limit_num', 'seismic_num', 'nameid')], id='nameid')
  bar_data$variable = factor(bar_data$variable, levels = c( 'seismic_num','limit_num'))
  dd <- p$data
  dd <- dd[dd$isTip,]
  lab <- dd$label[order(dd$y, decreasing = T)]
  score_data = score_data[lab, ]
  p2 = p+geom_tippoint(size=1)
  p2 = p2 +
    # new_scale_fill() +
    # geom_fruit(data=score_data, geom=geom_tile,
    #            mapping=aes(y=nameid, x=x1,  fill=wgd),
    #            color = "grey",
    #            offset = 0.05,size = 0.05, width=10,
    #            axis.params=list(
    #              axis="x",
    #              text = "WGD",
    #              text.angle = 45,
    #              text.size = 2,
    #              line.size = 0,
    #              vjust = 1,
    #              hjust= 1,
    #              inherit.aes=F,
    #            ))+
    # # scale_fill_npg()+
    # scale_fill_manual(values = c('WGD0'='#DEDEDE', 'WGD1'='#FFE53B', 'WGD2'='#FF2525'),drop = FALSE)+
    #
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=rearrange_score),
               color = "grey",
               offset = 0.05,size = 0.05, width=10,
               axis.params=list(
                 axis="x",
                 text = "RE",
                 text.angle = 45,
                 text.size = 2,
                 line.size = 0,
                 vjust = 1,
                 hjust= 1,
                 inherit.aes=F,
               ))+
    scale_fill_gradientn(colours = c('white', '#306DAA'))+
    #
    # new_scale_fill() +
    # geom_fruit(data=score_data, geom=geom_tile,
    #            mapping=aes(y=nameid, x=x1,  fill=limit_score),
    #            color = "grey",
    #            offset = 0.05,size = 0.05, width=10,
    #            axis.params=list(
    #              axis="x",
    #              text = "Limit",
    #              text.angle = 45,
    #              text.size = 2,
    #              line.size = 0,
    #              vjust = 1,
    #              hjust= 1,
    #              inherit.aes=F,
    #            ))+
    # scale_fill_gradientn(colours = c('white', '#08AEEA'))+
    #
    # new_scale_fill() +
    # geom_fruit(data=score_data, geom=geom_tile,
    #            mapping=aes(y=nameid, x=x1,  fill=seismic_score),
    #            color = "grey",
    #            offset = 0.05,size = 0.05, width=10,
    #            axis.params=list(
    #              axis="x",
    #              text = "Seismic",
    #              text.angle = 45,
    #              text.size = 2,
    #              line.size = 0,
    #              vjust = 1,
    #              hjust= 1,
    #              inherit.aes=F,
    #            ))+
    # scale_fill_gradientn(colours = c('white', '#B721FF'))+
    #
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=BFB),
               color = "grey",
               offset = 0.05,size = 0.05, width=10,
               axis.params=list(
                 axis="x",
                 text = "BFB",
                 text.angle = 45,
                 text.size = 2,
                 line.size = 0,
                 vjust = 1,
                 hjust= 1,
                 inherit.aes=F,
               ))+
    scale_fill_gradientn(colours = c('white', '#F4BA19'))+
    #
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=pde),
               color = "grey",
               offset = 0.05,size = 0.05, width=10,
               axis.params=list(
                 axis="x",
                 text = "CEE",
                 text.angle = 45,
                 text.size = 2,
                 line.size = 0,
                 vjust = 1,
                 hjust= 1,
                 inherit.aes=F,
               ))+
    scale_fill_gradientn(colours = c('white', '#8DC469'))+
    #
    new_scale_fill() +
    geom_fruit(mapping=aes(y=nameid, x=value, fill=variable),
               data=bar_data,
               geom=geom_bar,
               offset = 0.05,
               pwidth=0.3,
               size=0.01,
               color='gray',
               #fill='red',
               orientation="y",
               stat="identity",
               axis.params=list(
                 axis="x",
                 #text.angle = -45,
                 text.size = 2,
                 #line.size = 0,
                 vjust = 1,
                 inherit.aes=F,
               )
    )+
    scale_fill_manual(values = c('limit_num'='#08AEEA', 'seismic_num'='#B721FF'),drop = FALSE, name='chr_num')
  p2 = p2+    theme(legend.text = element_text(size=4),
                    legend.key.size = unit(2, 'mm'),
                    legend.title = element_text(size=6)
  )
  p2 = p2+scale_y_continuous(expand = expansion(mult=0.1))
  #p2+guides(fill=guide_legend(ncol =2),color=guide_legend(ncol =2))
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/进化树不同指标着色/hp_', prefix,'.pdf'),
         p2, width=160, height = 160, units='mm', dpi = 600)

  TNBC_list_with_score[[prefix]] = sctc
}

saveRDS(TNBC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all_with_score.rds')
TNBC_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all_with_score.rds')


#### 降维显示 #####
cna_gene = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_cna_gene.rds')
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')

same_gene = unique(gene_info$gene_name)
for(prefix in dna_file){
  same_gene = intersect(same_gene, colnames(cna_gene[[prefix]]))
}
all_cell_data = c()
all_cell_meta = c()
for(prefix in dna_file){
  print(prefix)
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = cna_gene[[prefix]][, same_gene]
  tmp_cna$clone = sctc$map_obj$cell_map_pos[rownames(tmp_cna), 'clone']
  tmp_cna = tmp_cna %>%
    group_by(clone) %>%
    summarise_all(list(mean)) %>%
    as.data.frame()
  rownames(tmp_cna) = tmp_cna[,1]
  tmp_cna = tmp_cna[,-1]
  all_cell_meta = rbind(all_cell_meta, data.frame(file=prefix, clone=rownames(tmp_cna)))
  all_cell_data = rbind(all_cell_data, tmp_cna)
}

aa = list(all_cell_data=all_cell_data, all_cell_meta=all_cell_meta)
saveRDS(aa, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_all_cell.rds')
aa = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_all_cell.rds')
all_cell_data = aa$all_cell_data
all_cell_meta = aa$all_cell_meta

library(Seurat)
filter_data = all_cell_data
rownames(filter_data) = paste0(all_cell_meta[,1], '-',all_cell_meta[,2])
all_cell_meta2 = c()
for(prefix in dna_file){
  sctc = TNBC_list_with_score[[prefix]]
  leaves_node = rownames(sctc$orig.data$all_node_data)
  leaves_node = leaves_node[!grepl('root|virtual', leaves_node)]
  tmp_meta = sctc$map_obj$cell_map_pos[leaves_node, ]
  tmp_meta$patient = prefix
  all_cell_meta2 = rbind(all_cell_meta2, tmp_meta)
}
all_cell_meta2 = all_cell_meta2 %>%
  dplyr::group_by(clone, patient) %>%
  dplyr::summarise(BFB=mean(BFB), rearrange_score=mean(rearrange_score), pde=mean(pde),limit_num=mean(limit_num), seismic_num=mean(seismic_num)) %>%
  as.data.frame()
rownames(all_cell_meta2) = paste0(all_cell_meta2[,2], '-',all_cell_meta2[,1])

####
srt = CreateSeuratObject(count=t(filter_data))
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
srt = process_srt(srt)
srt = AddMetaData(srt, all_cell_meta2)
srt[['cee_group']] = ifelse(srt$pde>median(srt$pde), 'High', 'Low')


# library(umap)
# filter_data = all_cell_data
# filter_meta = all_cell_meta
# var_gene = names(sort(apply(filter_data, 2, var), decreasing = T))
# var_gene = var_gene[!grepl('ENSG', var_gene)]

#var_clone = apply(filter_data, 1, var)
#filter_clone = names(var_clone[var_clone>0.05])
#filter_meta = filter_meta[var_clone>0.05, ]

# umap_res = umap(all_cell_data[, head(var_gene, 20000)])
# umap_res_df = umap_res$layout %>% as.data.frame()
# colnames(umap_res_df) = c('UMAP1', 'UMAP2')
# umap_res_df %>%
#   ggplot(aes(x=UMAP1, y=UMAP2)) +
#   geom_point(size=4)
#
# pca_res = prcomp(t(all_cell_data[, var_gene]))
# umap_res_df = pca_res$rotation[,c(3,5)] %>% as.data.frame()
# colnames(umap_res_df) = c('PC1', 'PC2')
#
#
# umap_res = umap(pca_res$rotation[,4:5] %>% as.data.frame())
# umap_res_df = umap_res$layout %>% as.data.frame()
# colnames(umap_res_df) = c('UMAP1', 'UMAP2')
#umap_res_df = as.data.frame(cbind(umap_res_df, all_cell_meta))
#umap_res_df %>%
#  ggplot(aes(x=UMAP1, y=UMAP2, color=file)) +
#  geom_point(size=4)+
#  scale_color_igv()

#colnames(umap_res_df) = c('PC1', 'PC2')
#rownames(umap_res_df) = paste0(filter_meta[,1], '-',filter_meta[,2])

#umap_res_df = as.data.frame(cbind(umap_res_df, all_cell_meta2[rownames(umap_res_df), ]))

umap_res_df = srt@reductions$umap@cell.embeddings %>% as.data.frame()
umap_res_df = as.data.frame(cbind(umap_res_df, srt@meta.data))
a=umap_res_df %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=patient)) +
  geom_point(size=2, shape=21, stroke=0.01)+
  scale_fill_igv()+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a
b=umap_res_df %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=pde)) +
  geom_point(size=2, shape=21, stroke=0.01)+
  scale_fill_gradientn(colours = c('white', '#8DC469'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a+b


#ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap_patient.pdf'),
#       a, width=76, height = 46, units='mm', dpi = 600)
#ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap_cee.pdf'),
#       b, width=50, height = 40, units='mm', dpi = 600)

#ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap_cee.pdf'),
#       a+b, width=120, height = 46, units='mm', dpi = 600)
#
c=umap_res_df %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=rearrange_score)) +
  geom_point(size=2, shape=21, stroke=0.01)+
  scale_fill_gradientn(colours = c('white', '#306DAA'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
c

d=umap_res_df %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=BFB)) +
  geom_point(size=2, shape=21, stroke=0.01)+
  scale_fill_gradientn(colours = c('white', '#F4BA19'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
d

a+b+c+d
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap.pdf'),
       a+b+c+d, width=130, height = 80, units='mm', dpi = 600)

#umap_res_df$cee_group = ifelse(umap_res_df$pde>median(umap_res_df$pde), 'High', 'Low')
e=umap_res_df %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=cee_group)) +
  geom_point(size=2, shape=21, stroke=0.01)+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
e
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap2.pdf'),
       a+e+c+d, width=130, height = 80, units='mm', dpi = 600)

ae = umap_res_df %>%
  mutate(patient=gsub('uber.', '',patient)) %>%
  mutate(patient=gsub('_cleaned.seg.txt', '',patient)) %>%
  mutate(patient=gsub('.seg.txt', '',patient)) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=patient, shape=cee_group, size=cee_group)) +
  geom_point(stroke=0.01)+
  #scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_shape_manual(values=c('High'=24, 'Low'=21))+
  scale_size_manual(values=c('High'=2, 'Low'=1))+
  scale_fill_igv()+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
ae
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_umap3.pdf'),
       cowplot::plot_grid(ae,c,d, nrow=2), width=120, height = 80, units='mm', dpi = 600)

a=umap_res_df %>%
  ggplot(aes(x=cee_group, y=rearrange_score, fill=cee_group)) +
  geom_boxplot(size=0.2, outlier.shape = NA, width=0.4)+
  geom_point(shape=21, size=1, stroke=0.1)+
  labs(y='Chromothripsis')+
  stat_compare_means(comparisons = list(c('High', 'Low')),label = 'p.signif')+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(0.2))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))

b=umap_res_df %>%
  ggplot(aes(x=cee_group, y=limit_num, fill=cee_group)) +
  geom_boxplot(size=0.2, outlier.shape = NA, width=0.4)+
  geom_point(shape=21, size=1, stroke=0.1)+
  labs(y='Gradual')+
  stat_compare_means(comparisons = list(c('High', 'Low')),label = 'p.signif')+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(0.2))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
c=umap_res_df %>%
  ggplot(aes(x=cee_group, y=seismic_num, fill=cee_group)) +
  geom_boxplot(size=0.2, outlier.shape = NA, width=0.4)+
  geom_point(shape=21, size=1, stroke=0.1)+
  labs(y='Seismic')+
  stat_compare_means(comparisons = list(c('High', 'Low')),label = 'p.signif')+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(0.2))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))

a+b+c

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_group_指标.pdf'),
       a+b+c, width=80, height = 40, units='mm', dpi = 600)


# 热图
pca_res = prcomp(t(all_cell_data))
#x = as.matrix(dist(all_cell_data))
#x[x==0] = NA
#x[is.na(x)] = min(x, na.rm = T)
#x = 1/x

x = cor(t(pca_res$rotation[,1:50]))
#x[is.na(x)] = min(x, na.rm = T)

rownames(x) = paste0(all_cell_meta[,1], '-',all_cell_meta[,2])
colnames(x) = paste0(all_cell_meta[,1], '-',all_cell_meta[,2])
width_in = 60 / 25.4
height_in = 60 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_PDE分组异质性热图.pdf', width=width_in, height=height_in)
Heatmap(x[rownames(umap_res_df), rownames(umap_res_df)],
        show_row_dend = F, show_column_dend = F,
        col=circlize::colorRamp2(c(-0.3,0,  0.3), c("blue","white", "red")),
        width = ncol(x)*unit(0.6, "mm"), # 宽度
        height = nrow(x)*unit(0.6, "mm"), # 高度
        show_row_names = F, show_column_names = F,
        row_split = factor(umap_res_df$cee_group, c('Low','High')),
        column_split = factor(umap_res_df$cee_group, c('High', 'Low')),
        )
dev.off()

hp = Heatmap(x[rownames(umap_res_df), rownames(umap_res_df)],
             show_row_dend = F, show_column_dend = F,
             col=circlize::colorRamp2(c(-0.3,0,  0.3), c("blue","white", "red")),
             width = ncol(x)*unit(0.6, "mm"), # 宽度
             height = nrow(x)*unit(0.6, "mm"), # 高度
             show_row_names = F, show_column_names = F,
             row_split = factor(umap_res_df$cee_group, c('Low','High')),
             column_split = factor(umap_res_df$cee_group, c('High', 'Low')),
)
hc_hp = row_dend(hp)
new_hc = merge(hc_hp$Low, hc_hp$High)
new_hc = as.hclust(new_hc)


### clone 层次聚类
library(ggtree)
library(ape)

all_cell_data2 = round(all_cell_data)
rownames(all_cell_data2) = paste0(all_cell_meta[,1], '-',all_cell_meta[,2])
#random_tree = rtree(ncol(aa), tip.label = sample(colnames(aa), ncol(aa)))
pyd = as.phyDat(as.matrix(all_cell_data2), type='USER', levels=0:max(all_cell_data2))
ham_dist = dist.hamming(pyd)
# MP_tree1 = optim.parsimony(random_tree, pyd)

hc = hclust(ham_dist)#dist(all_cell_data))
hc_tree = as.phylo(hc)
hc_tree = acctran(hc_tree, pyd)
hc_tree = as.phylo(new_hc)
gtree_df = fortify(hc_tree)
gtree_df[, c('patient', 'pde', 'rearrange_score')] = all_cell_meta2[gtree_df$label, c('patient', 'pde', 'rearrange_score')]
gtree_df$cee_group = ifelse(gtree_df$pde>median(gtree_df$pde, na.rm = T), 'High', 'Low')

gtree = ggtree(gtree_df,size=0.1)+
  geom_tippoint(mapping =aes(fill=patient, shape=cee_group), size=2, stroke=0.01)+
  scale_shape_manual(values=c('High'=24, 'Low'=21))+
  scale_fill_igv()+
  coord_flip()+
  scale_x_reverse()
gtree
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_tree_hp.pdf'),
       gtree, width=230, height = 30, units='mm', dpi = 600)


#### 统计高低两组病人种类分布 #####
tmp_meta= gtree_df %>% as.data.frame()
tmp_meta %>%
  na.omit() %>%
  group_by(cee_group) %>%
  summarise(length(unique(patient)))

# clone水平累计分布 #
# high
CNV_envents = filter_data
srt[['cnv_event']] = (rowSums((CNV_envents-2)!=0)/ncol(CNV_envents))[colnames(srt)]

tmp_func_share<-function(x){
  tmp_data = CNV_envents[x,,drop=F]
  tmp_res = (rowMeans(t(tmp_data))-apply(t(tmp_data), 1, min))==0
  return(sum(tmp_res)/ncol(tmp_data))
}
plot_df_high = srt@meta.data[srt@meta.data$cee_group=='High',]
plot_df_high$bc = rownames(plot_df_high)
plot_df_high = arrange(plot_df_high, cnv_event)
high_rate = c()
for(i in 1:nrow(plot_df_high)){
  print(i)
  bc = plot_df_high$bc[1:i]
  high_rate = rbind(high_rate, c(tmp_func_share(bc)))
}
high_rate = as.data.frame(high_rate)
colnames(high_rate) = c('share_CNA')
high_rate$group='High'
high_rate$cell_num = (1:nrow(high_rate))/nrow(high_rate)
# low
plot_df_low = srt@meta.data[srt@meta.data$cee_group=='Low',]
plot_df_low$bc = rownames(plot_df_low)
plot_df_low = arrange(plot_df_low, cnv_event)
low_rate = c()
for(i in 1:nrow(plot_df_low)){
  print(i)
  bc = plot_df_low$bc[1:i]
  low_rate = rbind(low_rate, c(tmp_func_share(bc)))
}
low_rate = as.data.frame(low_rate)
colnames(low_rate) = c('share_CNA')
low_rate$group='Low'
low_rate$cell_num = (1:nrow(low_rate))/nrow(low_rate)

plot_df2 = rbind(high_rate, low_rate)

a=reshape2::melt(plot_df2, id=c('group', 'cell_num')) %>%
  ggplot(aes(x=cell_num, y=value, color=group))+
  geom_line(size=0.5)+
  #geom_smooth(linewidth=1, se=F)+
  labs(x='Cumulative percentage of cells')+
  facet_wrap(~variable, scales = 'free_y', nrow=1)+
  labs(y='score')+
  scale_color_npg()+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone累计分布.pdf'), a,
       width=40, height=40, units='mm', dpi = 600, bg = 'transparent')


#### 计算gennome与指标相关性 ######
genom_index_plot_list = list()
genom_index_cor = c()
for(prefix in dna_file){
  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'clone')]
  meta2 = sctc$orig.data$cell_relation[rownames(meta),]
  # genome_cnv = rowSums((sctc$orig.data$all_node_data-2)!=0)[rownames(meta)] / ncol(sctc$orig.data$all_node_data)
  # 2
  genome_cnv = (meta2$Parent_gain_loc + meta2$Parent_loss_loc)/ ncol(sctc$orig.data$all_node_data)#(meta2$Root_gain_cn + meta2$Root_loss_loc) / (meta2$Root_gain_loc + meta2$Root_loss_loc)

  meta$genome_cnas = genome_cnv
  meta = meta[!grepl('root|virtual', rownames(meta)), ]

  plot_data = reshape2::melt(meta, id=c('clone', 'genome_cnas'))

  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  a=plot_data %>%
    ggplot(aes(x=value, y=genome_cnas))+
    geom_point(aes(color=clone), size=1)+
    stat_cor()+
    geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
    scale_color_manual(values = clone_col)+
    facet_wrap(~variable, scales = 'free_x')+
    labs(y='%genome_altered', title=new_prefix)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5))
  genom_index_plot_list[[new_prefix]] = a
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/指标与基因组变异相关性/',prefix,'.pdf'),
         a, width=160, height = 60, units='mm', dpi = 600)
  tmp_func = function(x,y){
    tmp_res = cor.test(x, y)
    return(cbind(tmp_res$estimate, tmp_res$p.value))
  }
  cor_res = plot_data %>%
    group_by(variable) %>%
    summarise(tmp_func(genome_cnas, value)) %>% as.data.frame()
  cor_res = as.data.frame(cbind(as.vector(cor_res[,1]), cor_res[,2]))
  cor_res$V4 = new_prefix
  genom_index_cor = rbind(genom_index_cor, cor_res)
}

p_order = c("KTN152", "KTN132", "t2" ,    "KTN126", "KTN206", "t7" ,    "t11"  ,  "KTN129", "t4" ,    "t8" ,    "KTN615", "t6" ,    "KTN302", "t5"  ,
            "t9",     "t10",    "t3" ,    "t1" ,    "KTN102", "t12")

p3=cowplot::plot_grid(plotlist = genom_index_plot_list[p_order], ncol=2)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/指标与基因组变异相关性/Allpatient.pdf'),
       p3, width=300, height = 600, units='mm', dpi = 600)

genom_index_cor = as.data.frame(genom_index_cor)
genom_index_cor$V2 = as.numeric(genom_index_cor$V2)
genom_index_cor$V3 = as.numeric(genom_index_cor$V3)
colnames(genom_index_cor) = c('index', 'cor', 'p_value', 'patient')
library(forcats)
genom_index_cor_plot =genom_index_cor %>%
  mutate(patient=fct_reorder(patient, cor)) %>%
  ggplot(aes(x=patient, y=cor, shape=index, color=index, size=-log10(p_value)))+
  geom_point()+
  scale_color_igv()+
  scale_size_continuous( range = c(1,3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )
genom_index_cor_plot
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/指标与基因组变异相关性/ALL.pdf'),
       genom_index_cor_plot, width=100, height = 60, units='mm', dpi = 600)


#### 计算不同指标比例分布 ######
library(ggrepel)
my_pie_plot <- function(x, layers, value,xwd=1){
  p = ggplot()
  p_num = 1
  x$my_id = 1:nrow(x)
  x_pos = xwd*(1:length(layers))
  for(l in layers){
    tmp_p = ggplot()+geom_bar(data=x,
                              mapping = aes_string(x=0, y=value, fill='my_id'),
                              stat='identity', position = 'fill')
    tmp_pied = layer_data(tmp_p)[, c('ymin','ymax')]
    tmp_pied[, 'tmp_group'] = apply(x[, layers[1:p_num], drop=F], 1, function(yy)paste0(yy,collapse = ';;;'))
    tmp_pied = tmp_pied %>%
      group_by(tmp_group) %>%
      summarise(ymin=min(ymin), ymax=max(ymax))
    tmp_pied[, l] = sapply(tmp_pied$tmp_group, function(yy)strsplit(yy, ';;;')[[1]][p_num])
    tmp_pied$label = round(tmp_pied$ymax - tmp_pied$ymin,2)
    tmp_pied$label_x = 0.5+x_pos[p_num]
    tmp_pied$label_y = (tmp_pied$ymax+tmp_pied$ymin)/2

    p = p+
      geom_rect(data=tmp_pied,
                mapping = aes_string(fill=l, ymax='ymax', ymin='ymin', xmax=1+x_pos[p_num], xmin=0+x_pos[p_num]), inherit.aes = F)+
      geom_text(data=tmp_pied,mapping = aes(label=label,  y=label_y, x=label_x), inherit.aes = F)
    p_num = p_num+1
  }

  p=p+
    theme(aspect.ratio=1)+
    coord_polar(theta="y")+
    theme_void()
  return(p)
}

prop_plot_list = list()
prop_data = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'clone','limit_num','seismic_num')]

  pie_1 = meta %>%
    group_by(clone) %>%
    summarise(limit_num_clone = sum(limit_num), seismic_num_clone=sum(seismic_num)) %>%
    filter(grepl('clone', clone)) %>%
    as.data.frame() %>%
    reshape2::melt(id='clone')
  colnames(pie_1) = c('clone', 'type', 'count')

  a=my_pie_plot(pie_1, c('type', 'clone'), 'count')+
    labs(title=new_prefix)+
    scale_fill_manual(values = c(setNames(pal_igv()(length(unique(pie_1$clone))),unique(pie_1$clone)),
                                 'limit_num_clone'='#DFEDD7','seismic_num_clone'='#C9E2F4'))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=4),
          legend.key.size = unit(2, 'mm'),
          legend.title = element_text(size=6))
  pie_1$patient = new_prefix
  pie_1$cell_num = nrow(sctc$map_obj$cell_map_pos)
  prop_data = rbind(prop_data, pie_1)
  prop_plot_list[[new_prefix]] = a
}

p=cowplot::plot_grid(plotlist = prop_plot_list[p_order], ncol=4)
p
chr_num_plot = prop_data %>%
  group_by(patient, type) %>%
  summarise(chr_num = sum(count), cell_num=unique(cell_num)) %>%
  mutate(chr_num = chr_num / cell_num) %>%
  as.data.frame() %>%
  mutate(patient = fct_reorder(patient, -chr_num, sum)) %>%
  ggplot(aes(x=patient, y=chr_num, fill=type))+
  geom_bar(stat='identity', color='black', size=0.2)+
  labs(y='#chr')+
  scale_fill_manual(values = c('limit_num_clone'='#08AEEA', 'seismic_num_clone'='#B721FF'))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),

  )

p_order = prop_data %>%
  group_by(patient, type) %>%
  summarise(chr_num = sum(count), cell_num=unique(cell_num)) %>%
  mutate(chr_num = chr_num / cell_num) %>%
  as.data.frame() %>%
  mutate(patient = fct_reorder(patient, -chr_num, sum))

genom_index_cor_plot =genom_index_cor %>%
  mutate(patient=factor(patient, levels(p_order$patient))) %>%
  ggplot(aes(x=patient, y=cor, shape=index, color=index, size=-log10(p_value)))+
  geom_point()+
  scale_color_igv()+
  scale_size_continuous( range = c(1,3))+
  lims(y=c(0,1))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )
a=cowplot::plot_grid(chr_num_plot,NULL,genom_index_cor_plot, ncol=1, align = 'vh',rel_heights = c(0.4,-0.12, 0.5) )

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/指标与基因组变异相关性_染色体构成.pdf'),
       a, width=120, height = 80, units='mm', dpi = 600)

# 细胞水平
# prop_data2 = c()
# for(prefix in dna_file){
#   print(prefix)
#   new_prefix = gsub('uber.', '', prefix)
#   new_prefix = gsub('.seg.txt', '', new_prefix)
#   new_prefix = gsub('_cleaned_uber', '', new_prefix)
#   new_prefix = gsub('_cleaned', '', new_prefix)
#
#   sctc = TNBC_list_with_score[[prefix]]
#   meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'clone','limit_num','seismic_num')]
#   meta$patient = new_prefix
#   meta$cell_num = nrow(sctc$map_obj$cell_map_pos)
#   prop_data2 = rbind(prop_data2, meta)
#   prop_data2 = rbind(prop_data2)
# }
# library(ggpubr)
# prop_data2[,c('limit_num', 'seismic_num', 'patient', 'cell_num')] %>%
#   reshape2::melt(id=c('patient', 'cell_num'), variable.name='type', value.name = "chr_num") %>%
#   mutate(patient = fct_reorder(patient, -chr_num, median)) %>%
#   ggboxplot(x='patient', y='chr_num', color='black',fill='type', add='jitter', shape=21)+
#   # geom_bar(stat='identity', color='black', size=0.2)+
#   stat_compare_means(aes(group=type), label = 'p.signif')+
#   labs(y='#chr_per_cell')+
#   scale_fill_manual(values = c('limit_num'='#DFEDD7','seismic_num'='#C9E2F4'))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1))


#### PDE与其他指标相关性 #####
pde_cor_plot_list = list()
pde_cor_data = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score','limit_score','seismic_score', 'BFB', 'pde', 'clone')]
  plot_data = reshape2::melt(meta, id=c('pde', 'clone'))
  a=plot_data %>%
    ggplot(aes(x=pde, y=value))+
    geom_point(aes(color=clone), size=1)+
    stat_cor()+
    geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
    scale_color_manual(values = clone_col)+
    facet_wrap(~variable, scales = 'free_y', nrow=1)+
    labs(y='', title = new_prefix)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_blank())

  tmp_func = function(x,y){
    tmp_res = cor.test(x, y)
    return(cbind(tmp_res$estimate, tmp_res$p.value))
  }
  cor_res = plot_data %>%
    group_by(variable) %>%
    summarise(tmp_func(pde, value)) %>% as.data.frame()
  cor_res = as.data.frame(cbind(as.vector(cor_res[,1]), cor_res[,2]))
  cor_res$V4 = new_prefix

  pde_cor_data = rbind(pde_cor_data, cor_res)
  pde_cor_plot_list[[new_prefix]] = a
}

p3=cowplot::plot_grid(plotlist = pde_cor_plot_list[p_order], ncol=2)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/指标与PDE相关性/Allpatient.pdf'),
       p3, width=360, height = 600, units='mm', dpi = 600)
pde_cor_data[is.na(pde_cor_data)] = 0
pde_cor_data$V2 = as.numeric(pde_cor_data$V2)
pde_cor_data$V3 = as.numeric(pde_cor_data$V3)
colnames(pde_cor_data) = c('type', 'cor', 'pvalue', 'patient')
pde_cor_data %>%
  mutate(type=factor(type, c('rearrange_score','limit_score','seismic_score', 'BFB'))) %>%
  ggplot(aes(x=type, y=patient, fill=cor, size=-log(pvalue)))+
  geom_point(shape=21, color='black', stroke=0.2)+
  scale_fill_gradientn(colours = c('white', 'red'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1))


tmp_cor = reshape2::dcast(pde_cor_data[, c(1,2,4)], formula = V1~V4, value.var='V2')
rownames(tmp_cor) = tmp_cor[,1]
tmp_cor = tmp_cor[,-1]
tmp_cor[is.na(tmp_cor)] = 0

#### 病人评分 箱式图 ######
library(forcats)
patient_index_score = c()
for(prefix in dna_file){
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'seismic_score', 'limit_score')]
  meta = meta[!grepl('root|virtual', rownames(meta)), ]
  meta = reshape2::melt(meta)
  colnames(meta) = c('index', 'score')
  meta$patient = new_prefix
  patient_index_score = rbind(patient_index_score, meta)
}

pde_order = patient_index_score %>%
  filter(index=='pde') %>%
  mutate(patient=fct_reorder(patient, score))
pde_order = levels(pde_order$patient)
a = patient_index_score %>%
  mutate(index=factor(index, c('pde', 'rearrange_score', 'BFB', 'limit_score','seismic_score'))) %>%
  mutate(patient=factor(patient, levels=pde_order)) %>%
  ggplot(aes(x=patient, y=score, fill=patient))+
  geom_boxplot(width=0.4,size=0.2, outlier.shape = NA)+
  geom_jitter(width=0.1, shape=21, size=0.1, stroke=0.1)+
  facet_wrap(~index, scales = 'free_y', ncol=1, strip.position="right")+
  scale_fill_igv()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/patient指标箱式图.pdf'),
       a, width=110, height = 120, units='mm', dpi = 600)

a = patient_index_score %>%
  filter(index=='pde') %>%
  mutate(patient=factor(patient, levels=pde_order)) %>%
  ggplot(aes(x=patient, y=score, fill= patient))+
  geom_boxplot(width=0.8,size=0.0, outlier.shape = NA)+
  geom_jitter(width=0.1, shape=21, size=0.01, stroke=0.1)+
  #facet_wrap(~index, scales = 'free_y', ncol=1, strip.position="right")+
  scale_fill_igv()+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        legend.position = 'none',
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/patient指标箱式图CEE.pdf'),
       a, width=70, height = 40, units='mm', dpi = 600)

#### 计算指标累计分布 ######
cum_prop_list = list()
cum_prop = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)


  sctc = TNBC_list_with_score[[prefix]]

  all_cell_cnv=sctc$orig.data$all_node_data
  cell_genome_var_prop = rowSums((all_cell_cnv-2)!=0)/ncol(all_cell_cnv)
  cell_genome_var_prop=cell_genome_var_prop[!grepl('root|virtual', names(cell_genome_var_prop))]
  plot_df = data.frame('genome_alt'=cell_genome_var_prop,
                       'pde' = sctc$map_obj$cell_map_pos[names(cell_genome_var_prop),'pde'])
  plot_df = plot_df %>% arrange(genome_alt)
  plot_df$pde_group = ifelse(plot_df$pde>median(plot_df$pde), 'High', 'low')
  # plot_df$pde_group = ifelse(rownames(plot_df)%in%High_patient, 'High', 'low')

  # high
  CNV_envents = all_cell_cnv
  tmp_func_share<-function(x){
    tmp_data = CNV_envents[x,,drop=F]
    tmp_res = (rowMeans(t(tmp_data))-apply(t(tmp_data), 1, min))==0
    return(sum(tmp_res)/ncol(tmp_data))
  }

  all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  cps_events = (all_rearrange_score_pvalue<pvalue_thr & all_rearrange_score>0)+0
  tmp_func_chromothripisis<-function(x){
    tmp_data = cps_events[x,,drop=F]
    tmp_data = sum(colSums(tmp_data)!=0) /ncol(cps_events)
    return(tmp_data)
  }

  BFB_pos_res = BFB_pos(sctc)
  rownames(BFB_pos_res) = BFB_pos_res[,1]
  BFB_pos_res = BFB_pos_res[,-1]
  tmp_func_BFB<-function(x){
    tmp_data = BFB_pos_res[x,,drop=F] == 'TRUE'
    tmp_data = sum(colSums(tmp_data)!=0) /ncol(tmp_data)
    return(tmp_data)
  }

  plot_df_high = plot_df[plot_df$pde_group=='High',]
  high_rate = c()
  for(i in 1:nrow(plot_df_high)){
    print(i)
    bc = rownames(plot_df_high)[1:i]
    high_rate = rbind(high_rate, c(tmp_func_share(bc), tmp_func_chromothripisis(bc), tmp_func_BFB(bc)))
  }
  high_rate = as.data.frame(high_rate)
  colnames(high_rate) = c('share_CNA', 'chromothripisis', 'BFB')
  high_rate$group='High'
  high_rate$cell_num = (1:nrow(high_rate))/nrow(high_rate)
  # low
  plot_df_low = plot_df[plot_df$pde_group=='low',]
  low_rate = c()
  for(i in 1:nrow(plot_df_low)){
    print(i)
    bc = rownames(plot_df_low)[1:i]
    low_rate = rbind(low_rate, c(tmp_func_share(bc), tmp_func_chromothripisis(bc), tmp_func_BFB(bc)))
  }
  low_rate = as.data.frame(low_rate)
  colnames(low_rate) = c('share_CNA', 'chromothripisis', 'BFB')
  low_rate$group='low'
  low_rate$cell_num = (1:nrow(low_rate))/nrow(low_rate)

  plot_df2 = rbind(high_rate, low_rate)
  plot_df2$patient = new_prefix
  cum_prop = rbind(cum_prop, plot_df2)

  a=reshape2::melt(plot_df2, id=c('group', 'cell_num', 'patient')) %>%
    ggplot(aes(x=cell_num, y=value, color=group))+
    geom_line(size=0.5)+
    #geom_smooth(linewidth=1, se=F)+
    labs(x='Cumulative percentage of cells')+
    facet_wrap(~variable, scales = 'free_y', nrow=1)+
    labs(y='score')+
    scale_color_npg()+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.title.y = element_blank()
          )
  a

  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/细胞累计分布/',prefix,'.pdf'), a,
         width=180, height=60, units='mm', dpi = 600, bg = 'transparent')
  a2=a+labs(title=new_prefix)+theme(plot.title = element_text(hjust = 0.5))
  cum_prop_list[[new_prefix]] = a2
}

p=cowplot::plot_grid(plotlist = cum_prop_list, ncol=2)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/细胞累计分布/Allpatient.pdf'),
       p, width=400, height = 600, units='mm', dpi = 600)

cum_prop = as.data.frame(cum_prop)

saveRDS(cum_prop, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/细胞累计分布prop.rds')
cum_prop = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/细胞累计分布prop.rds')
# 比较曲线下面积
library(pracma)
line_area = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  tmp_high = subset(cum_prop, patient==new_prefix&group=='High')
  high_area = c(trapz(tmp_high$cell_num, tmp_high$share_CNA),
                trapz(tmp_high$cell_num, tmp_high$chromothripisis),
                trapz(tmp_high$cell_num, tmp_high$BFB)
                )
  tmp_low = subset(cum_prop, patient==new_prefix&group=='low')
  low_area = c(trapz(tmp_low$cell_num, tmp_low$share_CNA),
                trapz(tmp_low$cell_num, tmp_low$chromothripisis),
                trapz(tmp_low$cell_num, tmp_low$BFB)
  )
  tmp_area = as.data.frame(rbind(high_area,low_area))
  colnames(tmp_area) = c('share_CNA', 'chromothripisis', 'BFB')
  tmp_area$group = c('High', 'Low')
  tmp_area$patient=prefix
  line_area = rbind(line_area, tmp_area)
}

library(ggbeeswarm)
a=line_area %>%
  reshape2::melt(id=c('group', 'patient')) %>%
  ggplot(aes(x=group, y=value, fill=group))+
  geom_boxplot(size=0.2, outlier.shape = NA)+
  geom_point(shape=21, size=2, stroke=0.2)+
  facet_wrap(~variable, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('High', 'Low')),label = 'p.signif', paired = T)+
  labs(y='AUC')+
  scale_fill_npg()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/细胞累计分布/AUC.pdf'),
       a, width=120, height = 60, units='mm', dpi = 600)

a=line_area %>%
  reshape2::melt(id=c('group', 'patient')) %>%
  filter(variable=='share_CNA') %>%
  ggplot(aes(x=group, y=value, fill=group))+
  geom_boxplot(size=0.1, outlier.shape = NA, width=0.6)+
  geom_point(shape=21, size=1, stroke=0.01)+
  stat_compare_means(comparisons = list(c('High', 'Low')), paired = T, size=2)+
  labs(y='AUC')+
  scale_fill_npg()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/细胞累计分布/AUC2.pdf'),
       a, width=26, height = 32, units='mm', dpi = 600)

#### PDE 与ITH 相关性 #####
all_cell_ITH_mat = c()
all_cell_ITH_plot_list = list()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]

  all_cell_cnv = sctc$orig.data$all_node_data
  cell_pos = sctc$map_obj$cell_map_pos
  cell_ITH_mat = c()
  for(clone in unique(cell_pos$clone)){
    if(grepl('clone', clone)){
      print(clone)
      tmp_clone_cells = rownames(cell_pos[cell_pos$clone==clone, ])

      for(ii in 1:(length(tmp_clone_cells)-1)){
        print(ii)
        i = tmp_clone_cells[ii]
        for(jj in (ii+1):length(tmp_clone_cells)){
          print(jj)
          j = tmp_clone_cells[jj]
          x = unlist(all_cell_cnv[i,,drop=T])
          y = unlist(all_cell_cnv[j,,drop=T])

          up_value = sum((x-mean(x))*(y-mean(y)))
          down_value = sqrt(sum((x-mean(x))^2)) * sqrt(sum((y-mean(y))^2))
          cell_ITH_mat = rbind(cell_ITH_mat, c(clone, i, j, 1-up_value/down_value))
        }
      }
    }
  }

  cell_ITH_mat = as.data.frame(cell_ITH_mat)
  colnames(cell_ITH_mat) = c('clone','c1', 'c2', 'cell_ITH')
  cell_ITH_mat$cell_ITH = as.numeric(cell_ITH_mat$cell_ITH)
  cell_ITH_mat$patient=prefix
  all_cell_ITH_mat = rbind(all_cell_ITH_mat, cell_ITH_mat)

  #saveRDS(cell_ITH_mat, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Huh7_clone_ITH.rds')
  #cell_ITH_mat = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Huh7_clone_ITH.rds')
  clone_ITH = cell_ITH_mat %>%
    group_by(clone, c1) %>%
    summarise(cell_ITH=median(cell_ITH)) %>% as.data.frame()
  clone_ITH$PDE = sctc$map_obj$cell_map_pos[clone_ITH$c1, 'pde']
  a = clone_ITH %>%
    ggplot(aes(x=PDE, y=cell_ITH))+
    geom_point(aes(color=clone), size=1)+
    geom_smooth(method = 'lm', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
    labs(x='PDE', y='cell_ITH')+
    #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
    scale_color_manual(values = clone_col)+
    # scale_fill_igv()+
    labs(title=new_prefix)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(1, 'mm'),
          legend.title = element_text(size=6),
          plot.title = element_text(hjust=0.5)
    )

  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE与ITH相关/',prefix,'.pdf'),
         a, width=70, height = 50, units='mm', dpi = 600)
  all_cell_ITH_plot_list[[new_prefix]] = a
}

saveRDS(all_cell_ITH_mat, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_Cell_ITH.rds')

p=cowplot::plot_grid(plotlist = all_cell_ITH_plot_list, ncol=4)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE与ITH相关/Allpatient.pdf'),
       p, width=300, height = 300, units='mm', dpi = 600)


# all_cell_ITH_mat$patient = all_cell_ITH_mat$c1
# all_cell_ITH_mat$patient = gsub('_.*', '', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('c.*', '', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('T1', 't1', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('T2', 't2', all_cell_ITH_mat$patient)
# saveRDS(all_cell_ITH_mat, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_Cell_ITH.rds')
#### 比较两组指标差异 #####
all_cell_ITH_mat = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_Cell_ITH.rds')

all_cell_ITH_mat = na.omit(all_cell_ITH_mat)
all_cell_ITH_mat = all_cell_ITH_mat %>%
  dplyr::group_by(c1, patient) %>%
  dplyr::summarise(cell_ITH=median(cell_ITH, na.rm = T)) %>% as.data.frame()

rownames(all_cell_ITH_mat) = all_cell_ITH_mat$c1
all_cell_ITH_mat$index='ITH'
all_cell_ITH_mat =all_cell_ITH_mat[,c('index', 'cell_ITH','patient')]
patient_index_score = c()
for(prefix in dna_file){
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'seismic_score', 'limit_score')]
  meta = meta[!grepl('root|virtual', rownames(meta)), ]
  # meta$bc = rownames(meta)
  meta = reshape2::melt(meta)
  colnames(meta) = c('index', 'score')
  meta$patient = new_prefix
  patient_index_score = rbind(patient_index_score, meta)
}

colnames(all_cell_ITH_mat) = colnames(patient_index_score)

patient_index_score = rbind(patient_index_score, all_cell_ITH_mat)
patient_index_score$patient = gsub('uber.', '', patient_index_score$patient)
patient_index_score$patient = gsub('.seg.txt', '', patient_index_score$patient)
patient_index_score$patient = gsub('_cleaned_uber', '', patient_index_score$patient)
patient_index_score$patient = gsub('_cleaned', '', patient_index_score$patient)

library(ggpubr)
### PDE 分组比较 ####
patient_index_score2 = patient_index_score[patient_index_score$index!='pde', ]
patient_index_score_pde = patient_index_score[patient_index_score$index=='pde', ]
patient_index_score_pde = patient_index_score_pde %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(pde_score=median(score)) %>%as.data.frame()

patient_index_score_pde$PDE_group = ifelse(patient_index_score_pde$pde_score>median(patient_index_score_pde$pde_score), 'PDE_high', 'PDE_low')

new_score = merge(patient_index_score_pde, patient_index_score2, by.x = 'patient', by.y='patient')

a= new_score %>%
  dplyr::group_by(patient,PDE_group,index) %>%
  dplyr::summarise(score=median(score)) %>%
  ggplot(aes(x=PDE_group, y=score, fill=PDE_group))+
  geom_boxplot(width=0.6, size=0.2, outlier.shape = NA)+
  geom_point(shape=21, size=2, stroke=0.2)+
  facet_wrap(~index, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('PDE_low', 'PDE_high')),label = 'p.format')+
  labs(y='Score')+
  scale_fill_manual(values=c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(mult = 0.2, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较指标.pdf'),
       a, width=150, height = 60, units='mm', dpi = 600)

a= new_score %>%
  dplyr::group_by(patient,PDE_group,index) %>%
  dplyr::summarise(score=median(score)) %>%
  filter(index=='ITH') %>%
  ggplot(aes(x=PDE_group, y=score, fill=PDE_group))+
  geom_boxplot(size=0.1, outlier.shape = NA, width=0.6)+
  geom_point(shape=21, size=1, stroke=0.01)+
  stat_compare_means(comparisons = list(c('PDE_high','PDE_low')),label = 'p.format',method.args = list(alternative = "greater"), size=2)+#,  method='wilcox.test')+
  labs(y='Score')+
  scale_fill_manual(values=c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(mult = 0.2, add = 0))+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较指标2.pdf'),
       a, width=26, height = 32, units='mm', dpi = 600)

###### PDE高低组病人chr num 比例 ####
prop_data = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'BFB', 'pde', 'clone','limit_num','seismic_num')]
  meta$patient = new_prefix
  meta$cell_num = nrow(sctc$map_obj$cell_map_pos)
  prop_data = rbind(prop_data, meta)
}
rownames(patient_index_score_pde) = patient_index_score_pde$patient
prop_data$PDE_group = patient_index_score_pde[prop_data$patient, 'PDE_group']


prop_data[,c('limit_num', 'seismic_num', 'patient', 'PDE_group', 'cell_num')] %>%
  reshape2::melt(id=c('patient', 'PDE_group', 'cell_num'), variable.name='type', value.name = "chr_num") %>%
  mutate(chr_num = chr_num / cell_num) %>%
  mutate(patient=fct_reorder(patient, -chr_num, sum)) %>%
  dplyr::group_by(patient, type, PDE_group) %>%
  dplyr::summarise(chr_num=sum(chr_num, na.rm = T)) %>%
  ggplot(aes(x=patient, y=chr_num, fill=type))+
  geom_bar(stat='identity', color='black', size=0.2, position = 'fill')+
  facet_wrap(~PDE_group, scales = 'free_x')+
  scale_fill_manual(values = c('limit_num'='#DFEDD7','seismic_num'='#C9E2F4'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1))
#2
prop_data[,c('limit_num', 'seismic_num', 'patient', 'PDE_group', 'cell_num')] %>%
  reshape2::melt(id=c('patient', 'PDE_group', 'cell_num'), variable.name='type', value.name = "chr_num") %>%
  mutate(chr_num = chr_num / cell_num) %>%
  mutate(patient=fct_reorder(patient, -chr_num, sum)) %>%
  dplyr::group_by(patient, type, PDE_group) %>%
  dplyr::summarise(chr_num=sum(chr_num, na.rm = T)) %>%
  ggboxplot(x='PDE_group', y='chr_num', fill = 'type', add = 'jitter', shape=21)+
  scale_fill_manual(values = c('limit_num'='#DFEDD7','seismic_num'='#C9E2F4'))+
  stat_compare_means(aes(group=type), label = 'p.signif')+
  theme(legend.position = 'right')

# 3
a=prop_data[,c('limit_num', 'seismic_num', 'patient', 'PDE_group', 'cell_num')] %>%
  reshape2::melt(id=c('patient', 'PDE_group', 'cell_num'), variable.name='type', value.name = "chr_num") %>%
  mutate(chr_num = chr_num / cell_num) %>%
  mutate(patient=fct_reorder(patient, -chr_num, sum)) %>%
  dplyr::group_by(patient, type, PDE_group) %>%
  dplyr::summarise(chr_num=sum(chr_num, na.rm = T)) %>%
  ggboxplot(x='type', y='chr_num', fill = 'PDE_group', add = 'jitter', shape=21,size=0.1, outlier.shape = NA, width=0.6,
            add.params=list(shape=21, size=1))+
  scale_fill_manual(values = c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  stat_compare_means(aes(group=PDE_group), size=2)+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较chr.pdf'),
       a, width=40, height = 32, units='mm', dpi = 600)


###### PDE高低组病人 两两CNA距离 ####
#High_patient = c("KTN152",  "KTN126", "KTN206","KTN129","KTN132","KTN615", "KTN302","t7" ,   "t9", "t4" )

gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')
# install.packages('valr')
library(valr)
cna_gene = list()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]

  # 位点转换成基因
  cna_loc = as.data.frame(t(sapply(colnames(tmp_cna), function(x)strsplit(x, '_')[[1]])))
  cna_loc$V2 = as.numeric(cna_loc$V2)
  cna_loc$V3 = c(as.numeric(cna_loc$V2[2:nrow(cna_loc)]), 0)
  cna_loc = cna_loc[cna_loc$V2<cna_loc$V3, ]
  cna_loc$cna_name = rownames(cna_loc)
  colnames(cna_loc) = c('chrom', 'start', 'end', 'cna_name')
  res = bed_intersect(cna_loc,gene_info)
  res = res%>%
    dplyr::mutate(percent_overlap = .overlap / (end.y-start.y)) %>%
    dplyr::filter(percent_overlap >= 1)

  # cna gene level
  cna_data_gene = t(tmp_cna)
  cna_data_gene = as.data.frame(cna_data_gene[res$cna_name.x, ])
  cna_data_gene$gene = res$gene_name.y
  cna_data_gene = cna_data_gene %>%
    group_by(gene) %>%
    summarise_all(list(mean)) %>% as.data.frame()

  rownames(cna_data_gene) = cna_data_gene[, 1]
  cna_data_gene = as.data.frame(t(cna_data_gene[,-1]))

  cna_gene[[prefix]] = cna_data_gene
}

saveRDS(cna_gene, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_cna_gene.rds')
cna_gene = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_cna_gene.rds')

same_gene = unique(gene_info$gene_name)
for(prefix in dna_file){
  same_gene = intersect(same_gene, colnames(cna_gene[[prefix]]))
}
patient_level_CNA = c()
for(prefix in dna_file){
  tmp_cna = cna_gene[[prefix]][, same_gene]
  patient_level_CNA = rbind(patient_level_CNA, round(colMeans(tmp_cna)))
}
new_prefix = gsub('uber.', '', dna_file)
new_prefix = gsub('.seg.txt', '', new_prefix)
new_prefix = gsub('_cleaned_uber', '', new_prefix)
new_prefix = gsub('_cleaned', '', new_prefix)
rownames(patient_level_CNA) = new_prefix
patient_index_score_pde = patient_index_score_pde %>% arrange(PDE_group)

cna_dist = 1/as.matrix(dist(patient_level_CNA))
# cna_dist = as.matrix(cor(t(patient_level_CNA)))

# cna_dist = log1p(cna_dist)
diag(cna_dist) = NA
# cna_dist = scale(cna_dist,scale = F, center = F)
width_in = 100 / 25.4
height_in = 80 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE分组异质性热图.pdf', width=width_in, height=height_in)
ht = Heatmap(cna_dist[patient_index_score_pde$patient, patient_index_score_pde$patient],
        cluster_rows = F, cluster_columns = F,
        column_split = patient_index_score_pde$PDE_group,
        row_split = patient_index_score_pde$PDE_group,
        width = ncol(cna_dist)*unit(2.5, "mm"), # 宽度
        height = nrow(cna_dist)*unit(2.5, "mm"), # 高度
        na_col = 'white',
        col = circlize::colorRamp2(c(0, 0.015/2, 0.015), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        column_title_gp = gpar(fontsize = 8),
        row_title_gp = gpar(fontsize = 8),
        #rect_gp = gpar(col='white'),
)
draw(ht)
dev.off()
#### 识别driver基因 01 ######
# 保存数据用与GISTIC
gistic_High = c()
gistic_Low = c()
tmp_name = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  tmp_name = c(tmp_name, new_prefix)
  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  sctc = TNBC_list_with_score[[prefix]]
  colnames(sctc$orig.data$all_node_data) = sapply(colnames(sctc$orig.data$all_node_data), function(x)strsplit(x,'\\.')[[1]][1])
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]
  tmp_cna = data.frame(cna=colMeans(tmp_cna))
  tmp_cna$chr = as.numeric(sapply(rownames(tmp_cna), function(x)strsplit(x, '_')[[1]][1]))
  tmp_cna$seg_start = as.numeric(sapply(rownames(tmp_cna), function(x)strsplit(x, '_')[[1]][2]))
  tmp_cna$seg_end = as.numeric(sapply(rownames(tmp_cna), function(x)strsplit(x, '_')[[1]][3]))
  #tmp_cna$sample
  tmp_cna$markers = 10
  tmp_cna$cna = log2(tmp_cna$cna) - 1
  tmp_cna$sample = new_prefix
  tmp_cna = tmp_cna[, c(6,2,3,4,5,1)]
  for(i in unique(tmp_cna$chr)){
    tmp_chr_data = tmp_cna[tmp_cna$chr==i,]
    tmp_cna[tmp_cna$chr==i, 'seg_end'] = tmp_chr_data$seg_start[c(2:nrow(tmp_chr_data), nrow(tmp_chr_data))]
  }
  if(tmp_group=='PDE_low'){
    gistic_Low = rbind(gistic_Low, tmp_cna)
  }else{
    gistic_High= rbind(gistic_High, tmp_cna)
  }
}
gistic_Low = gistic_Low[gistic_Low$seg_start!=gistic_Low$seg_end, ]
gistic_High = gistic_High[gistic_High$seg_start!=gistic_High$seg_end, ]

write.table(gistic_Low, paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low_seg.txt'),
            quote = F, row.names = F, col.names = F, sep='\t')
write.table(gistic_High, paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high_seg.txt'),
            quote = F, row.names = F, col.names = F, sep='\t')


# 加载gistic结果
high_gist = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/scores.gistic', sep='\t', header = T)
high_gist$group='High'
low_gist = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low/scores.gistic', sep='\t', header = T)
low_gist$group='Low'
amp_gist = rbind(high_gist[high_gist$Type=='Amp', ], low_gist[low_gist$Type=='Amp', ])
del_gist = rbind(high_gist[high_gist$Type=='Del', ], low_gist[low_gist$Type=='Del', ])

chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
chrinfo_chr_absstart$chromNum = gsub('chr0', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('chr', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('X', '23', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('Y', '24', chrinfo_chr_absstart$chromNum)

amp_gist = merge(amp_gist, chrinfo_chr_absstart, by.x='Chromosome', by.y = 'chromNum')
amp_gist$seg_start = amp_gist$Start+amp_gist$absStart
amp_gist$seg_end = amp_gist$End+amp_gist$absStart

del_gist = merge(del_gist, chrinfo_chr_absstart, by.x='Chromosome', by.y = 'chromNum')
del_gist$seg_start = del_gist$Start+del_gist$absStart
del_gist$seg_end = del_gist$End+del_gist$absStart


tmp_fun = function(x,y){
  (min(x)+max(y))/2
}
chr_center_pos = chrinfo %>%
  group_by(chromNum) %>%
  summarise(center = tmp_fun(absStart, absEnd))

chrinfo_other = list(
  chr_center_pos = as.numeric(chr_center_pos$center), # 染色体中心坐标
  chr_name = as.vector(chr_center_pos$chromNum), # 染色体名字
  chr_band = chrinfo[chrinfo$chromStart==0, 'absStart'] # 染色体边界线
)
p = plot_genome(chr_bg_color='white', show_x_axis = T)
#amp_gist_set$seg_end2 = amp_gist_set$seg_start[c(2:nrow(amp_gist_set),nrow(amp_gist_set))]
ids = setdiff(colnames(amp_gist), c('seg_start', 'seg_end'))
amp_gist_seg = reshape2::melt(amp_gist, id=c(ids))
del_gist_seg = reshape2::melt(del_gist, id=c(ids))



pp=p +
  geom_area(data=amp_gist_seg, mapping = aes(x=value, y=G.score, fill=group),color=NA,  alpha=0.6)+
  geom_area(data=del_gist_seg, mapping = aes(x=value, y=-G.score, fill=group),color=NA, alpha=0.6)+
  #geom_text_repel(data=signif_peak, mapping=aes(x=seg_start, y=2,label=Descriptor),nudge_y = 1)+
  scale_fill_manual(values = c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  geom_hline(yintercept = 0, color='black', linewidth=0.5, linetype='dashed')+
  geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
  labs('y'='Deletion<-GISTIC_score->Amplication')+
  #theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
#pp
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较基因组累计2.pdf'),
       pp, width=100, height = 30, units='mm', dpi = 600)

## 识别显著区域
signif_peak_high = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/all_lesions.conf_90.txt', sep='\t', header = T)
signif_peak_high = signif_peak_high[, 1:9]
signif_peak_high$group='High'
signif_peak_low = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low/all_lesions.conf_90.txt', sep='\t', header = T)
signif_peak_low = signif_peak_low[, 1:9]
signif_peak_low$group='Low'
signif_peak = rbind(signif_peak_high, signif_peak_low)
#signif_peak_low = signif_peak_low[signif_peak_low$q.values<0.05, ]
signif_peak = signif_peak[,c(1:3,6,10)]
signif_peak = cbind(signif_peak, t(sapply(signif_peak$Wide.Peak.Limits, function(x)strsplit(x, ':|-|\\(')[[1]][1:3])))
signif_peak$Unique.Name = ifelse(grepl('Amplification', signif_peak$Unique.Name), 'Amp', 'Del')
signif_peak = distinct(signif_peak)
signif_peak = as.data.frame(signif_peak)
signif_peak$chr = gsub('chr','', signif_peak$`1`)
signif_peak = merge(signif_peak, chrinfo_chr_absstart, by.x='chr', by.y = 'chromNum')
signif_peak$seg_start = as.numeric(signif_peak$`2`)+signif_peak$absStart
signif_peak$seg_end = as.numeric(signif_peak$`3`)+signif_peak$absStart
signif_peak$Gistic_diff = 0

for(i in 1:nrow(signif_peak)){
  print(i)
  tmp_chr = as.vector(signif_peak[i, 'chr'])
  tmp_start = as.numeric(signif_peak[i, '2'])
  tmp_end = as.numeric(signif_peak[i, '3'])
  tmp_band = as.vector(signif_peak[i, 'Descriptor'])
  tmp_type = as.vector(signif_peak[i, 'Unique.Name'])
  if(tmp_type=='Del'){
    # Del
    tmp_hg = high_gist[high_gist$Type=='Del', ]
    tmp_lg = low_gist[low_gist$Type=='Del', ]
    tmp_prop1 = subset(tmp_hg, Chromosome==tmp_chr&Start>=tmp_start&Start<=tmp_end)
    tmp_prop2 = subset(tmp_lg, Chromosome==tmp_chr&Start>=tmp_start&Start<=tmp_end)
    if(nrow(tmp_prop2)==0){
      diff_score = mean(tmp_prop1$G.score)
    }else{
      diff_score = mean(tmp_prop1$G.score) - mean(tmp_prop2$G.score)
    }
  }else{
    # Amp
    tmp_hg = high_gist[high_gist$Type=='Amp', ]
    tmp_lg = low_gist[low_gist$Type=='Amp', ]
    tmp_prop1 = subset(tmp_hg, Chromosome==tmp_chr&Start>=tmp_start&Start<=tmp_end)
    tmp_prop2 = subset(tmp_lg, Chromosome==tmp_chr&Start>=tmp_start&Start<=tmp_end)
    if(nrow(tmp_prop2)==0){
      diff_score = mean(tmp_prop1$G.score)
    }else{
      diff_score = mean(tmp_prop1$G.score) - mean(tmp_prop2$G.score)
    }
  }
  signif_peak[i, 'Gistic_diff'] = diff_score
}
signif_peak_label = signif_peak[signif_peak$Gistic_diff>1&signif_peak$q.values<0.05, ]
a=signif_peak %>%
  ggplot(aes(x=Gistic_diff, y=-log10(q.values)))+
  geom_point(aes(color=Unique.Name),size=2)+
  geom_label_repel(data=signif_peak_label, mapping=aes(label=Descriptor), fill=NA)+
  geom_hline(yintercept = -log10(0.05), linetype='dashed')+
  geom_vline(xintercept = 1, linetype='dashed')+
  labs(x='GISTIC_score_difference')+
  scale_color_manual(values = c('Amp'='red', 'Del'='blue'), name='')+
  theme_bw()+
  theme(panel.grid = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/Gistic_Peak筛选.pdf'),
       a, width=100, height = 80, units='mm', dpi = 600)


## 标记显著区域
signif_peak_label2 = signif_peak_label
signif_peak_label2$y=ifelse(signif_peak_label2$Unique.Name=='Amp',2,-2)

pp=p +
  geom_area(data=amp_gist_seg, mapping = aes(x=value, y=G.score, fill=group),color=NA,  alpha=0.6)+
  geom_area(data=del_gist_seg, mapping = aes(x=value, y=-G.score, fill=group),color=NA, alpha=0.6)+
  geom_text_repel(data=signif_peak_label2, mapping=aes(x=seg_start, y=y,label=Descriptor),nudge_y = 0.1)+
  scale_fill_manual(values = c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  geom_hline(yintercept = 0, color='black', linewidth=0.5, linetype='dashed')+
  geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
  labs('y'='Deletion<-GISTIC_score->Amplication')+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))

pp
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较基因组累计2label.pdf'),
       pp, width=260, height = 100, units='mm', dpi = 600)


#### 查看显著区域CNA及基因 ####
signif_peak_label
High_gene_amp = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/amp_genes.conf_90.txt', sep='\t', header = F)
gene_df = data.frame('gene'=High_gene_amp[5:nrow(High_gene_amp),2], 'band'='9p23')
gene_df = rbind(gene_df, data.frame('gene'=High_gene_amp[5:nrow(High_gene_amp),3], 'band'='17q25.3'))
gene_df$Type='Amp'
High_gene_del = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/del_genes.conf_90.txt', sep='\t', header = F)
gene_df = rbind(gene_df, data.frame('gene'=High_gene_del[5:nrow(High_gene_del),2], 'band'='10q23.31', 'Type'='Del'))
gene_df = gene_df[gene_df$gene!='', ]
#
high_mat = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/all_thresholded.by_genes.txt', sep='\t', header = T)
#high_mat = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/all_data_by_genes.txt', sep='\t', header = T)
high_mat = subset(high_mat, Gene.Symbol%in%gene_df$gene&Cytoband%in%gene_df$band)
rownames(high_mat) = high_mat$Gene.Symbol
low_mat = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low/all_thresholded.by_genes.txt', sep='\t', header = T)
#low_mat = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low/all_data_by_genes.txt', sep='\t', header = T)
low_mat = subset(low_mat, Gene.Symbol%in%gene_df$gene&Cytoband%in%gene_df$band)
rownames(low_mat) = low_mat$Gene.Symbol

driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因_Cosmic.csv', sep=',', header = T, check.names = F)
driver_gene = driver_gene$`Gene Symbol`
same_gene = intersect(driver_gene, high_mat$Gene.Symbol)
all_mat = cbind(high_mat[gene_df$gene,4:13], low_mat[gene_df$gene,4:13])

library(ComplexHeatmap)

ht = Heatmap(t(all_mat),
        cluster_rows = F, cluster_columns = F,
        #col=circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        left_annotation = rowAnnotation(df=c(rep('high', 10), rep('low', 10))),
        row_split = c(rep('high', 10), rep('low', 10)),
        top_annotation = HeatmapAnnotation(df=gene_df[,2:3]),
        rect_gp = gpar(color='white')
        )
ht

##
# 只看一个region
same_gene = intersect(driver_gene, high_mat2[grepl('17q25', high_mat2$Cytoband),'Gene.Symbol'])
2500171864
s = 2500171864/1000/1000
label_gene = data.frame('gene'=c('ASPSCR1'),
                        'start'=(2500171864+79934683)/1000/1000-s)

a=p +
  geom_area(data=amp_gist_seg, mapping = aes(x=value/1000/1000-s, y=G.score, color=group),fill=NA, size=1.2)+
  #geom_area(data=del_gist_seg, mapping = aes(x=value, y=-G.score, color=group),fill=NA, alpha=0.6)+
  #geom_text_repel(data=signif_peak_label2, mapping=aes(x=seg_start, y=y,label=Descriptor),nudge_y = 0.1)+
  scale_color_manual(values = c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  #geom_hline(yintercept = 0, color='black', linewidth=0.5, linetype='dashed')+
  #geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
  geom_text_repel(data=label_gene, mapping = aes(x=start, y=2.3, label=gene),nudge_y = 0.2)+
  lims(x=c(2524171864/1000/1000-s, 2581367074/1000/1000-s))+
  labs('y'='GISTIC_score', title = '17q', x='Genomic_pos(Mb)')+
  theme()
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/17qRegion.pdf'),
       a, width=100, height = 80, units='mm', dpi = 600)


# 箱式图
high_mat2 = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/high/all_data_by_genes.txt', sep='\t', header = T)
#high_mat2 = subset(high_mat2, Gene.Symbol%in%gene_df$gene&Cytoband%in%gene_df$band)
rownames(high_mat2) = high_mat2$Gene.Symbol
low_mat2 = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_gistic/low/all_data_by_genes.txt', sep='\t', header = T)
#low_mat2 = subset(low_mat2, Gene.Symbol%in%gene_df$gene&Cytoband%in%gene_df$band)
rownames(low_mat2) = low_mat2$Gene.Symbol

apply(high_mat2['CSNK1D',4:13], 1, function(x)sum(x>0)/10)
apply(low_mat2['CSNK1D',4:13], 1, function(x)sum(x>0)/10)
apply(high_mat2['CSNK1D',4:13], 1, function(x)sd(x))
apply(low_mat2['CSNK1D',4:13], 1, function(x)sd(x))

boxdata = data.frame('avg_amp'=c(unlist(high_mat2['ASPSCR1',4:13]),unlist(low_mat2['ASPSCR1',4:13])),
           'group' = c(rep('High', 10), rep('Low', 10)))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

boxdata2 = data_summary(boxdata, 'avg_amp', 'group')

a=boxdata2 %>%
  ggplot(aes(x=group, y=avg_amp, fill=group))+
  geom_bar(stat = 'identity', width=0.6)+
  scale_fill_manual(values = c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  labs(title='17q25:ASPSCR1')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组ASPSCR1.pdf'),
       a, width=80, height = 80, units='mm', dpi = 600)

### TCGA 分组ITH差异 ####
TCGA_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                      sep='\t', header = T, check.names = F)
colnames(TCGA_cna) = gsub('-01', '',colnames(TCGA_cna))
colnames(TCGA_cna) = gsub('-11', '',colnames(TCGA_cna))

signif_peak_label

del_10q23 = subset(gene_info,chrom=='10'&start>=89597154&end<=90565599)
amp_17q25 = subset(gene_info,chrom=='17'&start>=79501785&end<=81195210)
amp_9p23 = subset(gene_info,chrom=='9'&start>=12652132&end<=14381484)

del_10q23_data = na.omit(TCGA_cna[match(del_10q23$gene_name, TCGA_cna$Hugo_Symbol),])
amp_17q25_data = na.omit(TCGA_cna[match(amp_17q25$gene_name, TCGA_cna$Hugo_Symbol),])
amp_9p23_data = na.omit(TCGA_cna[match(amp_9p23$gene_name, TCGA_cna$Hugo_Symbol),])
# 分组
del_10q23_group = colMeans(del_10q23_data[,3:ncol(del_10q23_data)])
amp_17q25_group = colMeans(amp_17q25_data[,3:ncol(amp_17q25_data)])
amp_9p23_group = colMeans(del_10q23_data[,3:ncol(amp_9p23_data)])

tcga_ITH = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                      sep='\t', header = T)
tcga_ITH$patient = sapply(tcga_ITH$array, function(x)substr(x, 1,12))
tcga_ITH$del_10q23 = del_10q23_group[match(tcga_ITH$patient,names(del_10q23_group))]
tcga_ITH$amp_17q25 = amp_17q25_group[match(tcga_ITH$patient,names(amp_17q25_group))]
tcga_ITH$amp_9p23 = amp_9p23_group[match(tcga_ITH$patient,names(amp_9p23_group))]

tcga_ITH = na.omit(tcga_ITH)

a=tcga_ITH %>%
  filter(del_10q23<=0) %>%
  mutate(del_10q23=round(del_10q23)) %>%
  mutate(del_10q23=ifelse(del_10q23<0, 'Del', 'Neutral')) %>%
  mutate(del_10q23=factor(del_10q23, c('Del', 'Neutral'))) %>%
  ggplot(aes(x=del_10q23, y=Subclonal.genome.fraction, fill=del_10q23))+
  geom_boxplot(size=0.1, outlier.shape = NA, width=0.8)+
  stat_compare_means(comparisons = list(c('Del', 'Neutral')), label='p.format', size=2)+
  scale_fill_manual(values = c('Neutral'='gray', 'Del'='blue'))+
  theme_classic()+  theme(strip.background = element_blank(),
                          legend.position = 'none',
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size=6),
                          plot.title = element_text(hjust=0.5, size=8),
                          axis.line = element_line(linewidth = 0.15),
                          axis.ticks = element_line(linewidth = 0.15),
                          axis.ticks.length=unit(0.025, "cm"),
                          axis.text = element_text(size=4))
b=tcga_ITH %>%
  filter(amp_17q25>=0) %>%
  mutate(amp_17q25=round(amp_17q25)) %>%
  mutate(amp_17q25=ifelse(amp_17q25>0, 'Amp', 'Neutral')) %>%
  mutate(amp_17q25=factor(amp_17q25, c('Amp', 'Neutral'))) %>%
  ggplot(aes(x=amp_17q25, y=Subclonal.genome.fraction, fill=amp_17q25))+
  geom_boxplot(size=0.1, outlier.shape = NA, width=0.8)+
  stat_compare_means(comparisons = list(c('Amp', 'Neutral')), label='p.format', size=2)+
  scale_fill_manual(values = c('Neutral'='gray', 'Amp'='red'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
c=tcga_ITH %>%
  filter(amp_9p23>=0) %>%
  mutate(amp_9p23=round(amp_9p23)) %>%
  mutate(amp_9p23=ifelse(amp_9p23>0, 'Amp', 'Neutral')) %>%
  mutate(amp_9p23=factor(amp_9p23, c('Amp', 'Neutral'))) %>%
  ggplot(aes(x=amp_9p23, y=Subclonal.genome.fraction, fill=amp_9p23))+
  geom_boxplot(size=0.1, outlier.shape = NA, width=0.8)+
  stat_compare_means(comparisons = list(c('Amp', 'Neutral')), label='p.format', size=2)+
  scale_fill_manual(values = c('Neutral'='gray', 'Amp'='red'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))


#pp
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/TCGA分组ITH差异.pdf'),
       a+b+c, width=80, height = 40, units='mm', dpi = 600)

tcga_ITH %>%
  ggplot(aes(x=round(amp_9p23), y=Subclonal.genome.fraction, group=round(amp_9p23)))+
  geom_boxplot()


#### 三个区域gradual seismic覆盖比例 ####
pvalue_thr = 0.001
TNBC_list_with_coverage = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = TNBC_list_with_score[[prefix]]
  all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  # 识别gradual, seismic region #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5 & all_rearrange_score>0)] = NA
  all_region_seismic = list()
  for(i in rownames(all_rearrange_score2)){
    print(i)
    tmp_list = list()
    for(chr in 1:ncol(all_rearrange_score2)){
      if(!is.na(all_rearrange_score2[i, chr])){
        tmp_x = unlist(sctc$orig.data$all_node_data[i, grepl(paste0('^',chr,'_'), colnames(sctc$orig.data$all_node_data))])
        region_res = tryCatch({
          find_event_region(tmp_x, window_size=round(length(tmp_x)/10))
          },error = function(cond){return(NA)})
        if(!is.na(region_res)){
          tmp_list[[as.character(chr)]] = region_res
        }
      }
    }
    all_region_seismic[[i]] = tmp_list
  }

  TNBC_list_with_coverage[[prefix]] = all_region_seismic
}

TNBC_list_with_coverage_limit = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = TNBC_list_with_score[[prefix]]
  all_rearrange_score = read.table(glue('{scTrace_dir}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  # 识别gradual, seismic region #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5 & all_rearrange_score>0)] = NA
  all_region_seismic = list()
  for(i in rownames(all_rearrange_score2)){
    print(i)
    tmp_list = list()
    for(chr in 1:ncol(all_rearrange_score2)){
      if(!is.na(all_rearrange_score2[i, chr])){
        tmp_x = unlist(sctc$orig.data$all_node_data[i, grepl(paste0('^',chr,'_'), colnames(sctc$orig.data$all_node_data))])
        region_res = tryCatch({
          find_event_region(tmp_x, window_size=round(length(tmp_x)/10))
        },error = function(cond){return(NA)})
        if(!is.na(region_res)){
          tmp_list[[as.character(chr)]] = region_res
        }
      }
    }
    all_region_seismic[[i]] = tmp_list
  }

  TNBC_list_with_coverage_limit[[prefix]] = all_region_seismic
}

### 覆盖比例 ###
del_10q23 = subset(gene_info,chrom=='10'&start>=89597154&end<=90565599)
amp_17q25 = subset(gene_info,chrom=='17'&start>=79501785&end<=81195210)
amp_9p23 = subset(gene_info,chrom=='9'&start>=12652132&end<=14381484)


event_coverage_gradual = c()
for(prefix in dna_file){
  print(prefix)
  sctc = TNBC_list_with_score[[prefix]]
  tmp_data = sctc$orig.data$all_node_data
  tmp_all_region = TNBC_list_with_coverage_limit[[prefix]]
  for(r in names(tmp_all_region)){
    for(chr in names(tmp_all_region[[r]])){
      if(chr =='10') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^10_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          #num = nrow(chr_pos[chr_pos$V2>=89597154 & chr_pos$V2<=90565599,,drop=F]) / nrow(chr_pos)

          tmp_cove = chr_pos[chr_pos$V2>=89597154 & chr_pos$V2<=90565599,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(90565599 - 89597154)
          }
          event_coverage_gradual = rbind(event_coverage_gradual, c(prefix,'gradual', '10', r, num))
        }
      }else if(chr =='17') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^17_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          # num = nrow(chr_pos[chr_pos$V2>=79501785 & chr_pos$V2<=81195210,,drop=F])/ nrow(chr_pos)
          tmp_cove = chr_pos[chr_pos$V2>=79501785 & chr_pos$V2<=81195210,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(81195210 - 79501785)
          }

          event_coverage_gradual = rbind(event_coverage_gradual, c(prefix,'gradual', '17', r, num))
        }
      }else if(chr =='9') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^9_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          # num = nrow(chr_pos[chr_pos$V2>=12652132 & chr_pos$V2<=14381484,,drop=F])/ nrow(chr_pos)
          tmp_cove = chr_pos[chr_pos$V2>=12652132 & chr_pos$V2<=14381484,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(14381484 - 12652132)
          }

          event_coverage_gradual = rbind(event_coverage_gradual, c(prefix,'gradual', '9', r, num))
        }
      }
    }
  }
}
event_coverage_gradual = as.data.frame(event_coverage_gradual)
event_coverage_gradual$type='Gradual'
###
event_coverage_seismic = c()
for(prefix in dna_file){
  print(prefix)
  sctc = TNBC_list_with_score[[prefix]]
  tmp_data = sctc$orig.data$all_node_data
  tmp_all_region = TNBC_list_with_coverage[[prefix]]
  for(r in names(tmp_all_region)){
    for(chr in names(tmp_all_region[[r]])){
      if(chr =='10') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^10_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          # num = nrow(chr_pos[chr_pos$V2>=89597154 & chr_pos$V2<=90565599,,drop=F])/ nrow(chr_pos)

          tmp_cove = chr_pos[chr_pos$V2>=89597154 & chr_pos$V2<=90565599,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(90565599 - 89597154)
          }
          event_coverage_seismic = rbind(event_coverage_seismic, c(prefix,'seismic', '10', r, num))
        }
      }else if(chr =='17') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^17_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          #num = nrow(chr_pos[chr_pos$V2>=79501785 & chr_pos$V2<=81195210,,drop=F])/ nrow(chr_pos)
          tmp_cove = chr_pos[chr_pos$V2>=79501785 & chr_pos$V2<=81195210,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(81195210 - 79501785)
          }
          event_coverage_seismic = rbind(event_coverage_seismic, c(prefix,'seismic', '17', r, num))
        }
      }else if(chr =='9') {
        for(pos in tmp_all_region[[r]][[chr]]){
          chr_pos = grep(paste0('^9_'), colnames(tmp_data), value=T)
          chr_pos = chr_pos[pos[1]:pos[2]]
          chr_pos = t(sapply(chr_pos, function(x)strsplit(x, '_')[[1]]))
          chr_pos = as.data.frame(chr_pos)
          chr_pos[,1] = as.numeric(chr_pos[,1])
          chr_pos[,2] = as.numeric(chr_pos[,2] )
          chr_pos[,3] = as.numeric(chr_pos[,3] )
          # num = nrow(chr_pos[chr_pos$V2>=12652132 & chr_pos$V2<=14381484,,drop=F])/ nrow(chr_pos)
          tmp_cove = chr_pos[chr_pos$V2>=12652132 & chr_pos$V2<=14381484,,drop=F]
          if(nrow(tmp_cove)==0){
            num = 0
          }else{
            num = (max(tmp_cove$V2) - min(tmp_cove$V2) )/(14381484 - 12652132)
          }
          event_coverage_seismic = rbind(event_coverage_seismic, c(prefix,'seismic', '9', r, num))
        }
      }
    }
  }
}
event_coverage_seismic = as.data.frame(event_coverage_seismic)
event_coverage_seismic$type='Seismic'


event_coverage = as.data.frame(rbind(event_coverage_gradual, event_coverage_seismic))

colnames(event_coverage) = c('sample','error','chr','bc', 'coverage', 'type')
event_coverage$coverage = as.numeric(event_coverage$coverage)


event_coverage %>%
  #filter(coverage>0) %>%
  group_by(chr, type, sample) %>%
  summarise(coverage = mean(coverage)) %>%
  ggplot(aes(x=chr, y=coverage, fill=type)) +
  geom_boxplot()+
  stat_compare_means(aes(group=type),size=2, label = 'p.format')



df.summary <- event_coverage %>%
  #filter(coverage>0) %>%
  group_by(chr, type) %>%
  summarise(
    sd = sd(coverage, na.rm = TRUE),
    coverage = mean(coverage)
  )

a=ggplot()+
    stat_compare_means(data=event_coverage %>%
                         filter(coverage>0) ,
                       mapping = aes(x=chr, y=coverage, fill=type, group=type),
                       size=2, label = 'p.format', label.y=0.3)+
    geom_pointrange(data=df.summary, mapping = aes(x=chr,y=coverage, ymin = coverage-sd, ymax = coverage+sd, color=type),
                    position = position_dodge(0.5), size=0.5)+
    scale_color_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
    labs(y='Event coverage')+
    theme_classic()+
    theme(strip.background = element_blank(),
          #legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=6),
          plot.title = element_text(hjust=0.5, size=8),
          axis.line = element_line(linewidth = 0.15),
          axis.ticks = element_line(linewidth = 0.15),
          axis.ticks.length=unit(0.025, "cm"),
          axis.text = element_text(size=4))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_比较event覆盖度2.pdf',
         a, width=66, height = 36, units='mm', dpi = 600)

df.summary <- event_coverage %>%
  filter(coverage>0) %>%
  group_by(chr, type) %>%
  summarise(
    up75 = quantile(coverage, na.rm = TRUE,0.75),
    down25 = quantile(coverage, na.rm = TRUE,0.25),
    coverage = mean(coverage)
  )
a=ggplot()+
  stat_compare_means(data=event_coverage %>%
                       filter(coverage>0) ,
                     mapping = aes(x=chr, y=coverage, fill=type, group=type),
                     size=2, label = 'p.format', label.y=0.3)+
  geom_pointrange(data=df.summary, mapping = aes(x=chr,y=coverage, ymin = down25, ymax = up75, color=type),
                  position = position_dodge(0.5), size=0.5)+
  scale_color_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  labs(y='Event coverage')+
  theme_classic()+
  theme(strip.background = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
a
a=event_coverage %>%
  filter(coverage>0) %>%
  ggplot(aes(x=chr, y=coverage, color=type))+
  geom_boxplot(outlier.shape = NA, size=0.4)+
  geom_point(aes(fill=type),size=0.5, shape=21, color='black', stroke=0.01, position = position_jitterdodge(jitter.width=0.2))+
  #stat_compare_means(aes(group=type),size=2, label = 'p.format', label.y=0.3)+
  scale_color_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  scale_fill_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  scale_y_continuous(limits = c(0,0.2))+
  labs(y='Event coverage')+
  theme_classic()+
  theme(strip.background = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_比较event覆盖度3.pdf',
       a, width=66, height = 36, units='mm', dpi = 600)
aa=patient_index_score_pde

event_coverage$sample = gsub('_cleaned_uber.seg.txt', '', event_coverage$sample)
event_coverage$sample = gsub('uber.', '', event_coverage$sample)
event_coverage$sample = gsub('.seg.txt', '', event_coverage$sample)
event_coverage$group = patient_index_score_pde[event_coverage$sample, 'PDE_group']

aa=event_coverage %>%
  group_by(sample, chr, type, group, bc) %>%
  summarise(coverage=sum(coverage)) %>%
  # filter(coverage>0) %>%
  group_by(sample, chr, type, group) %>%
  summarise(prop=sum(coverage>0), all=n()) %>%
  mutate(prop=prop/all) %>%
  as.data.frame()
bb=c()
for(i in unique(aa$sample)){
  for(j in unique(aa$chr)){
    for(k in unique(aa$type)){
      for(m in unique(aa$group)){
        tmp = subset(aa, sample==i&chr==j&type==k&group==m)
        if(nrow(tmp)==0){
          bb = rbind(bb, c(i,j,k,m, 0))
        }else{
          bb = rbind(bb, tmp)
        }
      }
    }
  }

}
bb = as.data.frame(bb)
colnames(bb) = colnames(aa)
bb$prop = as.numeric(bb$prop)
a4 = bb %>%
  ggboxplot(x='chr', y='prop', color='type', add='jitter', add.params = list(size=0.2), size=0.5)+
  stat_compare_means(aes(group=type), label='p.format', size=2)+
  #ggplot(aes(x=chr, y=coverage, color=type))+
  #geom_boxplot(outlier.shape = NA, size=0.4)+
  #geom_point(aes(fill=type),size=4, shape=21, color='black', stroke=0.01, position = position_jitterdodge(jitter.width=0.2))+
  facet_wrap(~chr, scales = 'free')+
  scale_color_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  scale_fill_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  #scale_y_continuous(limits = c(0,0.2))+
  labs(y='Cell proportion')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
a4
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_比较event覆盖度5.pdf',
       a4, width=86, height = 36, units='mm', dpi = 600)


aa=event_coverage %>%
  group_by(sample, chr, type, group, bc) %>%
  summarise(coverage=sum(coverage)) %>%
  # filter(coverage>0) %>%
  group_by(sample, chr, type, group) %>%
  summarise(coverage=mean(coverage)) %>%
  as.data.frame()
bb=c()
for(i in unique(aa$sample)){
  for(j in unique(aa$chr)){
    for(k in unique(aa$type)){
      for(m in unique(aa$group)){
        tmp = subset(aa, sample==i&chr==j&type==k&group==m)
        if(nrow(tmp)==0){
          bb = rbind(bb, c(i,j,k,m, 0))
        }else{
          bb = rbind(bb, tmp)
        }
      }
    }
  }

}
bb = as.data.frame(bb)
colnames(bb) = colnames(aa)
bb$coverage = as.numeric(bb$coverage)
a4 = bb %>%
  ggboxplot(x='chr', y='coverage', color='type', add='jitter', add.params = list(size=0.2), size=0.5)+
  stat_compare_means(aes(group=type), label='p.format', size=2)+
  #ggplot(aes(x=chr, y=coverage, color=type))+
  #geom_boxplot(outlier.shape = NA, size=0.4)+
  #geom_point(aes(fill=type),size=4, shape=21, color='black', stroke=0.01, position = position_jitterdodge(jitter.width=0.2))+
  facet_wrap(~chr, scales = 'free')+
  scale_color_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  scale_fill_manual(values = c('Gradual'='#DFEDD7','Seismic'='#C9E2F4'))+
  #scale_y_continuous(limits = c(0,0.2))+
  labs(y='Event coverage')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=4))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/clone_比较event覆盖度4.pdf',
       a4, width=86, height = 36, units='mm', dpi = 600)
# library(circlize) # >= 0.4.10
# col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "gray", "red"))
# circos.clear()
# circos.heatmap(t(all_mat), split = c(rep('high', 10), rep('low', 10)), col = col_fun1,cell.border='black')
### TNBC 生存分析 #####
library(survival)
library(survminer)

metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

TCGA_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                      sep='\t', header = T, check.names = F)
colnames(TCGA_cna) = gsub('-01', '',colnames(TCGA_cna))
colnames(TCGA_cna) = gsub('-11', '',colnames(TCGA_cna))


all_pvalue = c()
#high_mat2[grepl('10q',high_mat2$Cytoband),'Gene.Symbol']
for(i in high_mat2[grepl('17q25', high_mat2$Cytoband),'Gene.Symbol']){
  print(i)
  select_gene = i#'PDCD6'
  metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                               sep='\t', header = T, check.names = F)
  rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
  same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

  # setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
  if(select_gene %in% metabric_cna$Hugo_Symbol){
    metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
    #metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']>0, 'Gain', 'other')
    metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
    x1 = surv_pvalue(km)
  }else{
    x1 = list('pval'=1)
  }

  #
  metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                               sep='\t', header = T, check.names = F)

  rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
  same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(TCGA_cna))

  metabric_clinic = metabric_clinic[same_sample, ]
  if(select_gene %in% TCGA_cna$Hugo_Symbol){
    metabric_clinic[, 'group'] = unlist(TCGA_cna[TCGA_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
    #metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Gain', 'other')

    metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
    x2=surv_pvalue(km)


    metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                                 sep='\t', header = T, check.names = F)

    rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
    same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(TCGA_cna))

    metabric_clinic = metabric_clinic[same_sample, ]
    metabric_clinic[, 'group'] = unlist(TCGA_cna[TCGA_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
    #metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Gain', 'other')

    metabric_clinic$PFS_STATUS = as.numeric(sapply(metabric_clinic$PFS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(PFS_MONTHS, PFS_STATUS)~group, data = metabric_clinic)
    x3=surv_pvalue(km)

  }else{
    x2 = list('pval'=1)
    x3 = list('pval'=1)
  }


  all_pvalue = rbind(all_pvalue, c(i, x1$pval, x2$pval, x3$pval))
}

all_pvalue = as.data.frame(all_pvalue)
all_pvalue[,2] = as.numeric(all_pvalue[,2])
all_pvalue[,3] = as.numeric(all_pvalue[,3])
all_pvalue[,4] = as.numeric(all_pvalue[,4])

all_pvalue %>% arrange(-V2, -V3)

###
select_gene = 'ASPSCR1'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)
rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

# setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']>0, 'Amp', 'Neutral')
metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
#metabric_clinic = metabric_clinic[metabric_clinic$THREEGENE=='ER+/HER2- Low Prolif', ]
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
#leg = as.data.frame(table(metabric_clinic[same_sample, 'group']))

ggsurvplot(km,
           pval=T,
           #pval.method=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           xlab="M",
           risk.table = T,
           ylab="OS(METABRIC)"
)
width_in = 120 / 25.4
height_in = 100 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/Metabric生存.pdf', width=width_in, height=height_in)
ggsurvplot(km,
           #pval=T,
           pval = "p-value=0.0089\nAmp=372\nNeutral=1461",
           #legend = "right", #将图例移动到下方
           legend.title = "CNA",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           surv.median.line = "hv",
           xlab="M",
           risk.table = F,
           ylab="OS(METABRIC)",
           palette=c('red', 'gray')
)
dev.off()





# TCGA
# phenotype_file <- read.table("/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
#                              header = T, sep = '\t', quote = "")
# phenotype_colnames <- as.data.frame(colnames(phenotype_file))
# #读取文件并去读列名（penotype的类型）
# table(phenotype_file$breast_carcinoma_estrogen_receptor_status)#查看雌激素受体(ER)的表达状态。
# table(phenotype_file$breast_carcinoma_progesterone_receptor_status)#查看孕激素受体（PR）状态
# table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status)#查看原癌基因Her-2状态
#
# phenotype_colnames <- colnames(phenotype_file)[ grep("receptor_status",colnames(phenotype_file))]
# eph <- phenotype_file[,phenotype_colnames[1:3]]
# tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
# tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]
#
# tnbc_sample = gsub('-01', '',tnbc_sample)
# tnbc_sample = gsub('-11', '',tnbc_sample)


#select_gene = 'TEX19'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                          sep='\t', header = T, check.names = F)
colnames(metabric_cna) = gsub('-01', '',colnames(metabric_cna))
colnames(metabric_cna) = gsub('-11', '',colnames(metabric_cna))

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)

rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

metabric_clinic = metabric_clinic[same_sample, ]
metabric_clinic[, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(1,0))
metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Amp', 'Neutral')

metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
# metabric_clinic = metabric_clinic[metabric_clinic$PATIENT_ID %in%tnbc_sample,]
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
#km <- survdiff(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           #pval=T,
           pval = "p-value=0.08\nAmp=357\nNeutral=504",
           #legend = "right", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           surv.median.line = "hv",
           xlab="M",
           risk.table = F,
           ylab="OS(TCGA)",
           #palette=c('red', 'gray')
)
width_in = 120 / 25.4
height_in = 100 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/TCGA生存.pdf', width=width_in, height=height_in)
ggsurvplot(km,
           #pval=T,
           pval = "p-value=0.26\nAmp=341\nNeutral=506",
           #legend = "right", #将图例移动到下方
           legend.title = "CNA",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           surv.median.line = "hv",
           xlab="M",
           risk.table = F,
           ylab="OS(TCGA)",
           palette=c('red', 'gray')
)
dev.off()






# PFS
metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)

rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

metabric_clinic = metabric_clinic[same_sample, ]
metabric_clinic[, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
#metabric_clinic = metabric_clinic %>% filter(group%in%c(-1,0))
#metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Gain', 'Neutral')

metabric_clinic$PFS_STATUS = as.numeric(sapply(metabric_clinic$PFS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(PFS_MONTHS, PFS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="PFS(TCGA)"
)


















# 绘制位点累计频率
rownames(patient_index_score_pde) = patient_index_score_pde$patient
driver_gene_rate2_gain = c()
driver_gene_rate2_loss = c()
same_loc = colnames(TNBC_list_with_score$uber.t12.seg.txt$orig.data$all_node_data)
tmp_name = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  tmp_name = c(tmp_name, new_prefix)
  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  # tmp_cna = cna_gene[[prefix]][, driver_gene$driver_gene]
  sctc = TNBC_list_with_score[[prefix]]
  colnames(sctc$orig.data$all_node_data) = sapply(colnames(sctc$orig.data$all_node_data), function(x)strsplit(x,'\\.')[[1]][1])
  tmp_cna = sctc$orig.data$all_node_data[,same_loc]
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]
  driver_gene_rate2_gain = rbind(driver_gene_rate2_gain, colMeans(tmp_cna>2))
  driver_gene_rate2_loss = rbind(driver_gene_rate2_loss, colMeans(tmp_cna<2))
}
rownames(driver_gene_rate2_gain) = tmp_name
rownames(driver_gene_rate2_loss) = tmp_name

high_gain_rate = colSums(driver_gene_rate2_gain[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_high', 'patient'], ])
low_gain_rate = colSums(driver_gene_rate2_gain[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_low', 'patient'], ])

plot_gain = data.frame(high=high_gain_rate, low=low_gain_rate)
plot_gain$index = colnames(driver_gene_rate2_gain)
plot_gain$chr = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][1]))
plot_gain$seg_start = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][2]))
plot_gain$seg_end = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][3]))

plot_gain$seg_end = plot_gain$seg_start
plot_gain$seg_start = c(0, plot_gain$seg_start[1:(nrow(plot_gain)-1)])

plot_gain$x = 1:nrow(plot_gain)
plot_gain %>%
  ggplot(aes(x=x))+
  geom_area(aes(y=high), fill='red', alpha=0.6)+
  geom_area(aes(y=low), fill='blue', alpha=0.6)

#
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
chrinfo_chr_absstart$chromNum = gsub('chr0', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('chr', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('X', '23', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('Y', '24', chrinfo_chr_absstart$chromNum)

plot_gain = merge(plot_gain, chrinfo_chr_absstart, by.x='chr', by.y = 'chromNum')
plot_gain$seg_start = plot_gain$seg_start+plot_gain$absStart
plot_gain$seg_end = plot_gain$seg_end+plot_gain$absStart



p = plot_genome(chr_bg_color='white', show_x_axis = T)

a1=p +
  #geom_hline(yintercept = 0, color='black')+
  #geom_area(data=plot_gain, aes(y=high), fill='red', alpha=0.6)+
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=high),fill='#E64B35FF', alpha=0.5)+
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=low),fill='#4DBBD5FF', alpha=0.5)+
  labs('y'='Amplication')+
  # geom_line(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC)))+
  #geom_text_repel(data=know_gene_data,
  #                mapping = aes(x=seg_start, y=7, label=gene),
  #                nudge_y = 1,
  #                #arrow = arrow(length = unit(0.02, "npc")),box.padding =1,
  #                max.overlaps = Inf,
#
  #)+
  #scale_fill_manual(values = c('Up'='red', 'Down'='blue'))+
  #scale_size(range = c(0.1,4))+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))

a1
#
high_loss_rate = colSums(driver_gene_rate2_loss[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_high', 'patient'], ])
low_loss_rate = colSums(driver_gene_rate2_loss[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_low', 'patient'], ])
plot_loss = data.frame(high=high_loss_rate, low=low_loss_rate)
plot_loss$index = colnames(driver_gene_rate2_loss)
plot_loss$chr = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][1]))
plot_loss$seg_start = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][2]))
plot_loss$seg_end = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][3]))
plot_loss$x = 1:nrow(plot_loss)
plot_loss$seg_end = plot_loss$seg_start
plot_loss$seg_start = c(0, plot_loss$seg_start[1:(nrow(plot_loss)-1)])

plot_loss = merge(plot_loss, chrinfo_chr_absstart, by.x='chr', by.y = 'chromNum')
plot_loss$seg_start = plot_loss$seg_start+plot_loss$absStart
plot_loss$seg_end = plot_loss$seg_end+plot_loss$absStart


a2=p +
  #geom_hline(yintercept = 0, color='black')+
  #geom_area(data=plot_gain, aes(y=high), fill='red', alpha=0.6)+
  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=high),fill='#E64B35FF', alpha=0.5)+
  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=low),fill='#4DBBD5FF', alpha=0.5)+
  # geom_line(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC)))+
  #geom_text_repel(data=know_gene_data,
  #                mapping = aes(x=seg_start, y=7, label=gene),
  #                nudge_y = 1,
  #                #arrow = arrow(length = unit(0.02, "npc")),box.padding =1,
  #                max.overlaps = Inf,
#
  #)+
  #scale_fill_manual(values = c('Up'='red', 'Down'='blue'))+
  #scale_size(range = c(0.1,4))+
  labs('y'='Deletion')+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))
pp = cowplot::plot_grid(a1,a2, ncol=1, align = 'vh')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较基因组累计.pdf'),
       pp, width=260, height = 100, units='mm', dpi = 600)

#
tmp_fun = function(x,y){
  (min(x)+max(y))/2
}
chr_center_pos = chrinfo %>%
  group_by(chromNum) %>%
  summarise(center = tmp_fun(absStart, absEnd))

chrinfo_other = list(
  chr_center_pos = as.numeric(chr_center_pos$center), # 染色体中心坐标
  chr_name = as.vector(chr_center_pos$chromNum), # 染色体名字
  chr_band = chrinfo[chrinfo$chromStart==0, 'absStart'] # 染色体边界线
)
p = plot_genome(chr_bg_color='white', show_x_axis = T, hlinewidth = 1)

pp=p +
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=high),fill='#E64B35FF', alpha=0.5)+
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=low),fill='#4DBBD5FF', alpha=0.5)+

  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=-high),fill='#E64B35FF', alpha=0.5)+
  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=-low),fill='#4DBBD5FF', alpha=0.5)+

  geom_hline(yintercept = 0, color='black', linewidth=0.5, linetype='dashed')+
  geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
  labs('y'='Deletion<------>Amplication')+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较基因组累计2.pdf'),
       pp, width=260, height = 100, units='mm', dpi = 600)

## 识别差异位点
chrinfo_band = chrinfo
#chrinfo_band$bandname = sapply(as.vector(chrinfo_band$bandname), function(x)strsplit(x, '\\.')[[1]][1])
chrinfo_band$bandname = sapply(as.vector(chrinfo_band$bandname), function(x)substr(x,1,1))

chrinfo_band = chrinfo_band %>%
  group_by(chromNum, bandname) %>%
  summarise(chromStart=min(chromStart), chromEnd=max(chromEnd),
            absStart = min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
tmp_prop = plot_gain
gain_FC = c()
for(i in 1:nrow(chrinfo_band)){
  print(i)
  tmp_chr = as.vector(chrinfo_band[i, 'chromNum'])
  tmp_chr = gsub('chrX', '23', tmp_chr)
  tmp_chr = gsub('chrY', '24', tmp_chr)
  tmp_chr = gsub('chr0', '', tmp_chr)
  tmp_chr = gsub('chr', '', tmp_chr)


  tmp_start = chrinfo_band[i, 'absStart']
  tmp_end = chrinfo_band[i, 'absEnd']
  tmp_band = as.vector(chrinfo_band[i, 'bandname'])
  tmp_prop2 = subset(tmp_prop, chr==tmp_chr&seg_start>=tmp_start&seg_start<=tmp_end)
  if(nrow(tmp_prop2)>0){
    tmp_fc = mean(tmp_prop2$high) / mean(tmp_prop2$low)
    if(all(tmp_prop2$high== mean(tmp_prop2$high)) | all(tmp_prop2$low== mean(tmp_prop2$low))){
      tmp_p = Inf
    }else{
      tmp_p = t.test(tmp_prop2$high, tmp_prop2$low)$p.value
    }
    gain_FC = rbind(gain_FC, c(tmp_chr, tmp_band, tmp_start, tmp_end, tmp_fc, tmp_p, mean(tmp_prop2$high), mean(tmp_prop2$low)))
  }
}

gain_FC = as.data.frame(gain_FC)
colnames(gain_FC) = c('chr', 'band', 'start', 'end', 'FC', 'pvalue', 'mean1', 'mean2')
gain_FC$FC = as.numeric(gain_FC$FC)
gain_FC$pvalue = as.numeric(gain_FC$pvalue)
gain_FC$mean1 = as.numeric(gain_FC$mean1)
gain_FC$mean2 = as.numeric(gain_FC$mean2)
gain_FC$group='Gain'
#
tmp_prop = plot_loss
loss_FC = c()
for(i in 1:nrow(chrinfo_band)){
  print(i)
  tmp_chr = as.vector(chrinfo_band[i, 'chromNum'])
  tmp_chr = gsub('chrX', '23', tmp_chr)
  tmp_chr = gsub('chrY', '24', tmp_chr)
  tmp_chr = gsub('chr0', '', tmp_chr)
  tmp_chr = gsub('chr', '', tmp_chr)


  tmp_start = chrinfo_band[i, 'absStart']
  tmp_end = chrinfo_band[i, 'absEnd']
  tmp_band = as.vector(chrinfo_band[i, 'bandname'])
  tmp_prop2 = subset(tmp_prop, chr==tmp_chr&seg_start>=tmp_start&seg_start<=tmp_end)
  if(nrow(tmp_prop2)>0){
    tmp_fc = mean(tmp_prop2$high) / mean(tmp_prop2$low)
    if(all(tmp_prop2$high== mean(tmp_prop2$high)) | all(tmp_prop2$low== mean(tmp_prop2$low))){
      tmp_p = Inf
    }else{
      tmp_p = t.test(tmp_prop2$high, tmp_prop2$low)$p.value
    }
    loss_FC = rbind(loss_FC, c(tmp_chr, tmp_band, tmp_start, tmp_end, tmp_fc, tmp_p,mean(tmp_prop2$high),mean(tmp_prop2$low)))
  }
}

loss_FC = as.data.frame(loss_FC)
colnames(loss_FC) = c('chr', 'band', 'start', 'end', 'FC', 'pvalue', 'mean1', 'mean2')
loss_FC$FC = as.numeric(loss_FC$FC)
loss_FC$pvalue = as.numeric(loss_FC$pvalue)
loss_FC$mean1 = as.numeric(loss_FC$mean1)
loss_FC$mean2 = as.numeric(loss_FC$mean2)
loss_FC$group='Loss'

all_FC = rbind(gain_FC, loss_FC)#%>%mutate(pvalue=p.adjust(pvalue))
#all_FC = all_FC %>% filter(FC>2, pvalue<0.001)


label_all_FC = all_FC %>% filter(mean1>2.243446, log2(FC)>1) %>% filter(chr!='23')%>% filter(chr!='24')
all_FC %>%
  filter(chr!='23')%>% filter(chr!='24')%>%
  mutate(pvalue=p.adjust(pvalue))%>%
  #filter(group=='Gain') %>%
  ggplot(aes(x=mean1, y=mean2, color=group))+
  geom_point(aes(size=abs(mean1-mean2)))+
  geom_hline(yintercept = 2.243446, color='black', linetype='dashed')+
  geom_vline(xintercept = 1, color='black', linetype='dashed')+
  geom_text_repel(data=label_all_FC, mapping = aes(label=paste0(chr,'_', band)), nudge_x=2)+
  geom_abline(slope = 1,intercept = -2)+
  #scale_size_continuous(limits=c(0,10))+
  scale_fill_igv()+
  theme_classic()
#
all_FC2 = merge(gain_FC, loss_FC, by = c('chr', 'band'))

label_all_FC = all_FC2 %>% filter(FC.x>1, FC.y<1) %>% filter(chr!='23')%>% filter(chr!='24')
label_all_FC = all_FC2 %>% filter(FC.x<1, FC.y>1) %>% filter(chr!='23')%>% filter(chr!='24')

all_FC2 %>%
  filter(chr!='23')%>% filter(chr!='24')%>%
  ggplot(aes(x=log2(FC.y), y=log2(FC.x)))+
  geom_point(aes(color=mean1.y-mean1.x, size=mean1.y-mean2.y))+
  geom_text_repel(data=label_all_FC, mapping=aes(label=paste0(chr, band)))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_gradientn(colors = c('blue', 'white','red'))+
  theme_classic()

all_FC2 %>%
  filter(chr!='23')%>% filter(chr!='24')%>%
  mutate(band=paste0(chr,band)) %>%
  mutate(band=fct_reorder(band, mean1.x-mean2.x-mean1.y+mean2.y)) %>%
  ggplot(aes(x=band))+
  geom_point(aes(y=mean1.x-mean2.x, size=(mean1.x+mean2.x)/2), color='red')+
  geom_point(aes(y=mean1.y-mean2.y, size=(mean1.y+mean2.y)/2), color='blue')+
  geom_segment(aes(x=band, xend=band,y=mean1.x-mean2.x, yend=mean1.y-mean2.y))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))


all_FC2 %>%
  filter(chr!='23')%>% filter(chr!='24')%>%
  mutate(band=paste0(chr,band)) %>%
  ggplot(aes(x=mean1.x-mean2.x, y=mean1.y-mean2.y, color=chr))+
  geom_point()+
  scale_color_igv()


# all_FC2 %>%
#   filter(mean1.x>4) %>%
#   ggplot(aes(x=FC.x, y=(FC.y)))+
#   geom_point()+
#   geom_hline(yintercept = c(0, log2(1.5), log2(0.5)))+
#   geom_vline(xintercept = c(0, log2(1.5)))+
#   geom_text_repel(aes(label=paste0(chr,'_', band)))

# all_FC2$pvalue.x = p.adjust(all_FC2$pvalue.x,method='BH')
# all_FC2$pvalue.y = p.adjust(all_FC2$pvalue.y,method='BH')

#label_all_FC = all_FC2 %>% filter(pvalue.x<0.05, pvalue.y>=0.05)
#label_all_FC = label_all_FC %>% filter(FC.x>2, FC.y<2)
#label_all_FC = label_all_FC %>% filter((mean1.x-mean2.x)>1,(mean1.y-mean2.y)<1,(mean1.y-mean2.y)>-1)
#label_all_FC = label_all_FC %>% filter(mean1.x>2)

# label_all_FC = all_FC2 %>% filter(FC.x>2, pvalue.x<0.001)
# label_all_FC = label_all_FC %>% filter(mean1.x>2)
# label_all_FC = label_all_FC %>% filter(FC.y<2,FC.y>0.5)
# label_all_FC = label_all_FC %>% filter(mean1.y<=2)
# label_group =label_all_FC[,c(1,2)]
# label_group$label_group = 'Signif'
# all_FC2 = left_join(all_FC2, label_group, by = c('chr', 'band'))

#all_FC2 %>%
#  #mutate(pvalue.x=p.adjust(pvalue.x))%>%
#  ggplot(aes(x=log2(FC.x), y=log2(FC.y)))+
#  geom_point(aes(size=mean1.x, fill=label_group), shape=21)+
#  geom_hline(yintercept = c(-1, 1), color='gray', linetype='dashed')+
#  geom_vline(xintercept = c(1), color='gray', linetype='dashed')+
#  geom_text_repel(data=label_all_FC, mapping = aes(label=paste0(chr,'_', band)),
#                  color='red', max.overlaps = Inf,
#                  arrow=arrow(angle = 30, length = unit(0.1, "inches")),
#                  nudge_x = 0, nudge_y = -3)+
#  scale_fill_igv()+
#  theme_classic()
#
#all_FC2 %>%
#  #mutate(pvalue.x=p.adjust(pvalue.x))%>%
#  ggplot(aes(x=mean1.x-mean2.x, y=mean1.y-mean2.y))+
#  geom_point(aes(size=mean1.x, fill=chr), shape=21)+
#  geom_hline(yintercept = c(-1, 1), color='gray', linetype='dashed')+
#  geom_vline(xintercept = c(1), color='gray', linetype='dashed')+
#  geom_text_repel(data=label_all_FC, mapping = aes(label=paste0(chr,'_', band)),
#                  color='red', max.overlaps = Inf,
#                  arrow=arrow(angle = 30, length = unit(0.1, "inches")),
#                  nudge_x = 0, nudge_y = -3)+
#  #scale_fill_gradientn(colours = c('white', 'red'))+
#  scale_size_continuous(limits=c(0,10))+
#  scale_fill_igv()+
#  theme_classic()

# 单独查看每个band区域
#chrinfo_band = chrinfo
#chrinfo_band$bandname = sapply(as.vector(chrinfo_band$bandname), function(x)strsplit(x, '\\.')[[1]][1])
# chrinfo_band = chrinfo_band %>%
#   group_by(chromNum, bandname) %>%
#   summarise(chromStart=min(chromStart), chromEnd=max(chromEnd),
#             absStart = min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')

driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因_Cosmic.csv', sep=',', header = T, check.names = F)
driver_gene = driver_gene$`Gene Symbol`

# label_all_FC =all_FC2 %>% filter((mean1.x-mean2.x)>4)

label_gene_all = c()
pg_list = list()
for(i in 1:nrow(label_all_FC)){
  #i=1
  tmp_chr = label_all_FC[i, 'chr']
  tmp_band=label_all_FC[i, 'band']
  tmp_start = label_all_FC[i, 'start.x']
  tmp_end = label_all_FC[i, 'end.x']
  tmp_prop = plot_loss
  # if(label_all_FC[i, 'group']=='Gain'){
  #   tmp_prop = plot_gain
  #   laby = 'Amplication'
  # }else{
  #   tmp_prop = plot_loss
  laby = 'Deletion'
  # }
  tmp_prop$seg_start = (tmp_prop$seg_start-tmp_prop$absStart)
  tmp_prop$seg_end = (tmp_prop$seg_end-tmp_prop$absStart)
  tmp_prop=tmp_prop[,c(1,5,6, 4,2,3,7,8,9)]
  colnames(tmp_prop)[1:3] = c('chrom', 'start', 'end')
  tmp_prop$chrom = as.character(tmp_prop$chrom)
  tmp_prop = bed_intersect(tmp_prop,gene_info) %>% as.data.frame()

  tmp_prop$start = tmp_prop$start.x / (1*1000*1000)
  ##### 绘制画布 ######
  tmp_prop_chr = tmp_prop[tmp_prop$chrom==tmp_chr,]
  tmp_chrinfo_band = chrinfo_band
  tmp_chrinfo_band$chromNum = gsub('chr0','', tmp_chrinfo_band$chromNum)
  tmp_chrinfo_band$chromNum = gsub('chr','', tmp_chrinfo_band$chromNum)
  tmp_chrinfo_band$chromNum = gsub('X','23', tmp_chrinfo_band$chromNum)
  tmp_chrinfo_band$chromNum = gsub('Y','24', tmp_chrinfo_band$chromNum)
  tmp_chrinfo_band = tmp_chrinfo_band %>%filter(chromNum==tmp_chr&bandname==tmp_band)
  tmp_chrinfo_band$absStart = as.numeric(tmp_chrinfo_band$absStart)
  tmp_chrinfo_band$absEnd = as.numeric(tmp_chrinfo_band$absEnd)
#
  # pg = plot_genome_chr(chr_num = tmp_chr,bg_color = 'white',
  #                      start = tmp_chrinfo_band$absStart,
  #                      end=tmp_chrinfo_band$absEnd)
  left_x = as.numeric(tmp_chrinfo_band$chromStart)/ (1*1000*1000)
  right_x = as.numeric(tmp_chrinfo_band$chromEnd)/ (1*1000*1000)
  label_gene =tmp_prop_chr
  label_gene = label_gene[match(driver_gene,label_gene$gene_name.y), ] %>% na.omit()
  label_gene$group = label_all_FC[i, 'group']
  label_gene_all = rbind(label_gene_all, label_gene)
  pg=ggplot()+
    #geom_vline(xintercept = tmp_chrinfo_band$absStart, color='gray', linewidth=0.1, linetype='dashed')+
    #geom_area(mapping = aes(x=c(as.numeric(tmp_start)/ (1*1000*1000), as.numeric(tmp_end)/ (1*1000*1000)),
    #                        y=c(max(tmp_prop_chr$high), max(tmp_prop_chr$high))), fill='gray')+
    geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=high.x),color='#E64B35FF')+
    geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=low.x),color='#4DBBD5FF')+
    geom_text_repel(data=label_gene, mapping = aes(x=start, y=high.x, label=gene_name.y), nudge_y = 1)+
    lims(x=c(left_x, right_x))+
    labs(title = paste0(tmp_chr,tmp_band), y=laby, x='Genomic_pos(Mb)')+
    #scale_x_continuous(breaks = c(left_x,right_x), labels = c(0, round(right_x-left_x)))
    theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))+
    theme_bw()+
    theme(panel.grid = element_blank())
  pg

  pg_list[[paste0(tmp_chr,tmp_band)]] = pg
}
cowplot::plot_grid(plotlist = pg_list, ncol=2)

# 只看一个region
i=4
tmp_chr = label_all_FC[i, 'chr']
tmp_band=label_all_FC[i, 'band']
tmp_start = label_all_FC[i, 'start']
tmp_end = label_all_FC[i, 'end']
tmp_prop = plot_gain
tmp_prop$seg_start = (tmp_prop$seg_start-tmp_prop$absStart)
tmp_prop$seg_end = (tmp_prop$seg_end-tmp_prop$absStart)
tmp_prop=tmp_prop[,c(1,5,6, 4,2,3,7,8,9)]
colnames(tmp_prop)[1:3] = c('chrom', 'start', 'end')
tmp_prop$chrom = as.character(tmp_prop$chrom)
tmp_prop = bed_intersect(tmp_prop,gene_info) %>% as.data.frame()

tmp_prop$start = tmp_prop$start.x / (1*1000*1000)
#
tmp_prop_chr = tmp_prop[tmp_prop$chrom==tmp_chr,]
tmp_chrinfo_band = chrinfo_band
tmp_chrinfo_band$chromNum = gsub('chr0','', tmp_chrinfo_band$chromNum)
tmp_chrinfo_band$chromNum = gsub('chr','', tmp_chrinfo_band$chromNum)
tmp_chrinfo_band$chromNum = gsub('X','23', tmp_chrinfo_band$chromNum)
tmp_chrinfo_band$chromNum = gsub('Y','24', tmp_chrinfo_band$chromNum)
tmp_chrinfo_band = tmp_chrinfo_band %>%filter(chromNum==tmp_chr&bandname==tmp_band)
tmp_chrinfo_band$absStart = as.numeric(tmp_chrinfo_band$absStart)
tmp_chrinfo_band$absEnd = as.numeric(tmp_chrinfo_band$absEnd)
#
# pg = plot_genome_chr(chr_num = tmp_chr,bg_color = 'white',
#                      start = tmp_chrinfo_band$absStart,
#                      end=tmp_chrinfo_band$absEnd)
left_x = as.numeric(tmp_chrinfo_band$chromStart)/ (1*1000*1000)
right_x = as.numeric(tmp_chrinfo_band$chromEnd)/ (1*1000*1000)
label_gene =tmp_prop_chr
label_gene = label_gene[match(driver_gene,label_gene$gene_name.y), ] %>% na.omit()
pg=ggplot()+
  #geom_vline(xintercept = tmp_chrinfo_band$absStart, color='gray', linewidth=0.1, linetype='dashed')+
  #geom_area(mapping = aes(x=c(as.numeric(tmp_start)/ (1*1000*1000), as.numeric(tmp_end)/ (1*1000*1000)),
  #                        y=c(max(tmp_prop_chr$high), max(tmp_prop_chr$high))), fill='gray')+
  geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=high.x),color='#E64B35FF')+
  geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=low.x),color='#4DBBD5FF')+
  geom_text_repel(data=label_gene, mapping = aes(x=start, y=high.x, label=gene_name.y), nudge_y = 1)+
  lims(x=c(left_x, right_x))+
  labs(title = paste0(tmp_chr,tmp_band), y=laby, x='Genomic_pos(Mb)')+
  #scale_x_continuous(breaks = c(left_x,right_x), labels = c(0, round(right_x-left_x)))
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))+
  theme_bw()+
  theme(panel.grid = element_blank())
pg



ggplot()+
  #geom_vline(xintercept = tmp_chrinfo_band$absStart, color='gray', linewidth=0.1, linetype='dashed')+
  #geom_area(mapping = aes(x=c(as.numeric(tmp_start)/ (1*1000*1000), as.numeric(tmp_end)/ (1*1000*1000)),
  #                        y=c(max(tmp_prop_chr$high), max(tmp_prop_chr$high))), fill='gray')+
  geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=high.x),color='#E64B35FF')+
  geom_line(data=tmp_prop_chr, mapping = aes(x=start, y=low.x),color='#4DBBD5FF')+
  geom_text_repel(data=label_gene, mapping = aes(x=start, y=high.x, label=gene_name.y), nudge_y = 1)+
  lims(x=c(0, 1.5))+
  labs(title = paste0(tmp_chr,tmp_band), y=laby, x='Genomic_pos(Mb)')+
  #scale_x_continuous(breaks = c(left_x,right_x), labels = c(0, round(right_x-left_x)))
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))+
  theme_bw()+
  theme(panel.grid = element_blank())


##### 基因表达与异质性相关性 ######
tcga_ITH = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                      sep='\t', header = T)

TCGA_bulk = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data//TCGA.BRCA.sampleMap_HiSeqV2',
                       header = T, row.names = 'sample', check.names = F)

same_sample = intersect(tcga_ITH$array, colnames(TCGA_bulk))

tcga_ITH = tcga_ITH[tcga_ITH$array%in%same_sample, ]
rownames(tcga_ITH) = tcga_ITH$array

tmp_gene = gene_info %>% filter(gene_info$chrom=='10',end>40200000)
#tmp_gene = tmp_prop_chr%>%filter(start>=150, start<=180)

tmp_gene = tmp_gene[!grepl('ENSG|LINC|RNA|RNU|SNO', tmp_gene$gene_name), ]

#keep_driver_gene = label_gene_all[label_gene_all$group=='Gain',]
keep_driver_gene = tmp_gene$gene_name#label_gene_all %>% filter(chrom==5)

TCGA_bulk = TCGA_bulk[intersect(rownames(TCGA_bulk), keep_driver_gene), same_sample]

exp_ITH_cor = cor(tcga_ITH$Subclonal.genome.fraction,  t(TCGA_bulk[, rownames(tcga_ITH)]))
exp_ITH = as.data.frame(cbind('ITH'=tcga_ITH$Subclonal.genome.fraction, t(TCGA_bulk[, rownames(tcga_ITH)])))

exp_ITH_cor = cor(na.omit(exp_ITH))
exp_ITH_cor = data.frame(cor_res = exp_ITH_cor[1,2:ncol(exp_ITH_cor)])
exp_ITH_cor$label = rownames(exp_ITH_cor)
exp_ITH_cor = exp_ITH_cor %>% arrange(cor_res)
exp_ITH_cor$x = 1:nrow(exp_ITH_cor)

#exp_ITH_cor$chr = label_gene_all[match(exp_ITH_cor$label, label_gene_all$gene_name.y), 'chrom']
exp_ITH_cor %>%
  ggplot(aes(x=x, y=cor_res))+
  geom_point(size=1)+
  #geom_text(aes(label=ifelse(cor_res>0.3, label, '')), nudge_x = -12)+
  scale_color_igv()+
  theme_classic()
## 箱线图验证拷贝数差异区域 #####


### TNBC 生存分析 #####
phenotype_file <- read.table("/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
                             header = T, sep = '\t', quote = "")
phenotype_colnames <- as.data.frame(colnames(phenotype_file))
#读取文件并去读列名（penotype的类型）
table(phenotype_file$breast_carcinoma_estrogen_receptor_status)#查看雌激素受体(ER)的表达状态。
table(phenotype_file$breast_carcinoma_progesterone_receptor_status)#查看孕激素受体（PR）状态
table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status)#查看原癌基因Her-2状态

phenotype_colnames <- colnames(phenotype_file)[ grep("receptor_status",colnames(phenotype_file))]
eph <- phenotype_file[,colnames_num[1:3]]
tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]


#
library(survival)
library(survminer)

metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

TCGA_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                          sep='\t', header = T, check.names = F)
colnames(TCGA_cna) = gsub('-01', '',colnames(TCGA_cna))
colnames(TCGA_cna) = gsub('-11', '',colnames(TCGA_cna))


all_pvalue = c()
for(i in rev(rownames(exp_ITH_cor))){
  print(i)
  select_gene = i#'PDCD6'
  metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                               sep='\t', header = T, check.names = F)
  rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
  same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

  # setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
  if(select_gene %in% metabric_cna$Hugo_Symbol){
    metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
    metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']<0, 'Gain', 'other')
    metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
    x1 = surv_pvalue(km)
  }else{
    x1 = list('pval'=1)
  }

  #
  metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                               sep='\t', header = T, check.names = F)

  rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
  same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(TCGA_cna))

  metabric_clinic = metabric_clinic[same_sample, ]
  if(select_gene %in% TCGA_cna$Hugo_Symbol){
    metabric_clinic[, 'group'] = unlist(TCGA_cna[TCGA_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
    metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']<0, 'Gain', 'other')

    metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
    x2=surv_pvalue(km)


    metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                                 sep='\t', header = T, check.names = F)

    rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
    same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(TCGA_cna))

    metabric_clinic = metabric_clinic[same_sample, ]
    metabric_clinic[, 'group'] = unlist(TCGA_cna[TCGA_cna$Hugo_Symbol==select_gene, same_sample])
    metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
    metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']<0, 'Gain', 'other')

    metabric_clinic$PFS_STATUS = as.numeric(sapply(metabric_clinic$PFS_STATUS, function(x)strsplit(x,':')[[1]][1]))
    km <- survfit(Surv(PFS_MONTHS, PFS_STATUS)~group, data = metabric_clinic)
    x3=surv_pvalue(km)

  }else{
    x2 = list('pval'=1)
    x3 = list('pval'=1)
  }


  all_pvalue = rbind(all_pvalue, c(i, x1$pval, x2$pval, x3$pval))
}
all_pvalue = as.data.frame(all_pvalue)
all_pvalue[,2] = as.numeric(all_pvalue[,2])
all_pvalue[,3] = as.numeric(all_pvalue[,3])
all_pvalue[,4] = as.numeric(all_pvalue[,4])

select_gene = 'PRDM2'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)
rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

# setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']<0, 'Loss', 'Neutral')
metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="OS(METABRIC)"
)


# TCGA
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                          sep='\t', header = T, check.names = F)
colnames(metabric_cna) = gsub('-01', '',colnames(metabric_cna))
colnames(metabric_cna) = gsub('-11', '',colnames(metabric_cna))

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)

rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

metabric_clinic = metabric_clinic[same_sample, ]
metabric_clinic[, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']<0, 'Loss', 'Neutral')

metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="OS(TCGA)"
)
# PFS
metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)

rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

metabric_clinic = metabric_clinic[same_sample, ]
metabric_clinic[, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,-1,-2))
metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']<0, 'Loss', 'Neutral')

metabric_clinic$PFS_STATUS = as.numeric(sapply(metabric_clinic$PFS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(PFS_MONTHS, PFS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "solid",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="PFS(TCGA)"
)


##### driver gene 统计 ######
# 统计基因数量
driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因_Cosmic.csv', sep=',', header = T, check.names = F)
breast_pos = grepl('breast',driver_gene$`Tumour Types(Germline)`)
breast_pos = breast_pos | grepl('breast',driver_gene$`Tumour Types(Germline)`)
driver_gene = driver_gene[breast_pos, 'Gene Symbol']

driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因.txt', sep='\t', header = T, check.names = F)
driver_gene = driver_gene[driver_gene$cancer_type_abbr=='BRCA', 'driver_gene']
driver_gene = strsplit(gsub(' ', '',driver_gene), ',')[[1]]

gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')
gene_info = gene_info[gene_info$gene_name%in%driver_gene, ]

driver_gene = gene_info$gene_name

rownames(patient_index_score_pde) = patient_index_score_pde$patient
driver_gene_rate = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  # tmp_cna = cna_gene[[prefix]][, driver_gene$driver_gene]
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]
  cna_loc = as.data.frame(t(sapply(colnames(tmp_cna), function(x)strsplit(x, '_')[[1]])))
  cna_loc$V2 = as.numeric(cna_loc$V2)
  cna_loc$V3 = c(as.numeric(cna_loc$V2[2:nrow(cna_loc)]), 0)
  cna_loc = cna_loc[cna_loc$V2<cna_loc$V3, ]
  cna_loc$cna_name = rownames(cna_loc)
  colnames(cna_loc) = c('chrom', 'start', 'end', 'cna_name')
  res = bed_intersect(cna_loc,gene_info) %>% as.data.frame()

  tmp_cna = as.data.frame(t(tmp_cna[, res$cna_name.x]))
  tmp_cna$gene = res$gene_name.y
  tmp_cna = tmp_cna %>%
    group_by(gene) %>%
    summarise_all(list(mean)) %>% as.data.frame()
  rownames(tmp_cna) = tmp_cna[,1]
  tmp_cna = tmp_cna[, -1]
  tmp_cna = round(tmp_cna)

  for(g in driver_gene){
    gain_num = sum(tmp_cna[g,]>2)/ncol(tmp_cna)
    loss_num = sum(tmp_cna[g,]<2)/ncol(tmp_cna)
    driver_gene_rate = rbind(driver_gene_rate, c(new_prefix, tmp_group, g, gain_num, loss_num))
  }
}
driver_gene_rate = as.data.frame(driver_gene_rate)
colnames(driver_gene_rate) = c('patient', 'PDE_group', 'gene', 'gain_num', 'loss_num')
driver_gene_rate$gain_num = as.numeric(driver_gene_rate$gain_num)
driver_gene_rate$loss_num = as.numeric(driver_gene_rate$loss_num)

a = driver_gene_rate %>%
  group_by(PDE_group, gene) %>%
  summarise(gain_num=median(gain_num), loss_num=median(loss_num)) %>%
  mutate(gain_num=ifelse(PDE_group=='PDE_low', -gain_num,gain_num),
         loss_num=ifelse(PDE_group=='PDE_low', -loss_num,loss_num)) %>%
  as.data.frame()%>%
  mutate(gene = fct_reorder(gene, abs(gain_num))) %>%
  reshape2::melt(id=c('PDE_group', 'gene'), variable.name='type', value.name = 'prop') %>%
  ggplot()+
  #geom_bar(aes(x=loss_num+gain_num, y=gene), stat='identity', fill='gray')+
  #geom_bar(aes(x=loss_num, y=gene), stat='identity', fill='#FFF5F5')+
  geom_bar(aes(x=prop, y=gene, fill=PDE_group), stat='identity',  color='black', width=0.8)+
  geom_text(aes(x=ifelse(prop>0, prop+0.1, prop-0.1), y=gene, label=abs(round(prop,2))))+
  scale_fill_manual(values=c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  # geom_vline(xintercept = 0, linewidth=2, color='white')+
  facet_wrap(~type)+
  theme_classic()
a

#### gene 与ITH 相关性 #####
tcga_ITH = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                      sep='\t', header = T)

TCGA_bulk = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data//TCGA.BRCA.sampleMap_HiSeqV2',
                       header = T, row.names = 'sample', check.names = F)

same_sample = intersect(tcga_ITH$array, colnames(TCGA_bulk))

tcga_ITH = tcga_ITH[tcga_ITH$array%in%same_sample, ]
rownames(tcga_ITH) = tcga_ITH$array

TCGA_bulk = TCGA_bulk[driver_gene, same_sample]

exp_ITH = as.data.frame(cbind('ITH'=tcga_ITH$Subclonal.genome.fraction, t(TCGA_bulk[, rownames(tcga_ITH)])))
exp_ITH = reshape2::melt(exp_ITH, id='ITH', variable.name='gene', value.name = 'expr')

exp_ITH %>%
  ggplot(aes(x=ITH, y=expr))+
  geom_point(size=0.5, color='gray')+
  facet_wrap(~gene, scales = 'free_y')+
  geom_smooth(method='lm', formula='y~x', se=F, color='black', linetype='dashed')+
  stat_cor()+
  theme_classic()


##### 识别driver 基因 ######
tcga_ITH = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                      sep='\t', header = T)

TCGA_bulk = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data//TCGA.BRCA.sampleMap_HiSeqV2',
                       header = T, row.names = 'sample', check.names = F)

same_sample = intersect(tcga_ITH$array, colnames(TCGA_bulk))

tcga_ITH = tcga_ITH[tcga_ITH$array%in%same_sample, ]
rownames(tcga_ITH) = tcga_ITH$array

driver_gene = rownames(TCGA_bulk)

exp_ITH = as.data.frame(cbind('ITH'=tcga_ITH$Subclonal.genome.fraction, t(TCGA_bulk[, rownames(tcga_ITH)])))

rownames(patient_index_score_pde) = patient_index_score_pde$patient
driver_gene_rate = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  # tmp_cna = cna_gene[[prefix]][, driver_gene$driver_gene]
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]
  cna_loc = as.data.frame(t(sapply(colnames(tmp_cna), function(x)strsplit(x, '_')[[1]])))
  cna_loc$V2 = as.numeric(cna_loc$V2)
  cna_loc$V3 = c(as.numeric(cna_loc$V2[2:nrow(cna_loc)]), 0)
  cna_loc = cna_loc[cna_loc$V2<cna_loc$V3, ]
  cna_loc$cna_name = rownames(cna_loc)
  colnames(cna_loc) = c('chrom', 'start', 'end', 'cna_name')
  res = bed_intersect(cna_loc,gene_info) %>% as.data.frame()

  tmp_cna = as.data.frame(t(tmp_cna[, res$cna_name.x]))
  tmp_cna$gene = res$gene_name.y
  tmp_cna = tmp_cna %>%
    group_by(gene) %>%
    summarise_all(list(mean)) %>% as.data.frame()
  rownames(tmp_cna) = tmp_cna[,1]
  tmp_cna = tmp_cna[, -1]
  tmp_cna = round(tmp_cna)

  for(g in driver_gene){
    gain_num = sum(tmp_cna[g,]>2)/ncol(tmp_cna)
    loss_num = sum(tmp_cna[g,]<2)/ncol(tmp_cna)
    driver_gene_rate = rbind(driver_gene_rate, c(new_prefix, tmp_group, g, gain_num, loss_num))
  }
}


##
driver_gene_rate2_gain = c()
driver_gene_rate2_loss = c()
tmp_name = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  tmp_name = c(tmp_name, new_prefix)
  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  # tmp_cna = cna_gene[[prefix]][, driver_gene$driver_gene]
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]

  driver_gene_rate2_gain = rbind(driver_gene_rate2_gain, colMeans(tmp_cna>2))
  driver_gene_rate2_loss = rbind(driver_gene_rate2_loss, colMeans(tmp_cna<2))
}
rownames(driver_gene_rate2_gain) = tmp_name
rownames(driver_gene_rate2_loss) = tmp_name

high_gain_rate = colSums(driver_gene_rate2_gain[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_high', 'patient'], ])
low_gain_rate = colSums(driver_gene_rate2_gain[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_low', 'patient'], ])
plot_gain = data.frame(high=high_gain_rate, low=low_gain_rate)
plot_gain$index = colnames(driver_gene_rate2_gain)
plot_gain$chr = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][1]))
plot_gain$seg_start = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][2]))
plot_gain$seg_end = as.numeric(sapply(plot_gain$index, function(x)strsplit(x, '_')[[1]][3]))
plot_gain$x = 1:nrow(plot_gain)
plot_gain %>%
  ggplot(aes(x=x))+
  geom_area(aes(y=high), fill='red', alpha=0.6)+
  geom_area(aes(y=low), fill='blue', alpha=0.6)

#
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
chrinfo_chr_absstart$chromNum = gsub('chr0', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('chr', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('X', '23', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('Y', '24', chrinfo_chr_absstart$chromNum)

plot_gain = merge(plot_gain, chrinfo_chr_absstart, by.x='chr', by.y = 'chromNum')
plot_gain$seg_start = plot_gain$seg_start+plot_gain$absStart
plot_gain$seg_end = plot_gain$seg_end+plot_gain$absStart



p = plot_genome(chr_bg_color='white', show_x_axis = T)

a1=p +
  #geom_hline(yintercept = 0, color='black')+
  #geom_area(data=plot_gain, aes(y=high), fill='red', alpha=0.6)+
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=high),fill='red', alpha=0.5)+
  geom_area(data=plot_gain, mapping = aes(x=seg_start, y=low),fill='blue', alpha=0.5)+
  # geom_line(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC)))+
  geom_text_repel(data=know_gene_data,
                  mapping = aes(x=seg_start, y=7, label=gene),
                  nudge_y = 1,
                  #arrow = arrow(length = unit(0.02, "npc")),box.padding =1,
                  max.overlaps = Inf,

  )+
  #scale_fill_manual(values = c('Up'='red', 'Down'='blue'))+
  #scale_size(range = c(0.1,4))+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))

#
high_loss_rate = colSums(driver_gene_rate2_loss[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_high', 'patient'], ])
low_loss_rate = colSums(driver_gene_rate2_loss[patient_index_score_pde[patient_index_score_pde$PDE_group=='PDE_low', 'patient'], ])
plot_loss = data.frame(high=high_loss_rate, low=low_loss_rate)
plot_loss$index = colnames(driver_gene_rate2_loss)
plot_loss$chr = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][1]))
plot_loss$seg_start = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][2]))
plot_loss$seg_end = as.numeric(sapply(plot_loss$index, function(x)strsplit(x, '_')[[1]][3]))
plot_loss$x = 1:nrow(plot_loss)
plot_loss = merge(plot_loss, chrinfo_chr_absstart, by.x='chr', by.y = 'chromNum')
plot_loss$seg_start = plot_loss$seg_start+plot_loss$absStart
plot_loss$seg_end = plot_loss$seg_end+plot_loss$absStart


a2=p +
  #geom_hline(yintercept = 0, color='black')+
  #geom_area(data=plot_gain, aes(y=high), fill='red', alpha=0.6)+
  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=high),fill='red', alpha=0.5)+
  geom_area(data=plot_loss, mapping = aes(x=seg_start, y=low),fill='blue', alpha=0.5)+
  # geom_line(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC)))+
  geom_text_repel(data=know_gene_data,
                  mapping = aes(x=seg_start, y=7, label=gene),
                  nudge_y = 1,
                  #arrow = arrow(length = unit(0.02, "npc")),box.padding =1,
                  max.overlaps = Inf,

  )+
  #scale_fill_manual(values = c('Up'='red', 'Down'='blue'))+
  #scale_size(range = c(0.1,4))+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))
cowplot::plot_grid(a1,a2, ncol=1, align = 'vh')
##



driver_gene_rate = as.data.frame(driver_gene_rate)
colnames(driver_gene_rate) = c('patient', 'PDE_group', 'gene', 'gain_num', 'loss_num')
driver_gene_rate$gain_num = as.numeric(driver_gene_rate$gain_num)
driver_gene_rate$loss_num = as.numeric(driver_gene_rate$loss_num)
saveRDS(driver_gene_rate, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_driver_gene_rate.rds')

old_driver_gene_rate = driver_gene_rate

driver_gene_rate = old_driver_gene_rate
driver_gene_rate = na.omit(driver_gene_rate)
gene_fc = c()
for(i in unique(driver_gene_rate$gene)){
  tmp = driver_gene_rate[driver_gene_rate$gene==i,]
  high = tmp[tmp$PDE_group=='PDE_high', ]
  low = tmp[tmp$PDE_group=='PDE_low', ]
  #a=(high$gain_num + high$loss_num)
  #b=low$gain_num+low$loss_num
  a=high$gain_num
  b=low$gain_num
  tmp_fc = mean(a) / mean(b)
  tmp_pvalue = t.test(a, b)


  gene_fc = rbind(gene_fc, c(i, tmp_fc, tmp_pvalue$p.value))
}

gene_fc = as.data.frame(gene_fc)
colnames(gene_fc) = c('gene', 'FC', 'pvalue')
gene_fc$FC = as.numeric(gene_fc$FC)
gene_fc$pvalue = as.numeric(gene_fc$pvalue)
gene_fc$adj_p = p.adjust(gene_fc$pvalue)


# cor
cor_res = c()
exp_ITH = na.omit(exp_ITH)
for(i in 1:ncol(exp_ITH)){
  cor_res = c(cor_res, cor(exp_ITH[,1], exp_ITH[,i]))
}
names(cor_res) = colnames(exp_ITH)
#
gene_fc$cor = cor_res[gene_fc$gene]

gene_fc_filter = gene_fc %>% filter(FC>1, pvalue<0.01,cor>0.2)

gene_fc_filter %>% arrange(-FC, -cor, pvalue)

gene_fc_filter %>%
  ggplot(aes(x=FC, y=cor, color=-log(pvalue)))+
  geom_point()+
  geom_text(aes(label=gene))+
  geom_vline(xintercept = 1, )+
  geom_hline(yintercept = 0.2)+
  scale_color_gradientn(colours = c('blue','white', 'red'))+
  theme_classic()

keep_gene = gene_fc_filter$gene


a = driver_gene_rate %>%
  filter(gene%in%keep_gene) %>%
  group_by(PDE_group, gene) %>%
  summarise(gain_num=median(gain_num), loss_num=median(loss_num)) %>%
  mutate(gain_num=ifelse(PDE_group=='PDE_low', -gain_num,gain_num),
         loss_num=ifelse(PDE_group=='PDE_low', -loss_num,loss_num)) %>%
  as.data.frame()%>%
  mutate(gene = fct_reorder(gene, abs(gain_num))) %>%
  reshape2::melt(id=c('PDE_group', 'gene'), variable.name='type', value.name = 'prop') %>%
  ggplot()+
  #geom_bar(aes(x=loss_num+gain_num, y=gene), stat='identity', fill='gray')+
  #geom_bar(aes(x=loss_num, y=gene), stat='identity', fill='#FFF5F5')+
  geom_bar(aes(x=prop, y=gene, fill=PDE_group), stat='identity',  color='black', width=0.8)+
  geom_text(aes(x=ifelse(prop>0, prop+0.1, prop-0.1), y=gene, label=abs(round(prop,2))))+
  scale_fill_manual(values=c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  # geom_vline(xintercept = 0, linewidth=2, color='white')+
  facet_wrap(~type)+
  theme_classic()
a

reshape2::melt(exp_ITH[, c('ITH', keep_gene)], id='ITH', variable.name='gene', value.name = 'expr') %>%
  ggplot(aes(x=ITH, y=expr))+
  geom_point(size=0.5, color='gray')+
  facet_wrap(~gene, scales = 'free')+
  geom_smooth(method='lm', formula='y~x', se=F, color='black', linetype='dashed')+
  stat_cor()+
  theme_classic()


know_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因_Cosmic.csv', sep=',', header = T, check.names = F)
#breast_pos = grepl('breast',know_gene$`Tumour Types(Germline)`)
#breast_pos = breast_pos | grepl('breast',know_gene$`Tumour Types(Germline)`)
#know_gene = know_gene[breast_pos, 'Gene Symbol']

know_gene = know_gene$`Gene Symbol`



know_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因.txt', sep='\t', header = T, check.names = F)
know_gene = know_gene[know_gene$cancer_type_abbr=='BRCA', 'driver_gene']
know_gene = strsplit(gsub(' ', '',know_gene), ',')[[1]]


intersect(know_gene, keep_gene)

p = plot_genome(chr_bg_color='white', show_x_axis = T)

gene_fc_data = gene_fc
gene_fc_data[, colnames(gene_info)] = gene_info[match(gene_fc_data$gene, gene_info$gene_name), ]
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
chrinfo_chr_absstart$chromNum = gsub('chr0', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('chr', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('X', '23', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('Y', '24', chrinfo_chr_absstart$chromNum)

gene_fc_data = merge(gene_fc_data, chrinfo_chr_absstart, by.x='chrom', by.y = 'chromNum')
gene_fc_data$seg_start = gene_fc_data$start+gene_fc_data$absStart
gene_fc_data$seg_end = gene_fc_data$end+gene_fc_data$absStart

gene_fc_data$color = ifelse(gene_fc_data$FC>1, 'Up', 'Down')

know_gene_data = gene_fc_data[gene_fc_data$gene%in%know_gene,]

pathway = read.csv('/Users/lab/wangxin/MEDALT/SCTMS/data/pathway_cell2018.csv', header = T)
pathway = subset(pathway, type%in%c('TSG', 'OG'))
keep_gene = pathway[pathway$type=='TSG','Gene']

#gene_fc_filter = gene_fc# %>% filter(FC>2, pvalue<0.05, cor>0)
#keep_gene = gene_fc_filter$gene
know_gene_data = gene_fc_data[gene_fc_data$gene%in%know_gene,]
know_gene_data = know_gene_data %>% filter( pvalue<0.05)
p +
  geom_hline(yintercept = 0, color='black')+
  # geom_segment(data=gene_fc_data, mapping = aes(x=seg_start, xend=seg_end,
  #                                               y=log2(FC), yend=log2(FC),
  #                                               color=color,
  #                                               linewidth=-log(pvalue)))+
  geom_point(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC),fill=color), size=0.8, shape=21, stroke=0)+
  # geom_line(data=gene_fc_data, mapping = aes(x=seg_start, y=log2(FC)))+
  geom_text_repel(data=know_gene_data,
                  mapping = aes(x=seg_start, y=log2(FC), label=gene),
                  arrow = arrow(length = unit(0.02, "npc")),box.padding =1,
                  max.overlaps = Inf,

                  )+
  scale_fill_manual(values = c('Up'='red', 'Down'='blue'))+
  scale_size(range = c(0.1,4))+
  theme(axis.text.x = element_text(angle = 50, hjust=1,vjust=1))

## 高斯混合 ###
library(mixtools)
logfc = log2(gene_fc_data[gene_fc_data$FC>1, 'FC'])

em <- normalmixEM(logfc)
plot(em, whichplots = 2)




### TNBC 生存分析 #####
phenotype_file <- read.table("/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
                             header = T, sep = '\t', quote = "")
phenotype_colnames <- as.data.frame(colnames(phenotype_file))
#读取文件并去读列名（penotype的类型）
table(phenotype_file$breast_carcinoma_estrogen_receptor_status)#查看雌激素受体(ER)的表达状态。
table(phenotype_file$breast_carcinoma_progesterone_receptor_status)#查看孕激素受体（PR）状态
table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status)#查看原癌基因Her-2状态

phenotype_colnames <- colnames(phenotype_file)[ grep("receptor_status",colnames(phenotype_file))]
eph <- phenotype_file[,colnames_num[1:3]]
tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]


#
library(survival)
library(survminer)
select_gene = 'RBBP8'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)
rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']<0, 'Gain', 'other')
metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "strata",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="OS"
)


# TCGA
select_gene = 'EXT1'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_cna.txt',
                          sep='\t', header = T, check.names = F)
colnames(metabric_cna) = gsub('-01', '',colnames(metabric_cna))
colnames(metabric_cna) = gsub('-11', '',colnames(metabric_cna))

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)

rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']!=0, 'Gain', 'other')
metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "Exp",#改变图例名称
           #legend.labs = c("High", "Low"),
           linetype = "strata",# 改变线条类型
           xlab="M",
           risk.table = F,
           ylab="OS"
)
