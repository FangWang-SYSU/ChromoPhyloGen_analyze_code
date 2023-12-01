source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
library(ggnewscale)
library(ggtreeExtra)
library(ggpubr)

dna_file = c('Fig3a_sub_dataset_P9T_20190409_1.txt',
             'Fig1c_sub_dataset_P19BT_20181009.txt',
             'Fig4a_sub_dataset_P9T_4_20200304.txt')
#
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sub/'
subNG_list = list()
for(prefix in dna_file){
  print(prefix)
  #min_cell_num=
  #aa = load_scTrace(scTrace_dir, prefix)
  sctc = scTraceClass(scTrace_dir,
                      prefix,
                      layout = 'slanted',
                      rate = 0.8, min_cell_num = 1, min_cell_num2 = 1)

  subNG_list[[prefix]] = sctc
}

saveRDS(subNG_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all.rds')
NG_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all.rds')

#### 计算score ######
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sub/'
NG_list_with_score = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = NG_list[[prefix]]
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

  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/进化树不同指标着色/tree_', prefix,'.pdf'),
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
  new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=rearrange_score),
               color = "grey",
               offset = 0.05,size = 0.05,width=2,
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
  new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=BFB),
               color = "grey",
               offset = 0.05,size = 0.05, width=2,
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
               offset = 0.05,size = 0.05, width=2,
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
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/进化树不同指标着色/hp_', prefix,'.pdf'),
         p2, width=160, height = 160, units='mm', dpi = 600)

  NG_list_with_score[[prefix]] = sctc
}
saveRDS(NG_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all_with_score.rds')
#saveRDS(NG_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all_with_score.rds')
NG_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all_with_score.rds')

dna_file = c('Fig1b_dataset_P9T_20190331.txt',
             'Fig1c_dataset_P19BT_20181009.txt',
             'Fig2c_dataset_P9T_20190409_2.txt',
             'Fig3a_dataset_P9T_20190409_1.txt',
             'Fig4a_dataset_P9T_4_20200304.txt',
             'Fig4b_dataset_photoconverted_M3_M4.txt'
)

scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG/'
NG_list = list()
for(prefix in dna_file){
  print(prefix)
  #min_cell_num=
  #aa = load_scTrace(scTrace_dir, prefix)
  if(prefix=='Fig4b_dataset_photoconverted_M3_M4.txt'){
    sctc = scTraceClass(scTrace_dir,
                        prefix,
                        layout = 'slanted',
                        rate = 0.8, min_cell_num = 1, min_cell_num2 = 1)
  }else{
    sctc = scTraceClass(scTrace_dir,
                        prefix,
                        layout = 'slanted')
  }


  NG_list[[prefix]] = sctc
}

saveRDS(NG_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sctc_all.rds')
NG_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sctc_all.rds')

#### 计算score ######
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG/'
NG_list_with_score = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = NG_list[[prefix]]
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

  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/进化树不同指标着色/tree_', prefix,'.pdf'),
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
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=rearrange_score),
               color = "grey",
               offset = 0.05,size = 0.05,width=2,
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
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=BFB),
               color = "grey",
               offset = 0.05,size = 0.05, width=2,
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
               offset = 0.05,size = 0.05, width=2,
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
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/进化树不同指标着色/hp_', prefix,'.pdf'),
         p2, width=160, height = 160, units='mm', dpi = 600)

  NG_list_with_score[[prefix]] = sctc
}
saveRDS(NG_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sctc_all_with_score.rds')

subNG_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all.rds')
NG_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sctc_all.rds')
NG_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sctc_all_with_score.rds')
subNG_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all_with_score.rds')


NG_list = c(subNG_list_with_score, NG_list_with_score)
##### 与其他方法比较细胞关系划分准确性 #####
library(ggplot2)
library(ggpubr)
library(treeio)
library(ape)
library(pROC)
library(phangorn)
library(ggtree)
library(phangorn)
library(dplyr)
library(ComplexHeatmap)

newick_to_phylo <- function(newick_dir){
  newick_file = file(newick_dir,open='r')
  newick_str = readLines(newick_file, n = 1)
  close(newick_file)
  newick_str = paste0('(', newick_str, ');')
  newick_str <- read.tree(text = newick_str)
  return(newick_str)
}
newick_to_hclust <- function(newick_dir){
  newick_str <- newick_to_phylo(newick_dir)

  dendrogram <- unroot(chronos(newick_str))
  dendrogram = as.hclust.phylo(dendrogram)
  return(dendrogram)
}
sitka_newick_to_phylo<-function(newick_dir){
  newick_file = file(newick_dir,open='r')
  newick_str = readLines(newick_file, n = 1)
  close(newick_file)
  newick_str <- read.tree(text = newick_str)
  return(newick_str)
}

medalt_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/'
sitka_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/sitka/'
medicc_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/medicc2/'

dna_file = c('Fig1b_dataset_P9T_20190331.txt',
             'Fig1c_dataset_P19BT_20181009.txt',
             'Fig2c_dataset_P9T_20190409_2.txt',
             'Fig3a_dataset_P9T_20190409_1.txt',
             'Fig4a_dataset_P9T_4_20200304.txt',
             'Fig4b_dataset_photoconverted_M3_M4.txt',
             'Fig3a_sub_dataset_P9T_20190409_1.txt',
             'Fig1c_sub_dataset_P19BT_20181009.txt',
             'Fig4a_sub_dataset_P9T_4_20200304.txt'
)

all_trees = c()
for(f in names(NG_list)){
  print(f)
  sctc_tree = NG_list[[f]]
  medalt_tree = sitka_newick_to_phylo(paste0(medalt_dir, 'medalt_',f, '_tree.txt'))
  # sitka tree
  sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, 'sitka_',f, '.newick'))
  sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))
  #medicc2
  if(file.exists(paste0(medicc_dir,'medicc2_',strsplit(f, '\\.')[[1]][1],'_final_tree.new'))){
    medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,'medicc2_',strsplit(f, '\\.')[[1]][1],'_final_tree.new'))
    medicc2_tree = drop.tip(medicc2_tree, 'diploid')
    medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'), function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
  }else{
    medicc2_tree = rtree(10)
  }
  #
  sim_data = sctc_tree$orig.data$all_node_data[!grepl('virtual|root',rownames(sctc_tree$orig.data$all_node_data)), ]
  sim_data = t(sim_data)
  edu_dist = dist(t(sim_data), method = "euclidean")
  #hclust_dist = hclust(edu_dist, method = "average")
  #cluster_tree = as.phylo.hclust(hclust_dist)
  # NJ
  pyd = as.phyDat(t(sim_data), type='USER', levels=c(0,1:max(sim_data)))
  mm = dist.hamming(pyd) # as.matrix(dist(t(sim_data), upper = T))
  NJ_tree<-NJ(mm)
  # random_tree
  set.seed(1234)
  random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
  # MP
  MP_tree = optim.parsimony(random_tree, pyd)
  #MP_tree = pratchet(pyd)
  # ML
  ML_tree = optim.pml(pml(NJ_tree, pyd))$tree

  all_trees[[f]] = list(sctc=sctc_tree$orig.data$tree,
                        medalt=medalt_tree,
                        sitka=sitka_tree,
                        medicc2=medicc2_tree,
                        NJ=NJ_tree,
                        MP=MP_tree,
                        ML=ML_tree
                        )
}

# tmp_tree = all_trees$Fig1c_sub_dataset_P19BT_20181009.txt$sctc
# ##### 单细胞进化树 rearrange着色 ####
# gtree_res = ggtree(tmp_tree, size=0.2)
# #gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
# a1=gtree_res +
#   geom_tippoint(aes(color=label))+
#   #scale_color_igv()+
#   #scale_color_gradientn(colours = c('yellow','red'))+
#   guides(color=guide_colourbar(title='Rearrange_score'))+
#   theme_tree(legend.position='right')+
#   #coord_flip()+
#   #scale_x_reverse()+
#   theme(legend.text = element_text(size=4),
#         legend.key.size = unit(2, 'mm'),
#         legend.title = element_text(size=6)
#   )
#
# tmp_tree = all_trees$Fig1c_sub_dataset_P19BT_20181009.txt$medalt
# ##### 单细胞进化树 rearrange着色 ####
# gtree_res = ggtree(tmp_tree, size=0.2)
# #gtree_res$data$group = CNA_mechnism[gtree_res$data$label, colorby]
# a2=gtree_res +
#   geom_tippoint(aes(color=label))+
#   #scale_color_igv()+
#   #scale_color_gradientn(colours = c('yellow','red'))+
#   guides(color=guide_colourbar(title='Rearrange_score'))+
#   theme_tree(legend.position='right')+
#   #coord_flip()+
#   #scale_x_reverse()+
#   theme(legend.text = element_text(size=4),
#         legend.key.size = unit(2, 'mm'),
#         legend.title = element_text(size=6)
#   )
#
# a1_data= fortify(all_trees$Fig1c_sub_dataset_P19BT_20181009.txt$sctc)
# a2_data= fortify(all_trees$Fig1c_sub_dataset_P19BT_20181009.txt$ML)
# a2_data$x <- a2_data$x + max(a1_data$x) + 1
#
# aa = bind_rows(a1_data, a2_data) %>%
#   filter(!is.na(label))
#
# a1 +
#   geom_tree(data = a2_data) +
#   geom_line(aes(x, y, group=label, color=node < 15), data=aa, alpha=.3)


### 统计成对细胞出现的次数 ####
# tip_dist <- function(tree, label=T){
#   x = fortify(tree)[, c('node', 'label')]
#   label_name = x$node
#   names(label_name) = x$label
#   r_label = data.frame(matrix(0, nrow=length(tree$tip.label),
#                               ncol = length(tree$tip.label)),
#                        row.names = tree$tip.label)
#   colnames(r_label) = tree$tip.label
#   for(i in tree$tip.label){
#     for(j in tree$tip.label){
#       d = nodepath(tree, from=label_name[i], to=label_name[j])
#       r_label[i,j] = length(d)
#     }
#   }
#   return(r_label)
# }
#

pheatmap::pheatmap(NG_list$Fig1b_dataset_P9T_20190331.txt$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('SC06.P9T.31032019_dedup', 'SC09.P9T.31032019_dedup')
c2=c('SC22.P9T.31032019_dedup','SC25.P9T.31032019_dedup', 'SC28.P9T.31032019_dedup')
all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$medalt)
a3_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sitka)
a4_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$NJ)
a6_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$MP)
a7_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$ML)

a1_data$group = ifelse(a1_data$label%in%c1,'a','b')
a2_data$group = ifelse(a2_data$label%in%c1,'a','b')
a3_data$group = ifelse(a3_data$label%in%c1,'a','b')
a4_data$group = ifelse(a4_data$label%in%c1,'a','b')
a5_data$group = ifelse(a5_data$label%in%c1,'a','b')
a6_data$group = ifelse(a6_data$label%in%c1,'a','b')
a7_data$group = ifelse(a7_data$label%in%c1,'a','b')


a1=ggtree(a1_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('SCTC')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a2=ggtree(a2_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('MEDALT')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a3=ggtree(a3_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('Sitka')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a4=ggtree(a4_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a5=ggtree(a5_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('NJ')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a6=ggtree(a6_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('MP')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
a7=ggtree(a7_data, size=0.2)+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+ggtitle('ML')+theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)

cowplot::plot_grid(a1, p2, nrow=1)



#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = ifelse(a1_data$label%in%c1,'a','b')
a2_data$group = ifelse(a2_data$label%in%c1,'a','b')
a3_data$group = ifelse(a3_data$label%in%c1,'a','b')
a4_data$group = ifelse(a4_data$label%in%c1,'a','b')
a5_data$group = ifelse(a5_data$label%in%c1,'a','b')
a6_data$group = ifelse(a6_data$label%in%c1,'a','b')
a7_data$group = ifelse(a7_data$label%in%c1,'a','b')

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa)


tmp_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
medalt_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$medalt
sitka_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$sitka
medicc2_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$medicc2
medicc2_tree$label = gsub('_NA_', '', medicc2_tree$label)
NJ_tree =  all_trees$Fig1b_dataset_P9T_20190331.txt$NJ
MP_tree = all_trees$Fig1b_dataset_P9T_20190331.txt$MP
ML_tree= all_trees$Fig1b_dataset_P9T_20190331.txt$ML




# heatmap
sctc = NG_list[['Fig1b_dataset_P9T_20190331.txt']]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/Fig1b_dataset_P9T_20190331.txt.pdf',
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)


########
file_name = 'Fig1c_sub_dataset_P19BT_20181009.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')

all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
a3_data= fortify(all_trees[[file_name]]$sitka)
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
a7_data= fortify(all_trees[[file_name]]$ML)
#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = 'a'#ifelse(a1_data$label%in%c1,'a','b')
a2_data$group = 'a'#ifelse(a2_data$label%in%c1,'a','b')
a3_data$group = 'a'#ifelse(a3_data$label%in%c1,'a','b')
a4_data$group = 'a'#ifelse(a4_data$label%in%c1,'a','b')
a5_data$group = 'a'#ifelse(a5_data$label%in%c1,'a','b')
a6_data$group = 'a'#ifelse(a6_data$label%in%c1,'a','b')
a7_data$group = 'a'#ifelse(a7_data$label%in%c1,'a','b')

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2, c2,c3,c4,c5))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2, c3, c4, c5),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa, size=0.5)


# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)


#######
file_name = 'Fig2c_dataset_P9T_20190409_2.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c('SC26.P9T.09042019_dedup','SC31.P9T.09042019_dedup','SC34.P9T.09042019_dedup','SC28.P9T.09042019_dedup','SC29.P9T.09042019_dedup','SC33.P9T.09042019_dedup', 'SC30.P9T.09042019_dedup')
c2 = c('SC32.P9T.09042019_dedup','SC35.P9T.09042019_dedup','SC36.P9T.09042019_dedup')
c3 = c('SC27.P9T.09042019_dedup','SC37.P9T.09042019_dedup')


all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp_drop = setdiff(all_trees[[file_name]]$sitka$tip.label, all_trees[[file_name]]$sctc$tip.label)
a3_data= fortify(drop.tip(all_trees[[file_name]]$sitka, tmp_drop))
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
a7_data= fortify(all_trees[[file_name]]$ML)
#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else{NA}})

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2, c3))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2, c3),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa, size=0.5)


# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)



#######
file_name = 'Fig3a_dataset_P9T_20190409_1.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('SC11.P9T.09042019_dedup','SC10.P9T.09042019_dedup','SC17.P9T.09042019_dedup')
c2=c('SC13.P9T.09042019_dedup','SC23.P9T.09042019_dedup')

all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp_drop = setdiff(all_trees[[file_name]]$sitka$tip.label, all_trees[[file_name]]$sctc$tip.label)
a3_data= fortify(drop.tip(all_trees[[file_name]]$sitka, tmp_drop))
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
a7_data= fortify(all_trees[[file_name]]$ML)
#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa, size=0.5)


# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)




######
file_name = 'Fig4a_dataset_P9T_4_20200304.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('X4.SC18.04032020_dedup','X4.SC03.04032020_dedup','X4.SC15.04032020_dedup')
c2=c('X4.SC14.04032020_dedup','X4.SC04.04032020_dedup')

all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp_drop = setdiff(all_trees[[file_name]]$sitka$tip.label, all_trees[[file_name]]$sctc$tip.label)
a3_data= fortify(drop.tip(all_trees[[file_name]]$sitka, tmp_drop))
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
a7_data= fortify(all_trees[[file_name]]$ML)
#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa, size=0.5)


# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)



######
file_name = 'Fig4b_dataset_photoconverted_M3_M4.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('M3.1.P9T.07092020_dedup','M3.2.P9T.07092020_dedup')
c2=c('M4.2.P9T.22052020_dedup','M4.1.P9T.22052020_dedup')

all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp_drop = setdiff(all_trees[[file_name]]$sitka$tip.label, all_trees[[file_name]]$sctc$tip.label)
a3_data= fortify(drop.tip(all_trees[[file_name]]$sitka, tmp_drop))
#a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data = a3_data
a4_data$label = gsub('_NA_', '', a4_data$label)
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
a7_data= fortify(all_trees[[file_name]]$ML)
#a1_data$x = (a1_data$x - min(a1_data$x)) / (max(a1_data$x) - min(a1_data$x))
a2_data$x = (a2_data$x - min(a2_data$x)) / (max(a2_data$x) - min(a2_data$x))
a3_data$x = (a3_data$x - min(a3_data$x)) / (max(a3_data$x) - min(a3_data$x))
a4_data$x = (a4_data$x - min(a4_data$x)) / (max(a4_data$x) - min(a4_data$x))
a5_data$x = (a5_data$x - min(a5_data$x)) / (max(a5_data$x) - min(a5_data$x))
a6_data$x = (a6_data$x - min(a6_data$x)) / (max(a6_data$x) - min(a6_data$x))

a2_data$x <- a2_data$x + max(a1_data$x) + 1
a3_data$x <- a3_data$x + max(a2_data$x) + 1
a4_data$x <- a4_data$x + max(a3_data$x) + 1
a5_data$x <- a5_data$x + max(a4_data$x) + 1
a6_data$x <- a6_data$x + max(a5_data$x) + 1
a7_data$x <- a7_data$x + max(a6_data$x) + 1

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

aa = bind_rows(a1_data, a2_data, a3_data, a4_data, a5_data,a6_data,a7_data) %>%
  filter(!is.na(label)) %>%
  filter(label %in% c(c1,c2))

#a1_data$branch.length = (a1_data$branch.length - min(a1_data$branch.length)) / (max(a1_data$branch.length) - min(a1_data$branch.length))
#a1_data$angle = (a1_data$angle - min(a1_data$angle)) / (max(a1_data$angle) - min(a1_data$angle))

p1 <- ggtree(a1_data, size=0.2) +
  geom_hilight(
    mapping=aes(subset = label %in% c(c1,c2),
                node = node,
                fill = as.factor(group)
    ), align="right"
  ) +
  labs(fill = "Group" )

pp = p1 +
  geom_tree(data = a2_data, size=0.2) +
  geom_tree(data = a3_data, size=0.2) +
  geom_tree(data = a4_data, size=0.2) +
  geom_tree(data = a5_data, size=0.2) +
  geom_tree(data = a6_data, size=0.2) +
  geom_tree(data = a7_data, size=0.2) +
  geom_line(aes(x, y, group=label, color=group), data=aa, size=0.5)

#######
min_subtree <-function(tree, tips){
  #tree = tmp_tree
  #tips = c1
  st = subtrees(tree)
  min_tree = tree
  min_nodes = length(min_tree$tip.label)
  for(i in st){
    diff_node = setdiff(tips, i$tip.label)
    #print(diff_node)
    if(length(diff_node)==0 & length(i$tip.label)<min_nodes){
      min_tree = i
      min_nodes = length(i$tip.label)
    }
  }
  return(min_tree)
}
all_tmp_dist = c()
file_name = 'Fig1b_dataset_P9T_20190331.txt'
c1=c('SC06.P9T.31032019_dedup', 'SC09.P9T.31032019_dedup')
c2=c('SC22.P9T.31032019_dedup','SC25.P9T.31032019_dedup', 'SC28.P9T.31032019_dedup')

for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
}
##
file_name = 'Fig1c_sub_dataset_P19BT_20181009.txt'
c1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
  mt = min_subtree(tmp_tree, c3)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
  mt = min_subtree(tmp_tree, c4)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
  mt = min_subtree(tmp_tree, c5)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))

}
##
file_name = 'Fig2c_dataset_P9T_20190409_2.txt'
c1 = c('SC26.P9T.09042019_dedup','SC31.P9T.09042019_dedup','SC34.P9T.09042019_dedup','SC28.P9T.09042019_dedup','SC29.P9T.09042019_dedup','SC33.P9T.09042019_dedup', 'SC30.P9T.09042019_dedup')
c2 = c('SC32.P9T.09042019_dedup','SC35.P9T.09042019_dedup','SC36.P9T.09042019_dedup')
c3 = c('SC27.P9T.09042019_dedup','SC37.P9T.09042019_dedup')
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
  mt = min_subtree(tmp_tree, c3)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))

}
##
file_name = 'Fig3a_dataset_P9T_20190409_1.txt'
c1=c('SC11.P9T.09042019_dedup','SC10.P9T.09042019_dedup','SC17.P9T.09042019_dedup')
c2=c('SC13.P9T.09042019_dedup','SC23.P9T.09042019_dedup')
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
}
##
file_name = 'Fig4a_dataset_P9T_4_20200304.txt'
c1=c('X4.SC18.04032020_dedup','X4.SC03.04032020_dedup','X4.SC15.04032020_dedup')
c2=c('X4.SC14.04032020_dedup','X4.SC04.04032020_dedup')
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
}
##
file_name = 'Fig4b_dataset_photoconverted_M3_M4.txt'
c1=c('M3.1.P9T.07092020_dedup','M3.2.P9T.07092020_dedup')
c2=c('M4.2.P9T.22052020_dedup','M4.1.P9T.22052020_dedup')
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  tmp_tree_dist = tip_dist(tmp_tree)
  mt = min_subtree(tmp_tree, c1)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f, file_name))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(mt$tip.label), f,file_name))
}

all_tmp_dist = as.data.frame(all_tmp_dist)
colnames(all_tmp_dist) = c('min_size', 'method', 'file')
all_tmp_dist$min_size = as.numeric(all_tmp_dist$min_size)
library(forcats)
all_tmp_dist %>%
  mutate(method=fct_reorder(method, min_size, median)) %>%
  mutate(group=ifelse(method=='sctc', 'sctc', 'other')) %>%
  ggplot(aes(x=method,  y=min_size, fill=group)) +
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_jitter(width = 0.2, size=0.5)+
  scale_fill_manual(values = c('sctc'='red', 'other'='gray')) +
  stat_compare_means(comparisons = list(c('sctc', 'MP')), label = 'p.format', method = 't.test')+
  theme_classic()+
  theme(legend.position = 'none')


######准确识别发生seismic模式的细胞 ####




# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL
)
gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)


###### 综合统计 #######



# Huh7细胞系热图
sctc = NG_list_with_score[['dataset_P9T_4_20200304_sub.txt']]
all_cell_cnv = sctc$orig.data$all_node_data
#Fn = sapply(rownames(sctc$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
#Fn = toupper(Fn)
#cna_meta = data.frame('Fn'=Fn, 'cell_name'=rownames(sctc$map_obj$cell_map_pos))
cna_meta = data.frame('clone'=sctc$map_obj$cell_map_pos$clone, 'cell_name'=rownames(sctc$map_obj$cell_map_pos))
rownames(cna_meta) = cna_meta$cell_name

col_meta = as.data.frame(t(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]])))
col_meta$V1 = as.numeric(col_meta$V1)
col_meta$V2 = as.numeric(col_meta$V2)
col_meta = col_meta %>% arrange(V1, V2)
col_chr =col_meta$V1

library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white",col="black",lwd=0.01),
                                               labels = 1:length(unique(col_chr)),
                                               labels_gp = gpar(cex = 0.6)))
cna_meta = cna_meta[, 'clone', drop=F]
cna_meta = cna_meta %>% arrange(clone)
cna_meta = cna_meta[grepl('clone', cna_meta$clone), , drop=F]

dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
cna_meta = cna_meta[lab,,drop=F]
# Fn_col = c('F1'='#6300FF', 'F2'='#C75127FF')
clone_col = setNames(pal_npg()(length(unique(cna_meta$clone))),unique(cna_meta$clone))
#clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')


left_anno = rowAnnotation(df=cna_meta,
                          col=list('clone' = clone_col),
                          width=unit(0.2, "mm"))
draw(left_anno)
cn_col = structure(c('blue','#91FF88', '#C6C6C6',  '#FFEB97',  '#FCCC00', '#ec9336', '#7d170e', 'darkred'),
                   names=c(0, 1,2,3,4,5,6, 7))
cn_col = structure(c('#2E6FAB','#6EA647', '#C6C6C6',  '#F8E596',  '#F4BA19', '#E37933', '#EC6A1A', '#B41D23'),
                   names=c(0, 1,2,3,4,5,6, 7))
width_in =  100/ 25.4
height_in = 70 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/cna热图.pdf', width=width_in, height=height_in)

ht = Heatmap(all_cell_cnv[rownames(cna_meta),rownames(col_meta)],
             cluster_rows = FALSE,cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = F,
             top_annotation = top_anno,
             left_annotation = left_anno,
             row_split = cna_meta$Fn,
             col = cn_col,
             #rect_gp = gpar(col = F),
             #rect_gp = gpar(col = 'white'),
             border=F,
             column_split = col_chr,
             column_gap = unit(0.5,  "mm"),
             heatmap_legend_param = list(
               title = "CNV",
               color_bar = "discrete",
               #at = c(0,1,2),
               direction = "horizontal",
               title_position = "leftcenter-rot"
               #legend_height = unit(3, "cm")
             ),
             row_title = NULL,
             column_title = NULL,
             use_raster=T,
             raster_quality = 5
)
draw(ht)
dev.off()

Heatmap(all_cell_cnv[rownames(cna_meta),rownames(col_meta)],
        cluster_rows = T,cluster_columns = FALSE,
        show_column_names = FALSE, show_row_names = F,
        top_annotation = top_anno,
        left_annotation = left_anno,
        #row_split = cna_meta$Fn,
        col = cn_col,
        #rect_gp = gpar(col = F),
        #rect_gp = gpar(col = 'white'),
        border=F,
        column_split = col_chr,
        column_gap = unit(0.5,  "mm"),
        heatmap_legend_param = list(
          title = "CNV",
          color_bar = "discrete",
          #at = c(0,1,2),
          direction = "horizontal",
          title_position = "leftcenter-rot"
          #legend_height = unit(3, "cm")
        ),
        row_title = NULL,
        column_title = NULL,
        #use_raster=T,
        #raster_quality = 20
)



















