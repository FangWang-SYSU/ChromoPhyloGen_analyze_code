source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
library(ggnewscale)
library(ggtreeExtra)
library(ggpubr)
library(ComplexHeatmap)
dna_file = c('Fig3a_sub_dataset_P9T_20190409_1.txt',
             'Fig1c_sub_dataset_P19BT_20181009.txt',
             'Fig4a_sub_dataset_P9T_4_20200304.txt')
#

# scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sub/'
# subNG_list = list()
# for(prefix in dna_file){
#   print(prefix)
#   #min_cell_num=
#   #aa = load_scTrace(scTrace_dir, prefix)
#   sctc = scTraceClass(scTrace_dir,
#                       prefix,
#                       layout = 'slanted',
#                       rate = 0.8, min_cell_num = 1, min_cell_num2 = 1)
#
#   subNG_list[[prefix]] = sctc
# }
# saveRDS(subNG_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all.rds')
NG_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all.rds')

#### 计算score ######
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_sub2/'
NG_list_with_score = list()
for(prefix in dna_file){
  print(prefix)
  #prefix='uber.t3.seg.txt'
  sctc = NG_list[[prefix]]
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
  pdf(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG2/{prefix}_cna热图.pdf'), width=width_in, height=height_in)

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
               use_raster=F,
               raster_quality = 5
  )
  draw(ht)
  dev.off()



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
  score_data = score_data[rownames(cna_meta), ]

  p2 = p+geom_tippoint(size=1)
  p2 = p2 +
    new_scale_fill() +
    geom_fruit(data=score_data, geom=geom_tile,
               mapping=aes(y=nameid, x=x1,  fill=wgd),
               color = "grey",
               offset = 0.05,size = 0.05, width=2,
               axis.params=list(
                 axis="x",
                 text = "WGD",
                 text.angle = 45,
                 text.size = 2,
                 line.size = 0,
               vjust = 1,
               hjust= 1,
               inherit.aes=F,
             ))+
    scale_fill_npg()+
    scale_fill_manual(values = c('WGD0'='#DEDEDE', 'WGD1'='#FFE53B', 'WGD2'='#FF2525'),drop = FALSE)+
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
    #
    # new_scale_fill() +
    # geom_fruit(data=score_data, geom=geom_tile,
    #            mapping=aes(y=nameid, x=x1,  fill=pde),
    #            color = "grey",
    #            offset = 0.05,size = 0.05, width=2,
    #            axis.params=list(
    #              axis="x",
    #              text = "CEE",
    #              text.angle = 45,
    #              text.size = 2,
    #              line.size = 0,
    #              vjust = 1,
    #              hjust= 1,
    #              inherit.aes=F,
    #            ))+
    # scale_fill_gradientn(colours = c('white', '#8DC469'))+
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
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG2/hp_', prefix,'.pdf'),
         p2, width=160, height = 160, units='mm', dpi = 600)

  NG_list_with_score[[prefix]] = sctc
}
# saveRDS(NG_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/subNG_sctc_all_with_score.rds')

###
cn_col = structure(c('#2E6FAB','#6EA647', '#C6C6C6',  '#F8E596',  '#F4BA19', '#E37933', '#EC6A1A', '#B41D23'),
                   names=c(0, 1,2,3,4,5,6, 7))
sctc = NG_list[['Fig3a_sub_dataset_P9T_20190409_1.txt']]
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG2/不同方法比较划相同祖先划分准确性/Fig1b_dataset_P9T_20190331.txt.pdf',
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(.12, 1), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=220, height = 60, units='mm', dpi = 600)


# ---------
sctc = NG_list_with_score[['Fig3a_sub_dataset_P9T_20190409_1.txt']]
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
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG2/cna热图.pdf', width=width_in, height=height_in)

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
             use_raster=F,
             raster_quality = 5
)
draw(ht)
dev.off()

# right meta
rownames(cna_meta)





