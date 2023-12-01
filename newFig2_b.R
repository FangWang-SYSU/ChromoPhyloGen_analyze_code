tmp_func <- function(tree, data, vec){
  pnodes = c()
  root_num = getRoot(tree)
  for(i in 1:nrow(data)){
    if(data$label[i] %in% vec){
      node_num = data$node[i]
      pnodes = c(pnodes, nodepath(tree, from=root_num,to=node_num))
    }
  }
  pnodes = unique(pnodes)
  return(pnodes)
}
cn_col = structure(c('#2E6FAB','#6EA647', 'white',  '#F8E596',  '#F4BA19', '#E37933', '#EC6A1A', '#B41D23'),
                   names=c(0, 1,2,3,4,5,6, 7))
########
file_name = 'Fig1c_sub_dataset_P19BT_20181009.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c()#c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c()#c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c()#c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c()#c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')

l1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup',
       'SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup',
       'SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup',
       'SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')

#all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
#tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
a3_data= fortify(all_trees[[file_name]]$sitka)
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
tmp =all_trees[[file_name]]$NJ
tmp$edge.length = tmp$edge.length+1
a5_data= fortify(all_trees[[file_name]]$NJ)
a6_data= fortify(all_trees[[file_name]]$MP)
tmp =all_trees[[file_name]]$ML
tmp$edge.length = tmp$edge.length+1
a7_data= fortify(tmp)

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, l1), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, l1), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, l1), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, l1), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, l1), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, l1), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, l1), 'a', 'b')})

a1_data$text = sapply(a1_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a2_data$text = sapply(a2_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a3_data$text = sapply(a3_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a4_data$text = sapply(a4_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a5_data$text = sapply(a5_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a6_data$text = sapply(a6_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})
a7_data$text = sapply(a7_data$label,function(x){if(x %in% c5){strsplit(x, '\\.')[[1]][1]}else if(x %in% c2){''}else{''}})

a1=ggtree(a1_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('b'='black', 'a'='gray'))+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 50, units='mm', dpi = 600)


dist_df = data.frame(dist=c(1,8,3, 7, 2,1, 2),
           method=c('SCTC', 'MEDALT', 'Sitka', 'MEDICC2', 'NJ', 'MP', 'ML'))
library(forcats)
a=dist_df %>%
  #mutate(method = fct_reorder(method, dist)) %>%
  mutate(method = factor(method, c('SCTC', 'MP', 'ML', 'NJ', 'Sitka', 'MEDICC2', 'MEDALT'))) %>%
  ggplot(aes(x=method, y=dist))+
  geom_point(size=2)+
  geom_text(aes(label=dist), color='white', size=1)+
  labs(y='Reciprocal cells distance')+
  scale_y_continuous(expand = expansion(0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_dist.pdf'),
       a,
       width=24, height = 26, units='mm', dpi = 600)


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
               border_gp = gpar(col = "gray", lwd=0.05),
               height = nrow(all_cell_cnv)*unit(1, "mm"), # 高度
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.3))),
               col = cn_col,
               column_split = as.numeric(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1])),
               row_split = factor(rownames(all_cell_cnv),rownames(all_cell_cnv)),
               column_gap = unit(0, "mm"),
               row_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               row_title = NULL,
               use_raster = F
)
width_in = 90 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

gb2 <- grid.grabExpr(draw(ht2))
space <- ggdraw()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'.pdf'),
       plot_grid(gb2,
                 plot_grid(space,pp, rel_heights = c(0.01, 0.9), ncol = 1),
                 rel_widths = c(.4, 0.6)),
       width=180, height = 100, units='mm', dpi = 600)


sctc = NG_list[[file_name]]
sctc$map_obj$cell_map_pos[, c('cell_name', 'rearrange_score', 'pde')] %>%
  reshape2::melt(id='cell_name') %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  ggplot(aes(x=value,y=cell_name, fill=value))+
  geom_bar(stat='identity')+
  facet_wrap(~variable, scales = 'free_x')

b1 = sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=rearrange_score,y=cell_name))+
  geom_segment(aes(x=0, xend=rearrange_score, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=rearrange_score,size=rearrange_score),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#306DAA'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('Chromothripsis')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
b1
b2=sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=pde,y=cell_name, fill=pde))+
  geom_segment(aes(x=0, xend=pde, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=pde,size=pde),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#8DC469'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('CEE')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

bb = cowplot::plot_grid(b1,b2,nrow=1, align = 'h')
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/score_',file_name,'.pdf'),
       bb,
       width=110, height = 40, units='mm', dpi = 600)

##################
########
file_name = 'Fig3a_dataset_P9T_20190409_1.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('SC11.P9T.09042019_dedup','SC10.P9T.09042019_dedup','SC17.P9T.09042019_dedup')
c2=c('SC13.P9T.09042019_dedup','SC23.P9T.09042019_dedup')
#all_tmp_dist = c()
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

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, c(c1,c2)), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, c(c1,c2)), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, c(c1,c2)), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, c(c1,c2)), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, c(c1,c2)), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, c(c1,c2)), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, c(c1,c2)), 'a', 'b')})

a1=ggtree(a1_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')

cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)

pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1),
       width=210, height = 40, units='mm', dpi = 600)

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
               border_gp = gpar(col = "black", width=0.1),
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               column_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               use_raster = F
)
width_in = 120 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

sctc = NG_list[[file_name]]
sctc$map_obj$cell_map_pos[, c('cell_name', 'rearrange_score', 'pde')] %>%
  reshape2::melt(id='cell_name') %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  ggplot(aes(x=value,y=cell_name, fill=value))+
  geom_bar(stat='identity')+
  facet_wrap(~variable, scales = 'free_x')

b1 = sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=rearrange_score,y=cell_name))+
  geom_segment(aes(x=0, xend=rearrange_score, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=rearrange_score,size=rearrange_score),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#306DAA'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('Chromothripsis')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
b1
b2=sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=pde,y=cell_name, fill=pde))+
  geom_segment(aes(x=0, xend=pde, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=pde,size=pde),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#8DC469'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('CEE')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

bb = cowplot::plot_grid(b1,b2,nrow=1, align = 'h')
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/score_',file_name,'.pdf'),
       bb,
       width=110, height = 40, units='mm', dpi = 600)



########
file_name = 'Fig4a_dataset_P9T_4_20200304.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('X4.SC18.04032020_dedup','X4.SC03.04032020_dedup','X4.SC15.04032020_dedup')
c2=c('X4.SC14.04032020_dedup','X4.SC04.04032020_dedup')
#all_tmp_dist = c()
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

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, c(c1,c2)), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, c(c1,c2)), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, c(c1,c2)), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, c(c1,c2)), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, c(c1,c2)), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, c(c1,c2)), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, c(c1,c2)), 'a', 'b')})

a1=ggtree(a1_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.4, mapping=aes(linetype=line, color=line))+geom_hilight(mapping=aes(subset = label %in% c(c1,c2),node = node,fill = as.factor(group)), align="right")+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')

cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)

pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1),
       width=210, height = 40, units='mm', dpi = 600)

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
               border_gp = gpar(col = "black", width=0.1),
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.6))),
               col = cn_col,
               column_split = sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1]),
               column_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               use_raster = F
)
width_in = 120 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

sctc = NG_list[[file_name]]
sctc$map_obj$cell_map_pos[, c('cell_name', 'rearrange_score', 'pde')] %>%
  reshape2::melt(id='cell_name') %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  ggplot(aes(x=value,y=cell_name, fill=value))+
  geom_bar(stat='identity')+
  facet_wrap(~variable, scales = 'free_x')

b1 = sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=rearrange_score,y=cell_name))+
  geom_segment(aes(x=0, xend=rearrange_score, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=rearrange_score,size=rearrange_score),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#306DAA'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('Chromothripsis')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
b1
b2=sctc$map_obj$cell_map_pos %>%
  filter(cell_name %in% lab) %>%
  mutate(cell_name=factor(cell_name, rev(lab))) %>%
  #mutate(cell_name = gsub('.P19BT.09102018_dedup','', cell_name)) %>%
  ggplot(aes(x=pde,y=cell_name, fill=pde))+
  geom_segment(aes(x=0, xend=pde, y=cell_name, yend=cell_name), linewidth=0.1)+
  geom_point(aes(fill=pde,size=pde),shape=21, stroke=0.1)+
  scale_fill_gradientn(colours = c('white', '#8DC469'))+
  scale_size_continuous(range=c(1,2))+
  ggtitle('CEE')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

bb = cowplot::plot_grid(b1,b2,nrow=1, align = 'h')
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/score_',file_name,'.pdf'),
       bb,
       width=110, height = 40, units='mm', dpi = 600)








