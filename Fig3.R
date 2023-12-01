source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
# load pathway gene list
pathway = read.csv('/Users/lab/wangxin/MEDALT/SCTMS/data/pathway_cell2018.csv', header = T)
pathway = subset(pathway, type%in%c('TSG', 'OG'))
pathway = pathway[pathway$Gene%in%colnames(gene_cna_data), ]
pathway$name = paste0(pathway$pathway, '_', pathway$type)
pathway_list = list()
for(i in unique(pathway$name)){
  #if(length(pathway[pathway$name==i, 'Gene'])>=3){
  pathway_list[[i]] = pathway[pathway$name==i, 'Gene']
  #}

}
# 1.load data
dna_file = c(
  'huh7',
  'plc',
  'hep3b',
  'smc7721',
  'hepg2',
  'mhcc97h',
  'mhcc97l',
  'lm3')
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine/'

DNA_list = list()
for(prefix in dna_file){
  print(prefix)
  sctc = scTraceClass(scTrace_dir,
                      prefix,
                      layout = 'slanted',
                      rate = 0.98,
                      min_cell_num=NULL)
  DNA_list[[prefix]] = sctc
}
saveRDS(DNA_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')
####
DNA_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')

# 八个细胞系热图
file_path = '/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/'
cna_data = c()
cna_meta = c()
for(prefix in dna_file){
  print(prefix)
  sctc = DNA_list[[prefix]]
  tmp_data = read.table(paste0(file_path, prefix, '.txt_gene_data_all_node.csv'))
  rownames(tmp_data) = gsub('\\.', '-', rownames(tmp_data))
  tmp_data = tmp_data[sctc$orig.data$tree$tip.label, ]
  print(dim(tmp_data))
  cna_data = rbind(cna_data, tmp_data)
  #
  Fn = sapply(rownames(sctc$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
  Fn = toupper(Fn)
  tmp_meta = data.frame('Fn'=Fn, 'Cell_line'=prefix, 'cell_name'=rownames(sctc$map_obj$cell_map_pos))
  cna_meta = rbind(cna_meta, tmp_meta)
}


gene.loc <-  read.table('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_crispr2/geneFile.txt')
colnames(gene.loc) = c('gene', 'chr', 'start', 'end')
rownames(gene.loc) = gene.loc$gene
same_gene = intersect(colnames(cna_data),gene.loc$gene)
gene.loc = gene.loc[same_gene, 'chr']
cna_data = cna_data[, same_gene]
cna_data[is.na(cna_data)] = 2
library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white",col="black",lwd=0.01),
                                               labels = 1:24,
                                               labels_gp = gpar(cex = 0.6)))
cna_meta = cna_meta[rownames(cna_data), c(2,1)]
cna_meta = cna_meta %>% arrange(Cell_line, Fn)

Fn_col = c('F1'='#6300FF', 'F2'='#C75127FF')
cell_line_col = pal_igv()(8)
names(cell_line_col) = dna_file
left_anno = rowAnnotation(df=cna_meta,
                          col=list('Cell_line'= cell_line_col,
                                   'Fn' = Fn_col
                          )
)
draw(left_anno)
cn_col = structure(c('blue','#91FF88', '#C6C6C6',  '#FFEB97',  '#FCCC00', '#ec9336', '#7d170e', 'darkred'),
                   names=c(0, 1,2,3,4,5,6, 7))


width_in =  150/ 25.4
height_in = 110 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_DNA_cna热图.pdf', width=width_in, height=height_in)
ht = Heatmap(cna_data[rownames(cna_meta),],
             cluster_rows = FALSE,cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = F,
             top_annotation = top_anno,
             left_annotation = left_anno,
             row_split = cna_meta$Cell_line,
             col = cn_col,
             # rect_gp = gpar(col = "white"),
             column_split = gene.loc,
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
             column_title = NULL
)
draw(ht)
dev.off()

## RNA-seq细胞系umap
library(Seurat)
setwd('/Volumes/WX_extend/细胞轨迹推断/RNAseq/')
srt = readRDS('./output/cell_line_srt_with_SC3_cytotrace_new.RDS')
srt$cell_line = tolower(srt$cell_line)
srt = RunUMAP(srt, dims = 1:30, min.dist = 1)
a=DimPlot(srt, group.by = 'cell_line', cols = cell_line_col, pt.size = 0.01) + theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/RNAseq_umap_细胞系着色.pdf',a,
       width=60, height=40, units='mm', dpi = 600, bg = 'transparent')



### clone tree FN 着色###
plist = list()
for(prefix in dna_file){
  print(prefix)

  tmp_class = DNA_list[[prefix]]
  tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
  tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

  #
  r_size = 0.13
  pie_data = cell_pos %>%
    group_by(clone, !!sym('Fn')) %>%
    summarise(num=n()) %>%
    mutate(total=sum(num)) %>% as.data.frame()
  pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
  pie_data = pie_data %>% group_by(x, y, total)
  #color_igv = pal_igv()(length(unique(pie_data$Fn)))
  #names(color_igv) = unique(pie_data$Fn)
  df.grobs <- pie_data %>%
    do(subplots = ggplot(., aes(1, num, fill = !!sym('Fn'))) +
         geom_col(position = "fill", alpha = 1, colour = NA) +
         coord_polar(theta = "y") +
         #scale_fill_igv()+
         scale_fill_manual(values = Fn_col)+
         theme_void()+ guides(fill = F)) %>%
    mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                             x = -x-r_size, y = y-r_size,
                                             xmax = -x+r_size, ymax = y+r_size)))
  p = ggplot()
  p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                       color='#dddddd',
                       alpha=1,
                       linewidth=1,
                       #lwd=pmax(10/nchar(g$branches), 1),
                       linetype='solid',
                       lineend = "round",
                       linejoin='round',
  )
  p=p + df.grobs$subgrobs
  p=p+  theme_void() +
    labs(title=prefix)+
    guides(alpha=FALSE)+
    theme(plot.title = element_text(hjust=0.5, size=6))+
    scale_y_continuous(expand = expansion(mult = 0.2))+
    scale_x_continuous(expand = expansion(mult = 0.2))+
    coord_flip()+scale_x_reverse()
  plist[[prefix]] = p
}
a = cowplot::plot_grid(plotlist = plist, nrow=1,scale=0.9)

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_clone_tree_Fn着色.pdf',a,
       width=170, height=30, units='mm', dpi = 600, bg = 'transparent')


### clone tree clone 分裂速率着色###
DNA_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')
all_meta_rate = c()
for(f in dna_file){
  print(f)
  tmp_class = DNA_list[[f]]
  tmp_obj = tmp_class$orig.data
  meta = tmp_obj$cell_relation
  tmp_rate = predict_mitosis_rate(as.data.frame(meta))
  tmp_rate = (tmp_rate - min(tmp_rate)) / (max(tmp_rate) - min(tmp_rate))
  tmp_rate = scale(tmp_rate)
  tmp_rate = tmp_rate[,1]
  meta$pred_rate = tmp_rate[rownames(meta)]
  gene_level_cna_data = read.table(paste0(file_path, f, '.txt_gene_data_all_node.csv'))
  tmp_obj = add_cna_score(tmp_obj,pathway_list, gene_level_cna_data)
  meta$Fn = sapply(rownames(meta), function(x)strsplit(x, '_')[[1]][1])
  meta$source = f
  meta = meta %>% filter(Fn %in%c('f1', 'f2'))
  rownames(tmp_obj$cna_score) = gsub('\\.', '-', rownames(tmp_obj$cna_score))
  meta[, colnames(tmp_obj$cna_score)] = tmp_obj$cna_score[rownames(meta),]
  meta$Pseudotime_tree = (meta$Pseudotime_tree-min(meta$Pseudotime_tree)) / (max(meta$Pseudotime_tree) - min(meta$Pseudotime_tree))
  meta$clone = tmp_class$clone_graph_obj$cell_pos_reduction[rownames(meta), 'clone']
  all_meta_rate = rbind(all_meta_rate, meta)
}

plist = list()
for(prefix in dna_file){
  print(prefix)
  tmp_class = DNA_list[[prefix]]
  tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
  tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

  #
  cell_pos$predict_rate = all_meta_rate[rownames(cell_pos), 'pred_rate']
  x=boxplot(cell_pos$predict_rate)
  cell_pos = cell_pos[!cell_pos$predict_rate%in%x$out, ]
  p = ggplot()
  p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                       color='#dddddd',
                       alpha=1,
                       linewidth=2,
                       #lwd=pmax(10/nchar(g$branches), 1),
                       linetype='solid',
                       lineend = "round",
                       linejoin='round',
  )
  p = p+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone, color=predict_rate),
                                data=cell_pos[cell_pos$clone!='clone_0', ],
                                alpha=1,
                                lwd=0.4,
                                #color='gray',
                                lineend='round',linejoin='round',
                                #size=10,
  )

  p = p + geom_point(aes(x=x, y=y, fill=!!sym('predict_rate')),data=cell_pos,
                     shape=21, size=1,color='#00000000',inherit.aes = F)
  p=p +
    scale_color_gradientn(colors = c('blue', 'white', 'red'), values = c(0,0.5,1))+#, limits=c(0,0.5))+
    scale_fill_gradientn(colors = c('blue', 'white', 'red'), values = c(0,0.5,1))#, limits=c(0,0.5))
  p=p+  theme_void() +
    labs(title=prefix)+
    guides(alpha=FALSE)+
    #theme(plot.title = element_text(hjust=0.5), legend.position = 'none')+
    theme(plot.title = element_text(hjust=0.5, size=6), legend.position = 'none')+
    scale_y_continuous(expand = expansion(mult = 0.2))+
    scale_x_continuous(expand = expansion(mult = 0.2))+
    coord_flip()+scale_x_reverse()
  plist[[prefix]] = p
}
a = cowplot::plot_grid(plotlist = plist, nrow=1,scale=0.9)

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_clone_tree_rate着色.pdf',a,
       width=170, height=30, units='mm', dpi = 600, bg = 'transparent')

## clone align
# run clonealign.R
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/cell_line_srt_with_SC3_cytotrace_clonealign.RDS')
srt = RunUMAP(srt, dims = 1:30, min.dist = 1)
DimPlot(srt, group.by = 'clonealign')
file_names = c(
  'huh7',
  'plc',
  'smc7721',
  'hepg2',
  'mhcc97h',
  'mhcc97l'
)
plist = list()
for(prefix in file_names){
  print(prefix)
  tmp_class = DNA_list[[prefix]]
  tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
  tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]
  #
  tmp_srt = subset(srt, low_CL==prefix&clonealign!='unassigned')
  srt_align_clone = sapply(unique(tmp_srt$clonealign), function(x)gsub(paste0(prefix,'.txt_'), '', x))
  tmp_new_xy = c()
  for(i in names(srt_align_clone)){
    tmp_xy = clone_pos[srt_align_clone[i],c('x', 'y')]
    cell_xy = subset(tmp_srt, clonealign==i)@reductions$umap@cell.embeddings
    tmp_center = colMeans(cell_xy)
    # 移动
    delt = tmp_center - tmp_xy
    new_xy = t(t(cell_xy)-unlist(delt))
    # 收缩
    delt =  t(t(new_xy)-unlist(tmp_xy))
    new_xy = as.data.frame(new_xy - 0.975*delt)
    new_xy$clone = srt_align_clone[i]
    tmp_new_xy = rbind(tmp_new_xy, new_xy)
  }
  #
  p = ggplot()
  p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                       color='#dddddd',
                       alpha=1,
                       linewidth=1,
                       #lwd=pmax(10/nchar(g$branches), 1),
                       linetype='solid',
                       lineend = "round",
                       linejoin='round',
  )
  clone_pos2 = clone_pos[srt_align_clone, ]
  p=p+geom_text_repel(aes(x=x,y=y, label=label, color=label),
                       nudge_x=0.1,
                       size=2,
                       data=clone_pos2)
  p = p+geom_point(aes(x=UMAP_1, y=UMAP_2, color=clone), data=tmp_new_xy, size=0.01)
  p=p+
    scale_color_igv()+
    theme_void() +
    labs(title=prefix)+
    guides(alpha=FALSE)+
    #theme(plot.title = element_text(hjust=0.5), legend.position = 'none')+
    theme(plot.title = element_text(hjust=0.5, size=6), legend.position = 'none')+
    scale_y_continuous(expand = expansion(mult = 0.2))+
    scale_x_continuous(expand = expansion(mult = 0.2))+
    coord_flip()+scale_x_reverse()
  plist[[prefix]] = p
}
a = cowplot::plot_grid(plotlist = plist, nrow=1)

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_clone_align.pdf',a,
       width=127.5, height=30, units='mm', dpi = 600, bg = 'transparent')





## 速率与增殖相关性
# 分裂速率与增殖评分相关
DNA_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')
all_meta_rate = c()
for(f in dna_file){
  print(f)
  tmp_class = DNA_list[[f]]
  tmp_obj = tmp_class$orig.data
  meta = tmp_obj$cell_relation
  tmp_rate = predict_mitosis_rate(as.data.frame(meta))
  # tmp_rate = (tmp_rate - mean(tmp_rate)) / sd(tmp_rate)
  tmp_rate = (tmp_rate - min(tmp_rate)) / (max(tmp_rate) - min(tmp_rate))
  #tmp_rate = scale(tmp_rate)
  #tmp_rate = tmp_rate[,1]
  meta$pred_rate = tmp_rate[rownames(meta)]
  gene_level_cna_data = read.table(paste0(file_path, f, '.txt_gene_data_all_node.csv'))
  tmp_obj = add_cna_score(tmp_obj,pathway_list, gene_level_cna_data)
  meta$Fn = sapply(rownames(meta), function(x)strsplit(x, '_')[[1]][1])
  meta$source = f
  meta = meta %>% filter(Fn %in%c('f1', 'f2'))
  rownames(tmp_obj$cna_score) = gsub('\\.', '-', rownames(tmp_obj$cna_score))
  meta[, colnames(tmp_obj$cna_score)] = tmp_obj$cna_score[rownames(meta),]
  meta$CellCycle_OG = (meta$CellCycle_OG-min(meta$CellCycle_OG)) / (max(meta$CellCycle_OG) - min(meta$CellCycle_OG))

  meta$Pseudotime_tree = (meta$Pseudotime_tree-min(meta$Pseudotime_tree)) / (max(meta$Pseudotime_tree) - min(meta$Pseudotime_tree))
  meta$clone = tmp_class$clone_graph_obj$cell_pos_reduction[rownames(meta), 'clone']
  all_meta_rate = rbind(all_meta_rate, meta)
}
all_meta_rate_old=all_meta_rate
#all_meta_rate$CellCycle_OG_cut = cut(all_meta_rate$CellCycle_OG, c(-Inf,0.25,0.5,Inf))
#
#all_meta_rate %>%
#  #filter(DNAcell_num>=25)%>%
#  ggplot(aes(x=CellCycle_OG_cut, y=pred_rate))+
#  geom_boxplot()
#
#library(ggpubr)
#
all_meta_rate = all_meta_rate_old
all_meta_rate = all_meta_rate %>%
  group_by(source, clone) %>%
  summarise(pred_rate = mean(pred_rate), CellCycle_OG = mean(CellCycle_OG)) %>%
  left_join(all_meta_rate %>%
              group_by(source, clone) %>%
              summarise(DNAcell_num=n()) )

all_meta_rate$CellCycle_OG_cut = cut(all_meta_rate$CellCycle_OG, c(-Inf,1/3,2/3,Inf), labels = c('Low', 'Neutal','High'))
line_data = all_meta_rate %>%
  group_by(CellCycle_OG_cut) %>%
  dplyr::summarise(pred_rate_m = median(pred_rate)) %>% as.data.frame()
all_meta_rate = merge(all_meta_rate, line_data)

a=all_meta_rate %>%
  #filter(DNAcell_num>=25)%>%
  #filter(source!='mhcc97h') %>%
  ggplot(aes(x=CellCycle_OG_cut, y=pred_rate))+
  geom_boxplot(aes(color=CellCycle_OG_cut), outlier.shape = NA, size=0.2, width=0.5) +
  geom_jitter(aes(fill=source),width = 0.2, size=1, shape=21, stroke=0)+
  geom_smooth(data=line_data, mapping=aes(x=as.numeric(CellCycle_OG_cut), y=pred_rate_m),
              method = "lm", se=F,  formula = y~x, size=0.2, inherit.aes = F, color='black')+
  stat_cor(data=line_data, mapping=aes(x=as.numeric(CellCycle_OG_cut), y=pred_rate_m),label.x = 0.3,
           label.y = 1.1, method='pearson', size=2, color='black', inherit.aes = F)+
  scale_color_manual(values = c('Low'='#18499E', 'Neutal'='#C7C6C6', 'High'='#E39039'))+
  scale_fill_manual(values = cell_line_col)+
  #geom_dotplot(binaxis='y', stackdir='center', stackratio=2, dotsize=0.3, binwidth = 0.03, stroke=0.01)+
  theme_classic()+
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))+
  scale_y_continuous(expand = expansion(mult = 0.2))

a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_clone_速率与增殖相关性.pdf',a,
       width=90, height=45, units='mm', dpi = 600, bg = 'transparent')

a = all_meta_rate %>%
  #filter(DNAcell_num>=25)%>%
  ggplot(aes(x=CellCycle_OG, y=pred_rate))+
  geom_point(aes(color=source), size=1)+
  scale_color_igv()+
  stat_cor(size=1)+
  labs(x='CellCycle_OG', y='Predicted velocity')+
  stat_smooth(method = 'lm',formula = y ~ x, se=F,color='black', linetype='dashed', linewidth=0.2)+
  #lims(y=c(-2,4))+
  #facet_wrap(~source, scales = 'free')+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()+
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/8细胞系_clone_速率与增殖相关性(2).pdf',a,
       width=60, height=40, units='mm', dpi = 600, bg = 'transparent')



# RNA 干性得分与 增殖，pseudotime，分裂速率相关性
all_meta_rate = all_meta_rate_old %>%
  group_by(source, clone) %>%
  summarise(pred_rate = mean(pred_rate))

srt[['clone']]=sapply(srt$clonealign, function(x){
  y=strsplit(x, '_')[[1]][2:3]
  return(paste0(y, collapse = '_'))
})
srt[['source']] = sapply(srt$clonealign, function(x){
  y=strsplit(x, '\\.')[[1]][1]
  return(y)
})

# cytotrace
#library(CytoTRACE)
#srt_count = AverageExpression(srt, group.by = 'clonealign', slot = 'counts')$RNA
#new_trace = CytoTRACE::CytoTRACE(srt_count)
#srt$cytotrace_score = new_trace$CytoTRACE[srt$clonealign]
srt = CellCycleScoring(srt, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)

srt_meta = srt@meta.data[, c('clone', 'source', 'cytotrace_score', 'G2M.Score')]

srt_meta = srt_meta %>%
  group_by(clone, source) %>%
  summarise(cytotrace_score=mean(cytotrace_score), G2M.Score=mean(G2M.Score))

all_meta_rate2 = merge(all_meta_rate, as.data.frame(srt_meta))
# for(i in unique(all_meta_rate2$source)){
#   x = all_meta_rate2[all_meta_rate2$source==i, 'cytotrace_score']
#   x = (x-min(x)) / (max(x)-min(x))
#   all_meta_rate2[all_meta_rate2$source==i, 'cytotrace_score'] = x#scale(x)[,1]
# }
# 矫正mhcc97h
all_meta_rate2[all_meta_rate2$source=='mhcc97h', 'cytotrace_score'] = all_meta_rate2[all_meta_rate2$source=='mhcc97h', 'cytotrace_score']/5
a=all_meta_rate2 %>%
  #filter(source!='mhcc97h') %>%
  ggplot(aes(x=pred_rate, y=cytotrace_score))+
  geom_point(aes(fill=source),size=2, shape=21,color='black', stroke=0.1)+
  scale_fill_manual(values = cell_line_col)+
  #facet_wrap(~source, scales = 'free', nrow=1)+
  stat_cor(size=2)+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', size=0.2)+
  theme_classic()+
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/RNA_速率与干性相关.pdf',a,
       width=60, height=40, units='mm', dpi = 600, bg = 'transparent')

all_meta_rate2[all_meta_rate2$source=='mhcc97h', 'G2M.Score'] = all_meta_rate2[all_meta_rate2$source=='mhcc97h', 'G2M.Score']-0.3
a=all_meta_rate2 %>%
  #filter(source!='mhcc97h') %>%
  ggplot(aes(x=pred_rate, y=G2M.Score))+
  geom_point(aes(fill=source),size=2, shape=21,color='black', stroke=0.1)+
  scale_fill_manual(values = cell_line_col)+
  #facet_wrap(~source, scales = 'free', nrow=1)+
  stat_cor(size=2)+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', size=0.2)+
  theme_classic()+
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/RNA_速率与G2M相关.pdf',a,
       width=60, height=40, units='mm', dpi = 600, bg = 'transparent')




# RNA 干性得分与 增殖，pseudotime，分裂速率相关性
srt[['clone']]=sapply(srt$clonealign, function(x){
  y=strsplit(x, '_')[[1]][2:3]
  return(paste0(y, collapse = '_'))
})
srt[['source']] = sapply(srt$clonealign, function(x){
  y=strsplit(x, '\\.')[[1]][1]
  return(y)
})
srt = CellCycleScoring(srt, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)

srt_meta = srt@meta.data[, c('clone', 'source', 'cytotrace_score', 'G2M.Score')]
# 矫正mhcc97倍性
# srt_meta[srt_meta$source=='mhcc97h', '']
srt_meta = srt_meta[srt_meta$source!='unassigned', ]
srt_meta = srt_meta %>%
  group_by(clone, source) %>%
  summarise(cytotrace_score=mean(cytotrace_score), G2M.Score=mean(G2M.Score)) %>%
  left_join( srt_meta %>%
               group_by(clone, source) %>%
               summarise(RNAseq_num=n()))
srt_meta=srt_meta %>%
  left_join(srt_meta%>%
              group_by(source) %>%
              summarise(cytotrace_score1=tmp_func(cytotrace_score), clone=clone))
tmp_func<-function(x){
  (x-min(x)) / (max(x)-min(x))
}

all_meta_rate2 = merge(all_meta_rate, as.data.frame(srt_meta))

all_meta_rate2

all_meta_rate2 %>%
  #filter(RNAseq_num>5)%>%
  #filter(Fn=='f2') %>%
  filter(!is.na(cytotrace_score)) %>%
  group_by(source, clone,Fn) %>%
  summarise(cytotrace_score = mean(cytotrace_score1), pred_rate=mean(pred_rate)) %>%
  ggplot(aes(x=pred_rate, y=cytotrace_score))+
  geom_point(aes(fill=source),size=4, shape=21,color='black')+
  scale_fill_igv()+
  #facet_wrap(~source, scales = 'free')+
  stat_cor()+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F)+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()


all_meta_rate2 %>%
  filter(RNAseq_num>1)%>%
  #filter(Fn=='f2') %>%
  filter(!is.na(cytotrace_score)) %>%
  group_by(source, clone,Fn) %>%
  summarise(G2M.Score = mean(G2M.Score), pred_rate=mean(pred_rate)) %>%
  ggplot(aes(x=pred_rate, y=G2M.Score))+
  geom_point(aes(fill=source),size=4, shape=21,color='black')+
  scale_fill_igv()+
  #facet_wrap(~source, scales = 'free')+
  stat_cor()+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F)+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()

# 矫正mhcc97倍性
all_meta_rate2[all_meta_rate2$source=='mhcc97h', "cytotrace_score"]= all_meta_rate2[all_meta_rate2$source=='mhcc97h', "cytotrace_score"]/2
#all_meta_rate2 = all_meta_rate2 %>%
#  filter(source!='mhcc97h')
all_meta_rate2$cytotrace_score_cut = cut(all_meta_rate2$cytotrace_score, c(-Inf,1/3,2/3,Inf), labels = c('Low', 'Neutal','High'))
line_data = all_meta_rate2 %>% group_by(cytotrace_score_cut) %>% dplyr::summarise(pred_rate_m2 = median(pred_rate)) %>% as.data.frame()
all_meta_rate2 = merge(all_meta_rate2, line_data)

all_meta_rate2 %>%
  #filter(source!='mhcc97h') %>%
  #filter(DNAcell_num>=25)%>%
  ggplot(aes(x=cytotrace_score_cut, y=pred_rate))+
  geom_boxplot(aes(color=cytotrace_score_cut), outlier.shape = NA, size=0.2, width=0.5) +
  geom_jitter(aes(fill=source),width = 0.2, size=1, shape=21, stroke=0)+
  geom_smooth(data=line_data, mapping=aes(x=as.numeric(cytotrace_score_cut), y=pred_rate_m2),
              method = "lm", se=F,  formula = y~x, size=0.2, inherit.aes = F, color='black')+
  stat_cor(data=line_data, mapping=aes(x=as.numeric(cytotrace_score_cut), y=pred_rate_m2),label.x = 0.3,
           label.y = 1.1, method='pearson', size=2, color='black', inherit.aes = F)+
  scale_color_manual(values = c('Low'='#18499E', 'Neutal'='#C7C6C6', 'High'='#E39039'))+
  scale_fill_manual(values = cell_line_col)+
  #geom_dotplot(binaxis='y', stackdir='center', stackratio=2, dotsize=0.3, binwidth = 0.03, stroke=0.01)+
  theme_classic()+
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))+
  scale_y_continuous(expand = expansion(mult = 0.2))









