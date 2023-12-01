#############
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')

library(Seurat)
library(ggsci)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
process_seurat <-function(obj, dim=30, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = row.names(obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,n.neighbors = 10L,dims = 1:dim,reduction = "pca", min.dist = 0.3, n.components=n.components)
  return(obj)
}
MyDimPlot <-function(obj, ...){
  umap = as.data.frame(obj@reductions$umap@cell.embeddings)
  DimPlot(obj, seed = seed, ...)+
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(), # 刻度不显示
      axis.text = element_blank(), # 刻度text不显示
      axis.title = element_blank()
    )+
    labs('title'='')+
    geom_segment(aes(x=min(umap$UMAP_1), y=min(umap$UMAP_2), xend=min(umap$UMAP_1)+2, yend=min(umap$UMAP_2)),
                 colour="black", size=0.5,arrow = arrow(length=unit(0.2,"cm")))+
    geom_segment(aes(x = min(umap$UMAP_1), y = min(umap$UMAP_2), xend = min(umap$UMAP_1), yend=min(umap$UMAP_2)+2),
                 colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"))) +
    annotate("text", x = min(umap$UMAP_1)+1, y = min(umap$UMAP_2)-0.25, label = "UMAP1",
             color="black",size = 4) +
    annotate("text", x = min(umap$UMAP_1)-0.25, y = min(umap$UMAP_2)+1, label = "UMAP2",
             color="black",size = 4, angle=90)
}

setwd('/Volumes/WX_extend/细胞轨迹推断/RNAseq/')

srt1 = readRDS('./output/786O_crispr_srt_filter.txt')
srt2 = readRDS('./output/A549_crispr_srt_filter.txt')
srt3 = readRDS('./output/Huh7_crispr_srt_filter.txt')
srt4 = readRDS('./output/U251_crispr_srt_filter.txt')

srt = merge(srt1, srt2)
srt = merge(srt, srt3)
srt = merge(srt, srt4)
srt = process_seurat(srt)

# DimPlot(srt, group.by = 'cell_line')
# srt = subset(srt, cell_line!='doublet')
cell_line_col = c('786O'='#2C86C8', 'A549'='#B7BC19', 'Huh7'='#E66327', 'U251'='#433691')
a = MyDimPlot(srt, group.by='cell_line', cols=cell_line_col, pt.size=2)+NoLegend()
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_umap.pdf',a,
       width=80, height=80, units='mm', dpi = 600, bg = 'transparent')

srt[['Map_gRNA']] = ifelse(!is.na(srt$feature_call), 'Map', 'Unmap')

map_gRNA_color = c('Map'='#9F9FFF', 'Unmap'='black')
MyDimPlot(srt, group.by='Map_gRNA', cols=map_gRNA_color, pt.size=1)


unmap_srt = subset(readRDS('./output/crispr_02.rds'), cell_line!='doublet')
unmap_srt@meta.data[is.na(unmap_srt@meta.data$num_features), 'num_features'] = 0
unmap_srt = subset(unmap_srt, num_features==0)
# 统计比例
srt2 = merge(srt, unmap_srt)
# srt@meta.data[is.na(srt@meta.data$num_features), 'num_features'] = 0
label = srt@meta.data %>%
  group_by(cell_line) %>%
  summarise(all = n()) %>%
  mutate(all=paste0('n=',all))
srt@meta.data %>%
  group_by(cell_line, num_features) %>%
  summarise(num = n()) %>%
  mutate(num_features = cut(num_features, breaks=c(-1, 0,1, Inf), labels=c('No guide','Single guide', 'Multiplet')))%>%
  ggplot(aes(x=cell_line, y=num, fill=num_features)) +
  geom_bar(stat='identity', position = 'fill',width=0.4)+
  geom_text(aes(x=cell_line, y=1.05, label=all), data=label, inherit.aes = F)+
  labs(y='Fraction of cells')+
  scale_fill_manual(values = c('No guide'='gray','Single guide'='#9F9FFF', 'Multiplet'='#9FD4FF'))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())


#####
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_tree/'
cell_line = c('Huh7', 'A549', '786O', 'U251')
Crispr_list = list()
for(prefix in cell_line){
  print(prefix)
  sctc = scTraceClass(scTrace_dir,
                      prefix,
                      layout = 'slanted',
                      rate = 0.9,
                      min_cell_num=NULL)
  Crispr_list[[prefix]] = sctc
}
saveRDS(Crispr_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_sctc.rds')
####

plist = list()
for(prefix in cell_line){
  print(prefix)
  tmp_class = Crispr_list[[prefix]]
  srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
  tmp_class$map_obj$cell_map_pos$gRNA = srt@meta.data[rownames(tmp_class$map_obj$cell_map_pos), 'feature_call']

  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

  #
  r_size = 0.04
  pie_data = cell_pos %>%
    group_by(clone, !!sym('gRNA')) %>%
    summarise(num=n()) %>%
    mutate(total=sum(num)) %>% as.data.frame()
  pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
  pie_data = pie_data %>% group_by(x, y, total)
  color_igv = pal_igv()(length(unique(pie_data$gRNA)))
  names(color_igv) = unique(pie_data$gRNA)
  df.grobs <- pie_data %>%
    do(subplots = ggplot(., aes(1, num, fill = !!sym('gRNA'))) +
         geom_col(position = "fill", alpha = 1, colour = "white") +
         coord_polar(theta = "y") +
         #scale_fill_igv()+
         scale_fill_manual(values = color_igv)+
         theme_void()+ guides(fill = F)) %>%
    mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                             x = x-r_size, y = y-r_size,
                                             xmax = x+r_size, ymax = y+r_size)))
  p = ggplot()
  p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                       color='#dddddd',
                       alpha=0.6,
                       linewidth=8,
                       #lwd=pmax(10/nchar(g$branches), 1),
                       linetype='solid',
                       lineend = "round",
                       linejoin='round',
  )
  p=p + df.grobs$subgrobs
  p=p+  theme_void() +
    labs(title=prefix)+
    guides(alpha=FALSE)+
    theme(plot.title = element_text(hjust=0.5))
  plist[[prefix]] = p
}
cowplot::plot_grid(plotlist = plist, ncol=2)

# Crispr 热图

scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_tree/'
prefix = 'U251'
sctc = Crispr_list[[prefix]]

cell_pos = sctc$map_obj$cell_map_pos
# load seurat rds
srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_filter_srt.RDS'))
srt[['clone_tree']] = 'other'
srt[['clone_tree']] = cell_pos[colnames(srt), 'clone']
library(Seurat)
DimPlot(srt, group.by = 'clone_tree', pt.size = 2) + ggsci::scale_color_igv()
DimPlot(srt,group.by = 'feature_call', pt.size = 2) +NoLegend()#+ ggsci::scale_color_igv()

###########
sctc$map_obj$cell_map_pos$gRNA = srt@meta.data[rownames(sctc$map_obj$cell_map_pos), 'feature_call']
# sctc$map_obj$cell_map_pos$gRNA[is.na(sctc$map_obj$cell_map_pos$gRNA )] = 'Virtual'

grna = na.omit(unique(sctc$map_obj$cell_map_pos$gRNA))
grna_num = paste0('gRNA-',1:length(grna))
names(grna_num) = grna
sctc$map_obj$cell_map_pos$gRNA = grna_num[sctc$map_obj$cell_map_pos$gRNA]
plot_cells(sctc, colorby='gRNA',shadow=TRUE, clone_line=TRUE,
           point.size=2, plot_virtual=F)+
  ggsci::scale_fill_igv()+NoLegend()


##########
gene.loc <-  read.table('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_crispr2/geneFile.txt')
colnames(gene.loc) = c('gene', 'chr', 'start', 'end')
rownames(gene.loc) = gene.loc$gene
gene.loc = gene.loc[colnames(sctc$orig.data$all_node_data), 'chr']
plot_cna_tree(sctc,
              colorby = 'gRNA',
              output = paste0(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/{prefix}_tree_heatmap.pdf')),
              column_split=sort(as.numeric(gsub('chr', '',gene.loc))),
              pt.size = 1,
              width=220,
              height=100,
              rel_widths=c(0.3,1),
              rel_heights = c(.08, 1),
              top_ann = T
)





