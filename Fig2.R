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
library(glue)
library(ggsci)
library(igraph)
library(ggstream)

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

get_clone <- function(tree, groups){
  st = subtrees(tree)
  res = c()
  for(i in st){
    res = rbind(res, c(length(i$tip.label), vegan::diversity(table(groups[i$tip.label]), index='shannon')))
  }
  return(res)
}
get_Label_dist <-   function(tree, labels){
  tips_num = 1:length(tree$tip.label)
  names(tips_num) = tree$tip.label
  gRNA_label_names = names(labels)
  tmp_res = c()
  for(grna in unique(labels)){
    print(grna)
    tmp_grna_cells = gRNA_label_names[labels==grna]
    tmp_path = c()
    for(i in 1:(length(tmp_grna_cells)-1)){
      for(j in (i+1):length(tmp_grna_cells)){
        pt = nodepath(tree, from = tips_num[tmp_grna_cells[i]],
                      to = tips_num[tmp_grna_cells[j]])
        tmp_path = c(tmp_path, length(pt))
      }
    }
    tmp_res = rbind(tmp_res, c(grna,mean(tmp_path)))
  }
  return(tmp_res)
}
get_clone_func <- function(tree, min_k=5){
  all_path = nodepath(tree)
  max_path = max(sapply(all_path, length))
  for(i in 1:max_path){
    tmp_pos = c()
    for(j in all_path){
      if(length(j)>=i){
        tmp_pos = c(tmp_pos, j[i])
      }else{
        tmp_pos = c(tmp_pos, j[length(j)])
      }
    }
    tmp_pos = unique(tmp_pos)
    if(length(tmp_pos)>=min_k){
      break
    }
  }

  rt_data = fortify(tree) %>% as.data.frame()

  true_clone = c()
  ki=1
  for(i in tmp_pos){
    if(i<=length(tree$tip.label)){
      lvs = rt_data[rt_data$node==i,'label']
    }else{
      lvs = extract.clade(tree,i)$tip.label
    }
    true_clone = rbind(true_clone, data.frame(clone=ki, name=lvs))
    ki = ki+1
  }
  true_clone = as.data.frame(true_clone)
  rownames(true_clone) = true_clone$name
  return(true_clone)
}

get_clone <- function(tree, label){
  res = c()
  for(i in 2:length(tree$tip.label)){
    tmp_cl = get_clone_func(tree, min_k=i)
    tmp_ix = pdfCluster::adj.rand.index(label, tmp_cl[names(label), 'clone'])
    res = rbind(res, c(i, tmp_ix))
  }
  return(res)
}


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
srt <- RunUMAP(srt,n.neighbors = 20,dims = 1:20,reduction = "pca",
               min.dist = 2, n.components=2)

cell_line_col = c('786O'='#2C86C8', 'A549'='#B7BC19', 'Huh7'='#E66327', 'U251'='#009944')
feature_call_color = pal_igv()(length(unique(srt$feature_call)))
names(feature_call_color) = unique(srt$feature_call)

a = MyDimPlot(srt, group.by='cell_line', cols=cell_line_col, pt.size=0.1)+NoLegend()
a = a$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=cell_line))+
  geom_point(shape=21, size=1.5,stroke=0.1, color='black')+
  scale_fill_manual(values = cell_line_col)+
  theme_void()+
  theme(legend.position = 'none')
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_umap_CL着色.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

# gRNA着色
a = MyDimPlot(srt, group.by='feature_call', cols=cell_line_col, pt.size=0.1)+NoLegend()
a = a$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=feature_call))+
  geom_point(shape=21, size=1.5,stroke=0.1, color='black')+
  scale_fill_manual(values = feature_call_color)+
  #scale_fill_igv()+
  theme_void()+
  theme(legend.position = 'none')
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_umap_gRNA着色.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

##### 加载tree
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
## DNA 拷贝数谱
cell_line = c('Huh7', 'A549', '786O', 'U251')
Crispr_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_sctc.rds')
same_gene = c()
for(prefix in cell_line){
  sctc = Crispr_list[[prefix]]
  tmp_gene = colnames(sctc$orig.data$all_node_data)
  same_gene = c(same_gene, tmp_gene)
}
same_gene = unique(same_gene)
cna_data = c()
cna_meta = c()
for(prefix in cell_line){
  sctc = Crispr_list[[prefix]]
  tmp_data = sctc$orig.data$all_node_data[sctc$orig.data$tree$tip.label, ]
  tmp_gene = setdiff(same_gene, colnames(tmp_data))
  tmp_data2 = matrix(2, nrow(tmp_data), ncol = length(tmp_gene))
  colnames(tmp_data2) = tmp_gene
  tmp_data = cbind(tmp_data, tmp_data2)
  print(dim(tmp_data))
  cna_data = rbind(cna_data, tmp_data[, same_gene])
  #
  tmp_meta = srt@meta.data[sctc$orig.data$tree$tip.label, c('feature_call','cell_line')]
  cna_meta = rbind(cna_meta, tmp_meta)
}

#
gene.loc <-  read.table('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_crispr2/geneFile.txt')
colnames(gene.loc) = c('gene', 'chr', 'start', 'end')
rownames(gene.loc) = gene.loc$gene
gene.loc = gene.loc[same_gene, 'chr']
library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white",col="black",lwd=0.01),
                                               labels = 1:22,
                                               labels_gp = gpar(cex = 0.6)))
cna_meta = cna_meta[, c(2,1)]
cna_meta = cna_meta %>% arrange(cell_line, feature_call)


lgd = Legend(labels =names(feature_call_color),
             title = "gRNA",
             legend_gp = gpar(fill = feature_call_color))
width_in =  180/ 25.4
height_in = 150 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/DNA_cna热图_legend.pdf', width=width_in, height=height_in)
draw(lgd)
dev.off()

left_anno = rowAnnotation(df=cna_meta,
                          col=list('cell_line'= cell_line_col,
                                   'feature_call' =feature_call_color
                          ),
                          show_legend=F
)
#draw(left_anno)
cn_col = structure(c('blue','#91FF88', '#C6C6C6',  '#FFEB97',  '#FCCC00', '#ec9336', '#7d170e', 'darkred'),
                   names=c(0, 1,2,3,4,5,6, 7))

width_in =  110/ 25.4
height_in = 150 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/DNA_cna热图.pdf', width=width_in, height=height_in)
ht = Heatmap(cna_data[rownames(cna_meta),],
             cluster_rows = FALSE,cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = F,
             top_annotation = top_anno,
             left_annotation = left_anno,
             row_split = cna_meta$cell_line,
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
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# clone tree
plist = list()
cell_line = c('786O','A549', 'Huh7',   'U251')
for(prefix in cell_line){
  print(prefix)
  tmp_class = Crispr_list[[prefix]]
  tmp_srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{prefix}_crispr_srt_filter.txt'))
  tmp_class$map_obj$cell_map_pos$gRNA = tmp_srt@meta.data[rownames(tmp_class$map_obj$cell_map_pos), 'feature_call']

  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

  #
  r_size = 0.09
  pie_data = cell_pos %>%
    group_by(clone, !!sym('gRNA')) %>%
    summarise(num=n()) %>%
    mutate(total=sum(num)) %>% as.data.frame()
  pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
  pie_data = pie_data %>% group_by(x, y, total)
  #color_igv = pal_igv()(length(unique(pie_data$gRNA)))
  #names(color_igv) = unique(pie_data$gRNA)
  df.grobs <- pie_data %>%
    do(subplots = ggplot(., aes(1, num, fill = !!sym('gRNA'))) +
         geom_col(position = "fill", alpha = 1, colour = NA) +
         coord_polar(theta = "y") +
         #scale_fill_igv()+
         scale_fill_manual(values = feature_call_color)+
         theme_void()+ guides(fill = F)) %>%
    mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                             x = x-r_size, y = y-r_size,
                                             xmax = x+r_size, ymax = y+r_size)))
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
    theme(plot.title = element_text(hjust=0.5))+
    scale_y_continuous(expand = expansion(mult = 0.2))+
    scale_x_continuous(expand = expansion(mult = 0.1))
  plist[[prefix]] = p
}
a = cowplot::plot_grid(plotlist = plist, ncol=1, align = 'h')

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_clone_tree.pdf',a,
       width=50, height=150, units='mm', dpi = 600, bg = 'transparent')



# 比较ARI
medalt_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medalt/Crispr/'
sitka_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sitka/Crispr/'
medicc_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medicc2/Crispr'
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_tree/'


tree_rds = list()
for(f in c('Huh7','A549','786O','U251')){
  #f= 'Huh7'
  print(f)
  # 加载scTrace结果
  scTrace_res = load_scTrace(scTrace_dir, f)
  scTrace_tree = scTrace_res$tree
  # 创建传统进化树
  tmp_data = scTrace_res$all_node_data[scTrace_res$tree$tip.label, ]
  #edu_dist = dist(tmp_data, method='euclidean')
  #pyd = as.phyDat(as.data.frame(t(tmp_data)), type='USER', levels=min(tmp_data):max(tmp_data))
  # NJ
  # NJ_tree<-NJ(edu_dist)
  mm = dist(tmp_data, method = "maximum")
  NJ_tree<-NJ(mm)
  # random_tree
  random_tree = rtree(nrow(tmp_data), tip.label = rownames(tmp_data))
  # MP
  pyd=as.phyDat(as.matrix(mm), type='USER', levels=0:max(tmp_data))
  MP_tree = optim.parsimony(random_tree, pyd)

  # ML
  ML_tree = optim.pml(pml(random_tree, pyd))$tree
  # sitka tree
  sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
  sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))

  # medalt
  medalt_tree = sitka_newick_to_phylo(paste0(medalt_dir, f, '_medalt_tree.txt'))

  ## medicc2
  medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
  medicc2_tree = drop.tip(medicc2_tree, 'diploid')
  medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'),
                                  function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))

  tree_rds[[f]] = list(
    scTrace_tree = scTrace_tree,
    NJ_tree = NJ_tree,
    MP_tree=MP_tree,
    ML_tree=ML_tree,
    random_tree=random_tree,
    sitka_tree = sitka_tree,
    medalt_tree=medalt_tree,
    medicc2_tree=medicc2_tree
  )
}
saveRDS(tree_rds, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_all_tree.rds')
tree_rds = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_all_tree.rds')
res = c()
for(f in c('Huh7','A549','786O','U251')){
  print(f)
  # load seurat rds
  srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{f}_filter_srt.RDS'))
  ###
  scTrace_tree = tree_rds[[f]][['scTrace_tree']]
  NJ_tree = tree_rds[[f]][['NJ_tree']]
  MP_tree = tree_rds[[f]][['MP_tree']]
  ML_tree = tree_rds[[f]][['ML_tree']]
  sitka_tree = tree_rds[[f]][['sitka_tree']]
  random_tree = tree_rds[[f]][['random_tree']]
  medalt_tree = tree_rds[[f]][['medalt_tree']]
  medicc2_tree = tree_rds[[f]][['medicc2_tree']]
  ##
  # 计算指标与crispr label
  gRNA_label = srt$feature_call

  scTrace_tree_s1 = get_clone(scTrace_tree, gRNA_label)
  NJ_tree_s1 = get_clone(NJ_tree, gRNA_label)
  MP_tree_s1 = get_clone(MP_tree, gRNA_label)
  ML_tree_s1 = get_clone(ML_tree, gRNA_label)
  sitka_tree_s1 = get_clone(sitka_tree, gRNA_label)
  medalt_tree_s1 = tryCatch({get_clone(medalt_tree, gRNA_label)},error = function(cond){return(c())})
  medicc2_tree_s1 = tryCatch({get_clone(medicc2_tree, gRNA_label)},error = function(cond){return(c())})

  plot_data = rbind(scTrace_tree_s1, NJ_tree_s1, MP_tree_s1, ML_tree_s1, sitka_tree_s1, medalt_tree_s1, medicc2_tree_s1)
  plot_data = cbind(plot_data, c(rep('scTrace', nrow(scTrace_tree_s1)),
                                 rep('NJ', nrow(NJ_tree_s1)),
                                 rep('MP', nrow(MP_tree_s1)),
                                 rep('ML', nrow(ML_tree_s1)),
                                 rep('sitka', nrow(sitka_tree_s1)),
                                 tryCatch({rep('medalt', nrow(medalt_tree_s1))},error = function(cond){return(c())}),
                                 tryCatch({rep('medicc2', nrow(medicc2_tree_s1))},error = function(cond){return(c())})

  ))
  plot_data = as.data.frame(plot_data)
  colnames(plot_data) = c('leaves_num', 'Score', 'Method')
  plot_data$Score = as.numeric(plot_data$Score)
  plot_data$leaves_num = as.numeric(plot_data$leaves_num)
  plot_data$Data = f
  res = rbind(res, plot_data)
}

res = as.data.frame(res)
res$Score[is.na(res$Score)] = 0
write.table(res,'/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_ARI_比较.txt')
res = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Crispr_ARI_比较.txt')
library(ggbeeswarm)
set.seed(123)
a=res %>%
  filter(leaves_num<=10) %>%
  group_by(Method, Data) %>%
  summarise(ARIindex = mean(Score)) %>%
  mutate(Method=factor(Method, c('scTrace', 'NJ', 'MP','medicc2', 'sitka', 'ML', 'medalt'))) %>%
  ggplot(aes_string(x='Method', y='ARIindex'))+
  geom_boxplot(outlier.shape = NA, width=0.8, lwd=0.1)+
  geom_quasirandom(aes(fill=Data),width = 0.1, shape=21, color='black', size=1, stroke=0.1)+
  #geom_jitter(aes(fill=Data),width = 0.1, shape=21, color='black', size=1, stroke=0.1)+
  scale_fill_manual(values = c('786O'='#2C86C8', 'A549'='#B7BC19', 'Huh7'='#E66327', 'U251'='#433691'))+
  theme_bw()+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format', size=2)+
  theme(panel.grid =element_blank(),
        axis.text.x = element_text(size=6, angle = 40, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=6),
        )
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_ARI比较.pdf',a,
       width=70, height=40, units='mm', dpi = 600, bg = 'transparent')


####
# 统计gRNA在clone里的分布
get_node_subtree <- function(tree, tip_names){
  tip_num = 1:length(tree$tip.label)
  names(tip_num) = tree$tip.label
  all_tip_path = nodepath(tree)
  max_depth = max(sapply(all_tip_path, function(x)length(x)))

  all_tip_path_df = c()
  for(i in all_tip_path){
    all_tip_path_df = rbind(all_tip_path_df, c(i, rep(NA, max_depth-length(i))))
  }
  #print(tip_names)
  #print(tip_num)
  #print(tip_num[tip_names])
  all_nodes_parent = all_tip_path_df[tip_num[tip_names],]
  parent_pos = apply(t(na.omit(t(all_nodes_parent))), 2, function(x)length(unique(x))==1)
  parent_pos = which(parent_pos==TRUE)
  parent_pos = parent_pos[length(parent_pos)]
  parent_node = all_nodes_parent[1,parent_pos]
  #print(parent_pos)
  #print(parent_node)
  # subtree
  subtree_df = all_tip_path_df[all_tip_path_df[,parent_pos]==parent_node, ]

  subtree_nodes = apply(subtree_df, 1, function(x){
    tmp = na.omit(x)
    return(tmp[length(tmp)])
  })
  subtree_nodes = unlist(subtree_nodes)
  return(names(tip_num[subtree_nodes]))
}

get_gRNA_score <-function(tree, gRNA){
  res = c()
  for(i in unique(gRNA)){
    sub_tree_leaves = get_node_subtree(tree, names(gRNA[gRNA==i]))
    #print(sub_tree_leaves)
    #print('------')
    score = vegan::diversity(table(gRNA[sub_tree_leaves]), index='shannon')
    res = rbind(res, c(score, i))
  }
  return(res)
}

res = c()
for(f in c('Huh7','A549','786O','U251')){
  print(f)
  # load seurat rds
  srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{f}_crispr_srt_filter.txt'))
  ###
  scTrace_tree = tree_rds[[f]][['scTrace_tree']]
  NJ_tree = tree_rds[[f]][['NJ_tree']]
  MP_tree = tree_rds[[f]][['MP_tree']]
  ML_tree = tree_rds[[f]][['ML_tree']]
  sitka_tree = tree_rds[[f]][['sitka_tree']]
  random_tree = tree_rds[[f]][['random_tree']]

  ##
  # 计算指标与crispr label
  gRNA_label = srt$feature_call

  medalt_tree = tree_rds[[f]][['medalt_tree']]
  medicc2_tree = tree_rds[[f]][['medicc2_tree']]
  ##
  # 计算指标与crispr label
  gRNA_label = srt$feature_call

  scTrace_tree_s1 = get_gRNA_score(scTrace_tree, gRNA_label)
  NJ_tree_s1 = get_gRNA_score(NJ_tree, gRNA_label)
  MP_tree_s1 = get_gRNA_score(MP_tree, gRNA_label)
  ML_tree_s1 = get_gRNA_score(ML_tree, gRNA_label)
  sitka_tree_s1 = get_gRNA_score(sitka_tree, gRNA_label)
  medalt_tree$tip.label = gsub('\\.', '-', medalt_tree$tip.label)
  medicc2_tree$tip.label = substr(medalt_tree$tip.label, 1, 18)
  medalt_tree_s1 = tryCatch({get_gRNA_score(medalt_tree, gRNA_label)},error = function(cond){return(c())})
  medicc2_tree_s1 = tryCatch({get_gRNA_score(medicc2_tree, gRNA_label)},error = function(cond){return(c())})

  plot_data = rbind(scTrace_tree_s1, NJ_tree_s1, MP_tree_s1, ML_tree_s1, sitka_tree_s1, medalt_tree_s1, medicc2_tree_s1)
  plot_data = cbind(plot_data, c(rep('scTrace', nrow(scTrace_tree_s1)),
                                 rep('NJ', nrow(NJ_tree_s1)),
                                 rep('MP', nrow(MP_tree_s1)),
                                 rep('ML', nrow(ML_tree_s1)),
                                 rep('sitka', nrow(sitka_tree_s1)),
                                 tryCatch({rep('medalt', nrow(medalt_tree_s1))},error = function(cond){return(c())}),
                                 tryCatch({rep('medicc2', nrow(medicc2_tree_s1))},error = function(cond){return(c())})

  ))
  plot_data = as.data.frame(plot_data)
  colnames(plot_data) = c('Score', 'gRNA', 'Method')
  plot_data$Score = as.numeric(plot_data$Score)
  plot_data$Data = f
  res = rbind(res, plot_data)
}

res %>%
  ggboxplot(x='Method', y='Score',  add = 'jitter')+
  #facet_wrap(~Data)+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format')


a=res %>%
  mutate(Method=factor(Method, c('scTrace', 'NJ', 'MP','medicc2', 'sitka', 'ML', 'medalt'))) %>%
  ggplot(aes_string(x='Method', y='Score'))+
  geom_boxplot(outlier.shape = NA, width=0.8, lwd=0.1)+
  #geom_quasirandom(aes(fill=gRNA), shape=21, color='black', size=1, stroke=0.1)+
  geom_jitter(aes(fill=gRNA),width = 0.3, shape=21, color='black', size=1, stroke=0.1)+
  scale_fill_manual(values = feature_call_color)+
  theme_bw()+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format', size=2)+
  theme(panel.grid =element_blank(),
        axis.text.x = element_text(size=6, angle = 40, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=6),
        legend.position = 'none'
  )
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/crispr_gRNA_diversity比较.pdf',a,
       width=45, height=40, units='mm', dpi = 600, bg = 'transparent')


