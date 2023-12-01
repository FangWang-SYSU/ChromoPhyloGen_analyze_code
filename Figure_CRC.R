source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
library(Seurat)
seurat_process <-function(obj, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:20)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:20,reduction = "pca", n.components=n.components)
  return(obj)
}
findmarker_volcano <- function(findallmarker_res, pvalue=0.01, top_marker=5, self_genes=NULL, logfc.threshold = 0.25){
  findallmarker_res$label = paste0('adjust P-val<',pvalue)
  findallmarker_res[findallmarker_res$p_val_adj>=pvalue, 'label'] = paste0('adjust P-val>=',pvalue)
  if(is.null(self_genes)){
    srt_markers_top5 = findallmarker_res %>% group_by(cluster) %>% top_n(n = top_marker, wt = avg_log2FC)
  }else{
    top_g = intersect(findallmarker_res$gene, self_genes)
    srt_markers_top5 = findallmarker_res[findallmarker_res$gene%in%top_g,]
  }
  #根据图p中log2FC区间确定背景柱长度：
  dfbar = findallmarker_res %>%
    group_by(cluster) %>%
    summarise(top_bar=max(avg_log2FC)) %>%
    merge(findallmarker_res %>%
            group_by(cluster) %>%
            summarise(bottom_bar=min(avg_log2FC)))
  # 根据log2FC确定上下色块范围
  #top_col_height = min(srt_markers[srt_markers$avg_log2FC>0, 'avg_log2FC'])
  #bottom_col_height = max(srt_markers[srt_markers$avg_log2FC<0, 'avg_log2FC'])
  #col_height = round(min(abs(bottom_col_height), abs(top_col_height)), 2)
  col_height = logfc.threshold *2
  ggplot()+
    geom_point(data = findallmarker_res[!findallmarker_res$gene%in%srt_markers_top5$gene, ],
               aes(x = cluster, y = avg_log2FC, color = label),
               size = 0.85,
               position=position_jitter(seed = 1994),
    )+ # 散点
    geom_point(data = srt_markers_top5,
               aes(x = cluster, y = avg_log2FC, color = label),
               color='blue',
               position=position_jitter(seed = 1994),
               size = 1)+ # top 点
    geom_col(data = dfbar,
             mapping = aes(x = cluster,y = top_bar),
             fill = "#dcdcdc",alpha = 0.4)+ # 背景柱
    geom_col(data = dfbar,
             mapping = aes(x = cluster,y = bottom_bar),
             fill = "#dcdcdc",alpha = 0.4)+ # 背景柱
    geom_tile(data = dfbar,
              aes(x=cluster,y=0,fill = cluster),
              height= col_height,
              color = "black",
              alpha = 0.6,
              show.legend = F)+ # 添加色块
    geom_text(data=dfbar,
              aes(x=cluster,y=0,label=cluster),
              size =6,
              color ="white")+ # 色块添加文字
    geom_text_repel(
      data=srt_markers_top5,
      aes(x=cluster,y=avg_log2FC,label=gene),
      force = 1.2,
      position = position_jitter(seed = 1994),
      arrow = arrow(length = unit(0.008, "npc"),
                    type = "open", ends = "last")
    )+ # top 基因添加注释
    scale_color_manual(name=NULL, values = c("red","black"))+ #散点颜色
    scale_fill_igv()+ #散点颜色
    theme_minimal()+
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y = element_line(color = "black",
                                 size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 15)
    )

}
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC/'

file_name = c('Pt01', 'Pt02', 'Pt03', 'Pt04', 'Pt05', 'Pt06', 'Pt07', 'Pt08', 'Pt09', 'Pt10',
              'Pt11', 'Pt13', 'Pt15', 'Pt16', 'Pt17', 'Pt18')

CRC_list = list()
for(prefix in file_name){
  print(prefix)
  sctc = scTraceClass(scTrace_dir,
                      prefix,
                      layout = 'slanted',
                      rate = 0.95,
                      min_cell_num=NULL)
  CRC_list[[prefix]] = sctc
}
saveRDS(CRC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc.rds')

crc_rate = c()
for(prefix in file_name){
  print(prefix)
  tmp_class = CRC_list[[prefix]]
  tmp_obj = tmp_class$orig.data
  meta = tmp_obj$cell_relation
  tmp_rate = predict_mitosis_rate(as.data.frame(meta))
  crc_rate = c(crc_rate, tmp_rate)
}
plot(density(crc_rate))
srt = readRDS('/Volumes/WX_extend/CRC项目/tumors_filter_by_cnv_man.rds')
srt[['predict_rate']] = crc_rate[colnames(srt)]

FeaturePlot(srt, features='predict_rate')+scale_color_gradientn(colors = c('blue', 'gray','red'))
DimPlot(srt, group.by = 'lesions')
# 将细胞分两组，识别差异基因功能
srt[['Cell_state']] = 'Mid'
srt@meta.data[srt$predict_rate>quantile(srt$predict_rate,0.75), 'Cell_state']= 'High_rate'
srt@meta.data[srt$predict_rate<quantile(srt$predict_rate,0.25), 'Cell_state']= 'Low_rate'

DimPlot(srt, group.by = 'Cell_state')
# 识别差异基因
Idents(srt) = srt$Cell_state

srt_marker_all = FindAllMarkers(srt,  max.cells.per.ident = 300)
srt_marker_all$cluster = factor(srt_marker_all$cluster, c('Low_rate', 'Mid',  'High_rate'))
findmarker_volcano(srt_marker_all,)


srt_marker = FindMarkers(srt, ident.1 = 'High_rate', ident.2 = 'Low_rate', max.cells.per.ident = 300, logfc.threshold = 0)
srt_marker_top = srt_marker[order(srt_marker$avg_log2FC), ]
srt_marker_top$group = 'ns'
srt_marker_top[srt_marker_top$avg_log2FC>0.1&srt_marker_top$p_val_adj<0.05, 'group'] = 'Up'
srt_marker_top[srt_marker_top$avg_log2FC<(-0.1)&srt_marker_top$p_val_adj<0.05, 'group'] = 'Down'
srt_marker_top %>%
  ggplot(aes(x=avg_log2FC, y=-log(p_val_adj), color=group))+
  geom_point(size=1)+
  geom_hline(yintercept = -log(0.05), linetype='dashed', color='gray')+
  geom_vline(xintercept = c(-0.25, 0.25), linetype='dashed', color='gray')+
  scale_color_manual(values = c('ns'='gray', 'Up'='darkred', 'Down'='darkblue'))+
  theme_bw()+
  theme(panel.grid = element_blank())

## 功能富集
hallmark_data = GSEABase::getGmt('/Users/lab/wangxin/Rprojects/cancer/h.all.v7.4.symbols.gmt')
hallmark_data = hallmark_data@.Data
hallmark2gene = c()
for(i in hallmark_data){
  hallmark2gene = rbind(hallmark2gene, data.frame(termid=i@setName, geneid=i@geneIds))
}
hallmark2gene = as.data.frame(hallmark2gene)
library(clusterProfiler)
tmp_gene = rownames(srt_marker_top[srt_marker_top$group=='Up', ])

# tmp_gene = srt_marker_all[srt_marker_all$cluster=='High_rate', ]
# tmp_gene = tmp_gene[tmp_gene$avg_log2FC>0.25&tmp_gene$p_val_adj<0.05, 'gene']
tmp_enrich = enricher(tmp_gene, TERM2GENE = hallmark2gene)
#
enrichplot::dotplot(tmp_enrich)
tmp_enrich = tmp_enrich@result[tmp_enrich@result$p.adjust<0.05, ]
tmp_hm = head(tmp_enrich$ID, 5)


####
plist = list()
for(prefix in file_name){
  print(prefix)
  tmp_class = CRC_list[[prefix]]

  sub_srt = subset(srt, patient==prefix)
  sub_srt[['clone_tree']] = 'other'
  sub_srt[['clone_tree']] = tmp_class$clone_graph_obj$cell_pos_reduction[colnames(sub_srt), 'clone']
  DimPlot(sub_srt, group.by = 'clone_tree') + ggsci::scale_color_igv()
  #####
  tmp_class$map_obj$cell_map_pos$lesions = sub_srt@meta.data[rownames(tmp_class$map_obj$cell_map_pos), 'lesions']
  tmp_class$map_obj$cell_map_pos$lesions[is.na(tmp_class$map_obj$cell_map_pos$lesions)] = 'Liver_meta'
  plot_cells(tmp_class, colorby='lesions',shadow=TRUE, clone_line=TRUE,
             point.size=1, plot_virtual=F)+
    ggsci::scale_fill_igv()

  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]
  cell_pos = cell_pos[!is.na(cell_pos$lesions), ]

  #
  r_size = 0.04
  pie_data = cell_pos %>%
    group_by(clone, !!sym('lesions')) %>%
    summarise(num=n()) %>%
    mutate(total=sum(num)) %>% as.data.frame()
  pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
  pie_data = pie_data %>% group_by(x, y, total)
  df.grobs <- pie_data %>%
    do(subplots = ggplot(., aes(1, num, fill = !!sym('lesions'))) +
         geom_col(position = "fill", alpha = 1, colour = "white") +
         coord_polar(theta = "y") +
         # scale_fill_igv()+
         scale_fill_manual(values = c('Liver_meta'='#BA6338FF', 'Colon_tumor'='#5050FFFF'))+
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
cowplot::plot_grid(plotlist = plist, ncol=6)

order_file = c('Pt06','Pt13')
cowplot::plot_grid(plotlist = plist[order_file], ncol=2)


prefix = 'Pt13'
sub_srt = subset(srt, patient==prefix)
#sub_srt = seurat_process(sub_srt)
#DimPlot(sub_srt, group.by = 'lesions') + ggsci::scale_color_igv()
sctc = CRC_list[[prefix]]
sctc$clone_graph_obj$cell_phylo_tree_collapse

sub_srt[['clone_tree']] = 'other'
sub_srt[['clone_tree']] = sctc$clone_graph_obj$cell_pos_reduction[colnames(sub_srt), 'clone']
DimPlot(sub_srt, group.by = 'clone_tree') + ggsci::scale_color_igv()
#####
sctc$map_obj$cell_map_pos$lesions = sub_srt@meta.data[rownames(sctc$map_obj$cell_map_pos), 'lesions']
sctc$map_obj$cell_map_pos$lesions[is.na(sctc$map_obj$cell_map_pos$lesions)] = 'Liver_meta'
plot_cells(sctc, colorby='lesions',shadow=TRUE, clone_line=TRUE,
           point.size=1, plot_virtual=F)+
  ggsci::scale_fill_igv()

cell_pos = sctc$map_obj$cell_map_pos
clone_pos = sctc$map_obj$clone_seg_data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label
cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]
cell_pos = cell_pos[!is.na(cell_pos$lesions), ]

p = ggplot()
p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                     color='#dddddd',
                     alpha=0.6,
                     linewidth=6,
                     #lwd=pmax(10/nchar(g$branches), 1),
                     linetype='solid',
                     lineend = "round",
                     linejoin='round',
)
p = p+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone),
                                data=cell_pos[cell_pos$clone!='clone_0', ],
                      alpha=0.8,
                      lwd=2,
                      color='gray',
                      lineend='round',linejoin='round',
                      #size=10,
)

p = p + geom_point(aes(x=x, y=y, fill=!!sym('lesions')),data=cell_pos,
                   shape=21, size=2,inherit.aes = F)
p +
  theme_void() +
  guides(alpha=FALSE)+ggsci::scale_fill_igv()


#
r_size = 0.04
pie_data = cell_pos %>%
  group_by(clone, !!sym('lesions')) %>%
  summarise(num=n()) %>%
  mutate(total=sum(num)) %>% as.data.frame()
pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
pie_data = pie_data %>% group_by(x, y, clone)
df.grobs <- pie_data %>%
  do(subplots = ggplot(., aes(1, num, fill = !!sym('lesions'))) +
       geom_col(position = "fill", alpha = 1, colour = "white") +
       coord_polar(theta = "y") +
       scale_fill_manual(values = c('Liver_meta'='#BA6338FF', 'Colon_tumor'='#5050FFFF'))+
       #scale_fill_igv()+
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
p+  theme_void() +
  guides(alpha=FALSE)

###

plot_cna_tree(sctc,
              colorby = 'lesions',
              output = paste0(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC/{prefix}_tree_heatmap.pdf')),
              #column_split=sort(as.numeric(gsub('chr', '',gene.loc))),
              pt.size = 0.1,
              width=220,
              height=100,
              rel_widths=c(0.3,1),
              rel_heights = c(.0, 1),
              top_ann = NULL
)



# RNA增殖评分
sub_srt = AddModuleScore(sub_srt, features = list(cc.genes$g2m.genes))
sub_srt = CellCycleScoring(sub_srt, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)

tmp_rate = predict_mitosis_rate(as.data.frame(sctc$orig.data$cell_relation))

cell_pos$predict_rate = tmp_rate[rownames(cell_pos)]



p = ggplot()
p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                     color='#dddddd',
                     alpha=0.6,
                     linewidth=6,
                     #lwd=pmax(10/nchar(g$branches), 1),
                     linetype='solid',
                     lineend = "round",
                     linejoin='round',
)
p = p+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone, color=predict_rate),
                              data=cell_pos[cell_pos$clone!='clone_0', ],
                              alpha=0.8,
                              lwd=2,
                              #color='gray',
                              lineend='round',linejoin='round',
                              #size=10,
)

p = p + geom_point(aes(x=x, y=y, fill=!!sym('predict_rate')),data=cell_pos,
                   shape=21, size=2,inherit.aes = F)
p +
  theme_void() +
  guides(alpha=FALSE)+
  scale_color_gradientn(colors = c('blue', 'white', 'red'), limits=c(0,0.5))+
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits=c(0,0.5))


