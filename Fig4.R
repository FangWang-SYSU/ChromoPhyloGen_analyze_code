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

srt = readRDS('/Volumes/WX_extend/CRC项目/tumors_filter_by_cnv_man.rds')

lesions_col = c('Colon_tumor'='#5050FFFF', 'Liver_meta'='#CE3D32FF')
a = DimPlot(srt, group.by = 'lesions', cols = lesions_col, pt.size = 0.01) + theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5))

a = a$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=lesions))+
  geom_point(shape=21, color='#00000000', size=0.001)+
  scale_fill_manual(values = lesions_col)+
  theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_umap_转移着色.pdf',a,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')
file_name = c('Pt01', 'Pt02', 'Pt03', 'Pt04', 'Pt05', 'Pt06', 'Pt07', 'Pt08', 'Pt09', 'Pt10',
              'Pt11', 'Pt13', 'Pt15', 'Pt16', 'Pt17', 'Pt18')
patient_col = pal_igv()(20)[3:20]
names(patient_col) = file_name

a = DimPlot(srt, group.by = 'patient', cols = patient_col, pt.size = 0.001) + theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5),
        legend.key.size = unit(0.5,'mm'))
a = a$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=patient))+
  geom_point(shape=21, color='#00000000', size=0.001)+
  scale_fill_manual(values = patient_col)+
  theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5),
        legend.key.size = unit(0.5,'mm'),
        legend.text = element_text(size=6))

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_umap_病人着色.pdf',a,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')



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



CRC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CRC_sctc.rds')

# Pt13 展示
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
                     linewidth=7,
                     #lwd=pmax(10/nchar(g$branches), 1),
                     linetype='solid',
                     lineend = "round",
                     linejoin='round',
)
# p = p+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone),
#                               data=cell_pos[cell_pos$clone!='clone_0', ],
#                               alpha=0.8,
#                               lwd=2,
#                               color='gray',
#                               lineend='round',linejoin='round',
#                               #size=10,
# )

p = p + geom_point(aes(x=x, y=y, fill=!!sym('lesions')),data=cell_pos,
                   shape=21, size=1.2,inherit.aes = F, color='white', stroke=0.1)
p=p +
  theme_void() +
  scale_y_continuous(expand = expansion(mult = 0.3))+
  scale_x_continuous(expand = expansion(mult = 0.3))+
  guides(alpha=FALSE)+ggsci::scale_fill_igv()




ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_Pt13_转移着色.pdf',p,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')


####
new_srt =readRDS('/Volumes/WX_extend/CRC项目/Epithelial_raw.rds')
sub_srt = subset(new_srt, patient==prefix)
sub_srt = AddModuleScore(sub_srt, features = list(cc.genes$g2m.genes))
ccgene_zzm = c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
               'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
               'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1')
sub_srt = AddModuleScore(sub_srt, features = list(ccgene_zzm), name='zzm')

sub_srt = CellCycleScoring(sub_srt, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)
tmp_rate = predict_mitosis_rate(as.data.frame(sctc$orig.data$cell_relation))

cell_pos$predict_rate = tmp_rate[rownames(cell_pos)]
cell_pos$CellCycle_score=sub_srt@meta.data[rownames(cell_pos), 'zzm1']



p_rate = ggplot()
p_rate = p_rate + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                     color='#dddddd',
                     alpha=0.6,
                     linewidth=7,
                     #lwd=pmax(10/nchar(g$branches), 1),
                     linetype='solid',
                     lineend = "round",
                     linejoin='round',
)
p_rate = p_rate+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone, color=predict_rate),
                              data=cell_pos[cell_pos$clone!='clone_0', ],
                              alpha=1,
                              lwd=1,
                              #color='gray',
                              lineend='round',linejoin='round',
                              #size=10,
)

p_rate = p_rate + geom_point(aes(x=x, y=y, fill=!!sym('predict_rate')),data=cell_pos,
                   shape=21, size=1.2,inherit.aes = F, color='#00000000', stroke=0.1)
p_rate = p_rate +
  theme_void() +
  guides(alpha=FALSE)+
  #scale_fill_gradient(  low = "#white", high = "#red")
  scale_color_gradientn(colors = c('blue', 'white', 'red'), limits=c(0,0.5),na.value = "blue")+
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits=c(0,0.5),na.value = "blue")+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_x_continuous(expand = expansion(mult = 0.1))
#
p_cc = ggplot()
p_cc = p_cc + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                               color='#dddddd',
                           alpha=0.6,
                           linewidth=7,
                               #lwd=pmax(10/nchar(g$branches), 1),
                               linetype='solid',
                               lineend = "round",
                               linejoin='round',
)
p_cc = p_cc+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone, color=CellCycle_score),
                                        data=cell_pos[cell_pos$clone!='clone_0', ],
                                        alpha=1,
                                        lwd=1,
                                        #color='gray',
                                        lineend='round',linejoin='round',
                                        #size=10,
)

p_cc = p_cc + geom_point(aes(x=x, y=y, fill=!!sym('CellCycle_score')),data=cell_pos,
                             shape=21, size=1.2,inherit.aes = F, color='#00000000')
p_cc = p_cc +
  theme_void() +
  #guides(alpha=FALSE)
  scale_color_gradientn(colors = c('blue', 'white', 'red'), limits=c(-0.4,0.4),na.value = "blue",)+
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits=c(-0.4,0.4),na.value = "blue",)+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_x_continuous(expand = expansion(mult = 0.1))

a=cowplot::plot_grid(p,p_rate,p_cc, nrow=1, align = 'h')

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_Pt13_rate着色.pdf',p_rate,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_Pt13_增殖着色.pdf',p_cc,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_Pt13.pdf',a,
       width=200, height=60, units='mm', dpi = 600, bg = 'transparent')

### 速率和增殖评分着色umap
ccgene_zzm = c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2',
               'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1',
               'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS', 'DEK', 'CCNB1', 'CCNE1')
srt = AddModuleScore(srt, features = list(ccgene_zzm), name='zzm')
srt = AddModuleScore(srt, features = list(cc.genes$g2m.genes))
srt = CellCycleScoring(srt, g2m.features = cc.genes$g2m.genes, s.features=cc.genes$s.genes)
srt[['clone_tree']] = 'other'
for(prefix in file_name){
  sctc=CRC_list[[prefix]]
  same_name = intersect(colnames(srt), rownames(sctc$map_obj$cell_map_pos))
  srt@meta.data[same_name, 'clone_tree'] = paste0(prefix,'_',sctc$map_obj$cell_map_pos[same_name, 'clone'])
}
srt$mean_exp = rowMeans(FetchData(srt, vars = ccgene_zzm))
# 全部细胞rate
crc_rate = c()
for(prefix in file_name){
  print(prefix)
  tmp_class = CRC_list[[prefix]]
  tmp_obj = tmp_class$orig.data
  meta = tmp_obj$cell_relation
  tmp_rate = predict_mitosis_rate(as.data.frame(meta))
  tmp_rate = (tmp_rate-min(tmp_rate)) / (max(tmp_rate) - min(tmp_rate))

  #tmp_rate = scale(tmp_rate)
  #tmp_rate = tmp_rate[,1]
  crc_rate = c(crc_rate, tmp_rate)
}

srt[['predict_rate']] = crc_rate[colnames(srt)]

tmp_func <- function(x,y){
  x[which.max(y)]
}
new_Fn = srt@meta.data %>%
  group_by(clone_tree, lesions) %>%
  summarise(n=n()) %>%
  na.omit() %>%
  group_by(clone_tree) %>%
  summarise(clone_lesions=tmp_func(lesions, n))
srt[['clone_lesions']] = NA
srt@meta.data[,'clone_lesions'] = new_Fn[match(srt$clone_tree, new_Fn$clone_tree), 'clone_lesions']

### 速率着色umap

a = FeaturePlot(srt, features='predict_rate')+scale_color_gradientn(colors = c('blue', 'gray','red'))
a = a$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=predict_rate))+
  geom_point(shape=21, color='#00000000', size=0.001)+
  scale_fill_gradientn(colours = c('blue', 'white','red'),values=c(0,0.75,1))+
  theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
#a
#a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_umap_rate着色.pdf',a,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')

### 增殖评分着色umap
b = FeaturePlot(srt, features='G2M.Score')+scale_color_gradientn(colors = c('blue', 'gray','red'))
b = b$data %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=G2M.Score))+
  geom_point(shape=21, color='#00000000', size=0.001)+
  #scale_fill_gradient2(low='blue', mid='white',high='red',midpoint = 0)+
  scale_fill_gradientn(colours = c('blue', 'white','red'),values=c(0,0.2,1))+
  theme_void()+
  theme(plot.title = element_text(size=8, hjust = 0.5),
        legend.key.size = unit(4,'mm'),
        legend.text = element_text(size=6))

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_umap_增殖着色.pdf',b,
       width=70, height=50, units='mm', dpi = 600, bg = 'transparent')


# 样本水平相关性
library(ggpubr)

#DimPlot(srt, group.by = 'clone_tree')
library(ggExtra)
a=srt@meta.data %>%
  group_by(clone_tree, lesions) %>%
  summarise('velocity'=mean(predict_rate), 'cc_score'=mean(G2M.Score)) %>%
  na.omit() %>%
  ggplot(aes(x=cc_score, y=velocity))+
  geom_point(aes(color=lesions), size=0.2)+
  scale_color_manual(values = lesions_col)+
  geom_smooth(method='lm',formula = 'y ~ x', se=F, color='black', linewidth=0.3)+
  stat_cor(size=2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6, angle = 40, hjust = 1, vjust = 1),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))

a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_速率与G2M相关.pdf',a,
       width=70, height=40, units='mm', dpi = 600, bg = 'transparent')

clone_exp = AverageExpression(srt, group.by = 'clone_tree', slot = 'counts')$RNA
library(CytoTRACE)
results <- CytoTRACE(clone_exp)
results$CytoTRACE

x = srt@meta.data %>%
  group_by(clone_tree, lesions) %>%
  summarise('velocity'=mean(predict_rate)) %>% as.data.frame()
x$cytotrace = results$CytoTRACE[x$clone_tree]
x = na.omit(x)
a=ggplot(x, aes(x=cytotrace, y=velocity))+
  geom_point(aes(color=lesions), size=0.2)+
  scale_color_manual(values = lesions_col)+
  geom_smooth(method='lm',formula = 'y ~ x', se=F, color='black', linewidth=0.3)+
  stat_cor(size=2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_速率与Cytotrace相关.pdf',a,
       width=80, height=50, units='mm', dpi = 600, bg = 'transparent')

# b1 = ggplot(x, aes(x=lesions, y=velocity))+
#   geom_boxplot()+
#   geom_jitter(width = 0.1)+
#   stat_compare_means(comparisons = list(c('Liver_meta', 'Colon_tumor')), label='p.format')+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text = element_text(size=6, angle = 40, hjust = 1, vjust = 1),
#         axis.title = element_text(size=6),
#         legend.key.size = unit(3,'mm'),
#         legend.text = element_text(size=6))
#
# ggplot(x, aes(x=lesions, y=cytotrace ))+
#   geom_boxplot()

# 将细胞分两组，识别差异基因功能
x$Clone_state = 'High_velocity'
x[x$velocity<median(x$velocity), 'Clone_state'] = 'Low_velocity'
a = x%>%
  ggplot(aes(x=velocity,color=Clone_state))+
  geom_density(bw=0.1)+
  geom_point(aes(y=0), size=0.2)+
  scale_color_npg()+
  labs('y'='Density')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_速率密度图.pdf',a,
       width=70, height=40, units='mm', dpi = 600, bg = 'transparent')


srt[['Cell_state']] = 'High_velocity'
srt@meta.data[srt$clone_tree%in%x[x$Clone_state=='Low_velocity', 'clone_tree'], 'Cell_state']= 'Low_velocity'

DimPlot(srt, group.by = 'Cell_state')

# 将细胞分两组，识别差异基因功能
srt@meta.data =srt@meta.data[!is.na(srt$predict_rate),]
srt = AddMetaData(srt, srt@meta.data[!is.na(srt$predict_rate),])
srt[['Cell_state']] = 'Mid'
srt@meta.data[srt$predict_rate>quantile(srt$predict_rate,0.75,na.rm = T), 'Cell_state']= 'High_velocity'
srt@meta.data[srt$predict_rate<quantile(srt$predict_rate,0.25,na.rm = T), 'Cell_state']= 'Low_velocity'

DimPlot(srt, group.by = 'Cell_state')
# 识别差异基因
Idents(srt) = srt$Cell_state

srt_marker = FindMarkers(srt, ident.1 = 'High_velocity', ident.2 = 'Low_velocity',
                         max.cells.per.ident = 300, logfc.threshold = 0)

srt_marker_top = srt_marker[order(srt_marker$avg_log2FC), ]
srt_marker_top$group = 'ns'
srt_marker_top[srt_marker_top$avg_log2FC>0.25&srt_marker_top$p_val_adj<0.01, 'group'] = 'Up'
srt_marker_top[srt_marker_top$avg_log2FC<(-0.25)&srt_marker_top$p_val_adj<0.01, 'group'] = 'Down'
srt_marker_top %>%
  ggplot(aes(x=avg_log2FC, y=-log(p_val_adj), color=group))+
  geom_point(size=1)+
  geom_hline(yintercept = -log(0.05), linetype='dashed', color='gray')+
  geom_vline(xintercept = c(-0.25, 0.25), linetype='dashed', color='gray')+
  scale_color_manual(values = c('ns'='gray', 'Up'='darkred', 'Down'='darkblue'))+
  theme_bw()+
  theme(panel.grid = element_blank())

srt_marker_top$group = 'ns'
srt_marker_top[srt_marker_top$avg_log2FC>0.25&srt_marker_top$p_val_adj<0.01, 'group'] = 'Up'
srt_marker_top[srt_marker_top$avg_log2FC<(-0.25)&srt_marker_top$p_val_adj<0.01, 'group'] = 'Down'
marker_label = srt_marker_top[srt_marker_top$group=='Up', ]
diff_pct = marker_label$pct.1 - marker_label$pct.2
marker_label = marker_label[diff_pct > quantile(diff_pct, 0.95),]
marker_label$gene = rownames(marker_label)
marker_label = marker_label%>% arrange(-avg_log2FC)
a = srt_marker_top %>%
  ggplot(aes(x=pct.2, y=pct.1, size=-log(p_val_adj), color=avg_log2FC))+
  geom_point()+
  geom_text_repel(aes(x=pct.2, y=pct.1,label=gene),data=marker_label[1:5,], size=2, inherit.aes = F, )+
  scale_color_gradient2(low='darkblue', mid='white', high='darkred',midpoint=0)+
  scale_size_continuous(range = c(0.01,0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_差异基因图.pdf',a,
       width=80, height=40, units='mm', dpi = 600, bg = 'transparent')


###
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
a = enrichplot::dotplot(tmp_enrich)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))+
  scale_x_continuous(expand = expansion(mult = 0.15))

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_功能富集.pdf',a,
       width=120, height=40, units='mm', dpi = 600, bg = 'transparent')

### depmap 验证 ###

effect_score = read.csv('/Users/lab/Downloads/CRISPR_gene_effect.csv', check.names = F, row.names = 'DepMap_ID')
up_gene_effect = effect_score[, grep(paste0(marker_label$gene, collapse = '|'), colnames(effect_score)), drop=F]

set.seed(1234)
random_gene = sample(rownames(srt_marker_top[srt_marker_top$group=='ns', ]), 50)
random_gene_effect = effect_score[, grep(paste0(random_gene, collapse = '|'), colnames(effect_score)), drop=F]

up_gene_effect = reshape2::melt(as.matrix(up_gene_effect))
up_gene_effect$source = 'Up'

random_gene_effect = reshape2::melt(as.matrix(random_gene_effect))
random_gene_effect$source = 'random'

gene_effect = rbind(up_gene_effect, random_gene_effect)

a=gene_effect %>%
  group_by(Var1,source) %>%
  summarise(value=mean(value, na.rm = T)) %>%
  ggplot(aes(x=source, y=value, fill=source))+
  geom_boxplot(size=0.1, outlier.size = 0.01, width=0.2)+
  scale_fill_manual(values = c('Up'='#CE3D32FF', 'random'='gray'))+
  stat_compare_means(comparisons = list(c('Up', 'random')), label='p.format', size=2)+
  scale_y_continuous(expand = expansion(mult = 0.15))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=6),
        #axis.title.x = element_blank(),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/CRC_depmap验证.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')









