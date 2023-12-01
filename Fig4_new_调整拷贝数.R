source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
library(Seurat)
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

# load pathway gene list
#############
bc_file = list.files('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/patient_cnv/')
scTrace_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/patient_cnv_infer4_res//'


BC_list = list()
for(prefix in bc_file){
  print(prefix)
  if(prefix %in% names(BC_list)){

  }else{
    if(prefix%in%c("CID4290A.txt","CID44991.txt","P30.txt")){
      sctc = scTraceClass(scTrace_dir,
                          prefix,
                          layout = 'slanted',
                          rate = 0.98,
                          min_cell_num=100,
                          min_cell_num2=100)
    }else if(prefix%in%c("CID4461.txt", 'P1.txt')){
    }else{
      sctc = scTraceClass(scTrace_dir,
                          prefix,
                          layout = 'slanted',
                          rate = 0.98,
                          min_cell_num=NULL,
                          min_cell_num2=NULL)
    }

    BC_list[[prefix]] = sctc
  }

}

saveRDS(BC_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_sctc_infer4_res.rds')
BC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_sctc_infer4_res.rds')


srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor.rds')

# srt[['cytotrace_score']] = results$CytoTRACE[colnames(srt)]

scTrace_dir3 = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/patient_cnv_infer4_res/'
all_mechnism = c()
BC_list_with_score = list()
for(prefix in names(BC_list)){
  print(prefix)
  # prefix='CID3586.txt'
  sctc = BC_list[[prefix]]
  CNA_mechnism = read.table(glue('{scTrace_dir3}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score = read.table(glue('{scTrace_dir3}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir3}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir3}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop[all_limit_prop<0] = 0
  pvalue_thr = 0.001
  CNA_mechnism$rearrange_score = rowMeans(all_rearrange_score)
  CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr)
  CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
  CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
  CNA_mechnism$BFB = CNA_mechnism$BFB / ncol(sctc$orig.data$all_node_data)
  #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)] = NA
  CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
  #
  all_rearrange_score2 = all_rearrange_score
  all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)] = NA
  CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)

  CNA_mechnism[is.na(CNA_mechnism)] = 0
  pde = predict_PDE(sctc)
  CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
  CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

  sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

  BC_list_with_score[[prefix]] = sctc

  CNA_mechnism$file = prefix
  all_mechnism = rbind(all_mechnism, CNA_mechnism)

}

same_bc = intersect(colnames(srt), rownames(all_mechnism))

srt = subset(srt, cells=same_bc)
srt = process_srt(srt,resolution = 0.1)
DimPlot(srt)
all_mechnism = all_mechnism[same_bc, ]
for(i in colnames(all_mechnism)){
  if(i %in% colnames(srt@meta.data)){
    srt[[i]] = NULL
  }
}

srt = AddMetaData(srt, all_mechnism)
saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor_with_score_infer4_res.rds')
saveRDS(BC_list_with_score, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_score_sctc_infer4_res.rds')

srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor_with_score_infer4_res.rds')
BC_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_score_sctc_infer4_res.rds')

srt_old = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor_with_score2.rds')
srt_new = srt

srt = srt_new#subset(srt_new, cells=colnames(srt_old))
srt = process_srt(srt,resolution = 0.1)
DimPlot(srt)
srt@meta.data[, 'cytotrace_score'] = srt_old@meta.data[colnames(srt), 'cytotrace_score']

###
# cytotrace
library(CytoTRACE)

srt[['cytotrace_score']] = 0#results$CytoTRACE[colnames(srt)]
srt[['tmp_group']] = paste0(srt$file, '_',srt$clone)
count_mat = AverageExpression(srt, group.by = 'tmp_group')$RNA
results <- CytoTRACE(count_mat)
srt@meta.data[, 'cytotrace_score'] = results$CytoTRACE[srt$tmp_group]


count_mat = FetchData(srt, vars = VariableFeatures(srt), slot = 'counts')

results <- CytoTRACE(t(count_mat))
srt@meta.data[, 'cytotrace_score'] = results$CytoTRACE[colnames(srt)]

FeaturePlot(srt, features = 'cytotrace_score') +
  scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))
##
#srt = subset(srt, cells=rownames(srt@meta.data[!is.na(srt$clone),]))
#srt = process_srt(srt)
#DimPlot(srt)

# 增加转录组异质性
gene_expression = FetchData(srt, vars = VariableFeatures(srt))
res = apply(gene_expression, 2, function(y){
  x = cut(y, breaks = 5, labels = F)
  x = as.numeric(as.vector(x))
  names(x) = names(y)
  return(x)
})
cell_ITH = apply(res, 1, function(a)vegan::diversity(a[a!=1]))
#cell_ITH = as.data.frame(cell_ITH)
srt[['Cell_ITH']] = cell_ITH[colnames(srt)]

saveRDS(srt, '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor_with_score_infer4_res.rds')
srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/tumor_with_score_infer4_res.rds')

# 添加亚型信息
srt_info = read.csv('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNVs_BreastCancer/TNBC_Clinical_info.csv', header = F)
# ER-/RP-/HER2+ 定义为 HER2 阳性（HER2阳性型）
# ER+/ RP+/-/HER2+/-定义为 ER 阳性（Luminal型）
# ER-/RP-/HER2- 定义为 TNBC

# #HER2阳性型
# ER-/PR-/HER2+
# HER2+
# #Luminal型
# ER+
# ER+/PR-/HER2-
# ER+/PR+ /HER2+
# ER+/PR+/HER2-
# #TNBC
# TNBC
srt_info[srt_info$V2 %in%c('ER-/PR- /HER2+', 'HER2+'), 'type'] = 'HER2+'
srt_info[srt_info$V2 %in%c('ER+', 'ER+/PR-/HER2-', 'ER+/PR+ /HER2+', 'ER+/PR+/HER2-'), 'type'] = 'Luminal'
srt_info[srt_info$V2 %in%c('TNBC'), 'type'] = 'TNBC'
rownames(srt_info) = srt_info$V1
srt[['CancerSub']] = srt_info[srt$patient, 'type']
# umap 着色
p1 = DimPlot(srt, group.by = 'patient')+scale_color_igv()
p2 = FeaturePlot(srt, features = 'cytotrace_score') +
  scale_colour_gradientn(colours = c("lightgrey","#3288BD" ,"#ABDDA4","#FEE08B", "#F46D43","#9E0142"))
p3 = FeaturePlot(srt, features = 'CellCycle1')+scale_colour_gradientn(colours = c("lightgrey","red"))
p4 = FeaturePlot(srt, features = 'limit_score')+scale_colour_gradientn(colours = c("lightgrey","red"))
p5 = FeaturePlot(srt, features = 'seismic_score')+scale_colour_gradientn(colours = c("lightgrey","red"))
p6 = FeaturePlot(srt, features = 'CPS_num')+scale_colour_gradientn(colours = c("lightgrey","red"))
p7 = FeaturePlot(srt, features = 'BFB')+scale_colour_gradientn(colours = c("lightgrey","red"))
p8 = FeaturePlot(srt, features = 'rearrange_score')+scale_colour_gradientn(colours = c("lightgrey","red"))
p9 = FeaturePlot(srt, features = 'pde')+scale_colour_gradientn(colours = c("lightgrey","red"))

a = cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=3)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/整体UMAP.pdf'),
       a, width=380, height = 380, units='mm', dpi = 600)

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/整体UMAP_patient.pdf'),
       p1, width=180, height = 120, units='mm', dpi = 600)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/整体UMAP_PDE.pdf'),
       p9, width=120, height = 100, units='mm', dpi = 600)
# 临床表型着色
#srt[['CancerSub']] = sapply(srt$subtype, function(x)ifelse(is.na(x),'Unknown', x))
p10 = DimPlot(srt, group.by = 'CancerSub')+scale_color_viridis_d()
p10
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/整体UMAP_Sub.pdf'),
       p10, width=120, height = 100, units='mm', dpi = 600)

length(table(srt$patient))
### PDE分组 ###
library(forcats)
srt$pde_group = ifelse(srt$pde > median(srt$pde,na.rm = T), 'High', 'Low')

patient_pde_order = srt@meta.data %>%
  group_by(patient, pde_group) %>%
  summarise(cell_num=n()) %>%
  left_join(srt@meta.data %>%
              group_by(patient) %>% summarise(all=n())) %>%
  mutate(prop=cell_num/all)
patient_pde_order = na.omit(patient_pde_order)

high_patient = patient_pde_order[patient_pde_order$pde_group=='High', 'patient']
other_patient = setdiff(patient_pde_order$patient, high_patient$patient)
patient_pde_order_low = patient_pde_order[patient_pde_order$pde_group=='Low', ] %>% arrange(prop) %>% as.data.frame()
patient_pde_order_low = patient_pde_order_low[patient_pde_order_low$patient%in%other_patient, ]
patient_pde_order = patient_pde_order[patient_pde_order$pde_group=='High', ] %>% arrange(prop) %>% as.data.frame()
patient_pde_order = rbind(patient_pde_order_low,patient_pde_order)

# clone 按照cee评分着色
a=srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, pde_group, clone) %>%
  summarise(cell_num=n(), pde=mean(pde)) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=patient, y=cell_num,fill=pde))+
  geom_bar(stat = 'identity', position = 'fill', width=0.9)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  #scale_fill_gradientn(values = )+
  scale_fill_gradientn(colours = c('#4DBBD5FF', '#E64B35FF'), limits=c(0,0.8), na.value = '#E64B35FF')+
  labs(y='cell%')+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        legend.position = 'right')+
  coord_flip()+
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
        legend.title = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/PDE高低分组病人比例分布2.pdf'),
       a, width=60, height = 100, units='mm', dpi = 600)
###
srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, pde_group) %>%
  summarise(cell_num=n()) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  mutate(pde_group=factor(pde_group, c('Low','High'))) %>%
  ggplot(aes(x=patient, y=cell_num,fill=pde_group))+
  geom_bar(stat = 'identity', position = 'fill', width=0.9)+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_y_continuous(breaks = c(0,0.5,1))+
  labs(y='cell%')+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        legend.position = 'left')+
  coord_flip()



a = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, pde_group) %>%
  summarise(cell_num=n()) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  mutate(pde_group=factor(pde_group, c('Low','High'))) %>%
  ggplot(aes(x=patient, y=cell_num,fill=pde_group))+
  geom_bar(stat = 'identity', position = 'fill', width=0.9)+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_y_continuous(breaks = c(0,0.5,1))+
  labs(y='cell%')+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        legend.position = 'left')+
  coord_flip()
a
library(ggridges)
gr = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=rearrange_score,y=patient, color=pde_group))+
  geom_density_ridges(fill=NA,scale = 0.8)

ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
# Pick the highest point. Could easily add quantiles or other features here.
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()


p_density = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=rearrange_score,y=patient, color=pde_group))+
  geom_density_ridges(fill=NA,scale = 0.8)+
  geom_segment(data = density_lines,
               aes(x = x, y = ymin, xend = x,
                   yend = ymin+0.8),#density*scale*iscale),
               color=density_lines$colour,
               linewidth=2,
               inherit.aes = F)+
  scale_color_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_x_continuous(breaks = c(0.0015), limits = c(0,0.003))+
  labs(x='RE')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
##
gr = srt@meta.data %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=BFB,y=patient, color=pde_group))+
  geom_density_ridges(fill=NA,scale = 0.8)

ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
# Pick the highest point. Could easily add quantiles or other features here.
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()

p_BFB = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, pde_group) %>%
  summarise(BFB=mean(BFB)) %>%
  mutate(pde_group=factor(pde_group, c('Low','High'))) %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=BFB,y=patient, fill=pde_group))+
  geom_bar(stat='identity',position = 'dodge')+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_x_continuous(breaks = c(0,0.015), limits = c(0,0.03))+
  labs(x='BFB')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
##

tmp_func <- function(x){
  # mean(x[x>0])
  sum(x>0) / length(x)
}
p_type = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, pde_group) %>%
  summarise(limit=tmp_func(limit_num), seismic=tmp_func(seismic_num), cell_num=n()) %>%
  reshape2::melt(id=c('patient', 'pde_group','cell_num'), variable.name='type', value.name = 'num') %>%
  mutate(patient=factor(patient, rev(patient_pde_order$patient))) %>%
  mutate(pde_group=factor(pde_group, c('Low','High'))) %>%
  ggplot(aes(x=num, y=pde_group, fill=type))+
  geom_bar(stat='identity', position = 'fill') +
  geom_text(aes(x=-0.1, y=pde_group, label=cell_num, color=pde_group), size=2)+
  scale_fill_manual(values=c('limit'='#DFEDD7','seismic'='#C9E2F4'))+
  scale_color_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  scale_x_continuous(expand = expansion(mult = 0.1, add = 0),breaks = c(0,0.5,1))+
  facet_wrap(~patient, ncol=1,strip.position='right')+
  labs(x='cell%')+
  theme_classic()+
  theme(strip.background = element_blank(),
        #strip.text = element_text(angle = 90),
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

pp=cowplot::plot_grid(a, p_density,p_BFB, p_type, align = 'h', rel_widths = c(0.4,0.15, 0.15, 0.3), nrow = 1)
pp
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/PDE高低分组病人比例分布.pdf'),
       pp, width=200, height = 200, units='mm', dpi = 600)

### PDE分组，比较不同指标 ###
VlnPlot(srt, features = c('nCount_RNA.x', 'nFeature_RNA.x', 'percent.mt'), pt.size = 0,ncol=1)
DimPlot(srt, label = T)+scale_color_igv()
srt[['new_group']] = paste0(srt$patient, '_', srt$clone)
srt@meta.data[, c('patient', 'pde_group',
                  'CellCycle1', 'cytotrace_score', 'rearrange_score',
                  'limit_score', 'seismic_score', 'BFB')] %>%
  reshape2::melt(id=c('patient', 'pde_group'), variable.name='index', value.name = 'score') %>%
  mutate(patient=factor(patient, patient_pde_order$patient)) %>%
  ggplot(aes(x=patient, y=score))+
  geom_boxplot()+
  facet_wrap(index~pde_group, scales = 'free_y')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
####
library(ggpubr)
tmp_meta = srt@meta.data[,c('new_group', 'pde_group', 'patient')] %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(new_group,pde_group) %>%
  summarise(n=n()) %>% as.data.frame()

tmp_order = reshape2::dcast(tmp_meta, new_group~pde_group)
tmp_order[is.na(tmp_order)] = 0
tmp_order$prop = tmp_order$High / (tmp_order$High+tmp_order$Low)
tmp_order$prop2 = tmp_order$Low / (tmp_order$High+tmp_order$Low)
tmp_order = tmp_order %>% arrange(prop)
tmp_meta %>%
  mutate(new_group = factor(new_group, tmp_order$new_group)) %>%
  ggplot(aes(x=new_group, y=n, fill=pde_group))+
  geom_bar(stat='identity', position='fill')+
  geom_hline(yintercept = 0.5)+
  theme(axis.text.x = element_blank())

clone_pde_prop = c()
tmp_meta_pde = c()
for(i in unique(tmp_meta$new_group)){
  tmp_aaa = tmp_meta[tmp_meta$new_group==i, ,drop=F] %>% arrange(-n)
  tmp_meta_pde = rbind(tmp_meta_pde, c(i, tmp_aaa[1,'pde_group']))
  if(nrow(tmp_aaa) == 2){
    tmp_fc = tmp_aaa[tmp_aaa$pde_group=='High', 'n'] / tmp_aaa[tmp_aaa$pde_group=='Low', 'n']
    clone_pde_prop = rbind(clone_pde_prop, c(i, tmp_fc))
  }else{
    if(tmp_aaa[1,'pde_group']=='High'){
      clone_pde_prop = rbind(clone_pde_prop, c(i, Inf))
    }else{
      clone_pde_prop = rbind(clone_pde_prop, c(i, 0))
    }
  }

}
tmp_meta_pde = as.data.frame(tmp_meta_pde)
rownames(tmp_meta_pde) = tmp_meta_pde[,1]
srt$clone_pde = tmp_meta_pde[srt$new_group, 2]
clone_cell_num =  srt@meta.data[,c('new_group', 'patient'),] %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(new_group) %>%
  summarise(n=n()) %>% as.data.frame()
rownames(clone_cell_num) = clone_cell_num$new_group
srt$clone_cell_num = clone_cell_num[srt$new_group, 2]

###
clone_pde_prop  = as.data.frame(clone_pde_prop)
clone_pde_prop$V2 = as.numeric(clone_pde_prop$V2)
clone_pde_prop = clone_pde_prop %>% arrange(V2)
clone_pde_prop$x = 1:nrow(clone_pde_prop)
clone_pde_prop[is.infinite(clone_pde_prop$V2), 'V2'] = NA
clone_pde_prop[is.na(clone_pde_prop$V2), 'V2'] = max(clone_pde_prop$V2, na.rm = T)

ggplot(clone_pde_prop, aes(x=x, y=V2))+
  geom_point()+
  geom_hline(yintercept = 1)
# 1
srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  #filter(!(seurat_clusters %in%c(17,18,31))) %>%
  group_by(patient, pde_group,file,  clone) %>%
  #summarise(score = mean(cytotrace_score)) %>%
  ggplot(aes(x=pde_group, y=cytotrace_score))+
  geom_violin(aes(fill=pde_group))+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  # scale_fill_npg()+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  stat_compare_means(comparisons = list(c('High', 'Low')), stat='p.signif')+
  labs(y='Cytotrace_score')+
  theme_classic()+
  theme(legend.position = 'none')

srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  #filter(clone_cell_num>5) %>%
  group_by(patient, clone_pde,file,  clone) %>%
  summarise(score = mean(cytotrace_score)) %>%
  ggplot(aes(x=clone_pde, y=score))+
  geom_violin(aes(fill=clone_pde))+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  geom_jitter()+
  # scale_fill_npg()+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  stat_compare_means(comparisons = list(c('High', 'Low')), stat='p.signif')+
  labs(y='Cytotrace_score')+
  theme_classic()+
  theme(legend.position = 'none')

a = srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(patient, clone_pde,file,  clone) %>%
  summarise(score = mean(cytotrace_score)) %>%
  ggplot(aes(y=clone_pde, x=score, fill=clone_pde))+
  geom_density_ridges(scale = 1,bandwidth=0.05 )+
  theme_classic()+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/PDE高低分组评分比较_clone(PDE).pdf'),
       a, width=100, height = 60, units='mm', dpi = 600)
#2
a=srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  #filter(clone_cell_num>100) %>%
  group_by(patient, clone_pde,file,  clone) %>%
  summarise(score = mean(CellCycle1)) %>%
  ggplot(aes(x=clone_pde, y=score))+
  # geom_violin(aes(fill=clone_pde))+
  geom_boxplot(aes(fill=clone_pde), width=0.5, outlier.shape = NA, size=0.2)+
  #geom_jitter(width=0.2)+
  # scale_fill_npg()+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  stat_compare_means(comparisons = list(c('High', 'Low')), label='p.signif')+
  labs(y='CellCycle_score')+
  scale_y_continuous(expand = expansion(mult=0.1))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.x = element_blank())
a
# 3
#gene_expression = FetchData(srt, vars = VariableFeatures(srt))
#res = apply(gene_expression, 2, function(y){
#   x = cut(y, breaks = 3, labels = F)
#   x = as.numeric(as.vector(x))
#   names(x) = names(y)
#   return(x)
# })
#cell_ITH = apply(res, 1, function(a)vegan::diversity(a[a!=1]))
# srt[['new_group']] = paste0(srt$patient, '_', srt$clone)
# group = 'new_group'
# gene_expression = FetchData(srt, vars = VariableFeatures(srt))
# gene_expression = GetAssayData(srt)
# gene_expression = t(as.matrix(gene_expression))
# # gene_expression = as.data.frame(srt@reductions$pca@cell.embeddings)
# gene_expression_center = as.data.frame(gene_expression)
# gene_expression_center[,group] = srt@meta.data[, group]
# gene_expression_center = gene_expression_center %>%
#   group_by(!!sym(group)) %>%
#   summarise_all(list(mean)) %>% as.data.frame()
# rownames(gene_expression_center) = gene_expression_center[,1]
# gene_expression_center = gene_expression_center[,-1]
#
# # euc_dist = as.matrix(dist(gene_expression, method = "euclidean"))
# # cell_ITH_dist = euc_dist[colnames(srt),colnames(srt)]
# # tmp_func <- function(x){
# #   aa=sd(x)/mean(x)
# #   #if(is.na(aa)){
# #   #  aa = 0
# #   #}
# #   return(aa)
# # }
# cell_ITH = c()
# tmp_meta = srt@meta.data[!is.na(srt@meta.data[, group]), ]
# for(i in unique(tmp_meta[,group])){
#   print(i)
#   tmp_bc = rownames(tmp_meta[tmp_meta[,group]==i,])
#   tmp_mat = gene_expression[tmp_bc, , drop=F]
#   tmp_center = unlist(gene_expression_center[i, ])
#   tmp_dist = apply(tmp_mat, 1, function(x){
#     sqrt(sum((x - tmp_center)^2))
#   })
#   cell_ITH_1 = rep(sd(tmp_dist)/mean(tmp_dist), length(tmp_dist))
#   names(cell_ITH_1) = names(tmp_dist)
#   # cell_ITH_1 = apply(cell_ITH_dist[tmp_bc, tmp_bc], 1, tmp_func)
#   cell_ITH = c(cell_ITH, cell_ITH_1)
# }


# 1.VariableFeatures genes expression
srt[['new_group']] = paste0(srt$patient, '_', srt$clone)
group = 'new_group'
gene_expression = FetchData(srt, vars = VariableFeatures(srt))

cell_ITH = c()
tmp_meta = srt@meta.data[!is.na(srt@meta.data[, group]), ]
for(i in unique(tmp_meta[,group])){
  print(i)
  tmp_bc = rownames(tmp_meta[tmp_meta[,group]==i,])
  tmp_mat = gene_expression[tmp_bc, , drop=F]
  #tmp_center = unlist(gene_expression_center[i, ])
  # tmp_dist = apply(tmp_mat, 1, function(x){
  #   sqrt(sum((x - tmp_center)^2))
  # })
  tmp_dist = as.matrix(dist(tmp_mat))
  #cell_ITH_1 = rep(sd(tmp_dist)/mean(tmp_dist), length(tmp_dist))
  # cell_ITH_1 = rep(sd(tmp_dist), length(tmp_dist))
  # names(cell_ITH_1) = names(tmp_dist)
  cell_ITH_1 = apply(tmp_dist, 1, function(x)sd(x)) # #*mean(x))
  # cell_ITH_1 = apply(cell_ITH_dist[tmp_bc, tmp_bc], 1, tmp_func)
  cell_ITH = c(cell_ITH, cell_ITH_1)
}
#cell_ITH = as.data.frame(cell_ITH)
srt[['Cell_ITH']] = cell_ITH[colnames(srt)]
b=srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  # filter(clone_cell_num>100) %>%
  group_by(patient, clone_pde,file,  clone) %>%
  summarise(score = median(Cell_ITH)) %>%
  ggplot(aes(x=clone_pde, y=score))+
  # geom_violin(aes(fill=clone_pde))+
  geom_boxplot(aes(fill=clone_pde), width=0.5, outlier.shape = NA, size=0.2)+
  #geom_jitter(width=0.2)+
  # scale_fill_npg()+
  scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
  stat_compare_means(comparisons = list(c('High', 'Low')), label='p.signif')+
  labs(y='Cell_ITH')+
  scale_y_continuous(expand = expansion(mult=0.1))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.x = element_blank())
a+b
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/PDE高低分组评分比较_clone.pdf'),
       a+b, width=80, height = 60, units='mm', dpi = 600)

# 4 指标在高低分布
a=srt@meta.data %>%
  filter(patient%in%patient_pde_order$patient) %>%
  #filter(clone_cell_num>100) %>%
  group_by(clone_pde, CancerSub, patient) %>%
  summarise(cell_num = n()) %>%
  group_by(clone_pde, CancerSub) %>%
  summarise(cell_num = n()) %>%

  ggplot(aes(x=clone_pde, y=cell_num, fill=CancerSub))+
  geom_bar(stat='identity', position = 'fill', width=0.6)+
  scale_fill_viridis_d()+
  theme_classic()+
  theme(axis.title.x = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/PDE高低分组比较亚型.pdf'),
       a, width=60, height = 60, units='mm', dpi = 600)


### monocle ######
library(monocle)
run_monocle2 <- function(data, metadata, ordering_genes=NULL,FormulaStr=NULL){
  pd <- new('AnnotatedDataFrame', data = metadata)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  #Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())

  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  if(is.null(ordering_genes)){
    ### detectGenes计算每一个细胞表达的基因个数.
    monocle_cds <- detectGenes(monocle_cds, min_expr = 0.01)
    ### 保留在10个细胞中表达的基因，当做表达基因.
    expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 1))
    disp_table <- dispersionTable(monocle_cds)
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01&dispersion_empirical>0.5)
    expressed_genes = intersect(expressed_genes, unsup_clustering_genes$gene_id)
    if(is.null(FormulaStr)){
      ordering_genes <- expressed_genes#row.names (subset(diff_test_res, qval < 0.05)) ##
    }else{
      diff_test_res <- differentialGeneTest(monocle_cds, fullModelFormulaStr = paste0("~", FormulaStr))
      ordering_genes <- row.names (subset(diff_test_res, qval < 0.05)) ##
    }

  }
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  plot_ordering_genes(monocle_cds)

  ## dimension reduciton
  monocle_cds <- reduceDimension(monocle_cds, method = 'DDRTree',max_compenents=1)
  ## ordering cells
  monocle_cds <- orderCells(monocle_cds)
  return(monocle_cds)
}
srt[['new_group']] = paste0(srt$patient, '_', srt$clone)
srt@meta.data$Cell_ITH[is.na(srt@meta.data$Cell_ITH)] = 0
tmp_meta = srt@meta.data[, c('new_group','Cell_ITH',
                             'CellCycle1', 'cytotrace_score', 'rearrange_score','pde',
                             'limit_score', 'seismic_score', 'BFB')] %>%
  group_by(new_group) %>%
  summarise_all(list(mean)) %>%
  as.data.frame()

tmp_meta$patient = sapply(tmp_meta$new_group, function(x)strsplit(x, '_')[[1]][1])
rownames(tmp_meta) = tmp_meta$new_group
tmp_meta$clone_group = tmp_meta_pde[tmp_meta$new_group, 2]

tmp_data = AverageExpression(srt, group.by = 'new_group')$RNA
#FC_cds = run_monocle2(data=as(as.matrix(tmp_data), 'sparseMatrix'),
#                      metadata= tmp_meta,
#                      #ordering_genes = signif_genes,
#                      FormulaStr = NULL) # 'Tumor_metastasis')
x = apply(tmp_data, 1, var)



FC_cds = run_monocle2(data=as(as.matrix(tmp_data), 'sparseMatrix'),
                      metadata= tmp_meta,
                      ordering_genes = VariableFeatures(srt),
                      FormulaStr = NULL)

#FC_cds2 = reduceDimension(FC_cds, reduction_method = 'DDRTree',ncenter= 500)

#plot_cell_trajectory(FC_cds2, color_by = "State", cell_size = 1)


saveRDS(FC_cds, paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_samples_monocle2.rds'))
FC_cds = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/BC_samples_monocle2.rds')

# pData(FC_cds)$CancerSub = srt_info[pData(FC_cds)$patient, 'type']
# 不同属性展示
#scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
a = plot_cell_trajectory(FC_cds, color_by = "pde", cell_size = 1)+theme(legend.position = 'right') +
  scale_color_gradientn(colours = c('#4DBBD5FF', '#E64B35FF'))#, limits=c(0.2,0.6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_CEE.pdf'),
       a, width=120, height = 80, units='mm', dpi = 600)
b = plot_cell_trajectory(FC_cds, color_by = "cytotrace_score", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
c = plot_cell_trajectory(FC_cds, color_by = "rearrange_score", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
d = plot_cell_trajectory(FC_cds, color_by = "Cell_ITH", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
e = plot_cell_trajectory(FC_cds, color_by = "CellCycle1", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
f = plot_cell_trajectory(FC_cds, color_by = "State", cell_size = 0.5)+theme(legend.position = 'right')+ scale_color_d3()
a+b+c+d+e+f


a = plot_cell_trajectory(FC_cds, color_by = "pde", cell_size = 1)+
  theme(legend.position = 'right')+
  scale_color_gradientn(colours = c('gray', 'red'), limits=c(0,0.3))
b = plot_cell_trajectory(FC_cds, color_by = "cytotrace_score", cell_size = 1)+theme(legend.position = 'right')+
  scale_color_gradientn(colours = c('gray', 'red'))
c = plot_cell_trajectory(FC_cds, color_by = "Cell_ITH", cell_size = 1)+theme(legend.position = 'right')+
  scale_color_gradientn(colours = c('gray', 'red'))

a+b+c
plot_cell_trajectory(FC_cds, color_by = "clone_group", cell_size = 1)

f = plot_cell_trajectory(FC_cds, color_by = "State", cell_size = 0.2,cell_link_size=0.3)+theme(legend.position = 'right')+ scale_color_d3()
f = f + theme_classic()+
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
f
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_state.pdf'),
       f, width=70, height = 50, units='mm', dpi = 600)

g=plot_cell_trajectory(FC_cds, color_by = "CancerSub", cell_size = 1)+theme(legend.position = 'right')+ scale_color_viridis_d()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_Sub.pdf'),
       g, width=120, height = 80, units='mm', dpi = 600)

g=plot_cell_trajectory(FC_cds, color_by = "patient", cell_size = 1)+theme(legend.position = 'right')+ scale_color_igv()
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_patient.pdf'),
       g, width=180, height = 80, units='mm', dpi = 600)

# 排序细胞，root=4
FC_cds = orderCells(FC_cds, root_state = 1)
plot_cell_trajectory(FC_cds, color_by = "Pseudotime", cell_size = 1)


### 随着拟时序指标变化情况 #####
cds_meta = pData(FC_cds)

bh5 = cds_meta[cds_meta$State!='2', ]
bh1 = cds_meta[cds_meta$State!='3', ]
# pre_bh = cds_meta[cds_meta$State%in%c('4', '2','3'), ]

bh5$branch = 'branch5'
bh1$branch = 'branch1'
# pre_bh$branch = 'pre_branch'


cds_branch = as.data.frame(rbind(bh5, bh1))
cds_branch$Cell_ITH = tmp_meta[cds_branch$new_group, 'Cell_ITH']
cds_branch$pde = tmp_meta[cds_branch$new_group, 'pde']
cds_branch$cytotrace_score = tmp_meta[cds_branch$new_group, 'cytotrace_score']

# cut_num = 3
# cds_branch$Pseudotime_norm = (cds_branch$Pseudotime-min(cds_branch$Pseudotime)) / (max(cds_branch$Pseudotime)- min(cds_branch$Pseudotime))
# # cds_branch$Pseudotime_norm = cut(cds_branch$Pseudotime_norm, breaks = cut_num, labels = (0:(cut_num-1))/(cut_num-1))
# cds_branch$Pseudotime_norm = cut(cds_branch$Pseudotime_norm,
#                                  breaks =3,
#                                  labels = c(0,0.5,1))
#
# #chr_prop_mat$pseudotime = as.numeric(as.vector(chr_prop_mat$pseudotime))

cut_num = 7
cds_branch$Pseudotime_norm = (cds_branch$Pseudotime-min(cds_branch$Pseudotime)) / (max(cds_branch$Pseudotime)- min(cds_branch$Pseudotime))
cds_branch$Pseudotime_norm = cut(cds_branch$Pseudotime_norm, breaks = cut_num, labels = (0:(cut_num-1))/(cut_num-1))

cds_branch$Pseudotime_norm = as.numeric(as.vector(cds_branch$Pseudotime_norm))


hline_pos = max(as.numeric(as.vector(cds_branch[intersect(bh5$new_group, bh1$new_group), 'Pseudotime_norm'])), na.rm = T)
hline_pos = 0.5
a = cds_branch %>%
  group_by(Pseudotime_norm,branch) %>%
  summarise(pde = mean(pde, na.rm=T),
            Cell_ITH = mean(Cell_ITH, na.rm=T),
            cytotrace = mean(cytotrace_score, na.rm=T),
            #limit_score = mean(limit_score, na.rm=T),
            #seismic_score = mean(seismic_score, na.rm=T)
            #rearrange_score = mean(rearrange_score, na.rm=T)
  ) %>%
  reshape2::melt(id=c('Pseudotime_norm', 'branch'), variable.name='type', value.name = 'score') %>%
  mutate(Pseudotime_norm = as.numeric(as.vector(Pseudotime_norm))) %>%
  #tidyr::fill(score, .direction ='up') %>%
  ggplot(aes(x=Pseudotime_norm, y=score, color=branch))+
  geom_smooth(method = 'loess', formula = y ~ x, se=F, span=1)+
  geom_vline(xintercept = hline_pos, linetype='dashed', linewidth=0.5)+
  scale_color_manual(values = c('branch1'='#1F77B4FF', 'branch5'='#9467BDFF'))+
  facet_wrap(~type, scales = 'free_y', ncol=1, strip.position = 'right')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank())+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
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
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/分支曲线_nomean.pdf'),
       a, width=70, height = 100, units='mm', dpi = 600)

# 统计比例
pie_data = pData(FC_cds)[, c('State', 'clone_group')]
pie_list = list()
pie_prop = c()
for(num in 1:5){
  tmp=pie_data[pie_data$State==num,] %>%
    na.omit() %>%
    group_by(clone_group) %>%
    dplyr::summarise(value=n()) %>%
    mutate(labs=scales::percent(value/sum(value)),
           prop = value/sum(value))

  a = tmp %>% ggpie(x='value', label='labs',
                    fill='clone_group',color = "white",
                    lab.pos = "in", lab.font = "white")+
    scale_fill_manual(values=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))+
    theme(legend.position = 'none')
  pie_list[[num]] = a
  tmp$State = num
  pie_prop = rbind(pie_prop, tmp)
}
a = cowplot::plot_grid(plotlist = pie_list, nrow=1)
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/State_pie.pdf'),
       a, width=200, height = 60, units='mm', dpi = 600)

pie_prop = as.data.frame(pie_prop)

a = pie_prop %>%
  filter(clone_group=='High') %>%
  mutate(State = factor(State, c(4,3,2,1,5))) %>%
  ggplot(aes(x=State, y=prop, group=1, fill='a'))+
  geom_line()+
  geom_point(color='#E64B35FF', size=4)+
  labs(y='High_clone_prop')+

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
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/State_High_prop.pdf'),
       a, width=70, height = 50, units='mm', dpi = 600)

#### 比较stat5和state1差异 #####
tmp_data = pData(FC_cds)[,c('State', 'pde', 'Cell_ITH','cytotrace_score','CellCycle1')]

tmp_data = cds_branch[,c('State', 'pde', 'Cell_ITH','cytotrace_score','CellCycle1', 'limit_score', 'seismic_score')]
a = tmp_data %>%
  filter(State%in%c('4', '5')) %>%
  reshape2::melt(id='State', variable.name='type', value.name = 'score') %>%
  ggplot(aes(x=State,y=score, color=State))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, size=0.1)+
  facet_wrap(~type, scales = 'free_y', ncol = 2, strip.position = 'right')+
  scale_color_manual(values = c('4'='#1F77B4FF', '5'='#9467BDFF'))+
  stat_compare_means(comparisons = list(c('5','4')) , label='p.format')+
  scale_y_continuous(expand = expansion(mult = 0.2))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/分支boxplot.pdf'),
       a, width=140, height = 100, units='mm', dpi = 600)

a = tmp_data %>%
  filter(State%in%c('1', '5')) %>%
  reshape2::melt(id='State', variable.name='type', value.name = 'score') %>%
  filter(type=='CellCycle1') %>%
  ggplot(aes(x=State,y=score, color=State))+
  geom_boxplot(outlier.shape = NA, width=0.5, size=0.2)+
  geom_jitter(width=0.1, size=0.1)+
  #facet_wrap(~type, scales = 'free_y', ncol = 2, strip.position = 'right')+
  scale_color_manual(values = c('1'='#1F77B4FF', '5'='#9467BDFF'))+
  stat_compare_means(comparisons = list(c('5','1')), size=2)+
  #scale_y_continuous(expand = expansion(mult = 0.2))+
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
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/分支boxplot_cc.pdf'),
       a, width=50, height = 50, units='mm', dpi = 600)

## 1,5 亚型构成
a=pData(FC_cds) %>%
  filter(State%in%c('1', '5')) %>%
  filter(patient%in%patient_pde_order$patient) %>%
  group_by(State, CancerSub, patient) %>%
  summarise(cell_num = n()) %>%
  group_by(State,CancerSub) %>%
  summarise(cell_num = n()) %>%

  group_by(State, CancerSub) %>%
  summarise(cell_num = n()) %>%
  ggplot(aes(x=State, y=cell_num, fill=CancerSub))+
  geom_bar(stat='identity', position = 'fill', width=0.6)+
  scale_fill_viridis_d()+
  theme_classic()+
  theme(axis.title.x = element_blank())
a

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/State1v5组比较亚型.pdf'),
       a, width=60, height = 60, units='mm', dpi = 600)
#### TNBC识别区域基因在state表达情况 ####
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')
del_10q23 = subset(gene_info,chrom=='10'&start>=89597154&end<=90565599)
amp_17q25 = subset(gene_info,chrom=='17'&start>=79501785&end<=81195210)
amp_9p23 = subset(gene_info,chrom=='9'&start>=12652132&end<=14381484)



my_genes <- FC_cds[row.names(subset(fData(FC_cds), gene_short_name %in% del_10q23$gene_name)),]
plot_genes_violin(my_genes, grouping = 'State', ncol=2)

my_genes <- FC_cds[row.names(subset(fData(FC_cds), gene_short_name %in% amp_17q25$gene_name)),]
plot_genes_violin(my_genes, grouping = 'State', ncol=4)

my_genes <- FC_cds[row.names(subset(fData(FC_cds), gene_short_name %in% amp_9p23$gene_name)),]
plot_genes_violin(my_genes, grouping = 'State', ncol=2)

my_genes = c(del_10q23$gene_name, amp_17q25$gene_name, amp_9p23$gene_name)
my_genes <- FC_cds[row.names(subset(fData(FC_cds), gene_short_name %in% my_genes)),]

cds_exprs = exprs(my_genes)
cds_exprs = Matrix::t(Matrix::t(cds_exprs)/sizeFactors(my_genes))
cds_exprs = reshape2::melt(as.matrix(cds_exprs))
colnames(cds_exprs) = c("f_id", "Cell", "expression")
cds_pData = pData(my_genes)
cds_exprs = merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
cds_exprs = merge(cds_exprs, gene_info, by.x = "f_id", by.y = "gene_name")


keep_gene = cds_exprs %>%
  filter(State%in%c(1,5)) %>%
  group_by(State, f_id) %>%
  summarise(expression=mean(expression)) %>%
  reshape2::dcast(f_id~State) %>%
  filter(`1`>`5`)

driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因_Cosmic.csv', sep=',', header = T, check.names = F)
driver_gene = driver_gene$`Gene Symbol`
same_gene = intersect(driver_gene, keep_gene$f_id)

cds_exprs %>%
  filter(State%in%c(1,5)) %>%
  filter(f_id%in%keep_gene$f_id) %>%
  ggplot(aes(x=State, y=expression))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('1', '5')))+
  facet_wrap(~f_id+chrom, scales = 'free', ncol=6)

a=cds_exprs %>%
  filter(State%in%c(1,5)) %>%
  filter(f_id%in%same_gene) %>%
  ggplot(aes(x=State, y=expression, color=State))+
  geom_boxplot(outlier.shape = NA, width=0.5, size=0.2)+
  geom_jitter(width=0.1, size=0.1)+
  #facet_wrap(~type, scales = 'free_y', ncol = 2, strip.position = 'right')+
  scale_color_manual(values = c('1'='#1F77B4FF', '5'='#9467BDFF'))+
  stat_compare_means(comparisons = list(c('5','1')), size=2)+
  #scale_y_continuous(expand = expansion(mult = 0.2))+
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
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/State1v5组比较PTEN.pdf'),
       a, width=50, height = 50, units='mm', dpi = 600)


#### 比较stat5和state1染色体碎裂构成 ####
tmp_data = pData(FC_cds)
tmp_srt_meta = srt@meta.data
tmp_srt_meta$State = tmp_data[match(tmp_srt_meta$new_group, tmp_data$new_group),'State']

# 发生事件的细胞占比
a = tmp_srt_meta %>%
  filter(State%in%c('1', '5')) %>%
  group_by(State) %>%
  #summarise(CPS_num=mean(CPS_num), limit_num=mean(limit_num), seismic_num=mean(seismic_num)) %>%
  summarise(CPS_num=sum(CPS_num>0)/length(CPS_num),
            limit_num=sum(limit_num>0)/length(limit_num),
            seismic_num=sum(seismic_num>0)/length(seismic_num)) %>%

  reshape2::melt(id='State', value.name = 'cell_prop', variable.name='type') %>%
  ggplot(aes(x=State, y=cell_prop, fill=State))+
  geom_bar(stat='identity', position = 'dodge', width=0.5)+
  facet_wrap(~type, scales = 'free_y')+
  scale_fill_manual(values = c('1'='#1F77B4FF', '5'='#9467BDFF'))+
  theme_classic()+
  theme(strip.background = element_blank())
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/1v5_chr_num(01).pdf'),
       a, width=110, height = 60, units='mm', dpi = 600)

# seismic细胞，cytotrace评分更低
#tmp_srt_meta %>%
#  filter(State%in%c('1')) %>%
#  mutate(seismic_group = ifelse(seismic_num>0, 'S', 'NS')) %>%
#  group_by(new_group)+
#  ggplot(aes(x=seismic_group, y=cytotrace_score)) +
#  geom_boxplot()

box_d = tmp_srt_meta %>% filter(State%in%c('1'))
line_d = tmp_srt_meta %>% filter(State%in%c('1')) %>%
  group_by(seismic_num) %>%
  summarise(cytotrace_score = mean(cytotrace_score))
a = box_d %>%
  ggplot(aes(x=seismic_num, y=cytotrace_score, group=seismic_num, fill=as.factor(seismic_num))) +
  geom_boxplot(width=0.5,size=0.2)+
  geom_smooth(data=line_d, mapping = aes(x=seismic_num, y=cytotrace_score),
              method = 'lm', formula = 'y ~ x',se=F, inherit.aes = F, color='black')+
  stat_cor(data=line_d, mapping = aes(x=seismic_num, y=cytotrace_score),inherit.aes = F, color='red')+
  scale_fill_brewer(palette="Blues", name='Seismic')+
  theme_classic()
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/Seismic_cyto_cor.pdf'),
       a, width=80, height = 60, units='mm', dpi = 600)

box_d = tmp_srt_meta %>% filter(State%in%c('5'))
line_d = tmp_srt_meta %>% filter(State%in%c('5')) %>%
  group_by(limit_num) %>%
  summarise(cytotrace_score = mean(cytotrace_score))
a = box_d %>%
  ggplot(aes(x=limit_num, y=cytotrace_score, group=limit_num)) +
  geom_boxplot(width=0.5)+
  geom_smooth(data=line_d, mapping = aes(x=limit_num, y=cytotrace_score),
              method = 'lm', formula = 'y ~ x',se=F, inherit.aes = F)+
  stat_cor(data=line_d, mapping = aes(x=limit_num, y=cytotrace_score),inherit.aes = F, color='red')+
  theme_classic()
a

## 1,5两组基因组变异事件差异 ###
all_cell_events = c()
for(prefix in names(BC_list)){
  print(prefix)
  sctc = BC_list[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('virtual|root', rownames(tmp_cna)),]
  tmp_cna = (tmp_cna!=2)
  tmp_cna = data.frame(event_num = rowSums(tmp_cna), all_pos=ncol(tmp_cna))
  all_cell_events = rbind(all_cell_events, tmp_cna)
}

all_cell_events = as.data.frame(all_cell_events)
all_cell_events$new_group = srt@meta.data[rownames(all_cell_events), 'new_group']
tmp_data = pData(FC_cds)
all_cell_events$State = tmp_data[match(all_cell_events$new_group, tmp_data$new_group),'State']

a=all_cell_events %>%
  filter(State%in%c('1', '5')) %>%
  mutate(event_num=event_num/all_pos) %>%
  group_by(State, new_group) %>%
  summarise(mean_events = mean(event_num)) %>%
  ggplot(aes(x=State, y=mean_events, color=State))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, size=0.1)+
  scale_color_manual(values = c('1'='#1F77B4FF', '5'='#9467BDFF'))+
  stat_compare_means(comparisons = list(c('5','1')) , label='p.signif')+
  scale_y_continuous(expand = expansion(mult = 0.2))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/State1v5_Events.pdf'),
       a, width=60, height = 60, units='mm', dpi = 600)

#### 推断Limit和seismic的先后关系 ####
# 染色体
cps_mat = c()
limit_mat = c()
seismic_mat = c()
names_meta = c()
for(prefix in names(BC_list)){
  print(prefix)
  sctc = BC_list[[prefix]]
  CNA_mechnism = read.table(glue('{scTrace_dir3}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score = read.table(glue('{scTrace_dir3}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir3}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir3}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  pvalue_thr = 0.001
  CPS_num = all_rearrange_score_pvalue<pvalue_thr& all_rearrange_score>0
  limit_num = all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5& all_rearrange_score>0
  seismic_num = all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5& all_rearrange_score>0

  tmp_meta = srt@meta.data[srt$file==prefix, ]
  for(cl in unique(tmp_meta$new_group)){
    print(cl)
    tmp_bc = rownames(tmp_meta[tmp_meta$new_group==cl, ])
    cps_mat = rbind(cps_mat, colSums(CPS_num[tmp_bc,,drop=F])/nrow(CPS_num[tmp_bc,,drop=F]))
    limit_mat = rbind(limit_mat, colSums(limit_num[tmp_bc,,drop=F])/nrow(limit_num[tmp_bc,,drop=F]))
    seismic_mat = rbind(seismic_mat, colSums(seismic_num[tmp_bc,,drop=F])/nrow(seismic_num[tmp_bc,,drop=F]))
    names_meta = rbind(names_meta, c(prefix, cl))
  }
}
rownames(seismic_mat) = names_meta[,2]
rownames(limit_mat) = names_meta[,2]

#seismic_mat2 = apply(seismic_mat, 2, scale)
#limit_mat2 = apply(limit_mat, 2, scale)
limit1 = reshape2::melt(limit_mat) %>% mutate(type='Limit')
seismic1 = reshape2::melt(seismic_mat) %>% mutate(type='Seismic')


chr_prop_mat = as.data.frame(rbind(limit1, seismic1))
colnames(chr_prop_mat) = c('new_group', 'chr', 'prop', 'type')
chr_prop_mat$new_group = as.vector(chr_prop_mat$new_group)

tmp_cds_meta = pData(FC_cds)
branch1_bc = rownames(tmp_cds_meta[tmp_cds_meta$State!=5,])
branch5_bc = rownames(tmp_cds_meta[tmp_cds_meta$State!=1,])
branch1_mat = chr_prop_mat[chr_prop_mat$new_group%in%branch1_bc, ]
branch5_mat = chr_prop_mat[chr_prop_mat$new_group%in%branch5_bc, ]
branch1_mat$branch = 'branch1'
branch5_mat$branch = 'branch5'

chr_prop_mat = rbind(branch1_mat, branch5_mat)
# chr_prop_mat$pseudotime = as.numeric(as.vector(cds_branch[chr_prop_mat$new_group,'Pseudotime_norm']))
chr_prop_mat$pseudotime = as.numeric(as.vector(cds_branch[chr_prop_mat$new_group,'Pseudotime']))

cut_num = 6
chr_prop_mat$pseudotime = (chr_prop_mat$pseudotime-min(chr_prop_mat$pseudotime)) / (max(chr_prop_mat$pseudotime)- min(chr_prop_mat$pseudotime))
chr_prop_mat$pseudotime = cut(chr_prop_mat$pseudotime, breaks = cut_num, labels = (0:(cut_num-1))/(cut_num-1))
chr_prop_mat$pseudotime = as.numeric(as.vector(chr_prop_mat$pseudotime))

#cut_num = 6
#cds_branch$Pseudotime_norm = (cds_branch$Pseudotime-min(cds_branch$Pseudotime)) / (max(cds_branch$Pseudotime)- min(cds_branch$Pseudotime))
#cds_branch$Pseudotime_norm = cut(cds_branch$Pseudotime_norm, breaks = cut_num, labels = (0:(cut_num-1))/(cut_num-1))

# a=chr_prop_mat %>%filter(type=='Seismic', chr=='X1', pseudotime==0.1)
# pData(FC_cds)$test = 'aa'
# pData(FC_cds)[intersect(branch1_bc, branch5_bc), 'test']='select'
# plot_cell_trajectory(FC_cds, color_by = "test", cell_size = 1)+theme(legend.position = 'right')+ scale_color_igv()


a=chr_prop_mat%>%#[chr_prop_mat$new_group%in%rownames(cds_branch[cds_branch$State%in%c(1,5),]),] %>%
  mutate(chr=gsub('X', 'chr', chr)) %>%
  group_by(chr, pseudotime, branch,type) %>%
  summarise(prop = mean(prop)) %>%
  filter(type=='Seismic') %>% as.data.frame() %>%
  mutate(chr=factor(chr, paste0('chr', 1:22))) %>%
  ggplot(aes(x=pseudotime, y=prop, color=branch))+
  #geom_line(size=1)+
  #geom_smooth(method = 'glm', formula = y~poly(x,2),se=F, size=1)+
  geom_smooth(method = 'loess',se=F, size=1)+
  facet_grid(~chr, scales = 'free_y')+
  scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1))+
  theme_bw()+
  scale_color_manual(values = c('branch1'='#1F77B4FF', 'branch5'='#9467BDFF'))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing=unit(0, "lines"),
        plot.margin = unit(c(0,0.1,0,0.1), "lines"))
a
# cor_res = c()
# for(i in 1:ncol(limit_mat)){
#   cor_res = rbind(cor_res, cor(limit_mat[,i], seismic_mat[,i]))# method='spearman'))
# }
# cor_res = as.data.frame(cor_res)
# cor_res$chr = colnames(limit_mat)
# cor_res$x=1:nrow(cor_res)
# ggplot(cor_res, aes(x=chr, y=V1))+
#   geom_point()
#
# #
# cds_meta = pData(FC_cds)
# cds_meta$chr_prop = limit_mat[rownames(cds_meta), 'X3']
#
# cds_meta%>%
#   ggplot(aes(x=Pseudotime, y=chr_prop))+
#   #geom_point(aes(color=State))+
#   geom_smooth(method = 'loess', formula = y ~ x, se=F)+
#   scale_color_d3()


### 具体发生的染色体
scTrace_dir3 = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/BC_res/patient_cnv_infer3/'
cps_mat = c()
limit_mat = c()
seismic_mat = c()
names_meta = c()
for(prefix in names(BC_list)){
  print(prefix)
  sctc = BC_list[[prefix]]
  CNA_mechnism = read.table(glue('{scTrace_dir3}/{prefix}mode.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score = read.table(glue('{scTrace_dir3}/{prefix}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{scTrace_dir3}/{prefix}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{scTrace_dir3}/{prefix}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  pvalue_thr = 0.001
  CPS_num = all_rearrange_score_pvalue<pvalue_thr& all_rearrange_score>0
  limit_num = all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5& all_rearrange_score>0
  seismic_num = all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5& all_rearrange_score>0

  tmp_meta = srt@meta.data[srt$file==prefix, ]
  for(cl in unique(tmp_meta$new_group)){
    print(cl)
    tmp_bc = rownames(tmp_meta[tmp_meta$new_group==cl, ])

    cps_mat = rbind(cps_mat, colSums(CPS_num[tmp_bc,,drop=F])/nrow(CPS_num[tmp_bc,,drop=F]))
    limit_mat = rbind(limit_mat, colSums(limit_num[tmp_bc,,drop=F])/nrow(limit_num[tmp_bc,,drop=F]))
    seismic_mat = rbind(seismic_mat, colSums(seismic_num[tmp_bc,,drop=F])/nrow(seismic_num[tmp_bc,,drop=F]))
    names_meta = rbind(names_meta, c(prefix, cl))
  }

}
rownames(seismic_mat) = names_meta[,2]
rownames(limit_mat) = names_meta[,2]

names_meta = as.data.frame(names_meta)
rownames(names_meta) = names_meta[,2]
# state1, state5
meta_1_5 = tmp_srt_meta %>% filter(State%in%c('1','5'))

ht_mat = seismic_mat[unique(meta_1_5$new_group), ]
ht_meta = names_meta[unique(meta_1_5$new_group), ]
ht_meta$State = pData(FC_cds)[ht_meta$V2, 'State']
ht_meta = ht_meta %>% arrange(State)
ht_mat = ht_mat[ht_meta$V2, ]

top_data = cbind('group'=ht_meta[,3], as.data.frame(ht_mat)) %>%
  as.data.frame() %>%
  group_by(group) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(top_data) = top_data[,1]
top_data = top_data[,-1]

b = top_data %>%
  as.matrix() %>%
  reshape2::melt() %>%
  mutate(Var2=gsub('X', 'chr', Var2)) %>%
  mutate(Var1=ifelse(Var1=='1', 'branch1','branch5')) %>%
  mutate(Var2=factor(Var2, paste0('chr', 1:22))) %>%
  ggplot(aes(x=Var1, y=value, fill=Var1))+
  geom_bar(stat='identity', width=0.6)+
  facet_grid(~Var2, scales = 'free_y')+
  theme_bw()+
  scale_fill_manual(values = c('branch1'='#1F77B4FF', 'branch5'='#9467BDFF'))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.switch.pad.wrap = unit(0, 'pt'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,0.1,0,0.1), "lines"),
        panel.spacing=unit(0, "lines"),
        #plot.background = element_rect(fill='red')
  )
#plot.margin = margin(0, 0, 0, 0, "cm")
b
c=cowplot::plot_grid(b, a, align = 'v', ncol=1, rel_heights = c(0.4,0.6))

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/分支曲线染色体.pdf'),
       c, width=400, height = 60, units='mm', dpi = 600)

#ggarrange(b, a, ncol=1, align = 'v')

top_ann = HeatmapAnnotation(Chr = anno_barplot(t(top_data),
                                               gp = gpar(fill =c('1'='#1F77B4FF', '5'='#9467BDFF')),
                                               beside=T,
                                               height = unit(2, "cm")))
colnames(ht_mat) = paste0('chr',1:22)

ht_meta$patient = sapply(ht_meta$V2, function(x)strsplit(x,'_')[[1]][1])
patient_col = setNames(pal_igv()(length(unique(ht_meta$patient))), unique(ht_meta$patient))
Heatmap(ht_mat[, paste0('chr',1:22)],
        cluster_rows = F, cluster_columns = F,
        show_row_names = F,
        row_split = ht_meta[,3],
        column_split = factor(paste0('chr',1:22),paste0('chr',1:22)),
        col=circlize::colorRamp2(c(0,  0.25), c("white", "#08AEEA")),
        left_annotation = rowAnnotation(df=as.data.frame(ht_meta[,3:4,drop=F]),
                                        col=list('State'=c('1'='#1F77B4FF', '5'='#9467BDFF'),
                                                 'patient'=patient_col)),
        #top_annotation = top_ann,
        column_title = NULL,
        row_title = NULL,
        show_row_dend = F
)

# 统计8号染色体seismic病人分布
chr8_seismic = reshape2::melt(as.matrix(ht_mat)) ## data.frame(prop=ht_mat[, 'chr8'])
chr8_seismic[, c('State', 'patient')] = ht_meta[chr8_seismic$Var1, c('State', 'patient')]
chr8_seismic %>%
  filter(Var2=='chr8') %>%
  group_by(Var2,State, patient) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(x=State, y=value))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  stat_compare_means(comparisons = list(c('1', '5')), label='p.format')+
  theme_classic()

chr8_seismic %>%
  filter(Var2=='chr8') %>%
  group_by(Var2,State, patient) %>%
  summarise(value=mean(value)) %>%
  group_by(Var2,State, patient) %>%
  summarise(value=ifelse(value>0, 'Seismic', 'Non')) %>%
  ggplot(aes(x=State, fill=value))+
  geom_bar(stat='count', position = 'fill')+
  scale_fill_manual(values = c('Seismic'='#E64B35FF', 'Non'='gray'))+
  #stat_compare_means(comparisons = list(c('1', '5')), label='p.format')+
  theme_classic()
#coord_polar(theta = 'y')


# limit
ht_mat = limit_mat[unique(meta_1_5$new_group), ]
ht_meta = names_meta[unique(meta_1_5$new_group), ]
ht_meta$State = pData(FC_cds)[ht_meta$V2, 'State']
ht_meta = ht_meta %>% arrange(State)
ht_mat = ht_mat[ht_meta$V2, ]
top_data = cbind('group'=ht_meta[,3], as.data.frame(ht_mat)) %>%
  as.data.frame() %>%
  group_by(group) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(top_data) = top_data[,1]
top_data = top_data[,-1]
top_ann = HeatmapAnnotation(Chr = anno_barplot(t(top_data),
                                               gp = gpar(fill =c('1'='#1F77B4FF', '5'='#9467BDFF')),
                                               beside=T,
                                               height = unit(2, "cm")))
colnames(ht_mat) = paste0('chr',1:22)
Heatmap(ht_mat[, paste0('chr',1:22)],
        cluster_rows = F, cluster_columns = F,
        show_row_names = F,
        row_split = ht_meta[,3],
        column_split = factor(paste0('chr',1:22),paste0('chr',1:22)),
        col=circlize::colorRamp2(c(0,  1), c("white", "#08AEEA")),
        left_annotation = rowAnnotation(df=as.data.frame(ht_meta[,3,drop=F]),
                                        col=list('State'=c('1'='#1F77B4FF', '5'='#9467BDFF') )),
        top_annotation = top_ann,
        column_title = NULL,
        row_title = NULL,
        show_row_dend = F
)
#### 17号染色体展示#####
# 选择clone state1
sort(seismic_mat[rownames(pData(FC_cds)[pData(FC_cds)$State=='1', ]), 'X17'])
select_clone = 'P4_clone_3'
tmp_file = 'P4.txt'
tmp_clone= 'clone_3'
tmp_bc = rownames(srt@meta.data[srt$new_group==select_clone, ])
tmp_sctc = BC_list[[tmp_file]]
tmp_sctc_cna = tmp_sctc$orig.data$all_node_data

tmp_sctc_cna = tmp_sctc_cna[tmp_bc, grepl('^17_', colnames(tmp_sctc_cna))]
pheatmap::pheatmap(tmp_sctc_cna, cluster_cols = F, cluster_rows = F)
plot(unlist(tmp_sctc_cna[7,]))

bc_cna = as.data.frame(t(tmp_sctc_cna[7,]))
colnames(bc_cna) = 'CNA'
bc_cna$chr = sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[1]])
bc_cna$start = sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[2]])
bc_cna$end = sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[3]])
bc_cna$start = as.numeric(bc_cna$start) / (1*1000 * 1000)
bc_cna$end = as.numeric(bc_cna$end) / (1*1000 * 1000)

a = bc_cna %>%
  ggplot()+
  #geom_segment(aes(x=start,xend=end, y=CNA, yend=CNA), color='darkred', linewidth=2)+
  geom_point(aes(x=start, y=CNA), size=0.2, color='red')+
  labs(x='Genomic_position(MB)', title='Chr17')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold'))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/Seismic_展示chr17.pdf'),
       a, width=100, height = 60, units='mm', dpi = 600)

### 17号染色体拷贝数倍性 ###
all_cna = c()
names_meta = c()
for(prefix in names(BC_list)){
  print(prefix)
  sctc = BC_list[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_meta = srt@meta.data[srt$file==prefix, ]
  for(cl in unique(tmp_meta$new_group)){
    print(cl)
    tmp_bc = rownames(tmp_meta[tmp_meta$new_group==cl, ])
    tmp_cna_bc = tmp_cna[tmp_bc, ]
    all_cna = rbind(all_cna, colMeans(tmp_cna_bc))
    names_meta = rbind(names_meta, c(prefix, cl))
  }
}
rownames(all_cna) = names_meta[,2]

all_cna_mat = all_cna[unique(meta_1_5$new_group), grepl('^17_', colnames(all_cna))]

order_meta = arrange(distinct(meta_1_5[,c('State', 'new_group')]), State)

Heatmap(all_cna_mat[order_meta$new_group, ],
        cluster_rows = F, cluster_columns = F,
        show_row_names = F, show_column_names = F,
        left_annotation = rowAnnotation(df=order_meta[, 'State', drop=F]),
        row_split = order_meta[, 'State']
)
all_cna_mat = as.data.frame(all_cna_mat[order_meta$new_group, ])
all_cna_mat[, 'group'] = order_meta$State
all_cna_mat[is.na(all_cna_mat)] = 2
mean_ploidy = all_cna_mat %>%
  group_by(group) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(mean_ploidy) = mean_ploidy[,1]
mean_ploidy = mean_ploidy[, -1]
mean_ploidy = as.data.frame(t(mean_ploidy))
colnames(mean_ploidy) = c('State1', 'State5')
mean_ploidy$chr = sapply(rownames(mean_ploidy), function(x)strsplit(x,'_')[[1]][[1]])
mean_ploidy$start = as.numeric(sapply(rownames(mean_ploidy), function(x)strsplit(x,'_')[[1]][[2]]))
mean_ploidy$end = as.numeric(sapply(rownames(mean_ploidy), function(x)strsplit(x,'_')[[1]][[3]]))


#### 识别扩增基因 #####
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')
gene_info = gene_info[gene_info$chrom=='17', ]

bc_cna = as.data.frame(t(tmp_sctc_cna[7,]))
colnames(bc_cna) = 'CNA'
bc_cna$chr = sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[1]])
bc_cna$start = as.numeric(sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[2]]))
bc_cna$end = as.numeric(sapply(rownames(bc_cna), function(x)strsplit(x,'_')[[1]][[3]]))
bc_cna = bc_cna[, c(2,3,4,1)]
colnames(bc_cna) = c('chrom', 'start', 'end', 'cna_name')
bc_cna$chrom = as.character(bc_cna$chrom)
bc_cna_gene = valr::bed_intersect(bc_cna,gene_info) %>% as.data.frame()

label_gene = c("KRT17", "CAVIN1",'KRT16',"NFE2L1","NMT1","RPS6KB1")
label_data = bc_cna_gene[bc_cna_gene$gene_name.y%in%label_gene, ]
bc_cna_gene %>%
  ggplot()+
  #geom_segment(aes(x=start,xend=end, y=CNA, yend=CNA), color='darkred', linewidth=2)+
  geom_point(aes(x=start.x, y=cna_name.x), size=0.2, color='red')+
  geom_text_repel(data=label_data, mapping=aes(x=start.x, y=cna_name.x,label=gene_name.y), nudge_y = 0.5)+
  labs(x='Genomic_position(MB)', title='Chr17')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold'))


intersect(top_gene, bc_cna_gene$gene_name.y)
label_gene = c("KRT17", "CAVIN1",'KRT16',"NFE2L1","NMT1","RPS6KB1")
label_data2 = merge(mean_ploidy, label_data, by.x='start', by.y='start.x')

a=mean_ploidy %>%
  ggplot(aes(x=start / (1*1000 * 1000)))+
  geom_line(aes(y=State1), color='#1F77B4',size=1.2)+
  geom_line(aes(y=State5), color='#9467BD',size=1.2)+
  geom_text_repel(data=label_data2,
                  mapping=aes(x=start / (1*1000 * 1000), y=State1,label=gene_name.y),
                  max.overlaps = Inf, nudge_y = 0.5)+
  labs(x='Genomic_position(MB)',y='mean_CNA', title='Chr17')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold'))

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/Seismic_展示chr17_mean.pdf'),
       a, width=100, height = 60, units='mm', dpi = 600)

#### 差异Hallmark ####
hallmark_data = GSEABase::getGmt('/Users/lab/wangxin/Rprojects/cancer/h.all.v7.4.symbols.gmt')
hallmark_data = hallmark_data@.Data
hallmark_list = list()
for(i in hallmark_data){
  hallmark_list[[i@setName]] = i@geneIds
}
srt_with_score = AddModuleScore(srt, features = hallmark_list, name = names(hallmark_list))
srt_with_score = srt_with_score@meta.data
colnames(srt_with_score)[(ncol(srt_with_score)-length(hallmark_list)+1): ncol(srt_with_score)] = names(hallmark_list)
srt_with_score = srt_with_score[,names(hallmark_list)]

srt_with_score$new_group = srt$new_group
srt_with_score_clone = srt_with_score %>% group_by(new_group) %>% summarise_all(list(mean)) %>% as.data.frame()
rownames(srt_with_score_clone)  =srt_with_score_clone[,1]
srt_with_score_clone = srt_with_score_clone[, -1]

cds_meta = pData(FC_cds)

fc_hall_res = c()
for(i in unique(cds_meta$State)){
  print(i)
  tmp_bc1 = rownames(cds_meta[cds_meta$State==i,])
  tmp_bc2 = rownames(cds_meta[cds_meta$State!=i,])
  for(h in colnames(srt_with_score_clone)){
    tmp_fc = mean(srt_with_score_clone[tmp_bc1, h]) / mean(srt_with_score_clone[tmp_bc2, h])
    tmp_pv = t.test(srt_with_score_clone[tmp_bc1, h], srt_with_score_clone[tmp_bc2, h])$p.value
    fc_hall_res = rbind(fc_hall_res, c(i,h,tmp_fc,tmp_pv))
  }
}
fc_hall_res = as.data.frame(fc_hall_res)
colnames(fc_hall_res) = c('State', 'HALLMARK', 'FC', 'pvalue')
fc_hall_res$FC = as.numeric(fc_hall_res$FC)
fc_hall_res$pvalue = as.numeric(fc_hall_res$pvalue)

top_mark = fc_hall_res %>% group_by(State) %>% top_n(4, wt =FC) %>% as.data.frame()

pheatmap::pheatmap(srt_with_score_clone[,unique(top_mark$HALLMARK)])
tmp_score = srt_with_score_clone[,unique(top_mark$HALLMARK)]
tmp_score$State = cds_meta[rownames(tmp_score), 'State']
tmp_score = tmp_score %>% group_by(State) %>% summarise_all(list(mean)) %>% as.data.frame()
rownames(tmp_score)  =tmp_score[,1]
tmp_score = tmp_score[, -1]

plot_data = reshape2::melt(as.matrix(tmp_score))

aa = t(tmp_score)[c(1,5,6,7,10,11,12,13,14,15), c(4,3,2,1,5)]
width_in =  90/ 25.4
height_in = 50 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/state功能热图.pdf', width=width_in, height=height_in)

pheatmap::pheatmap(aa[c(7,6,5,4,10,9,3,1,2,8), ],
                   cluster_cols = F, cluster_rows = F,
                   fontsize=6,
                   border_color = 'white')
dev.off()


plot_data %>%
  ggplot(aes(x=Var1, y=Var2, size=value, color=value))+
  geom_point()+
  scale_color_gradientn(colours = c('white', 'red'))+
  theme_classic()

pData(FC_cds)$HALLMARK_score = srt_with_score_clone[rownames(pData(FC_cds)), 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
g=plot_cell_trajectory(FC_cds, color_by = "HALLMARK_score", cell_size = 1)+theme(legend.position = 'right')+
  scale_color_gradientn(colours = c('gray','red'))
g

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_Hallmark.pdf'),
       g, width=120, height = 80, units='mm', dpi = 600)











### 热图差异基因功能 ####
BEAM_res <- BEAM(FC_cds, branch_point = 2, cores = 1, branch_states = c('1','5'))
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
top_gene = row.names(subset(BEAM_res,qval < 1e-4))

pgbh = plot_genes_branched_heatmap(FC_cds[top_gene,],branch_point = 2, num_clusters = 3, return_heatmap = TRUE)
order_mat = pgbh$annotation_row
order_mat = order_mat[order(order_mat$Cluster),, drop=F]
# 绘图
keep_gene = intersect(top_gene, bc_cna_gene$gene_name.y)#c('ALOX5', 'CTSE', 'CXCL14', 'DDIT4', 'MGLL', 'PDZK1IP1', 'SLC20A1')
mark_gene = keep_gene#intersect(c2_gene, pc_df_keman_gene)
mark_pos = match(mark_gene, rownames(order_mat))
ht = Heatmap(pgbh$heatmap_matrix[rownames(order_mat), ],
             cluster_rows = T, cluster_columns = F,
             show_row_names = F, show_column_names = F,
             top_annotation = HeatmapAnnotation(df = pgbh$annotation_col, col=list('Cell Type'=c('Cell fate 1'="#1F77B4FF",
                                                                                                 'Cell fate 2'="#9467BDFF",
                                                                                                 'Pre-branch'="#CDDEB7FF"))),
             left_annotation = rowAnnotation(df = order_mat,
                                             col=list('Cluster'=c('1'="#1F77B4FF", '2'="#FF7F0EFF",
                                                                  '3'="#2CA02CFF"))),
             right_annotation = HeatmapAnnotation(a=anno_mark(at = mark_pos, labels = mark_gene,
                                                              which = "row", side='right',
                                                              link_width = unit(8, "mm"),
                                                              padding = unit(0.5, "mm"),
                                                              link_height = unit(5, "mm"),
                                                              labels_gp = gpar(cex=2)),
                                                  which = "row"),
             row_split = order_mat$Cluster,
             column_split = c(rep(1,100), rep(2,100)),
             row_title = NULL, column_title = NULL
)
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_DF.pdf'),width = 8,height = 8)
draw(ht, padding = unit(c(20, 10, 10, 10), "mm"))
dev.off()

# 功能注释
hallmark_data = GSEABase::getGmt('/Users/lab/wangxin/Rprojects/cancer/h.all.v7.4.symbols.gmt')
hallmark_data = hallmark_data@.Data
hallmark2gene = c()
for(i in hallmark_data){
  hallmark2gene = rbind(hallmark2gene, data.frame(termid=i@setName, geneid=i@geneIds))
}
hallmark2gene = as.data.frame(hallmark2gene)
library(clusterProfiler)
library(org.Hs.eg.db)
hm_res = list()
for(i in unique(order_mat$Cluster)){
  tmp_gene = rownames(order_mat[order_mat$Cluster==i,,drop=F])
  if(length(tmp_gene)<5){
    hm_res[[i]] = ''
  }else{
    tmp_enrich = enricher(tmp_gene, TERM2GENE = hallmark2gene)@result

    tmp_enrich = tmp_enrich[tmp_enrich$p.adjust<0.05, ]
    tmp_hm = head(tmp_enrich$ID, 5)
    #hm_res = rbind(hm_res, data.frame(name=tmp_hm, cluster=i))
    if(length(tmp_hm)<=5){
      test1 = bitr(tmp_gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
      tmp_go = enrichGO(
        gene=unique(test1$ENTREZID),
        keyType="ENTREZID",
        OrgDb=org.Hs.eg.db,
        ont="ALL",
        pAdjustMethod="BH",
        pvalueCutoff=0.01,
        qvalueCutoff=0.05,
        readable=TRUE
      )
      tmp_go = tmp_go@result %>% group_by(ONTOLOGY) %>% top_n(2, wt=-p.adjust)

      #tmp_kegg = enrichKEGG(gene =unique(test1$ENTREZID),
      #                      organism = 'hsa',
      #                      pvalueCutoff = 0.05,
      #                      pAdjustMethod = "BH",
      #                      qvalueCutoff = 0.05)
      #tmp_kegg = tmp_kegg@result  %>% top_n(2, wt=-p.adjust)

      tmp_hm = c(tmp_hm,tmp_go$Description)#, tmp_kegg$Description)
    }
    hm_res[[i]] = tmp_hm
  }

}
# 文字注释
library(ComplexHeatmap)
text = lapply(hm_res, function(x) {
  data.frame(x, col = 'black', fontsize = 12)
})

right_anno = rowAnnotation(textbox = anno_textbox(order_mat$Cluster, hm_res,
                                                  background_gp = gpar(fill = "NA", col = "#AAAAAA")))

ht = Heatmap(pgbh$heatmap_matrix[rownames(order_mat), ],
             cluster_rows = T, cluster_columns = F,
             show_row_names = F, show_column_names = F,
             top_annotation = HeatmapAnnotation(df = pgbh$annotation_col, col=list('Cell Type'=c('Cell fate 1'="#1F77B4FF",
                                                                                                 'Cell fate 2'="#9467BDFF",
                                                                                                 'Pre-branch'="#CDDEB7FF"))),
             left_annotation = rowAnnotation(df = order_mat,
                                             col=list('Cluster'=c('1'="#1F77B4FF", '2'="#FF7F0EFF",
                                                                  '3'="#2CA02CFF"))),
             right_annotation = right_anno,
             row_split = order_mat$Cluster,
             column_split = c(rep(1,100), rep(2,100)),
             row_title = NULL, column_title = NULL
)

pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/monocle_DF_mark.pdf'),width = 15,height = 12)
draw(ht, padding = unit(c(20, 10, 10, 10), "mm"))
dev.off()

## 显著基因功能##
tmp_gene = c("KRT17", "CAVIN1",'KRT16',"NFE2L1","NMT1","RPS6KB1")
tmp_enrich = enricher(tmp_gene, TERM2GENE = hallmark2gene)@result
tmp_enrich = tmp_enrich[tmp_enrich$pvalue<0.05, ]
tmp_enrich
test1 = bitr(tmp_gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
tmp_go = enrichGO(
  gene=unique(test1$ENTREZID),
  keyType="ENTREZID",
  OrgDb=org.Hs.eg.db,
  ont="ALL",
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  qvalueCutoff=0.05,
  readable=TRUE
)
dotplot(tmp_go)
tmp_go@result

#1
plot_cell_trajectory(FC_cds,  markers = c("KRT17", "CAVIN1",'GAA','KRT16',"NFE2L1","NMT1","RPS6KB1"),
                     use_color_gradient = TRUE)+scale_color_gradientn(colours = c('gray', 'red'))
#5
plot_cell_trajectory(FC_cds,  markers = c("CYBC1", "TUBG2",'GAA','MAPT'), use_color_gradient = TRUE)

### CAVIN1 比较 ####
gene_cna = data.frame(CAVIN1=all_cna_mat[, '17_40553769_40565472'])
gene_cna$new_group = rownames(all_cna_mat)
gene_cna$State = order_meta[match(gene_cna$new_group, order_meta$new_group), 'State']
a=gene_cna %>%
  ggplot(aes(x=State, y=CAVIN1, fill=State))+
  geom_boxplot(outlier.shape = NA, width=0.6, size=0.2)+
  geom_jitter(width=0.2, size=0.1)+
  labs(x='State',y='CAVIN1_CNA')+
  scale_fill_manual(values = c('1'='#1F77B4FF', '5'='#9467BDFF'))+
  stat_compare_means(comparisons = list(c('1', '5')), label = 'p.signif')+
  scale_y_continuous(expand = expansion(mult = 0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold'))
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/CAVIN1_CNA比较.pdf'),
       a, width=60, height = 60, units='mm', dpi = 600)

#### 生存分析 #######
library(survival)
library(survminer)


select_gene = 'PTRF'#'CAVIN1'
metabric_cna = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_cna.txt',
                          sep='\t', header = T, check.names = F)

metabric_clinic = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TCGA_BC/brca_metabric/data_clinical_patient.txt',
                             sep='\t', header = T, check.names = F)
rownames(metabric_clinic) = metabric_clinic$PATIENT_ID
same_sample = intersect(metabric_clinic$PATIENT_ID, colnames(metabric_cna))

# setdiff(metabric_clinic$PATIENT_ID, colnames(metabric_cna))
metabric_clinic[same_sample, 'group'] = unlist(metabric_cna[metabric_cna$Hugo_Symbol==select_gene, same_sample])
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
metabric_clinic[same_sample, 'group'] = ifelse(metabric_clinic[same_sample, 'group']>0, 'Gain', 'Neutral')
metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "CNA",#改变图例名称
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
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Gain', 'Neutral')

metabric_clinic$OS_STATUS = as.numeric(sapply(metabric_clinic$OS_STATUS, function(x)strsplit(x,':')[[1]][1]))
km <- survfit(Surv(OS_MONTHS, OS_STATUS)~group, data = metabric_clinic)
surv_pvalue(km)
ggsurvplot(km,
           pval=T,
           egend = "bottom", #将图例移动到下方
           legend.title = "CNA",#改变图例名称
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
metabric_clinic = metabric_clinic %>% filter(group%in%c(0,1))
metabric_clinic[, 'group'] = ifelse(metabric_clinic[, 'group']>0, 'Gain', 'Neutral')

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












library(colorRamps)
get_hp_data = function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2",
                        num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
                        add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
                        norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
                        trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
                        cores = 1)
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) ==
      FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  return(heatmap_matrix)
}
### 绘制热图 ###
pgbh=get_hp_data(FC_cds[sig_gene_names,],
                 num_clusters = 5,
                 cores = 1,
                 show_rownames = T, return_heatmap = T)
pgbh2 = pheatmap::pheatmap(pgbh, cluster_cols = F)
cluster = cutree(pgbh2$tree_row, 5)
cluster = as.data.frame(cluster) %>% arrange(cluster)
new_order = c()
for(i in 1:5){
  new_order = c(new_order, names(cluster[cluster==i]))
}

mark_gene = c('RAP2B', 'PLCH1', 'GMPS', 'PFN2')#intersect(c2_gene, pc_df_keman_gene)

order_mat = pgbh$annotation_row
order_mat = order_mat[order(order_mat$Cluster),, drop=F]

mark_pos = match(mark_gene, rownames(order_mat))
ht = Heatmap(pgbh$heatmap_matrix[rownames(order_mat), ],
             cluster_rows = T, cluster_columns = F,
             show_row_names = F, show_column_names = F,
             left_annotation = rowAnnotation(df = order_mat,
                                             col=list('Cluster'=c('1'="#1F77B4FF", '2'="#FF7F0EFF",
                                                                  '3'="#2CA02CFF", '4'="#D62728FF",'5'='red'))),
             right_annotation = HeatmapAnnotation(a=anno_mark(at = mark_pos, labels = mark_gene,
                                                              which = "row", side='right',
                                                              link_width = unit(8, "mm"),
                                                              padding = unit(0.5, "mm"),
                                                              link_height = unit(5, "mm"),
                                                              labels_gp = gpar(cex=2)),
                                                  which = "row"),
             row_split = order_mat$Cluster,
             column_split = c(rep(1,100), rep(2,100)),
             row_title = NULL, column_title = NULL
)
pdf(paste0('../output/', '样本水平trajectory', '.pdf'),width = 8,height = 8)
draw(ht, padding = unit(c(20, 10, 10, 10), "mm"))
dev.off()


#####
gene_name = c('RAP2B', 'PLCH1', 'GMPS', 'PFN2')
plot_genes_branched_pseudotime(FC_cds[gene_name,],
                               branch_point = 1,
                               color_by = "pde",
                               ncol = 1)+ scale_color_gradientn(colours = c('blue','white', 'red'))

tmp_data = pData(FC_cds)[,c('new_group','Pseudotime', 'State', 'cytotrace_score','pde', 'Cell_ITH')]

tmp_gene_exp = AverageExpression(srt, group.by = 'new_group')$RNA %>% log1p()

tmp_data$RAP2B = tmp_gene_exp['RAP2B', rownames(tmp_data)]
tmp_data$Pseudotime = (tmp_data$Pseudotime - min(tmp_data$Pseudotime)) / (max(tmp_data$Pseudotime)-min(tmp_data$Pseudotime))
tmp_data = tmp_data %>% arrange(Pseudotime)

tmp_data %>%
  ggplot(aes(x=Pseudotime,y=RAP2B, color=State))+
  geom_point(outlier.shape = NA)

tmp_data %>%
  ggplot(aes(x=State,y=RAP2B, color=State))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)+
  scale_color_igv()+
  stat_compare_means(comparisons = list(c('5','2'),c('5','3'),c('5','1'),c('5','4'),
                                        c('1','2'),c('1','3'),c('1','4')) , label='p.signif')+
  theme_classic()

tmp_data %>%
  ggplot(aes(x=State,y=Cell_ITH, color=State))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)+
  scale_color_igv()+
  stat_compare_means(comparisons = list(c('5','2'),c('5','3'),c('5','1'),c('5','4'),
                                        c('1','2'),c('1','3'),c('1','4')) , label='p.signif')+
  theme_classic()

##### hallmark #############
enrich_with_marker <- function(FindAllMarkers_res, ref_list, topn=20, scale=TRUE){
  toplist <- FindAllMarkers_res %>% group_by(cluster) %>% top_n(topn, wt = avg_log2FC)
  cell_types = names(ref_list)
  cluster = sort(unique(toplist$cluster))
  inter_mat = matrix(0, ncol = length(cluster), nrow = length(cell_types))

  for (i in 1:length(cell_types)) {
    for (j in 1:length(cluster)) {
      pm1 = ref_list[[i]]
      pm2 = toplist[toplist$cluster==cluster[j], "gene", drop=TRUE]
      inter_mat[i,j] = length(intersect(pm1, pm2))
    }
  }
  rownames(inter_mat) = c(cell_types)
  inter_mat <- inter_mat[!rowSums(inter_mat)==0,]


  sapply(rownames(inter_mat), function(x)length(ref_list[[x]]))

  norma = sapply(rownames(inter_mat), function(x)length(ref_list[[x]])) / rowSums(inter_mat)
  enri = apply(inter_mat,2,function(x)x/norma)
  enri = apply(enri, 2, function(x) x/sum(x))
  if(scale){
    enri = t(scale(t(enri)))
  }
  rownames(enri) = rownames(inter_mat)
  colnames(enri) = cluster
  enri[which(is.nan(enri))] = 0

  pheatmap::pheatmap(enri,
                     border_color = 'white')
  return(enri)
}

hallmark_data = GSEABase::getGmt('/Users/lab/wangxin/Rprojects/cancer/h.all.v7.4.symbols.gmt')
hallmark_data = hallmark_data@.Data
hallmark_list = list()
for(i in hallmark_data){
  hallmark_list[[i@setName]] = i@geneIds
}

Idents(srt) = srt$pde_group
marker_res = FindAllMarkers(srt,  max.cells.per.ident = 300)
x = enrich_with_marker(marker_res, hallmark_list, topn = 20,scale=F)
x[x<0.05] = NA
x = x[rowSums(is.na(x))<2,]
pheatmap::pheatmap(x, cluster_rows = F, cluster_cols = F,border_color = 'white', na_col = 'gray')

library(clusterProfiler)
library(org.Hs.eg.db)
tmp_marker = marker_res[marker_res$avg_log2FC>0&marker_res$p_val_adj<0.05,]
test1 = bitr(rownames(tmp_marker)[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
tmp_go = enrichGO(
  gene=unique(test1$ENTREZID),
  keyType="ENTREZID",
  OrgDb=org.Hs.eg.db,
  ont="ALL",
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  qvalueCutoff=0.05,
  readable=TRUE
)

a = dotplot(tmp_go)
a


##### hallmark #############
library(ComplexHeatmap)
srt_with_score = AddModuleScore(srt, features = hallmark_list, name = names(hallmark_list))
srt_with_score = srt_with_score@meta.data
colnames(srt_with_score)[(ncol(srt_with_score)-length(hallmark_list)+1): ncol(srt_with_score)] = names(hallmark_list)

srt_with_score = srt_with_score[,names(hallmark_list)]

srt_with_score = t(srt_with_score)

top_meta = srt@meta.data[, 'pde_group', drop=F] %>% arrange(pde_group)

ht = Heatmap(srt_with_score[, rownames(top_meta)],
             cluster_rows = T, cluster_columns = T,
             show_column_names = F,
             top_annotation = HeatmapAnnotation(df=top_meta,
                                                col = list('pde_group'=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))),
             column_split = top_meta$pde_group
)
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/hallmark.pdf'),width = 8,height = 8)
draw(ht, padding = unit(c(20, 10, 10, 10), "mm"))
dev.off()

srt_with_score = t(srt_with_score)

srt_with_score_new = apply(srt_with_score, 2, function(x)(x-min(x))/(max(x)-min(x)))


high_bc = grep('NA',rownames(srt@meta.data[srt$pde_group=='High',]), value = T, invert = T)
low_bc = grep('NA',rownames(srt@meta.data[srt$pde_group=='Low',]), value = T, invert = T)
hallmark_fc = c()
for(i in 1:ncol(srt_with_score_new)){
  print(i)
  tmp_h = srt_with_score_new[high_bc, i]
  tmp_l = srt_with_score_new[low_bc, i]
  x = t.test(tmp_h, tmp_l)
  tmp_fc = mean(tmp_h) / mean(tmp_l)
  hallmark_fc = rbind(hallmark_fc, c(names(hallmark_list)[i], tmp_fc, x$p.value))
}

hallmark_fc = as.data.frame(hallmark_fc)
colnames(hallmark_fc) = c('HALLMARK', 'FC', 'pValue')
hallmark_fc$FC = as.numeric(hallmark_fc$FC)
hallmark_fc$pValue = as.numeric(hallmark_fc$pValue)

#hallmark_fc = hallmark_fc[hallmark_fc$pValue<0.01, ]
#hallmark_fc = hallmark_fc[abs(log2(hallmark_fc$FC))>0.1, ]
label_data = hallmark_fc[log2(hallmark_fc$FC)>0.1 & hallmark_fc$pValue<0.01,]
hallmark_fc %>%
  ggplot(aes(x=log2(FC), y=-log(pValue)))+
  geom_point()+
  geom_text_repel(aes(label=HALLMARK), data=label_data)+
  theme_bw()


ht = Heatmap(t(srt_with_score_new)[, rownames(top_meta)],
             cluster_rows = T, cluster_columns = F,
             show_column_names = F,
             top_annotation = HeatmapAnnotation(df=top_meta,
                                                col = list('pde_group'=c('High'='#E64B35FF', 'Low'='#4DBBD5FF'))),
             column_split = top_meta$pde_group
)
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图BC_调整拷贝数/hallmark.pdf'),width = 8,height = 8)
draw(ht, padding = unit(c(20, 10, 10, 10), "mm"))
dev.off()
