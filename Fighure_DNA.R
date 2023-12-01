source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')


scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine/'
prefix = 'huh7'
sctc = scTraceClass(scTrace_dir,
                    prefix,
                    layout = 'slanted',
                    rate = 0.8,
                    min_cell_num=150)

###########
sctc$map_obj$cell_map_pos$Fn = sapply(rownames(sctc$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
plot_cells(sctc, colorby='Fn',shadow=TRUE, clone_line=TRUE,
           point.size=2, plot_virtual=F)+
  ggsci::scale_fill_igv()

plot_cells(sctc, colorby='Fn',
           type = 'pie',
           shadow=TRUE, clone_line=TRUE,
           point.size=2, plot_virtual=F,
           r_size= 0.02)+
  ggsci::scale_fill_igv()


##########
gene.loc <-  read.table('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_crispr2/geneFile.txt')
colnames(gene.loc) = c('gene', 'chr', 'start', 'end')
rownames(gene.loc) = gene.loc$gene
gene.loc = gene.loc[colnames(sctc$orig.data$all_node_data), 'chr']

plot_cna_tree(sctc,
              colorby = 'Fn',
              output = paste0(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/DNAseq/{prefix}_tree_heatmap.pdf')),
              #column_split=sort(as.numeric(gsub('chr', '',gene.loc))),
              pt.size = 0.1,
              width=220,
              height=100,
              rel_widths=c(0.3,1),
              rel_heights = c(.0, 1)
)



#############
#####
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
                      rate = 0.9,
                      min_cell_num=NULL)
  DNA_list[[prefix]] = sctc
}
saveRDS(DNA_list, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc.rds')



####

plist = list()
for(prefix in dna_file){
  print(prefix)
  tmp_class = DNA_list[[prefix]]
  tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])

  cell_pos = tmp_class$map_obj$cell_map_pos
  clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
  rownames(clone_pos) = clone_pos$label
  cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

  #
  r_size = 0.04
  pie_data = cell_pos %>%
    group_by(clone, !!sym('Fn')) %>%
    summarise(num=n()) %>%
    mutate(total=sum(num)) %>% as.data.frame()
  pie_data[, c('x', 'y')] = clone_pos[pie_data$clone,c('x', 'y')]
  pie_data = pie_data %>% group_by(x, y, total)
  color_igv = pal_igv()(length(unique(pie_data$Fn)))
  names(color_igv) = unique(pie_data$Fn)
  df.grobs <- pie_data %>%
    do(subplots = ggplot(., aes(1, num, fill = !!sym('Fn'))) +
         geom_col(position = "fill", alpha = 1, colour = NA) +
         coord_polar(theta = "y") +
         #scale_fill_igv()+
         scale_fill_manual(values = color_igv)+
         theme_void()+ guides(fill = F)) %>%
    mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                             x = -x-r_size, y = y-r_size,
                                             xmax = -x+r_size, ymax = y+r_size)))
  p = ggplot()
  p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                       color='#dddddd',
                       alpha=0.6,
                       linewidth=4,
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
    coord_flip()+scale_x_reverse()

  plist[[prefix]] = p
}
cowplot::plot_grid(plotlist = plist, nrow=1)






##############
# F1,F2比较----pseudotime
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

all_meta = c()
for(f in dna_file){
  print(f)
  tmp_obj = load_scTrace(scTrace_dir, f)
  meta = tmp_obj$cell_relation
  meta$Fn = sapply(rownames(meta), function(x)strsplit(x, '_')[[1]][1])
  meta$source = f
  meta = meta %>% filter(Fn %in%c('f1', 'f2'))
  meta$Pseudotime_tree = (meta$Pseudotime_tree-min(meta$Pseudotime_tree)) / (max(meta$Pseudotime_tree) - min(meta$Pseudotime_tree))
  all_meta = rbind(all_meta, meta)
}
all_meta = as.data.frame(all_meta)
library(ggpubr)
all_meta %>%
  filter(Fn %in%c('f1', 'f2')) %>%
  group_by(source, Fn) %>%
  summarise(mean_Pseudotime = median(Pseudotime_tree)) %>%
  ggplot(aes(x=Fn, y=mean_Pseudotime))+
  geom_boxplot()+
  geom_jitter(aes(fill=source),width=0.2, size=4,color='black', shape=21)+
  scale_fill_igv()+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()

all_meta %>%
  ggplot(aes(x=Pseudotime_tree, fill=Fn))+
  geom_density(alpha=0.5)+
  scale_fill_nejm()
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)

all_meta %>%
  ggplot(aes(x=Fn, y=Pseudotime_tree))+
  geom_boxplot(fill='blue', alpha=0.1)+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = FALSE)+
  theme_classic()


# score
file_path = '/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/'
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

all_meta_score = c()
for(f in dna_file){
  print(f)
  tmp_obj = load_scTrace(scTrace_dir, f)
  meta = tmp_obj$cell_relation
  gene_level_cna_data = read.table(paste0(file_path, f, '.txt_gene_data_all_node.csv'))
  tmp_obj = add_cna_score(tmp_obj,pathway_list, gene_level_cna_data)
  meta$Fn = sapply(rownames(meta), function(x)strsplit(x, '_')[[1]][1])
  meta$source = f
  meta = meta %>% filter(Fn %in%c('f1', 'f2'))
  rownames(tmp_obj$cna_score) = gsub('\\.', '-', rownames(tmp_obj$cna_score))
  meta[, colnames(tmp_obj$cna_score)] = tmp_obj$cna_score[rownames(meta),]
  meta$Pseudotime_tree = (meta$Pseudotime_tree-min(meta$Pseudotime_tree)) / (max(meta$Pseudotime_tree) - min(meta$Pseudotime_tree))
  all_meta_score = rbind(all_meta_score, meta)
}

# F1,F2 cell cycle 评分比较

all_meta_score %>%
  ggplot(aes(x=Fn, y=CellCycle_OG))+
  geom_boxplot(fill='blue', alpha=0.1)+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = FALSE)+
  theme_classic()
all_meta_score %>%
  filter(Fn %in%c('f1', 'f2')) %>%
  group_by(source, Fn) %>%
  summarise(CellCycle_OG = median(CellCycle_OG)) %>%
  ggplot(aes(x=Fn, y=CellCycle_OG))+
  geom_boxplot()+
  geom_jitter(aes(fill=source),width=0.2, size=4,color='black', shape=21)+
  scale_fill_igv()+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()
# F1,F2 增殖评分与分裂速度相关性

DNA_sctc = list()
for(f in dna_file){
  DNA_sctc[[f]] = scTraceClass(scTrace_dir, f,
                               layout = 'slanted',
                               rate=0.8,
                               min_cell_num=NULL)
}

saveRDS(DNA_sctc, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNAseq_scClass.rds')
all_meta_rate = c()
for(f in dna_file){
  print(f)
  tmp_class = DNA_sctc[[f]]
  tmp_obj = tmp_class$orig.data
  meta = tmp_obj$cell_relation
  tmp_rate = predict_mitosis_rate(as.data.frame(meta))
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

tmp_func <- function(x,y){
  x[which.max(y)]
}
new_Fn = all_meta_rate %>%
  group_by(Fn, source, clone) %>%
  summarise(n=n()) %>%
  group_by(clone, source) %>%
  summarise(clone_fn=tmp_func(Fn, n))
all_meta_rate = left_join(all_meta_rate, new_Fn)

# pseudotime
a=all_meta_rate %>%
  group_by(source, clone_fn,clone) %>%
  summarise(Pseudotime_tree = mean(Pseudotime_tree)) %>%
  ggplot(aes(x=clone_fn, y=Pseudotime_tree))+
  geom_boxplot()+
  geom_jitter(aes(fill=source),width=0.2, size=4,color='black', shape=21)+
  scale_fill_igv()+
  facet_wrap(~source, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format')+
  theme_classic()

# rate
b=all_meta_rate %>%
  group_by(source, clone_fn,clone) %>%
  summarise(pred_rate = mean(pred_rate)) %>%
  ggplot(aes(x=clone_fn, y=pred_rate))+
  geom_boxplot()+
  geom_jitter(aes(fill=source),width=0.2, size=4,color='black', shape=21)+
  scale_fill_igv()+
  facet_wrap(~source, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format')+
  theme_classic()

# score
c=all_meta_rate %>%
  group_by(source, clone_fn,clone) %>%
  summarise(CellCycle_OG = mean(CellCycle_OG)) %>%
  ggplot(aes(x=clone_fn, y=CellCycle_OG))+
  geom_boxplot()+
  geom_jitter(aes(fill=source),width=0.2, size=4,color='black', shape=21)+
  scale_fill_igv()+
  facet_wrap(~source, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format')+
  theme_classic()
cowplot::plot_grid(a,b,c,ncol=1)

# 分裂速率与增殖评分相关
all_meta_rate %>%
  group_by(source, clone,clone_fn) %>%
  summarise(pred_rate = median(pred_rate), score = mean(CellCycle_OG)) %>%
  ggplot(aes(x=score, y=pred_rate))+
  geom_point(aes(color=source, shape=clone_fn), size=2)+
  scale_color_igv()+
  facet_wrap(~source, scales = 'free')+
  stat_cor()+
  labs(x='CellCycle_OG', y='')+
  stat_smooth(method = 'lm', formula = 'y ~ x')+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()

# pseudotime与增殖评分相关
all_meta_rate %>%
  group_by(source, clone,clone_fn) %>%
  summarise(Pseudotime_tree = mean(Pseudotime_tree), score = mean(CellCycle_TSG)) %>%
  ggplot(aes(x=score, y=Pseudotime_tree))+
  geom_point(aes(color=source, shape=clone_fn), size=2)+
  scale_color_igv()+
  stat_cor()+
  labs(x='CellCycle_OG', y='Pseudotime')+
  #facet_wrap(~source, scales = 'free')+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linewidth=0.5)+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()


all_meta_rate %>%
  group_by(source, clone,clone_fn) %>%
  summarise(pred_rate = mean(pred_rate), score = mean(CellCycle_TSG)) %>%
  ggplot(aes(x=score, y=pred_rate))+
  geom_point(aes(color=source, shape=clone_fn),width=0.2, size=2)+
  scale_color_igv()+
  #facet_wrap(~source, scales = 'free')+
  stat_cor()+
  stat_smooth(method = 'lm', formula = 'y ~ x')+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()



# RNA 干性得分与 增殖，pseudotime，分裂速率相关性
srt[['clone']]=sapply(srt$clonealign, function(x){
  y=strsplit(x, '_')[[1]][2:3]
  return(paste0(y, collapse = '_'))
})
srt[['source']] = sapply(srt$clonealign, function(x){
  y=strsplit(x, '\\.')[[1]][1]
  return(y)
})

all_meta_rate2 = left_join(all_meta_rate, srt@meta.data[, c('clone', 'source', 'cytotrace_score')])


all_meta_rate2 %>%
  filter(clone_fn=='f2') %>%
  filter(!is.na(cytotrace_score)) %>%
  group_by(source, clone,clone_fn) %>%
  summarise(cytotrace_score = mean(cytotrace_score), Pseudotime_tree=mean(Pseudotime_tree)) %>%
  ggplot(aes(x=Pseudotime_tree, y=cytotrace_score))+
  geom_point(aes(fill=source),size=4, shape=21,color='black')+
  scale_fill_igv()+
  #facet_wrap(~source, scales = 'free')+
  stat_cor()+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F)+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()


all_meta_rate2 %>%
  filter(clone_fn=='f2') %>%
  filter(!is.na(cytotrace_score)) %>%
  group_by(source, clone,clone_fn) %>%
  summarise(cytotrace_score = mean(cytotrace_score), pred_rate=mean(pred_rate)) %>%
  ggplot(aes(x=pred_rate, y=cytotrace_score))+
  geom_point(aes(fill=source),size=4, shape=21,color='black')+
  scale_fill_igv()+
  #facet_wrap(~source, scales = 'free')+
  stat_cor()+
  stat_smooth(method = 'lm', formula = 'y ~ x', se=F)+
  #stat_compare_means(comparisons = list(c('f1', 'f2')), label='p.format', paired = TRUE)+
  theme_classic()






