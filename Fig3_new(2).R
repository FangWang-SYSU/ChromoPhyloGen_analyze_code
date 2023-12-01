source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
draw_label_theme <- function(label, theme = NULL, element = "text", ...) {
  if (is.null(theme)) {
    theme <- ggplot2::theme_get()
  }
  if (!element %in% names(theme)) {
    stop("Element must be a valid ggplot theme element name")
  }

  elements <- ggplot2::calc_element(element, theme)

  cowplot::draw_label(label,
                      fontfamily = elements$family,
                      fontface = elements$face,
                      colour = elements$color,
                      size = elements$size,
                      ...
  )
}

dna_file = c(
  'uber.t12.seg.txt',
  'uber.t11.seg.txt',
  'uber.t10.seg.txt',
  'uber.t9.seg.txt',
  'uber.t8.seg.txt',
  'uber.t7.seg.txt',
  'uber.t6.seg.txt',
  'uber.t5.seg.txt',
  'uber.t4.seg.txt',
  'uber.t3.seg.txt',
  'uber.t2.seg.txt',
  'uber.t1.seg.txt',
  'KTN615_cleaned_uber.seg.txt',
  'KTN302_cleaned_uber.seg.txt',
  'KTN206_cleaned_uber.seg.txt',
  'KTN152_cleaned_uber.seg.txt',
  'KTN132_cleaned_uber.seg.txt',
  'KTN129_cleaned_uber.seg.txt',
  'KTN126_cleaned_uber.seg.txt',
  'KTN102_cleaned_uber.seg.txt')

TNBC_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all.rds')

patient_model = c()
model_plot_list = list()
for(prefix in dna_file){
  print(prefix)
  sctc = TNBC_list[[prefix]]
  sctc_cna = sctc$orig.data$all_node_data
  sctc_cna = sctc_cna[!grepl('root|virtual',rownames(sctc_cna)),]

  sctc_cna_event = sctc_cna-2
  events = sort(rowSums(sctc_cna_event!=0)) / ncol(sctc_cna)
  events_df = data.frame(events=events,step=1:length(events))
  model1 = lm(events~step, events_df)
  # reg = lm(events ~ step * I(step>= 12), events_df)
  x_range = 1:length(events)
  errs = sapply(x_range, function(BREAK) {
    mean(lm(events ~  I(step >= BREAK), data = events_df)$residuals^2)
  })
  # plot(x_range, errs)
  x_min = x_range[which.min(errs)]
  # axis(side = 1L, at = x_min)
  # abline(v = x_min, col = 'red')
  model_step = lm(events ~  I(step>= x_min), events_df)

  patient_model = rbind(patient_model, c(summary(model1)$adj.r.squared, AIC(model1), BIC(model1),
                                         summary(model_step)$adj.r.squared, AIC(model_step), BIC(model_step),
                                         x_min
                                         ))
  events_df$linear_fit =  predict(model1,events_df)
  events_df$step_fit =  predict(model_step,events_df)
  events_df$cut_error =  errs
  events_df$group = ifelse(events_df$step>x_min, 'Aneu', 'other')

  p1 = ggplot(events_df, aes(x=step))+
    geom_point(aes(y=events, color=group))+
    geom_line(aes(y=linear_fit), color='black')+
    labs(title = paste0('R^2: ', round(summary(model1)$adj.r.squared, 2),'; ',
                        'BIC: ',round(BIC(model1), 2)#, '; ',
                        #'AIC: ', round(AIC(model1),2)
                        ))+
    scale_color_igv()+
    theme_bw()+
    theme(panel.grid = element_blank(), plot.title = element_text(size=10), legend.position = 'none')

  p2 = ggplot(events_df, aes(x=step))+
    geom_point(aes(y=events, color=group))+
    geom_line(aes(y=step_fit), color='black')+
    labs(title = paste0('R^2: ', round(summary(model_step)$adj.r.squared, 2),'; ',
                        'BIC: ',round(BIC(model_step), 2)#, '; ',
                        #'AIC: ', round(AIC(model_step),2)
                        )
         )+
    scale_color_igv()+
    theme_bw()+
    theme(panel.grid = element_blank(), plot.title = element_text(size=10), legend.position = 'none')

  p3 = ggplot(events_df, aes(x=step))+
    geom_point(aes(y=cut_error, color=group))+
    geom_vline(xintercept = x_min, linetype='dashed', color='gray')+
    labs(title = 'Cut_num')+
    scale_color_igv()+
    theme_bw()+
    theme(panel.grid = element_blank(), plot.title = element_text(size=10), legend.position = 'none')

  p = cowplot::plot_grid(p3,p2,p1, nrow=1)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)


  title <- ggdraw() +
    draw_label_theme(new_prefix,
                      element = "plot.title",
                    hjust = 0.5, vjust = 0)
  p2=cowplot::plot_grid(title,p, ncol=1, rel_heights = c(0.1, 1))

  model_plot_list[[new_prefix]] = p2
  ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PG分组/分组', prefix,'.pdf'),
         p, width=240, height = 60, units='mm', dpi = 600)
}

patient_model = as.data.frame(patient_model)
colnames(patient_model) = c('Linear_R2', 'Linear_AIC', 'Linear_BIC',
                            'Step_R2', 'Step_AIC', 'Step_BIC', 'cut_num')

patient_model$patient = names(TNBC_list)
patient_model$patient = gsub('uber.', '', patient_model$patient)
patient_model$patient = gsub('.seg.txt', '', patient_model$patient)
patient_model$patient = gsub('_cleaned_uber', '', patient_model$patient)
patient_model$patient = gsub('_cleaned', '', patient_model$patient)


a = patient_model %>%
  mutate(patient=fct_reorder(patient, Step_R2-Linear_R2))%>%
  ggplot(aes(x=patient))+
  geom_segment(mapping = aes(x=patient, xend=patient, y=Linear_R2, yend=Step_R2), color='gray', linewidth=0.5)+
  geom_point(aes(y=Linear_R2, fill=Linear_BIC, shape='G_model'),  stroke=0.1, size=2)+
  geom_point(aes(y=Step_R2, fill=Step_BIC, shape='P_model'), stroke=0.1, size=2)+
  scale_shape_manual(values = c('G_model'=21, 'P_model'=23))+
  scale_fill_gradientn(colors=c('#FF2525', '#FFE53B'))+

  labs(y='adj_R^2')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/分组.pdf'),
       a, width=100, height = 60, units='mm', dpi = 600)

# patient order
p_sort = patient_model %>%
  mutate(patient=fct_reorder(patient, Step_R2-Linear_R2))
p_sort = levels(p_sort$patient)

p3=cowplot::plot_grid(plotlist = model_plot_list[p_sort], ncol=2)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PG分组/All分组model.pdf'),
       p3, width=300, height = 600, units='mm', dpi = 600)

# Ggroup & Pgroup
Pgroup = p_sort[14:20]
Ggroup = p_sort[1:14]

model_meta = data.frame(patient=p_sort, model=c(rep('Gradual', 14), rep('Punctuated', 6)),
                        row.names = p_sort)
#### 比较两组指标差异 #####
TNBC_list_with_score = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_sctc_all_with_score.rds')
all_cell_ITH_mat = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_Cell_ITH.rds')
# all_cell_ITH_mat$patient = all_cell_ITH_mat$c1
# all_cell_ITH_mat$patient = gsub('_.*', '', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('c.*', '', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('T1', 't1', all_cell_ITH_mat$patient)
# all_cell_ITH_mat$patient = gsub('T2', 't2', all_cell_ITH_mat$patient)
# saveRDS(all_cell_ITH_mat, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/TNBC_Cell_ITH.rds')
all_cell_ITH_mat = na.omit(all_cell_ITH_mat)
all_cell_ITH_mat = all_cell_ITH_mat %>%
  group_by(c1, patient) %>%
  summarise(cell_ITH=median(cell_ITH, na.rm = T)) %>% as.data.frame()

rownames(all_cell_ITH_mat) = all_cell_ITH_mat$c1
all_cell_ITH_mat$index='ITH'
all_cell_ITH_mat =all_cell_ITH_mat[,c('index', 'cell_ITH','patient')]
patient_index_score = c()
for(prefix in dna_file){
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  sctc = TNBC_list_with_score[[prefix]]
  meta = sctc$map_obj$cell_map_pos[, c('rearrange_score', 'new_BFB', 'pde')]
  meta = meta[!grepl('root|virtual', rownames(meta)), ]
  # meta$bc = rownames(meta)
  meta = reshape2::melt(meta)
  colnames(meta) = c('index', 'score')
  meta$patient = new_prefix
  patient_index_score = rbind(patient_index_score, meta)
}

colnames(all_cell_ITH_mat) = colnames(patient_index_score)
patient_index_score = rbind(patient_index_score, all_cell_ITH_mat)
patient_index_score$Model = model_meta[patient_index_score$patient, 'model']

library(ggpubr)
a=patient_index_score %>%
  group_by(patient,Model,index) %>%
  summarise(score=median(score)) %>%
  ggplot(aes(x=Model, y=score, fill=Model))+
  geom_boxplot(width=0.6, size=0.2, outlier.shape = NA)+
  geom_point(shape=21, size=2, stroke=0.2)+
  facet_wrap(~index, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('Gradual', 'Punctuated')),label = 'p.format')+
  labs(y='Score')+
  scale_fill_manual(values=c('Gradual'='#00913A', 'Punctuated'='#601986'))+
  scale_y_continuous(expand = expansion(mult = 0.2, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
        legend.position = 'none')
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/分组比较指标.pdf'),
       a, width=140, height = 60, units='mm', dpi = 600)


a=patient_index_score %>%
  group_by(patient,Model,index) %>%
  summarise(score=median(score)) %>%
  ggplot(aes(x=Model, y=score))+
  geom_boxplot(width=0.6, size=0.2, outlier.shape = NA)+
  geom_jitter(aes(fill=patient),width=0.2,shape=21, size=4, stroke=0.2)+
  facet_wrap(~index, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('Gradual', 'Punctuated')),label = 'p.format')+
  labs(y='Score')+
  scale_fill_igv()+
  scale_y_continuous(expand = expansion(mult = 0.2, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1))
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/分组比较指标(病人着色).pdf'),
       a, width=240, height = 160, units='mm', dpi = 600)


### PDE 分组比较 ####
patient_index_score2 = patient_index_score[patient_index_score$index!='pde', ]
patient_index_score_pde = patient_index_score[patient_index_score$index=='pde', ]
patient_index_score_pde = patient_index_score_pde %>%
  group_by(patient) %>%
  summarise(pde_score=mean(score)) %>%as.data.frame()
patient_index_score_pde$PDE_group = ifelse(patient_index_score_pde$pde_score>median(patient_index_score_pde$pde_score), 'PDE_high', 'PDE_low')

new_score = merge(patient_index_score_pde, patient_index_score2, by.x = 'patient', by.y='patient')

a= new_score %>%
  group_by(patient,PDE_group,index) %>%
  summarise(score=median(score)) %>%
  ggplot(aes(x=PDE_group, y=score, fill=PDE_group))+
  geom_boxplot(width=0.6, size=0.2, outlier.shape = NA)+
  geom_point(shape=21, size=2, stroke=0.2)+
  facet_wrap(~index, scales = 'free_y', nrow=1)+
  stat_compare_means(comparisons = list(c('PDE_low', 'PDE_high')),label = 'p.format')+
  labs(y='Score')+
  scale_fill_manual(values=c('PDE_high'='#E64B35FF', 'PDE_low'='#4DBBD5FF'))+
  scale_y_continuous(expand = expansion(mult = 0.2, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图TNBC/PDE高低分组比较指标.pdf'),
       a, width=110, height = 60, units='mm', dpi = 600)

###### PDE高低组病人 两两CNA距离 ####
gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = gene_info[, c('seqnames', 'start', 'end', 'gene_name')]
colnames(gene_info) = c('chrom', 'start', 'end', 'gene_name')
# install.packages('valr')
library(valr)
cna_gene = list()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)
  sctc = TNBC_list_with_score[[prefix]]
  tmp_cna = sctc$orig.data$all_node_data
  tmp_cna = tmp_cna[!grepl('root|virtual', rownames(tmp_cna)), ]

  # 位点转换成基因
  cna_loc = as.data.frame(t(sapply(colnames(tmp_cna), function(x)strsplit(x, '_')[[1]])))
  cna_loc$V2 = as.numeric(cna_loc$V2)
  cna_loc$V3 = c(as.numeric(cna_loc$V2[2:nrow(cna_loc)]), 0)
  cna_loc = cna_loc[cna_loc$V2<cna_loc$V3, ]
  cna_loc$cna_name = rownames(cna_loc)
  colnames(cna_loc) = c('chrom', 'start', 'end', 'cna_name')
  res = bed_intersect(cna_loc,gene_info)
  res = res%>%
    dplyr::mutate(percent_overlap = .overlap / (end.y-start.y)) %>%
    dplyr::filter(percent_overlap >= 1)

  # cna gene level
  cna_data_gene = t(tmp_cna)
  cna_data_gene = as.data.frame(cna_data_gene[res$cna_name.x, ])
  cna_data_gene$gene = res$gene_name.y
  cna_data_gene = cna_data_gene %>%
    group_by(gene) %>%
    summarise_all(list(mean)) %>% as.data.frame()

  rownames(cna_data_gene) = cna_data_gene[, 1]
  cna_data_gene = as.data.frame(t(cna_data_gene[,-1]))

  cna_gene[[prefix]] = cna_data_gene
}

same_gene = unique(gene_info$gene_name)
for(prefix in dna_file){
  same_gene = intersect(same_gene, colnames(cna_gene[[prefix]]))
}
patient_level_CNA = c()
for(prefix in dna_file){
  tmp_cna = cna_gene[[prefix]][, same_gene]
  patient_level_CNA = rbind(patient_level_CNA, round(colMeans(tmp_cna)))
}
new_prefix = gsub('uber.', '', dna_file)
new_prefix = gsub('.seg.txt', '', new_prefix)
new_prefix = gsub('_cleaned_uber', '', new_prefix)
new_prefix = gsub('_cleaned', '', new_prefix)
rownames(patient_level_CNA) = new_prefix
patient_index_score_pde = patient_index_score_pde %>% arrange(PDE_group)

cna_dist = 1/as.matrix(dist(patient_level_CNA))
#cna_dist = as.matrix(cor(t(patient_level_CNA)))

# cna_dist = log1p(cna_dist)
diag(cna_dist) = NA
# cna_dist = scale(cna_dist,scale = F, center = F)
Heatmap(cna_dist[patient_index_score_pde$patient, patient_index_score_pde$patient],
        cluster_rows = F, cluster_columns = F,
        column_split = patient_index_score_pde$PDE_group,
        row_split = patient_index_score_pde$PDE_group,
        )

same_gene_rate = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  tmp_cna = cna_gene[[prefix]][, same_gene]
  for(g in same_gene){
    gain_num = sum(tmp_cna[,g]>2)/nrow(tmp_cna)
    loss_num = sum(tmp_cna[,g]<2)/nrow(tmp_cna)
    same_gene_rate = rbind(same_gene_rate, c(new_prefix, tmp_group, g, gain_num, loss_num))
  }
}

# 统计基因数量
driver_gene = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/CNV_驱动基因.txt', sep='\t', header = T)
driver_gene = driver_gene[driver_gene$cancer_type_abbr=='BRCA', 'driver_gene']
driver_gene = strsplit(driver_gene, ', ')[[1]]
driver_gene = intersect(driver_gene, same_gene)
driver_gene = as.data.frame(driver_gene)

rownames(patient_index_score_pde) = patient_index_score_pde$patient

driver_gene_rate = c()
for(prefix in dna_file){
  print(prefix)
  new_prefix = gsub('uber.', '', prefix)
  new_prefix = gsub('.seg.txt', '', new_prefix)
  new_prefix = gsub('_cleaned_uber', '', new_prefix)
  new_prefix = gsub('_cleaned', '', new_prefix)

  tmp_group = patient_index_score_pde[new_prefix, 'PDE_group']
  tmp_cna = cna_gene[[prefix]][, driver_gene$driver_gene]
  for(g in driver_gene$driver_gene){
    gain_num = sum(tmp_cna[,g]>2)/nrow(tmp_cna)
    loss_num = sum(tmp_cna[,g]<2)/nrow(tmp_cna)
    driver_gene_rate = rbind(driver_gene_rate, c(new_prefix, tmp_group, g, gain_num, loss_num))
  }
}
driver_gene_rate = as.data.frame(driver_gene_rate)
colnames(driver_gene_rate) = c('patient', 'PDE_group', 'gene', 'gain_num', 'loss_num')
driver_gene_rate$gain_num = as.numeric(driver_gene_rate$gain_num)
driver_gene_rate$loss_num = as.numeric(driver_gene_rate$loss_num)

driver_gene_rate %>%
  group_by(gene, PDE_group) %>%
  mutate(gain_num=ifelse(PDE_group=='PDE_low', -gain_num,gain_num),
         loss_num=ifelse(PDE_group=='PDE_low', -loss_num,loss_num))%>%
  summarise(gain_num=mean(gain_num), loss_num=mean(loss_num)) %>%
  as.data.frame()%>%
  #mutate(gene = fct_reorder(gene, abs(gain_num+loss_num))) %>%
  mutate(gene=factor(gene, gene_fc$gene)) %>%
  ggplot()+
    geom_bar(aes(x=loss_num+gain_num, y=gene), stat='identity', fill='#FF3C55')+
    geom_bar(aes(x=loss_num, y=gene), stat='identity', fill='#2B86C5')+
  geom_vline(xintercept = 0, linewidth=2, color='white')+
  theme_classic()

gene_fc = c()
for(i in unique(driver_gene_rate$gene)){
  tmp = driver_gene_rate[driver_gene_rate$gene==i,]
  high = tmp[tmp$PDE_group=='PDE_high', ]
  low = tmp[tmp$PDE_group=='PDE_low', ]
  gene_fc = rbind(gene_fc, c(i, mean(high$gain_num) / mean(low$gain_num)))
}
gene_fc = as.data.frame(gene_fc)
colnames(gene_fc) = c('gene', 'FC')
gene_fc$FC = as.data.frame(gene_fc$FC)
gene_fc = gene_fc %>% arrange(FC)

#### gene 与ITH 相关性 #####
tcga_ITH = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                           sep='\t', header = T)

TCGA_bulk = read.table('/Volumes/WX_extend/甲状腺PTC/数据分析/TCGA_bulk/TCGA_THCA_rnaseq_FPKM.tsv', header = T, row.names = 'sample')
TCGA_bulk = log2(TCGA_bulk+1)
TCGA_bulk = TCGA_bulk[keep_gene,]










