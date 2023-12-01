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
################################################################################
file_name = 'huh7'
file_name2 = 'huh7.txt'
scTrace_out = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine/'
# 1.load data
mydata_rds= list()
for(file_name in c('huh7.txt',
                     'plc.txt',
                     'hep3b.txt',
                     'smc7721.txt',
                     'hepg2.txt',
                     'mhcc97h.txt',
                     'mhcc97l.txt',
                     'lm3.txt',
                     'integrate_lm3.txt')){
  prefix = strsplit(file_name, '\\.')[[1]][1]
  sctc_obj = scTraceClass(scTrace_out, prefix, min_cell_num = 100)
  gene_cna_data = read.table(paste0('/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/', file_name,'_gene_data.csv'))
  rownames(gene_cna_data) = gsub('\\.', '-', rownames(gene_cna_data))
  sctc_obj = add_cna_score(sctc_obj, pathway_list, gene_cna_data)
  sctc_obj = predict_mitosis_rate(sctc_obj)
  sctc_obj$cell_pos_obj$cell_map_pos$Fn = sapply(sctc_obj$cell_pos_obj$cell_map_pos$cell_name, function(x)strsplit(x, '_')[[1]][1])
  mydata_rds[[file_name]]=sctc_obj
}

all_CL_res = c()
for(file_name in c('huh7.txt',
                   'plc.txt',
                   'hep3b.txt',
                   'smc7721.txt',
                   'hepg2.txt',
                   'mhcc97h.txt',
                   'mhcc97l.txt',
                   'lm3.txt',
                   'integrate_lm3.txt')){
  tmp_obj = mydata_rds[[file_name]]
  tmp_cell_data = tmp_obj$cell_pos_obj$cell_map_pos
  tmp_cell_data$CL = strsplit(file_name, '\\.')[[1]][1]
  all_CL_res = rbind(all_CL_res, tmp_cell_data)
}

## F1，F2 clone水平比例
all_CL_res = as.data.frame(all_CL_res)

all_CL_res %>%
  group_by(Fn, CL) %>%
  summarise(CellCycle_OG = mean(CellCycle_OG)) %>%
  ggboxplot(x='Fn', y='CellCycle_OG', add='jitter')+
  stat_compare_means(comparisons = list(c('f1', 'f2')))
  
all_CL_res %>%
  group_by(Fn,CL) %>%
  summarise(pred_rate = mean(pred_rate)) %>%
  ggboxplot(x='Fn', y='pred_rate', add='jitter', color='CL')+
  stat_compare_means(comparisons = list(c('f1', 'f2')))

library(ggpubr)
data = all_CL_res %>%
  group_by(Fn, CL) %>%
  summarise(CellCycle_OG = mean(pred_rate)) %>% as.data.frame() %>%
  na.omit()
compare = list(c('f1', 'f2'))
pd1 = as.data.frame(cbind(data[data$Fn=='f1', 3], data[data$Fn=='f2', 3]))
colnames(pd1) = c('f1', 'f2')
ggpaired(pd1, cond1 = "f1", cond2 = "f2",
             color = "condition", palette = "jco",title = 'Cluster_1') +
  stat_compare_means(comparisons = compare, method = "wilcox.test",size=6,paired = TRUE,
                     label = "p.format",angle=0)



plot_cells(mydata_rds[['huh7.txt']]$cell_pos_obj,  type='pie',colorby='Fn',shadow=TRUE, clone_line=TRUE,
           point.size=4, plot_virtual=F)+
  ggsci::scale_fill_igv()


packdata = sctms_obj$packdata
packdata[, 'from'] = sctms_obj$node_info[packdata$name, 'Class']
pie_cols = setdiff(unique(packdata$from), NA)
pie_cols = setdiff(pie_cols, 'virtual')
pie_data = packdata %>%
  group_by(group, from) %>%
  summarise(num = n()) %>% 
  filter(from %in% pie_cols)
pie_data = reshape2::dcast(pie_data, group~from)
pie_data[is.na(pie_data)] = 0
pie_data %>%
  mutate(f1=f1/(f2+f1), f2=f2/(f2+f1)) %>%
  ggplot(aes(x=f2,y=f1, color=group))+
  geom_point()+
  scale_color_igv()


plot_cells(sctc_obj$cell_pos_obj, colorby='clone',shadow=TRUE, clone_line=TRUE,
           point.size=4, plot_virtual=F)+
  ggsci::scale_fill_igv()



plot_cna_tree(sctc_obj$oirg.data$tree, sctc_obj$oirg.data$all_node_data,
              colorby = 'Fn',
              cell_info=sctc_obj$cell_pos_obj$cell_map_pos,
              output = paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/Huh7_tree_heatmap_CL.pdf')
)


# load gene level cna data

gene_cna_data = read.table(paste0('/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/', file_name,'_gene_data_all_node.csv'))
rownames(gene_cna_data) = gsub('\\.', '-', rownames(gene_cna_data))


# 4.add score
sctms_obj = add_cna_score(sctms_obj, gene_cna_data, features = pathway_list, add_clone=TRUE)

# 5.show tree
sctms_obj$node_info$Class = sapply(rownames(sctms_obj$node_info), function(x)strsplit(x, '_')[[1]][1])
sctms_obj$node_info$Class[rownames(sctms_obj$node_info)=='root'] = 'root'
sctms_obj$node_info$Clone = rownames(sctms_obj$node_info)

clone_color = ggsci::pal_igv()(nrow(sctms_obj$clone_info))
sctms_obj$clone_info$color = clone_color
names(clone_color) = sctms_obj$clone_info$clone
a=show_tree(sctms_obj, 
            display_type = 'clone',
            add_arrow = FALSE,
            color.by = 'clone',
            size.by = 'scale_size',
            show_clone_name = TRUE,
) +scale_color_manual(values = clone_color)+
  ggexpand(.1,  side = "hv")+
  theme(legend.position = 'non')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_clone_tree.pdf'), 
       a, width=120, height=120, units='mm', dpi = 600, bg = 'transparent')
#+ scale_color_gradientn(colours=c('white', 'red'))
b=show_tree(sctms_obj,
            display_type = 'piescatter',
            add_arrow = FALSE,
            color.by = 'Class',
            show_virtual=FALSE
) +scale_fill_manual(values = c('f1'='#60BD68FF', 'f2'='#F15854FF'))
b 
#c = cowplot::plot_grid(a,b, ncol = 1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_clone_tree_pie.pdf'), 
       b, width=120, height=120, units='mm', dpi = 600, bg = 'transparent')
# 7 branch analyze
get_children <- function(edges, node){
  res = c(node)
  num = 1
  while(num<=length(res)){
    curr_node = res[num]
    num = num + 1
    tmp_child = edges[edges[,1]==curr_node, 2]
    res = c(res, tmp_child)
  }
  return(res)
}
#################################################################################
newick_to_phylo <- function(newick_dir){
  newick_file = file(newick_dir,open='r')
  newick_str = readLines(newick_file, n = 1)
  close(newick_file)
  newick_str = paste0('(', newick_str, ');')
  newick_str <- read.tree(text = newick_str)
  return(newick_str)
}
# pseudo analyze
pseudo_res = sctms_obj$pseudo
rownames(pseudo_res) = pseudo_res$name
pseudo_res = pseudo_res[, c(-1,-2)]
cna_data = sctms_obj$all_nodes_data
rownames(cna_data) = cna_data[,1]
cna_data = cna_data[, -1]

#file_name = 'huh7.txt'
tree_dir = paste0('/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_keyCNA/' , file_name , '_infer_tree.txt')
tree_res = newick_to_phylo(tree_dir)

gtree_res = ggtree(tree_res, size=0.1)
gtree_data = gtree_res$data %>% as.data.frame()
leaves_data = cna_data[gtree_data[gtree_data$isTip==TRUE,'label'], ]

cn_col = c('blue','#91FF88', '#C6C6C6',  '#FFEB97',  '#FCCC00', '#ec9336', '#7d170e', 'darkred')
names(cn_col) = c('0', '1','2','3','4','5','6', '7')
gh_data = apply(leaves_data, 2, as.character)
rownames(gh_data) = rownames(leaves_data)

gtree_res$data$Fn = sapply(gtree_res$data$label, function(x)strsplit(x, '_')[[1]][1])
gtree_res = gtree_res+
  #theme_tree2() +
  geom_tippoint(aes(color=Fn), size=0.05)+
  scale_color_manual(values = c('f1'='#60BD68FF', 'f2'='#F15854FF'))

a = gheatmap(gtree_res, gh_data, width=5, offset=5, colnames=F,color=NA) +
  #scale_x_ggtree() +
  scale_fill_manual(values = cn_col) + 
  scale_y_continuous(expand=c(0, 1)) + 
  # xlab("Time") +
  theme(legend.text=element_text(size=8), 
        legend.key.height=unit(.2, "cm"),
        legend.key.width=unit(.3, "cm"), 
        #legend.position=c(.13, y=.945),
        #axis.text.x=element_text(size=10, angle = 90,vjust = 0), 
        axis.text.x=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
        #axis.title.x = element_text(size=12)
  )
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_tree_heatmap.pdf'), 
       a, width=120, height=60, units='mm', dpi = 600, bg = 'transparent')


## F1，F2 clone水平比例
packdata = sctms_obj$packdata
packdata[, 'from'] = sctms_obj$node_info[packdata$name, 'Class']
pie_cols = setdiff(unique(packdata$from), NA)
pie_cols = setdiff(pie_cols, 'virtual')
pie_data = packdata %>%
  group_by(group, from) %>%
  summarise(num = n()) %>% 
  filter(from %in% pie_cols)
pie_data = reshape2::dcast(pie_data, group~from)
pie_data[is.na(pie_data)] = 0
pie_data %>%
  mutate(f1=f1/(f2+f1), f2=f2/(f2+f1)) %>%
  ggplot(aes(x=f2,y=f1, color=group))+
  geom_point()+
  scale_color_igv()


# clone分支 水平
clone_pseudo_res = pseudo_res
clone_pseudo_res$clone = sctms_obj$node_info[rownames(clone_pseudo_res), 'clone']
branch1 = get_children(sctms_obj$edges, 'clone_1')
branch2 = c('clone_0', get_children(sctms_obj$edges, 'clone_2'))
branch3 =  get_children(sctms_obj$edges, 'clone_3')
clone_pseudo_res$branch = sapply(clone_pseudo_res$clone, function(x){
  if(x %in% branch1){
    return('branch1')
  }else if(x %in% branch2){
    return('branch2')
  }else if(x %in% branch3){
    return('branch3')
  }else{
    return('root')
  }
})
clone_pseudo_res = clone_pseudo_res[-1, ]

na_row = rowSums(clone_pseudo_res=='none') >0
clone_pseudo_res[!na_row, 'Mitosis_loc'] =  as.numeric(clone_pseudo_res$split_ad_loc[!na_row])#as.numeric(clone_pseudo_res$split_dd_loc[!na_row]) +
clone_pseudo_res[!na_row, 'Mitosis_rate'] = clone_pseudo_res[!na_row, 'Mitosis_loc'] / as.numeric(clone_pseudo_res$split_time[!na_row])

clone_pseudo_res$Root_loc =   as.numeric(clone_pseudo_res$Root_gain_loc) +   as.numeric(clone_pseudo_res$Root_loss_loc) 
clone_pseudo_res$Parent_loc = as.numeric(clone_pseudo_res$Parent_gain_loc) + as.numeric(clone_pseudo_res$Parent_loss_loc) 
clone_pseudo_res$Root_cn =    as.numeric(clone_pseudo_res$Root_gain_cn) +    as.numeric(clone_pseudo_res$Root_loss_cn) 
clone_pseudo_res$Parent_cn =  as.numeric(clone_pseudo_res$Parent_gain_cn) +  as.numeric(clone_pseudo_res$Parent_loss_cn) 


clone_pseudo_res_melt = reshape2::melt(clone_pseudo_res, id=c('clone', 'branch'))
clone_pseudo_res_melt = clone_pseudo_res_melt[clone_pseudo_res_melt$value!='none', ]
clone_pseudo_res_melt$value = as.numeric(clone_pseudo_res_melt$value)

clone_pseudo_res_melt$variable = gsub('split', 'Mitosis', clone_pseudo_res_melt$variable)

a = clone_pseudo_res_melt %>%
  na.omit() %>%
  filter(variable!='Mitosis_time2') %>%
  filter(variable %in%c("Root_loc", "Root_cn", "Parent_loc", "Parent_cn", 
                        'Mitosis_loc', 'Mitosis_rate','pseudotime','Mitosis_time')) %>%
  mutate(variable = factor(variable, levels = c("Root_loc", "Root_cn", "Parent_loc", "Parent_cn", 
                                                'Mitosis_loc', 'Mitosis_rate','pseudotime','Mitosis_time'))) %>%
  ggplot(aes(x=branch, y=value, color=branch))+
  geom_boxplot()+
  geom_jitter(width=0.1)+
  ggsci::scale_color_igv()+
  facet_wrap(~variable, scales = 'free_y', ncol=4)+
  stat_compare_means(comparisons = list(c('branch1', 'branch2'),
                                        c('branch1', 'branch3'),
                                        c('branch2', 'branch3')), label='p.format')+
  theme_bw()
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch.pdf'), 
       a, width=180, height=100, units='mm', dpi = 600, bg = 'transparent')

# pseudotime
a = clone_pseudo_res_melt %>%
  na.omit() %>%
  filter(variable=='pseudotime') %>%
  ggplot(aes(x=branch, y=value, color=branch))+
  geom_boxplot(width=0.5, outlier.size = 0.1,size=0.5)+
  #geom_jitter(width=0.01)+
  ggsci::scale_color_igv()+
  stat_compare_means(comparisons = list(c('branch1', 'branch2'),
                                        c('branch1', 'branch3'),
                                        c('branch2', 'branch3')), label='signif',size=2)+
  labs(y='Pseudotime')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
        axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
        axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1),
  )
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_pseudotime.pdf'), 
       a, width=40, height=60, units='mm', dpi = 600, bg = 'transparent')

# pseudotime
a = clone_pseudo_res_melt %>%
  na.omit() %>%
  filter(variable=='Mitosis_loc') %>%
  ggplot(aes(x=branch, y=value, color=branch))+
  geom_boxplot(width=0.5, outlier.size = 0.1,size=0.5)+
  #geom_jitter(width=0.01)+
  ggsci::scale_color_igv()+
  stat_compare_means(comparisons = list(c('branch1', 'branch2'),
                                        c('branch1', 'branch3'),
                                        c('branch2', 'branch3')), label='signif',size=2)+
  labs(y='Mitosis loc')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
        axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
        axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1),
  )
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_mitosis.pdf'), 
       a, width=40, height=60, units='mm', dpi = 600, bg = 'transparent')

# 评分比较
cna_score = sctms_obj$cna_score
cna_score$clone = sctms_obj$node_info[rownames(cna_score), 'clone']
cna_score$branch = sapply(cna_score$clone, function(x){
  if(x %in% branch1){
    return('branch1')
  }else if(x %in% branch2){
    return('branch2')
  }else if(x %in% branch3){
    return('branch3')
  }else{
    return('root')
  }
})

tmp_data = reshape2::melt(cna_score, id=c('clone', 'branch')) %>%
  filter(clone!='root') %>%
  filter(variable %in% grep('_OG',variable, value = T))

require(gridExtra)
require(dplyr)

a=tmp_data %>%
  group_by(branch,variable)%>%
  mutate(mean_score=median(value)) %>%
  group_by(variable) %>% 
  do(gg = {ggplot(., aes(branch, value, fill = mean_score)) + 
      #geom_boxplot() + 
      geom_violin()+
      #geom_jitter(width=0.1,size=0.1)+
      facet_wrap(~variable)+
      scale_fill_gradientn(colours = c('yellow', 'red'))+
      stat_compare_means(comparisons = list(c('branch1', 'branch2'),
                                            c('branch1', 'branch3'),
                                            c('branch2', 'branch3')), label='p.format',
                         size=3, label.y.npc='bottom')+
      theme_bw()+
      theme(legend.position = 'none')
  }) %>% 
  .$gg %>% arrangeGrob(grobs = ., nrow = 3) %>% grid.arrange()
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_score.pdf'), 
       a, width=180, height=200, units='mm', dpi = 600, bg = 'transparent')
# HIPPO
a=tmp_data %>%
  filter(variable=='HIPPO_OG')%>%
  group_by(branch,variable)%>%
  mutate(mean_score=mean(value)) %>%
  group_by(variable) %>% 
  ggplot(., aes(branch, value, fill = mean_score)) + 
  #geom_boxplot() + 
  geom_violin(size=0.1)+
  #geom_jitter(width=0.1,size=0.1)+
  scale_fill_gradientn(colours = c('yellow', 'red'))+
  stat_compare_means(comparisons = list(c('branch1', 'branch2'),
                                        c('branch1', 'branch3'),
                                        c('branch2', 'branch3')), label='p.format',
                     size=2, label.y.npc='bottom')+
  theme_classic()+
  labs(y='HIPPO_OG')+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
        axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
        axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1),
  )
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_HIPPO.pdf'), 
       a, width=40, height=60, units='mm', dpi = 600, bg = 'transparent')


# branch analyze
key_events = read.table(paste0('/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_keyCNA/',file_name,'_events.txt'),
                        sep = ',', header = T)

tmp_event_data = subset(key_events, clone_name1=='clone_0'& clone_name2=='clone_1')
colors = c('1'='#91FF88', '2'='#C6C6C6', '3'='#FFEB97', '4'='#FCCC00', '5'='#ec9336', '6'='#7d170e', '7'='darkred')

pheatmap::pheatmap(t(tmp_event_data[,c('clone1', 'clone2')]),
                   cluster_rows = F, cluster_cols = F,
                   show_colnames = F,
                   cellwidth = 6, cellheight = 6,
                   border_color = '#dddddd',
                   color = colors,
                   annotation_col = tmp_event_data[, 'type', drop=F],
                   annotation_colors = list(type=c('copy'="#9467BDFF", 'aneu'="#D62728FF", 'rearrange'="#2CA02CFF"))
)



library(reshape2)
score_res = sctms_obj$cna_score
pseudo_res = sctms_obj$pseudo
rownames(pseudo_res) = pseudo_res$name
score_res$mitosis1 = pseudo_res[rownames(score_res),]$split_ad_loc 
score_res$mitosis2 = pseudo_res[rownames(score_res),]$split_dd_loc 
score_res$split_time = pseudo_res[rownames(score_res),]$split_time

score_res$clone = sctms_obj$node_info[rownames(score_res), 'clone']
score_res$branch = sapply(score_res$clone, function(x){
  if(x %in% branch1){
    return('branch1')
  }else if(x %in% branch2){
    return('branch2')
  }else if(x %in% branch3){
    return('branch3')
  }else{
    return('root')
  }
})
score_res = melt(score_res, id=c('mitosis1','mitosis2', 'split_time','clone', 'branch'))
score_res$from=file_name
#x = score_res %>%
#  filter(mitosis1!='none') %>%
#  mutate(mitosis1 = as.numeric(mitosis1),
#         mitosis2 = as.numeric(mitosis2),
#         mitosis = mitosis1+mitosis2)%>%
#  group_by(branch, variable) %>%
#  summarise(score=mean(value), mitosis=mean(mitosis)) %>% as.data.frame()
#x$from = file_name
write.table(score_res, paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_score.txt'), quote = F)











## 分裂速率预测
model_dataset = pseudo_res[-1, ]
na_pos = rowSums(model_dataset=='none')==0

train_data = data.frame('root_loc'=as.numeric(model_dataset$Root_gain_loc)+as.numeric(model_dataset$Root_loss_loc),
                        'root_cn' = as.numeric(model_dataset$Root_gain_cn)+as.numeric(model_dataset$Root_loss_cn),
                        'parent_loc' = as.numeric(model_dataset$Parent_gain_loc)+as.numeric(model_dataset$Parent_loss_loc),
                        'parent_cn' = as.numeric(model_dataset$Parent_gain_cn)+as.numeric(model_dataset$Parent_loss_cn),
                        'pseudotime' = as.numeric(model_dataset$pseudotime),
                        'self_time' = as.numeric(model_dataset$split_time2)
)


train_d = train_data[na_pos, ]
train_label = model_dataset[na_pos, c('split_dd_loc', 'split_ad_loc', 'split_time')]
train_label = apply(train_label, 2, as.numeric)
train_label = (train_label[,1]+train_label[,2]) / train_label[,3]
train_d$label = train_label
train_d[is.na(train_d)] = 0
rownames(train_d) = rownames(model_dataset[na_pos, 1:8])


# test
test_data = train_data[!na_pos, ]
rownames(test_data) = rownames(model_dataset[!na_pos, ])

# 线性模型
# fm = label~Root_gain_loc+Root_loss_loc+Root_gain_cn+Root_loss_cn+Parent_gain_loc+Parent_loss_loc+Parent_gain_cn+Parent_loss_cn
fm = label~root_loc+root_cn+parent_loc+parent_cn+pseudotime#+self_time
#fm = label~root_loc+parent_loc+pseudotime

# fm = label~root_cn+pseudotime
#lm_model=lm(data = train_d, formula = fm)
#summary(lm_model)

glm_model=glm(data = train_d, formula = fm)
summary(glm_model)
glm_pred = predict(glm_model, train_d)
#plot(train_d$label, glm_pred)
#abline(a=0, b=1)

glm_pred = as.data.frame(predict(glm_model, test_data))
colnames(glm_pred) = 'pred_rate'
glm_pred$Fn = sapply(rownames(glm_pred), function(x)strsplit(x,'_')[[1]][1])
glm_pred$clone = sctms_obj$node_info[rownames(glm_pred), 'clone']
glm_pred$branch = sapply(glm_pred$clone, function(x){
  if(x %in% branch1){
    return('branch1')
  }else if(x %in% branch2){
    return('branch2')
  }else if(x %in% branch3){
    return('branch3')
  }else{
    return('root')
  }
})

#glm_pred %>%
#  ggplot(aes(x=Fn, y=pred_rate))+
#  geom_boxplot()+
#  geom_jitter(width = 0.1)
a=glm_pred %>%
  ggplot(aes(x=branch, y=pred_rate))+
  geom_boxplot()+
  geom_jitter(width = 0.1)
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_branch_rate_pred.pdf'), 
       a, width=100, height=80, units='mm', dpi = 600, bg = 'transparent')

clone_rate = glm_pred %>% group_by(clone) %>% summarise(clone_rate = mean(pred_rate)) %>% as.data.frame()
rownames(clone_rate) = clone_rate$clone
clone_rate$branch = sapply(clone_rate$clone, function(x){
  if(x %in% branch1){
    return('branch1')
  }else if(x %in% branch2){
    return('branch2')
  }else if(x %in% branch3){
    return('branch3')
  }else{
    return('root')
  }
})
clone_rate %>%
  ggplot(aes(x=branch, y=clone_rate))+
  geom_boxplot()+
  geom_jitter(width = 0.1)

sctms_obj$clone_info$clone_rate = clone_rate[sctms_obj$clone_info$clone, 'clone_rate']
a=show_tree(sctms_obj, 
            display_type = 'clone',
            add_arrow = FALSE,
            color.by = 'clone_rate',
            size.by = 'scale_size',
            show_clone_name = FALSE,
) +scale_color_gradientn(colours = c('yellow', 'red'))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/mydata/',file_name,'_clone_rate_pred.pdf'), 
       a, width=160, height=100, units='mm', dpi = 600, bg = 'transparent')

a=show_tree(sctms_obj, 
            display_type = 'points',
            add_arrow = FALSE,
            #color.by = 'clone_rate',
            #size.by = 'scale_size',
            show_clone_name = FALSE,
) +scale_color_gradientn(colours = c('yellow', 'red'))
a

show_tree(sctms_obj,
          display_type = 'points',
          color.by = 'Class',
          alpha = 0.6,
          size=1,
          add_arrow = FALSE,
          add_mark = F
)









