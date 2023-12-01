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

########
medalt_dir = '/Volumes/WX_extend/BioSoftware/MEDALT/Crispr/'
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
  edu_dist = dist(tmp_data, method='euclidean')
  pyd = as.phyDat(as.data.frame(t(tmp_data)), type='USER', levels=min(tmp_data):max(tmp_data))
  # NJ
  NJ_tree<-NJ(edu_dist)
  # random_tree
  random_tree = rtree(nrow(tmp_data), tip.label = rownames(tmp_data))
  # MP
  MP_tree = optim.parsimony(random_tree, pyd)
  # ML
  ML_tree = optim.pml(pml(random_tree, pyd))$tree
  # sitka tree
  sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
  sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))

  # medalt
  medalt_tree = tryCatch({
    sitka_newick_to_phylo(paste0(medalt_dir, f, '/CNV.tree.txt'))
    },
    error = function(cond){
      return(NULL)
    })

  ## medicc2
  medicc2_tree = tryCatch({
    medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
    medicc2_tree = drop.tip(medicc2_tree, 'diploid')
    medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'),
                                    function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
  },
  error = function(cond){
    return(NULL)
  })

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

saveRDS(tree_rds, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/all_tree.rds')
tree_rds = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/all_tree.rds')
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
#saveRDS(res, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/all_tree_score.rds')

set.seed(123)
res %>%
  filter(leaves_num<20) %>%
  group_by(Method, Data) %>%
  summarise(ARIindex = mean(Score)) %>%
  mutate(Method=factor(Method, c('scTrace', 'NJ', 'MP', 'sitka', 'ML'))) %>%
  ggplot(aes_string(x='Method', y='ARIindex'))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(fill=Data),width = 0.2, shape=21, color='black', size=4)+
  scale_fill_manual(values = c('786O'='#2C86C8', 'A549'='#B7BC19', 'Huh7'='#E66327', 'U251'='#433691'))+
  theme_classic()+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format')


# res$leaves_num_cut = cut(res$leaves_num, breaks = seq(0, 100, 20))
ggplot(res, aes(x=leaves_num, y=Score, color=Method))+
  geom_line()+
  facet_wrap(~Data, scales = 'free')+
  scale_color_igv()+
  theme_classic()


ggplot(res, aes(x=Score, color=Method))+
  geom_density()+
  facet_wrap(~Data, scales = 'free')+
  scale_color_igv()+
  theme_classic()


res %>%
  #group_by(Method, Data) %>%
  #summarise(Dist = mean(Dist)) %>%
  ggboxplot(x='Method', y='Diversity', color='Method', add = 'jitter')+
  facet_wrap(~Data)+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format')

res$leaves_num_cut = as.vector(cut(res$leaves_num, breaks = seq(0, 2000, 100)))
res %>%
  group_by(leaves_num_cut,Method,Data) %>%
  group_by(mean_d = mean(Diversity), sd_d=sd(Diversity))


ggplot(res, aes(x=leaves_num_cut, y=Diversity, color=Method))+
  geom_violin()+
  scale_color_igv()


# 计算每个gRNA在tree的分布距离
res = c()
for(f in c('Huh7','A549')){
  # load seurat rds
  srt = readRDS(glue('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/{f}_filter_srt.RDS'))
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

  scTrace_tree_s1 = get_Label_dist(scTrace_tree, gRNA_label)
  NJ_tree_s1 = get_Label_dist(NJ_tree, gRNA_label)
  MP_tree_s1 = get_Label_dist(MP_tree, gRNA_label)
  ML_tree_s1 = get_Label_dist(ML_tree, gRNA_label)
  sitka_tree_s1 = get_Label_dist(sitka_tree, gRNA_label)

  plot_data = rbind(scTrace_tree_s1, NJ_tree_s1, MP_tree_s1, ML_tree_s1, sitka_tree_s1)
  plot_data = cbind(plot_data, c(rep('scTrace', nrow(scTrace_tree_s1)),
                                 rep('NJ', nrow(NJ_tree_s1)),
                                 rep('MP', nrow(MP_tree_s1)),
                                 rep('ML', nrow(ML_tree_s1)),
                                 rep('sitka', nrow(sitka_tree_s1))
  ))
  plot_data = as.data.frame(plot_data)
  colnames(plot_data) = c('gRNA', 'Dist', 'Method')
  plot_data$Dist = as.numeric(plot_data$Dist)
  plot_data$Data = f
  res = rbind(res, plot_data)
}

res %>%
  #group_by(Method, Data) %>%
  #summarise(Dist = mean(Dist)) %>%
  ggboxplot(x='Method', y='Dist', color='Method', add = 'jitter')+
  facet_wrap(~Data)+
  stat_compare_means(comparisons = list(c('scTrace', 'NJ')), label='p.format')




# plot_data$leaves_num_cut = cut(plot_data$leaves_num, breaks = seq(0, 100, 20))
# ggplot(plot_data, aes(x=leaves_num_cut, y=Diversity, color=Method))+
#   geom_boxplot()+
#   scale_color_igv()
#

# boxplot(scTrace_tree_s1, NJ_tree_s1, MP_tree_s1, ML_tree_s1)
# tree <- bird.families$tree
# # Assign tip labels to groups
# groups <- gRNA_label[scTrace_tree$tip.label]
# tip_colors <- ggsci::pal_igv()(length(unique(groups)))
# names(tip_colors) <- unique(groups)
# tip_colors <- tip_colors[groups]

# # Plot tree with colored tip labels
# plot(scTrace_tree,tip.color=tip_colors,  main="Example Tree")
# # Calculate dissimilarity matrix
# dist_matrix <- cophenetic(MP_tree)
# dist_matrix = reshape2::melt(dist_matrix)
# dist_matrix$True_label = gRNA_label[dist_matrix$Var1] == gRNA_label[dist_matrix$Var2]

# pROC::plot.roc(dist_matrix$True_label, as.numeric(dist_matrix$value),
#                print.auc=TRUE,
#                auc.polygon=TRUE,
#                grid=c(0.1, 0.2),
#                grid.col=c("green", "red"),
#                max.auc.polygon=TRUE,
#                auc.polygon.col="skyblue",
#                #print.thres=TRUE,
#                main='KNN train set ROC')
# # dist_matrix <- dist_matrix[lower.tri(dist_matrix)]

# # Perform Permanova test
#



ggplot(res, aes(x=leaves_num, y=Diversity, color=Method))+
  geom_line()+
  facet_wrap(~Data, scales = 'free')+
  scale_color_igv()+
  theme_classic()



# library(paletteer)
# paletteer_d("colorblindr::OkabeIto")
library(showtext)
# a = font.files()
# a[grep('Times', a$file),]

font_add('Times New Roman', '/System/Library/Fonts/Supplemental/Times New Roman.ttf')
showtext_auto()
#a = t(apply(RF_data, 1, function(x) x[1:8] / max(x[1:8],na.rm = TRUE))) %>%
a = t(apply(RF_data, 1, function(x) x[1:8])) %>%
  reshape2::melt() %>%
  filter(Var2!='hclust') %>%
  mutate(type=sapply(Var2, function(x)ifelse(x=='scCNAT', 'scCNAT', 'other')))%>%
  mutate(variable = factor(Var2, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  #subset(variable%in%c('scCNAT',  'NJ')) %>%
  ggplot(aes(x=variable, y=value, color=type))+
  geom_boxplot(outlier.shape=NA, width=0.4)+
  #geom_jitter(width = 0.06,size=0.1)+
  scale_color_manual(values = c('scCNAT'='#D55E00FF', 'other'='#999999FF'))+
  stat_compare_means(comparisons = list(c('scCNAT',  'MEDICC2')),
                     label='p.format', label.y = 200,size=3, family='Times New Roman')+
  labs(title='Different methods', y='Robinson-Foulds distance')+
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

RF_data %>%
  reshape2::melt(id='cell_num') %>%
  filter(variable!='hclust') %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%
  #subset(variable%in%c('scCNAT',  'NJ')) %>%
  ggplot(aes(x=as.factor(cell_num), y=value, fill=variable))+
  geom_boxplot(outlier.shape=NA, width=0.4)+
  #geom_jitter(width = 0.06,size=0.1)+
  scale_fill_igv()+
  labs(title='Different methods', y='Robinson-Foulds distance')+
  theme_classic()+
  theme(#legend.position = 'none',
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
    axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
    axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
    axis.line = element_line(size=0.1),
    axis.ticks = element_line(size=0.1),
  )

RF_data %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('MP','NJ','hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=1, fill='white', shape=21) +
  labs(title='Different methods', y='Robinson-Foulds distance')+
  scale_color_igv()+
  theme_classic()+
  theme(#legend.position = 'none',
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
    axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
    axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
    axis.line = element_line(size=0.1),
    axis.ticks = element_line(size=0.1),
  )

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图//RFDist.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')
library(ggsci)
Quartet_res %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('MP','NJ','hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=1, fill='white', shape=21) +
  labs(title='Different methods', y='Quartet distance')+
  scale_color_igv()+
  theme_classic()+
  theme(#legend.position = 'none',
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.title = element_text(family='Times New Roman', face = 'bold', size = 8, color = 'black'),
    axis.text.y = element_text(family='Times New Roman',  size = 6, color = 'black'),
    axis.text.x = element_text(family='Times New Roman', size = 8, color = 'black', angle = 30, hjust=1, vjust = 1),
    axis.line = element_line(size=0.1),
    axis.ticks = element_line(size=0.1),
  )



#subset(variable%in%c('scCNAT',  'NJ')) %>%
ggplot(aes(x=cell_num, y=value, color=variable))+
  geom_boxplot(outlier.shape=NA, width=0.4)+
  #geom_jitter(width = 0.06,size=0.1)+
  scale_color_manual(values = c('scCNAT'='#D55E00FF', 'other'='#999999FF'))+
  #stat_compare_means(comparisons = list(c('scCNAT',  'MEDICC2')),
  #                   label='p.format', label.y = 200,size=3, family='Times New Roman')+
  labs(title='Different methods', y='Robinson-Foulds distance')+
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





RF_data %>%
  reshape2::melt(id.vars= 'cell_num') %>%
  mutate(cells = cut(cell_num, breaks = c(0,100,200,300,500,1000))) %>%
  mutate(variable = factor(variable, levels=c('scCNAT','MP','sitka', 'NJ',  'hclust', 'MEDALT','ML'))) %>%
  #subset(variable%in%c('scCNAT',  'NJ')) %>%
  ggboxplot(x='cells', y='value',color = 'variable', add='jitter', palette = 'npg')+
  stat_compare_means(comparisons = list(c('scCNAT',  'NJ')), label='p.format')+
  labs(x='Cell number', y='Distance')+
  theme(legend.position = 'right')



# SPR_data %>%
#   reshape2::melt(id.vars= 'cell_num') %>%
#   mutate(cells = cut(cell_num, breaks = c(0,100,200,300,500,1000))) %>%
#   # subset(variable%in%c('scCNAT',  'NJ')) %>%
#   #mutate(variable = factor(variable, levels=c('scCNAT', 'NJ',  'hclust','MEDALT', 'ML'))) %>%
#   ggboxplot(x='cells', y='value',color = 'variable', add='jitter', palette = 'npg')+
#   #stat_compare_means(aes(group=variable), label='p.signif') +
#   labs(x='Cell number', y='Distance')+
#   theme(legend.position = 'right')
##########
data = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_tree_F1F2/F1F2_data.txt'
sim_data = read.table(data, header = T)
sim_data = t(sim_data[,3:ncol(sim_data)])
cn_col = structure(c('#FFFFFF',
                     '#91FF88',
                     '#C6C6C6',
                     '#FFEB97',
                     '#FCCC00',
                     '#ec9336',
                     '#7d170e'
), names=c(0,1,2,3,4,5,6))

Heatmap(as.matrix(sim_data),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        #column_split = chr_cluster,
        heatmap_legend_param = list(
          title = 'CN',
          color_bar = "discrete",
          title_position = "leftcenter-rot", # 图例标题位置
          legend_height = unit(3, "cm") #图例长度
        ),
        #top_annotation = top_color,
        #left_annotation = left_ann,
        column_title = "Genomic Region",
        column_title_side = c("bottom"),
        col=cn_col)
##########




# 进化树
pyd = as.phyDat(sim_data, type='USER', levels=c(0,1,2,3,4,5,6))
mm = dist(sim_data) # as.matrix(dist(t(sim_data), upper = T))
NJ_tree<-NJ(mm) # NJ
plot(NJ_tree)
# random tree
random_tree = rtree(ncol(sim_data), tip.label = colnames(sim_data))
# MP

MP_tree = optim.parsimony(random.addition(pyd), pyd)
plot(MP_tree)
# ML
ML_tree = optim.pml(pml(random_tree, pyd))$tree
plot(ML_tree)
# rt
rt = newick_to_phylo(ref_tree)
###
TBRDist::TBRDist(unroot(rt), unroot(it), exact = TRUE)

TreeDist::JaccardRobinsonFoulds(unroot(rt), unroot(it))

treedist(rt, MP_tree)
treedist(rt, ML_tree)
treedist(rt, NJ_tree)
treedist(rt, it)


RF.dist(rt, optim.parsimony(NJ_tree, pyd))

RF.dist(rt, optim.parsimony(random_tree, pyd))
RF.dist(rt, ML_tree)
RF.dist(rt, optim.pml(pml(upgma(mm), pyd))$tree)

RF.dist(rt, optim.pml(pml(random_tree, pyd))$tree)

RF.dist(rt, optim.pml(pml(upgma(mm), pyd))$tree)


treedist(rt, MP_tree)
treedist(rt, it)

treedist(rt, ML_tree_fit$tree)

SPR.dist(unroot(rt), unroot(MP_tree))
#########
##########
SPR.dist(rt, it)

ML_tree<-pml(upgma(mm), pyd)
ML_tree_fit = optim.pml(pml(upgma(mm), pyd))
as.phylo(ML_tree_fit)



ggtree(ML_tree_fit$tree)



#########
x = rtree(6)
y = rtree(6)

y$tip.label = paste0('y_',y$tip.label)

plot(x+y)
RF.dist(x, y)



rt = newick_to_phylo(ref_tree)
it = newick_to_phylo(infer_tree)
cluster_tree = as.phylo.hclust(hclust_dist)
# Robinson-Foulds
RF.dist(rt, it)
RF.dist(rt, cluster_tree)

treedist(rt, it)
sprdist(rt, it)
SPR.dist(rt, it)

treedist(rt, cluster_tree)
sprdist(rt, cluster_tree)
SPR.dist(rt, cluster_tree)


###########
data = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_tree_100/SimID_1654561438.txt'
ref_tree = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_tree_100/SimID_1654561438.newick'
infer_tree = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_tree_100/SimID_1654561438_infer_SClone_tree.txt'
ref_den = newick_to_hclust(ref_tree)
infer_den = newick_to_hclust(infer_tree)

# cluster
sim_data = read.table(data, header = T)
sim_data = sim_data[,3:ncol(sim_data)]
edu_dist = dist(t(sim_data), method = "euclidean")
hclust_dist = hclust(edu_dist, method = "complete")
pheatmap::pheatmap(t(sim_data)[hclust_dist$order,],
                   show_rownames = F, cluster_cols = F, cluster_rows = T)
##

infer_auc = c()
for(k in 2:80){
  ref_sub = cutree(ref_den, k=k)
  infer_sub = cutree(infer_den, k=k)
  hclust_sub = cutree(hclust_dist, k=k)
  # rect.hclust(infer_den, k = 4, border = 1:4)
  auc1 = as.numeric(multiclass.roc(ref_sub,infer_sub[names(ref_sub)])$auc)
  auc2 = as.numeric(multiclass.roc(ref_sub,hclust_sub[names(ref_sub)])$auc)
  infer_auc=rbind(infer_auc, c(auc1,auc2))
}

boxplot(infer_auc)

infer_auc = as.data.frame(infer_auc)
colnames(infer_auc) = c('scCNAT', 'hclust')
infer_auc$cluster_num = 2:80

ggplot(infer_auc) +
  geom_line(mapping = aes(x=cluster_num, y=scCNAT, color='scCNAT'))+
  geom_line(mapping = aes(x=cluster_num, y=hclust, color='hclust')) +
  scale_color_manual(values=c('scCNAT'='red', 'hclust'='blue'))

