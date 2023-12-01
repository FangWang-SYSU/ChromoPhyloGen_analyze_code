library(ggplot2)
library(ggpubr)
library(treeio)
library(ape)
library(pROC)
library(phangorn)
library(phangorn)
library(dplyr)
library(ComplexHeatmap)
library(jsonlite)
library(igraph)
library(ggtree)
library(vegan)
library(pdfCluster)
library(ggsci)

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
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}
get_clone <- function(tree, min_k=5){
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
clone_accuracy <- function(ref_tree, infer_tree, clone_num=5){
  rt_clone = get_clone(ref_tree, min_k = clone_num)
  it_clone = get_clone(infer_tree, min_k = clone_num)
  rt_clone$it_clone = it_clone[rownames(rt_clone), 'clone']
  tmp_func = function(x){
    diversity(table(x))
  }
  res = rt_clone %>%
    group_by(clone) %>%
    summarise(d = tmp_func(it_clone))

  return(mean(res$d))
}
clone_ARI <- function(ref_tree, infer_tree, clone_num=NULL){
  if(is.null(clone_num)){
    cn = 10
    tmp_all_res = c()
    for(i in 2:cn){
      rt_clone = get_clone(ref_tree, min_k = i)
      it_clone = get_clone(infer_tree, min_k = i)
      rt_clone$it_clone = it_clone[rownames(rt_clone), 'clone']
      tmp_all_res = c(tmp_all_res, adj.rand.index(rt_clone$clone, rt_clone$it_clone))
    }
    res = mean(tmp_all_res)
  }else{
    rt_clone = get_clone(ref_tree, min_k = clone_num)
    it_clone = get_clone(infer_tree, min_k = clone_num)
    rt_clone$it_clone = it_clone[rownames(rt_clone), 'clone']
    res = adj.rand.index(rt_clone$clone, rt_clone$it_clone)
  }
  return(res)
}
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
####### 1.运行软件 #######




####### 2.比较RF距离 #######
library(Quartet)
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num/'
scTrace_dir='/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sim_for_RFdist/'
medalt_dir = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_for_RFdist/'
sitka_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_sitka/'
medicc_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_medicc2/'
sim_100_dir_res = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_res/'
fl = sapply(list.files(sim_100_dir), function(x)substring(x,1,16))
fl = unique(fl)
res = c()
Quartet_res = c()
for(f in fl){
  print(f)
  tmp_infer = paste0(scTrace_dir, f, 'cell_tree.newick')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    tmp_data = paste0(sim_100_dir, f, '.txt')
    rt = newick_to_phylo(tmp_ref)
    #it = newick_to_phylo(tmp_infer)
    scTrace_res = load_scTrace(scTrace_dir, f)
    it = scTrace_res$tree

    medalt_tree = sitka_newick_to_phylo(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    # sitka tree
    sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
    sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))
    #medicc2
    if(file.exists(paste0(medicc_dir,f,'/',f,'_final_tree.new'))){
      medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
      medicc2_tree = drop.tip(medicc2_tree, 'diploid')
      medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'), function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
    }else{
      medicc2_tree = rtree(10)
    }

    #
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]
    edu_dist = dist(t(sim_data), method = "manhattan")
    hclust_dist = hclust(edu_dist, method = "average")
    cluster_tree = as.phylo.hclust(hclust_dist)
    # NJ
    # pyd = as.phyDat(sim_data, type='USER', levels=0:max(sim_data))
    mm = dist(t(sim_data), method = "maximum") #dist(t(sim_data)) # dist.hamming(pyd) # as.matrix(dist(t(sim_data), upper = T))
    # NJ_tree<-NJ(edu_dist)
    NJ_tree<-NJ(mm)
    # random_tree
    random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
    # MP
    pyd=as.phyDat(as.matrix(mm), type='USER', levels=0:max(sim_data))
    MP_tree = optim.parsimony(random_tree, pyd)
    # ML
    ML_tree = optim.pml(pml(random_tree, pyd))$tree
    rf_dist = c(          TreeDist::RobinsonFoulds(rt, it, similarity = F,normalize = T),
                          TreeDist::RobinsonFoulds(rt, cluster_tree, similarity = F,normalize = T),
                          TreeDist::RobinsonFoulds(rt, NJ_tree, similarity = F,normalize = T),
                          TreeDist::RobinsonFoulds(rt, ML_tree, similarity = F,normalize = T),
                          TreeDist::RobinsonFoulds(rt, MP_tree, similarity = F,normalize = T),
                          tryCatch({TreeDist::RobinsonFoulds(rt, medalt_tree, similarity = F,normalize = T)},error = function(cond){return(NA)}),
                          tryCatch({TreeDist::RobinsonFoulds(rt, sitka_tree, similarity = F,normalize = T)},error = function(cond){return(NA)}),
                          tryCatch({TreeDist::RobinsonFoulds(rt, medicc2_tree, similarity = F,normalize = T)},error = function(cond){return(NA)})
    )
    reducedX = collapse.singles(rt)
    reducedX <- keep.tip(reducedX, reducedX$tip.label)

    #it <- keep.tip(it, reducedX$tip.label)

    Qt_res = c(tryCatch({QuartetDivergence(QuartetStatus(rt, it, nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(cluster_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(NJ_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(ML_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(MP_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(medalt_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(sitka_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)}),
               tryCatch({QuartetDivergence(QuartetStatus(trees=reducedX, cf=keep.tip(medicc2_tree, reducedX$tip.label), nTip = TRUE),similarity = TRUE)},error = function(cond){return(NA)})
    )
    Quartet_res = rbind(Quartet_res, c(Qt_res, ncol(sim_data)))
    res = rbind(res, c(rf_dist, ncol(sim_data)))
  }
}

res = as.data.frame(res)
Quartet_res = as.data.frame(Quartet_res)
colnames(Quartet_res) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
write.table(Quartet_res, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Quartet_res.txt')

# res = cbind(res[1:5], NA, res[6:12])
# RF_data2 = RF_data
RF_data = res
colnames(RF_data) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
write.table(RF_data, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/RFdist_res.txt')
RF_data = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/RFdist_res.txt')
# library(paletteer)
# paletteer_d("colorblindr::OkabeIto")
library(showtext)
# a = font.files()
# a[grep('Times', a$file),]

font_add('Times New Roman', '/System/Library/Fonts/Supplemental/Times New Roman.ttf')
showtext_auto()
a = RF_data %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=0.1, fill='white', shape=21) +
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
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_仿真/Fig1_RF距离.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')

Quartet_res = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Quartet_res.txt')
a = Quartet_res %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=0.1, fill='white', shape=21) +
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
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_仿真/Fig1_Quartet距离.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')


#### 仿真数据，拟时序准确性，秩相关 ####
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num//'
medalt_dir = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_for_RFdist/'
sitka_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_sitka/'
sim_100_dir_res = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_key/'
medicc_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_medicc2/'
fl = sapply(list.files(sim_100_dir_res), function(x)substring(x,1,16))
fl = unique(fl)
pseudo_cor = c()
for(f in fl){
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    tmp_data = paste0(sim_100_dir, f, '.txt')
    rt = newick_to_phylo(tmp_ref)
    it = newick_to_phylo(tmp_infer)
    medalt_tree = sitka_newick_to_phylo(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
    sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))
    #medicc2
    if(file.exists(paste0(medicc_dir,f,'/',f,'_final_tree.new'))){
      medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
      medicc2_tree = drop.tip(medicc2_tree, 'diploid')
      medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'), function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
    }else{
      medicc2_tree = NA
    }
    #
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]
    edu_dist = dist(t(sim_data), method = "euclidean")
    hclust_dist = hclust(edu_dist, method = "complete")
    cluster_tree = as.phylo.hclust(hclust_dist)
    # NJ
    pyd = as.phyDat(sim_data, type='USER', levels=c(0,1,2,3,4,5,6))
    #mm = dist(t(sim_data)) # dist.hamming(pyd) # as.matrix(dist(t(sim_data), upper = T))
    NJ_tree<-NJ(edu_dist)
    # random_tree
    random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
    # MP
    MP_tree = optim.parsimony(random_tree, pyd)
    # ML
    ML_tree = optim.pml(pml(random_tree, pyd))$tree
    #
    rt = collapse.singles(rt)
    rt <- keep.tip(rt, rt$tip.label)
    cluster_tree = keep.tip(cluster_tree, rt$tip.label)
    NJ_tree = keep.tip(NJ_tree, rt$tip.label)
    ML_tree = keep.tip(ML_tree, rt$tip.label)
    medalt_tree = keep.tip(medalt_tree, rt$tip.label)
    medalt_tree = keep.tip(medalt_tree, rt$tip.label)
    sitka_tree = keep.tip(sitka_tree, rt$tip.label)

    true_pseudo = sapply(nodepath(rt), length)
    sim_key_res  = load_SCTMS(sim_100_dir_res, f)
    #sim_key_res = build_clone_tree(sim_key_res)
    rownames(sim_key_res$pseudo) = sim_key_res$pseudo$name



    cor_df = data.frame(ref_pseudo=true_pseudo,
                        infer_pseudo=sim_key_res$pseudo[rt$tip.label, 'pseudotime'],
                        cluster_pseudo = sapply(nodepath(cluster_tree), length),
                        NJ_pseudo = sapply(nodepath(NJ_tree), length),
                        ML_pseudo = sapply(nodepath(ML_tree), length),
                        MP_pseudo = sapply(nodepath(MP_tree), length),
                        medalt_pseudo = sapply(nodepath(medalt_tree), length),
                        sitka_pseudo = sapply(nodepath(sitka_tree), length),
                        medicc2_pseudo = tryCatch({sapply(nodepath(medicc2_tree), length)}, error = function(cond){return(NA)})
    )

    cor_df2 = cor_df %>%
      group_by(ref_pseudo) %>%
      summarise_all(list(median))

    #cor_df2 = cor_df %>%
    #  group_by(infer_pseudo) %>%
    #  summarise_all(list(mean))
    #
    res = c(
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$ref_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$cluster_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$NJ_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$ML_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$MP_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$medalt_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$sitka_pseudo)$auc}, error = function(cond){return(NA)}),
      tryCatch({multiclass.roc(cor_df2$infer_pseudo, cor_df2$medicc2_pseudo)$auc}, error = function(cond){return(NA)})
    )
    # res = cor(cor_df2$ref_pseudo, cor_df2$infer_pseudo, method = 'spearman')#cor(true_pseudo, sim_key_res$pseudo[rt$tip.label, 'pseudotime'], method = 'pearson')
    # res = cor(cor_df2$ref_pseudo, cor_df2$infer_pseudo, method = 'pearson')#cor(true_pseudo, sim_key_res$pseudo[rt$tip.label, 'pseudotime'], method = 'pearson')
    pseudo_cor = rbind(pseudo_cor, c(res, dim(sim_data)[2]))
  }
}
pseudo_cor = as.data.frame(pseudo_cor)
colnames(pseudo_cor) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
#rownames(pseudo_cor) = fl
#pseudo_cor[rownames(use_name),]
# colnames(pseudo_cor) = c('correlation', 'cell_num')
write.table(pseudo_cor, '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/pseudo_cor_res.txt', quote = F)
pseudo_cor = read.table('/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/pseudo_cor_res.txt')
a = pseudo_cor %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=1, fill='white', shape=21) +
  labs(title='Different methods', y='MultiAUC: Pseudotime')+
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
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_仿真/pseudotime_correlation.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')


# 评估仿真数据整体的不等倍分离比率与推断tree的不等倍分离比率
# 仿真数据，评估模型预测不等倍分离的准确性
find_triple <- function(tree){
  lvs = tree$tip.label
  all_path = nodepath(tree)
  tmp_res = c()
  for(i in all_path){
    tmp_res = rbind(tmp_res, i[(length(i)-1):length(i)])
  }
  tmp_res = as.data.frame(tmp_res)
  colnames(tmp_res) = c('parent', 'child')
  triple_res = c()
  for(i in unique(tmp_res$parent)){
    x = tmp_res[tmp_res$parent==i, 'child']
    if(length(x)==2){
      triple_res = rbind(triple_res, c(i,x))
    }
  }
  triple_res = as.data.frame(triple_res)
  colnames(triple_res) = c('parent', 'child1', 'child2')
  return(triple_res)
}
aneu_rate <- function(triple_res, tips, data){
  ad_locs = c()
  for(i in 1:nrow(triple_res)){
    tmp_child1 = as.vector(triple_res[i,'child1'])
    tmp_child2 = as.vector(triple_res[i,'child2'])
    tmp_child1 = tips[tmp_child1]
    tmp_child2 = tips[tmp_child2]
    tmp_data = t(data[, c(tmp_child1, tmp_child2)])
    locs = (((tmp_data[1,] + tmp_data[2,]) %%2) == 0) & (tmp_data[1,] != tmp_data[2,])
    ad_locs = c(ad_locs, sum(locs))
  }
  return(ad_locs)
}

aneu_rate_res = c()
for(f in fl){
  thr = 0
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    tmp_data = paste0(sim_100_dir, f, '.txt')
    rt = newick_to_phylo(tmp_ref)
    it = newick_to_phylo(tmp_infer)
    medalt_tree = sitka_newick_to_phylo(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
    sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))
    #medicc2
    if(file.exists(paste0(medicc_dir,f,'/',f,'_final_tree.new'))){
      medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
      medicc2_tree = drop.tip(medicc2_tree, 'diploid')
      medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'), function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
    }else{
      medicc2_tree = NA
    }
    #
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]
    edu_dist = dist(t(sim_data), method = "euclidean")
    hclust_dist = hclust(edu_dist, method = "complete")
    cluster_tree = as.phylo.hclust(hclust_dist)
    # NJ
    pyd = as.phyDat(sim_data, type='USER', levels=c(0,1,2,3,4,5,6))
    #mm = dist(t(sim_data)) # dist.hamming(pyd) # as.matrix(dist(t(sim_data), upper = T))
    NJ_tree<-NJ(edu_dist)
    # random_tree
    random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
    # MP
    MP_tree = optim.parsimony(random_tree, pyd)
    # ML
    ML_tree = optim.pml(pml(random_tree, pyd))$tree

    rt_rate = aneu_rate(find_triple(rt), rt$tip.label, sim_data)
    it_rate = aneu_rate(find_triple(it), it$tip.label, sim_data)
    cluster_tree_rate = aneu_rate(find_triple(cluster_tree), cluster_tree$tip.label, sim_data)
    NJ_tree_rate = aneu_rate(find_triple(NJ_tree), NJ_tree$tip.label, sim_data)

    tmp_func <- function(tt){
      rr = aneu_rate(find_triple(tt), tt$tip.label, sim_data)
      sum(rr>thr) / length(rr)
    }

    res = c(
      tryCatch({tmp_func(rt)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(it)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(NJ_tree)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(ML_tree)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(MP_tree)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(medalt_tree)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(sitka_tree)}, error = function(cond){return(NA)}),
      tryCatch({tmp_func(medicc2_tree)}, error = function(cond){return(NA)})
    )
    aneu_rate_res = rbind(aneu_rate_res, c(res, ncol(sim_data)))
  }
}
aneu_rate_res = as.data.frame(aneu_rate_res)
colnames(aneu_rate_res) = c('ref_aneu_rate', 'infer_aneu_rate')

colnames(aneu_rate_res) = c('ref', 'scCNAT', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
write.table(aneu_rate_res, '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/aneu_rate_res.txt', quote = F)
aneu_rate_res = read.table('/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/aneu_rate_res.txt')

reshape2::melt(aneu_rate_res,id=c('ref', 'cell_num')) %>%
  ggplot(aes(x=ref, y=value, color=as.factor(cell_num)))+
  geom_point()+
  stat_cor()+
  facet_wrap(~variable, ncol=3)+
  coord_fixed()+
  theme_bw()+
  scale_color_igv()+
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

cnum = unique(aneu_rate_res$cell_num)
cor_data = data.frame(matrix(0, nrow=7, ncol=length(cnum)))
rownames(cor_data) = c('scCNAT', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2')
colnames(cor_data) = cnum
cor_pvalue = cor_data
for(i in rownames(cor_data)){
  for(j in colnames(cor_data)){
    tmp = aneu_rate_res[aneu_rate_res$cell_num==j, i]
    pos = is.na(tmp)
    tmp = na.omit(tmp)
    tmp2 = aneu_rate_res[aneu_rate_res$cell_num==j, 'ref'][!pos]
    if(length(tmp)!=0){
      tryCatch({
        cor_p = cor.test(tmp2, tmp)
        cor_data[i, j] = cor_p$estimate
        cor_pvalue[i, j] = cor_p$p.value
      }, error = function(cond){return(NA)})
    }
  }
}

cor_all = reshape2::melt(as.matrix(cor_data))
colnames(cor_all) = c('Var1', 'Var2', 'cor')
cor_all = merge(cor_all, reshape2::melt(as.matrix(cor_pvalue)))
cor_all %>%
  #filter(!Var1%in%c('MP', 'NJ')) %>%
  ggplot()+
  geom_line(aes(x=Var2, y=cor, color=Var1))+
  labs(title='Different methods', y='correlation aneu rate')+
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

pheatmap::pheatmap(cor_data,
                   cluster_rows = F, cluster_cols = F,
                   display_numbers = round(cor_pvalue, 2)
)

aneu_rate_res %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('hclust')) %>%
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


ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/aneu_rate_correlation.pdf',
       a, width=60, height=60, units='mm', dpi = 600, bg = 'transparent')




# 预测分裂速率的准确性
pred_cor = c()
all_auc_res =c()
for(f in fl){
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    it = newick_to_phylo(tmp_infer)
    tmp_data = paste0(sim_100_dir, f, '.txt')
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]
    sim_key_res  = load_SCTMS(sim_100_dir_res, f)
    #sim_key_res = build_clone_tree(sim_key_res)
    rownames(sim_key_res$pseudo) = sim_key_res$pseudo$name

    it_triple = find_triple(it)
    tmp_node_label = fortify(it) %>% as.data.frame()
    rownames(tmp_node_label) = tmp_node_label$node

    set.seed(1994)
    label_triple = sample(1:nrow(it_triple), as.integer(nrow(it_triple)/5))
    label_triple = it_triple[label_triple, ]
    drop_tips = c(label_triple[,2], label_triple[,3])
    new_it = drop.tip(it, drop_tips, trim.internal=FALSE)
    groth_leaves = new_it$node.label
    groth_nodes = tmp_node_label[label_triple[,1], 'label']

    #old_leaves = it$tip.label
    new_it = drop.tip(new_it, new_it$tip.label, trim.internal=FALSE)
    new_leaves = new_it$tip.label
    train_nodes = setdiff(new_it$node.label, 'virtual_1')
    # 训练模型
    fm = label~root_loc+root_cn+parent_loc+parent_cn+pseudotime#+self_time
    fm = label~root_loc+parent_loc+pseudotime#+self_time
    model_dataset = sim_key_res$pseudo[-1, ]
    #na_pos = rowSums(model_dataset=='none')==0
    train_data = data.frame('root_loc'=as.numeric(model_dataset$Root_gain_loc)+as.numeric(model_dataset$Root_loss_loc),
                            #'root_cn' = as.numeric(model_dataset$Root_gain_cn)+as.numeric(model_dataset$Root_loss_cn),
                            'parent_loc' = as.numeric(model_dataset$Parent_gain_loc)+as.numeric(model_dataset$Parent_loss_loc),
                            #'parent_cn' = as.numeric(model_dataset$Parent_gain_cn)+as.numeric(model_dataset$Parent_loss_cn),
                            'pseudotime' = as.numeric(model_dataset$pseudotime)#,
                            #'self_time' = as.numeric(model_dataset$split_time2)
    )
    rownames(train_data) = rownames(model_dataset)
    train_d = train_data[train_nodes, ]
    all_label = model_dataset[c(train_nodes, new_leaves), c('split_dd_loc', 'split_ad_loc', 'split_time')]
    all_label = apply(all_label, 2, as.numeric)
    all_label = (all_label[,1]+all_label[,2]) / (all_label[,3]+1)
    # all_label = sapply(all_label, function(x)ifelse(x>median(all_label), 1, 0))
    names(all_label) = c(train_nodes, new_leaves)

    train_d$label = all_label[rownames(train_d)]
    test_data = train_data[new_leaves, ]
    true_test_label = all_label[new_leaves]
    true_test_label = sapply(names(true_test_label), function(x)ifelse(x%in%groth_nodes,1,0))
    glm_model=glm(data = train_d, formula = fm)
    glm.pred <- predict(glm_model, test_data)
    # glm_model=loess(data = train_d, formula = fm)
    # glm.model <- glm(formula = fm, family = binomial(link = "logit"), data = train_d)
    # glm.pred <- predict(glm.model, test_data, type="response")
    # logit.pred <- factor(glm.pred > .5, levels=c(FALSE, TRUE), labels=c("Growth", "Non-growth"))

    #table(logit.pred, true_test_label)
    #
    #pROC::roc(true_test_label, as.numeric(glm.pred))$auc
    #
    #pROC::plot.roc(true_test_label, as.numeric(glm.pred),
    #               print.auc=TRUE,
    #               auc.polygon=TRUE,
    #               grid=c(0.1, 0.2),
    #               grid.col=c("green", "red"),
    #               max.auc.polygon=TRUE,
    #               auc.polygon.col="skyblue",
    #               #print.thres=TRUE,
    #               main='logistic regression train set ROC')

    glm_pred = predict(glm_model, test_data)
    # pred_cor = c(pred_cor, cor(true_test_label, glm_pred))

    pred_cor = rbind(pred_cor, data.frame('true_rate'=true_test_label, 'pred_rate'=glm_pred))


    min_cn = model_dataset[c(train_nodes, new_leaves), c('split_dd_loc', 'split_ad_loc')]
    min_cn = as.numeric(min_cn$split_dd_loc) + as.numeric(min_cn$split_ad_loc)
    min_cn = setdiff(min_cn, 0)
    min_cn = min(min_cn)
    ts = 1
    auc_res = c()
    napos = is.na(glm_pred)
    while(ts<20){
      pred_cn = glm_pred[!napos] * ts
      predict_label = (pred_cn>min_cn)+0
      if(length(unique(predict_label))==2){
        tmp_roc = pROC::roc(predict_label, true_test_label[!napos])$auc
      }else{
        tmp_roc = NA
      }
      auc_res = rbind(auc_res, c(ts, tmp_roc))
      ts = ts + 1
    }
    all_auc_res = rbind(all_auc_res, auc_res)
  }
}
pred_cor = as.data.frame(pred_cor)
pred_cor_stat = pred_cor
pred_cor_stat = pred_cor_stat[pred_cor_stat$true_rate <= 8,]
pred_cor_stat = pred_cor_stat[pred_cor_stat$true_rate >= 1,]
pred_cor_stat = pred_cor_stat[pred_cor_stat$pred_rate <= 20,]
pred_cor_stat = pred_cor_stat[pred_cor_stat$pred_rate >= 0,]
a = pred_cor %>%
  na.omit() %>%
  ggplot(aes(x=true_rate, y=pred_rate)) +
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(data=pred_cor_stat, mapping = aes(x=true_rate, y=pred_rate),
              formula = 'y~x', method='lm', se=F, color='red', inherit.aes = T)+
  stat_cor(data=pred_cor_stat, mapping = aes(x=true_rate, y=pred_rate),
           color='red', size=4, label.y = 1, label.x=12)+
  #lims(x=c(0,20), y=c(0,20))+
  labs(x='True mitosis rate', y='Predict mitosis rate', title='Mitosis rate correlation')+
  theme_bw()+
  theme(legend.position = 'none',
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=10))#+
#coord_equal()
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Mitosis_rate_correlation.pdf',
       a, width=120, height=120, units='mm', dpi = 600, bg = 'transparent')

#
all_auc_res = as.data.frame(all_auc_res)
colnames(all_auc_res) = c('time', 'roc')

a = all_auc_res %>%
  ggplot(aes(x=as.factor(time), y=roc))+
  geom_boxplot(fill='lightblue', alpha=0.6, width=0.5)+
  labs(x='Pseudotime', y='AUC', title='Mitosis rate AUC')+
  theme_bw()+
  theme(legend.position = 'none',
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=10))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Mitosis_rate_auc.pdf',
       a, width=120, height=60, units='mm', dpi = 600, bg = 'transparent')





# 预测分裂速率的准确性
predict_mitosis_rate <- function(obj){
  pseudo_res = obj$pseudo
  model_dataset = pseudo_res[-1, ]
  na_pos = rowSums(model_dataset=='none')==0
  train_data = data.frame('root_loc'=as.numeric(model_dataset$Root_gain_loc)+as.numeric(model_dataset$Root_loss_loc),
                          'root_cn' = as.numeric(model_dataset$Root_gain_cn)+as.numeric(model_dataset$Root_loss_cn),
                          #'parent_loc' = as.numeric(model_dataset$Parent_gain_loc)+as.numeric(model_dataset$Parent_loss_loc),
                          #'parent_cn' = as.numeric(model_dataset$Parent_gain_cn)+as.numeric(model_dataset$Parent_loss_cn),
                          'pseudotime' = as.numeric(model_dataset$pseudotime)
                          #'self_time' = as.numeric(model_dataset$split_time2)
  )
  train_d = train_data[na_pos, ]
  train_label = model_dataset[na_pos, c('split_dd_loc', 'split_ad_loc', 'split_time')]
  train_label = apply(train_label, 2, as.numeric)
  train_label[train_label[,3]==0,3] = 1
  train_label = (train_label[,1]+train_label[,2]) / train_label[,3]
  train_d$label = train_label
  train_d[is.na(train_d)] = 0
  rownames(train_d) = rownames(model_dataset[na_pos, 1:8])
  outliers = boxplot.stats(train_d$label)$out
  train_d = train_d[setdiff(rownames(train_d), names(outliers)),]
  # test
  test_data = train_data[!na_pos, ]
  rownames(test_data) = rownames(model_dataset[!na_pos, ])

  #fm = label~root_loc+root_cn+parent_loc+parent_cn+pseudotime+self_time
  fm = label~root_loc+root_cn+pseudotime
  glm_model=glm(data = train_d, formula = fm, family = gaussian())
  glm_pred = predict(glm_model, test_data)
  return(list(model=glm_model, pred_res=glm_pred))
}

calculate_mitosis_rate <- function(obj){
  pseudo_res = obj$pseudo
  model_dataset = pseudo_res[-1, ]
  na_pos = rowSums(model_dataset=='none')==0
  train_d = train_data[na_pos, ]
  train_label = model_dataset[na_pos, c('split_dd_loc', 'split_ad_loc', 'split_time')]
  train_label = apply(train_label, 2, as.numeric)
  train_label = (train_label[,1]+train_label[,2]) / train_label[,3]
  train_label[is.na(train_label)] = 0
  names(train_label) = rownames(model_dataset[na_pos, ])

  outliers = boxplot.stats(train_label)$out
  train_label = train_label[setdiff(names(train_label), names(outliers))]
  return(train_label)
}

pred_cor = c()
all_auc_res =c()
sample_tree_dir="/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_res_sample/"
for(f in fl){
  if(file.exists(paste0(sample_tree_dir, f, '_infer_tree.txt'))){

    sim_key_res  = load_SCTMS(sim_100_dir_res, f)
    sim_key_res = build_clone_tree(sim_key_res)
    rownames(sim_key_res$pseudo) = sim_key_res$pseudo$name
    #sample_tree = paste0(sample_tree_dir, f, '_infer_tree.txt')
    sim_key_res_sample  = load_SCTMS(sample_tree_dir, f)
    sim_key_res_sample = build_clone_tree(sim_key_res_sample)
    rownames(sim_key_res_sample$pseudo) = sim_key_res_sample$pseudo$name

    true_rate = calculate_mitosis_rate(sim_key_res)
    sample_rate = predict_mitosis_rate(sim_key_res_sample)
    names(sample_rate$pred_res) = gsub('new_', '',names(sample_rate$pred_res))
    share_name = intersect(names(true_rate), names(sample_rate$pred_res))
    if(length(share_name)!=0){
      true_rate = true_rate[share_name]
      sample_rate = sample_rate$pred_res[share_name]
      #cor(true_rate, sample_rate)
      #plot(true_rate, sample_rate)
      rt = newick_to_phylo(paste0(sim_100_dir_res,f,'_infer_tree.txt'))
      predict_label = c()
      for(i in share_name){
        tmp_lv = extract.clade(rt, i)
        max_l = max(sapply(nodepath(tmp_lv), length))
        if(max_l>2){
          predict_label = c(predict_label, 1)
        }else{
          predict_label = c(predict_label, 0)
        }
      }
      if(length(unique(predict_label))==2){
        tmp_roc = pROC::roc(predict_label, sample_rate)$auc
      }else{
        tmp_roc = NA
      }
      # pred_cor = c(pred_cor, cor(true_test_label, glm_pred))
      pred_cor = rbind(pred_cor, data.frame('true_rate'=true_rate, 'pred_rate'=sample_rate))
      #
      #min_cn = sim_key_res$pseudo[names(true_rate), c('split_dd_loc', 'split_ad_loc')]
      #min_cn = as.numeric(min_cn$split_dd_loc) + as.numeric(min_cn$split_ad_loc)
      #min_cn = setdiff(min_cn, 0)
      #min_cn = min(min_cn)
      #ts = 1
      #auc_res = c()
      #napos = is.na(sample_rate)
      #while(ts<20){
      #  pred_cn = sample_rate[!napos] * ts
      #  predict_label = (pred_cn>min_cn)+0
      #  if(length(unique(predict_label))==2){
      #    tmp_roc = pROC::roc(predict_label, true_rate[!napos])$auc
      #  }else{
      #    tmp_roc = NA
      #  }
      #  auc_res = rbind(auc_res, c(ts, tmp_roc))
      #  ts = ts + 1
      #}
      all_auc_res = rbind(all_auc_res, tmp_roc)
    }

  }
}
pred_cor = as.data.frame(pred_cor)
pred_cor_stat = pred_cor
pred_cor_stat = pred_cor_stat[pred_cor_stat$true_rate > 0,]
pred_cor_stat = pred_cor_stat[pred_cor_stat$pred_rate > 0,]
a = pred_cor %>%
  na.omit() %>%
  ggplot(aes(x=true_rate, y=pred_rate)) +
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(data=pred_cor_stat, mapping = aes(x=true_rate, y=pred_rate),
              formula = 'y~x', method='lm', se=F, color='red', inherit.aes = T)+
  stat_cor(data=pred_cor_stat, mapping = aes(x=true_rate, y=pred_rate),
           color='red', size=4#, label.y = 1, label.x=12
  )+
  lims(x=c(0,20), y=c(0,20))+
  labs(x='True mitosis rate', y='Predict mitosis rate', title='Mitosis rate correlation')+
  theme_bw()+
  theme(legend.position = 'none',
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=10))+
  coord_equal()
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Mitosis_rate_correlation.pdf',
       a, width=120, height=120, units='mm', dpi = 600, bg = 'transparent')

#
boxplot(all_auc_res)
all_auc_res = as.data.frame(all_auc_res)
colnames(all_auc_res) = c('time', 'roc')

a = all_auc_res %>%
  ggplot(aes(x=as.factor(time), y=roc))+
  geom_boxplot(fill='lightblue', alpha=0.6, width=0.5)+
  labs(x='Pseudotime', y='AUC', title='Mitosis rate AUC')+
  theme_bw()+
  theme(legend.position = 'none',
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=10))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Mitosis_rate_auc.pdf',
       a, width=120, height=60, units='mm', dpi = 600, bg = 'transparent')







####### 3.比较ARIindex #######
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num//'
scTrace_dir='/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sim_for_RFdist/'
medalt_dir = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_for_RFdist/'
sitka_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_sitka/'
sim_100_dir_res = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_key/'
medicc_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_medicc2/'
fl = sapply(list.files(sim_100_dir_res), function(x)substring(x,1,16))
fl = unique(fl)
res = c()
ari_res = c()
auc_res = c()
for(f in fl){
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    tmp_data = paste0(sim_100_dir, f, '.txt')
    rt = newick_to_phylo(tmp_ref)
    scTrace_res = load_scTrace(scTrace_dir, f)
    it = scTrace_res$tree
    medalt_tree = sitka_newick_to_phylo(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    sitka_tree = sitka_newick_to_phylo(paste0(sitka_dir, f, '_sitka.newick'))
    sitka_tree$tip.label = sapply(sitka_tree$tip.label, function(x)substr(x, 6, nchar(x)))
    #medicc2
    if(file.exists(paste0(medicc_dir,f,'/',f,'_final_tree.new'))){
      medicc2_tree = sitka_newick_to_phylo(paste0(medicc_dir,f,'/',f,'_final_tree.new'))
      medicc2_tree = drop.tip(medicc2_tree, 'diploid')
      medicc2_tree$tip.label = sapply(strsplit(medicc2_tree$tip.label, '_'), function(x)paste0(x[1], '_', x[2], '_', x[3], '_'))
    }else{
      medicc2_tree = NA
    }
    #
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]
    edu_dist = dist(t(sim_data), method = "euclidean")
    hclust_dist = hclust(edu_dist, method = "complete")
    cluster_tree = as.phylo.hclust(hclust_dist)

    # NJ_tree<-NJ(edu_dist)
    mm = dist(t(sim_data), method = "maximum")
    NJ_tree<-NJ(mm)
    # random_tree
    random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
    # MP
    pyd=as.phyDat(as.matrix(mm), type='USER', levels=0:max(sim_data))
    MP_tree = optim.parsimony(random_tree, pyd)


    # MEDALT
    # medalt_tree <- ML_tree# read.tree(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    ###
    cn = length(rt$tip.label)/4#length(unique(sim_key_res$node_info[rownames(rt_clone), 'clone'])) + 10
    #cn = length(unique(sim_key_res$node_info[rownames(rt_clone), 'clone']))
    rt_clone = get_clone(rt, min_k = cn)
    it_clone = get_clone(it, min_k = cn)

    #sim_key_res  = load_SCTMS(sim_100_dir_res, f)
    #sim_key_res = build_clone_tree(sim_key_res)
    #rt_clone$it_clone =  sim_key_res$node_info[rownames(rt_clone), 'clone']
    rt_clone$it_clone = it_clone[rt_clone$name, 'clone']


    tip_dist <- function(tree, label=T){
      x = fortify(tree)[, c('node', 'label')]
      label_name = x$node
      names(label_name) = x$label
      r_label = data.frame(matrix(0, nrow=length(tree$tip.label),
                                  ncol = length(tree$tip.label)),
                           row.names = tree$tip.label)
      colnames(r_label) = tree$tip.label
      for(i in tree$tip.label){
        for(j in tree$tip.label){
          d = nodepath(tree, from=label_name[i], to=label_name[j])
          if(label){
            if(length(d)==3){
              r_label[i,j] = 1
            }
          }else{
            r_label[i,j] = length(d)
          }
        }
      }
      return(r_label)
    }
    ref_label = tip_dist(rt, label = T)
    # lapply(st, function(x))
    # rt$tip.label

    #pROC::plot.roc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(NJ_score)))

    #x = pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(NJ_score)))

    tmp_func = function(x){
      diversity(table(x))
    }
    rt_clone_res = rt_clone %>%
      group_by(clone) %>%
      summarise(d = tmp_func(it_clone))
    auc_dist = c(
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(rt, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(cluster_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(NJ_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(ML_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(MP_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(medalt_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(sitka_tree, label = F))))},error = function(cond){return(NA)}),
      tryCatch({pROC::auc(as.vector(as.matrix(ref_label)), as.vector(as.matrix(tip_dist(medicc2_tree, label = F))))},error = function(cond){return(NA)})
    )
    ###
    #cn = length(unique(sim_key_res$node_info[rownames(rt_clone), 'clone']))# + 10
    rf_dist = c(mean(rt_clone_res$d),
                clone_accuracy(rt, cluster_tree, clone_num = cn),
                clone_accuracy(rt, NJ_tree, clone_num = cn),
                clone_accuracy(rt, ML_tree, clone_num = cn),
                tryCatch({clone_accuracy(rt, MP_tree, clone_num = cn)},error = function(cond){return(NA)}),
                clone_accuracy(rt, medalt_tree, clone_num = cn),
                clone_accuracy(rt, sitka_tree, clone_num = cn),
                tryCatch({clone_accuracy(rt, medicc2_tree, clone_num = cn)},error = function(cond){return(NA)})
    )
    #cn = NULL
    ARI_dist = c(adj.rand.index(rt_clone$clone, rt_clone$it_clone),
                 clone_ARI(rt, cluster_tree, clone_num = cn),
                 clone_ARI(rt, NJ_tree, clone_num = cn),
                 clone_ARI(rt, ML_tree, clone_num = cn),
                 clone_ARI(rt, MP_tree, clone_num = cn),
                 clone_ARI(rt, medalt_tree, clone_num = cn),
                 clone_ARI(rt, sitka_tree, clone_num = cn),
                 tryCatch({clone_ARI(rt, medicc2_tree, clone_num = cn)},error = function(cond){return(NA)})
    )
    res = rbind(res, c(rf_dist, ncol(sim_data)))
    ari_res = rbind(ari_res, c(ARI_dist, ncol(sim_data)))
    auc_res = rbind(auc_res, c(auc_dist, ncol(sim_data)))
  }
}
# res = as.data.frame(res)
# res = cbind(res[1:5], NA, res[6:12])
RF_data = ari_res
colnames(RF_data) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
write.table(RF_data, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ARIindex_res.txt')
RF_data = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/ARIindex_res.txt')

a = RF_data %>%
  reshape2::melt(id='cell_num') %>%
  #filter(!variable%in%c('MP','NJ','hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=1, fill='white', shape=21) +
  labs(title='Different methods', y='ARI index')+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/clone_ARI.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')



####### 4.比较AUC #######
RF_data = auc_res
colnames(RF_data) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
write.table(RF_data, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/AUCindex_res.txt')
RF_data = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/AUCindex_res.txt')

a=RF_data %>%
  reshape2::melt(id='cell_num') %>%
  filter(!variable%in%c('hclust')) %>%
  mutate(variable = factor(variable, levels=c( 'scCNAT','MEDICC2','MP','sitka', 'NJ', 'MEDALT','ML'))) %>%

  group_by(cell_num, value) %>%
  ggplot(aes(x = cell_num,y = value, group=variable,color=variable)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width=0.1)+
  stat_summary(fun.y="mean",geom="line") +
  stat_summary(fun.y="mean",geom="point",size=1, fill='white', shape=21) +
  labs(title='Different methods', y='AUC (ref, real)')+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/clone_AUC.pdf',
       a, width=80, height=50, units='mm', dpi = 600, bg = 'transparent')

####### 5.比较预测速率的相关性 #######
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
find_triple <- function(tree){
  lvs = tree$tip.label
  all_path = nodepath(tree)
  tmp_res = c()
  for(i in all_path){
    tmp_res = rbind(tmp_res, i[(length(i)-1):length(i)])
  }
  tmp_res = as.data.frame(tmp_res)
  colnames(tmp_res) = c('parent', 'child')
  triple_res = c()
  for(i in unique(tmp_res$parent)){
    x = tmp_res[tmp_res$parent==i, 'child']
    if(length(x)==2){
      triple_res = rbind(triple_res, c(i,x))
    }
  }
  triple_res = as.data.frame(triple_res)
  colnames(triple_res) = c('parent', 'child1', 'child2')
  return(triple_res)
}
aneu_rate <- function(triple_res, tips, data){
  ad_locs = c()
  for(i in 1:nrow(triple_res)){
    tmp_child1 = as.vector(triple_res[i,'child1'])
    tmp_child2 = as.vector(triple_res[i,'child2'])
    tmp_child1 = tips[tmp_child1]
    tmp_child2 = tips[tmp_child2]
    tmp_data = t(data[, c(tmp_child1, tmp_child2)])
    locs = (((tmp_data[1,] + tmp_data[2,]) %%2) == 0) & (tmp_data[1,] != tmp_data[2,])
    ad_locs = c(ad_locs, sum(locs))
  }
  return(ad_locs)
}
pred_cor = c()
all_auc_res =c()
for(f in fl){
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    it = newick_to_phylo(tmp_infer)
    rt = newick_to_phylo(tmp_ref)
    scTrace_res = load_scTrace(scTrace_dir, f)

    tmp_data = paste0(sim_100_dir, f, '.txt')
    sim_data = read.table(tmp_data, header = T)
    sim_data = sim_data[,3:ncol(sim_data)]

    meta = scTrace_res$cell_relation
    # tmp_rate = predict_mitosis_rate(as.data.frame(meta))

    it_triple = find_triple(it)
    tmp_node_label = fortify(it) %>% as.data.frame()
    rownames(tmp_node_label) = tmp_node_label$node

    tmp_res_auc = c()
    for(iter in 1:100){
      for(k in 2:10){
        if(k>=nrow(it_triple)){
          break
        }
        # set.seed(1994)
        label_triple = sample(1:nrow(it_triple), k)
        label_triple = it_triple[label_triple, ]
        drop_tips = c(label_triple[,2], label_triple[,3])
        new_it = drop.tip(it, drop_tips, trim.internal=FALSE)
        new_it_nodes = c(new_it$node.label, new_it$tip.label)
        groth_nodes = setdiff(new_it$tip.label, it$tip.label)
        x = as.data.frame(meta[new_it_nodes,])
        x[new_it$tip.label, 'Mitosis_dd_loc']=NA
        new_rate = predict_mitosis_rate(x)

        true_test_label = sapply(names(new_rate), function(x)ifelse(x%in%groth_nodes,1,0))
        if(length(unique(true_test_label))==1){
          break
        }
        auc_res = pROC::roc(true_test_label, as.numeric(new_rate))$auc

        tmp_res_auc = rbind(tmp_res_auc, c(k,as.numeric(auc_res), length(it$tip.label), iter))

        #true_rate = meta[groth_nodes, c('Mitosis_dd_loc', 'Mitosis_ad_loc', 'Mitosis_time'), drop=F]
        #true_rate = log1p((true_rate[,1]+true_rate[,2])) / true_rate[,3]
        #pred_cor = rbind(pred_cor, c(k, cor(new_rate[groth_nodes], true_rate), length(it$tip.label)))
      }
    }
    tmp_res_auc = as.data.frame(tmp_res_auc)
    tmp_res_auc = tmp_res_auc %>%
      group_by(V4, V3) %>%
      summarise(V2=max(V2)) %>% as.data.frame()
    all_auc_res = rbind(all_auc_res, tmp_res_auc)#c(k,as.numeric(auc_res), length(it$tip.label)))
  }
}
all_auc_res = as.data.frame(all_auc_res)
colnames(all_auc_res) = c('iter', 'cell_num','auc')

write.table(all_auc_res, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/AUCvelocity_res.txt')
all_auc_res = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/AUCvelocity_res.txt')


a = all_auc_res %>%
  #filter(k>26) %>%
  #group_by(cell_num,k) %>%
  #summarise(auc=min(auc)) %>%
  ggplot(aes(x=as.factor(cell_num), y=auc)) +
  geom_boxplot(size=0.1, outlier.size = 0.01, width=0.5)+
  geom_hline(yintercept = 0.5, color='gray', linetype='dashed')+
  labs(x='Cell number', y='AUC of velocity')+
  theme_bw()+
  theme(axis.text = element_text(size=6),
        panel.grid = element_blank(),
        axis.title = element_text(size=6),
        legend.key.size = unit(3,'mm'),
        legend.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/clone_AUC_rate.pdf',
       a, width=50, height=40, units='mm', dpi = 600, bg = 'transparent')


pred_cor = as.data.frame(pred_cor)
colnames(pred_cor) = c('k', 'roc', 'cell_num')
pred_cor %>%
  ggplot(aes(x=as.factor(cell_num), y=roc,fill=as.factor(k))) +
  geom_boxplot()+
  geom_hline(yintercept = 0.5, color='gray', linetype='dashed')+
  labs(x='Cell number', y='AUC', title='Mitosis velocity')+
  theme_bw()+
  theme(legend.position = 'none',
        #axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5, size=10))
####### 6.pseudotime与叶子节点深度相关性 #######


####### 7.统计运行时间和内存 #######


