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

########



# https://ms609.github.io/TreeDist/
#sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_time///'
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num//'

# sim_100_dir = '/Users/lab/wangxin/MEDALT/scCNAT/simulation_res/sim_tree_100/'
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
  tmp_infer = paste0(sim_100_dir_res, f, '_infer_tree.txt')
  if(file.exists(tmp_infer)){
    tmp_ref = paste0(sim_100_dir,f,'.newick')
    tmp_data = paste0(sim_100_dir, f, '.txt')
    rt = newick_to_phylo(tmp_ref)
    it = newick_to_phylo(tmp_infer)

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
    pyd = as.phyDat(sim_data, type='USER', levels=c(0,1,2,3,4,5,6))
    #mm = dist(t(sim_data)) # dist.hamming(pyd) # as.matrix(dist(t(sim_data), upper = T))
    NJ_tree<-NJ(edu_dist)
    # random_tree
    random_tree = rtree(ncol(sim_data), tip.label = sample(colnames(sim_data), ncol(sim_data)))
    # MP
    MP_tree = optim.parsimony(random_tree, pyd)
    # ML
    ML_tree = optim.pml(pml(random_tree, pyd))$tree

    # break
    # #
    # t1 <-ggtree(rt)
    # t2 <- ggtree(it)
    # d1 <- t1$data
    # d2 <- t2$data
    # d1$tree <-'t1'
    # d2$tree <-'t2'
    # d2$x <- max(d2$x) - d2$x + max(d1$x) +  max(d1$x)*0.3
    # pp <- t1 + geom_tree(data=d2)
    # dd <- bind_rows(d1, d2) %>%
    #   filter(isTip == TRUE)
    # dd1 <- as.data.frame(dd)
    # pp + geom_line(aes(x, y, group=label), data=dd1)

    # MEDALT
    # medalt_tree <- ML_tree# read.tree(paste0(sim_100_dir, f, '_medalt_tree.txt'))
    ###
    # RF.dist(rt, it)
    # TreeDist::RobinsonFoulds(rt, it)
    # TreeDist::RobinsonFoulds(rt, it, similarity = F,normalize = T)
    # TreeDist::RobinsonFoulds(rt, cluster_tree, similarity = F,normalize = T)
    # TreeDist::RobinsonFoulds(rt, NJ_tree, similarity = F,normalize = T)
    # TreeDist::RobinsonFoulds(rt, ML_tree, similarity = F,normalize = T)
    # TreeDist::RobinsonFoulds(rt, MP_tree, similarity = F,normalize = T)

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
# res = cbind(res[1:5], NA, res[6:12])
# RF_data2 = RF_data
RF_data = res
colnames(RF_data) = c('scCNAT', 'hclust', 'NJ', 'ML', 'MP', 'MEDALT', 'sitka', 'MEDICC2','cell_num')
#rownames(RF_data) = fl
#RF_data  = RF_data[RF_data$scCNAT<150, ]

write.table(RF_data, '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/RFdist_res.txt', quote = F)


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

