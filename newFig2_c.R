
####### Fig1c ############

file_name = 'Fig1c_sub_dataset_P19BT_20181009.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')

other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c1,c2,c3,c4))
other_c = grep('virtual|root', other_c, invert = T, value = T)
l1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup',
       'SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup',
       'SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup',
       'SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')

#all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
#tmp_tree$edge.length=NULL
#tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')
tmp_tree = all_trees[[file_name]]$medalt
#tmp_tree$edge.length=NULL
a2_data= fortify(tmp_tree)
tmp_tree = all_trees[[file_name]]$sitka
#tmp_tree$edge.length=NULL
a3_data= fortify(tmp_tree)
tmp_tree = all_trees[[file_name]]$medicc2
#tmp_tree$edge.length=NULL
a4_data= fortify(tmp_tree)
a4_data$label = gsub('_NA_', '', a4_data$label)
tmp_tree = all_trees[[file_name]]$NJ
tmp_tree$edge.length = tmp_tree$edge.length+0.1
a5_data= fortify(tmp_tree)
a6_data= fortify(all_trees[[file_name]]$MP)
tmp_tree = all_trees[[file_name]]$ML
tmp_tree$edge.length = tmp_tree$edge.length+0.1
a7_data= fortify(tmp_tree)

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'a'}else if(x %in% c3){'b'}else if(x %in% c4){'c'}else if(x %in% other_c){'d'}else{'virtual'}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, l1), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, l1), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, l1), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, l1), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, l1), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, l1), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, l1), 'a', 'b')})

a1_data$text = sapply(a1_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a2_data$text = sapply(a2_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a3_data$text = sapply(a3_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a4_data$text = sapply(a4_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a5_data$text = sapply(a5_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a6_data$text = sapply(a6_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})
a7_data$text = sapply(a7_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c5){''}else{''}})

cols = c('a'="#F0E685FF", 'b'="#466983FF",'c'= "#BA6338FF",'d'="#5DB1DDFF",'e'="#802268FF", 'virtual'='gray')
a1=ggtree(a1_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, mapping=aes(color=group), size=0.2, branch.length='none')+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 50, units='mm', dpi = 600)

a1_data$text = sapply(a1_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a2_data$text = sapply(a2_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a3_data$text = sapply(a3_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a4_data$text = sapply(a4_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a5_data$text = sapply(a5_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a6_data$text = sapply(a6_data$label,function(x){strsplit(x, '\\.')[[1]][1]})
a7_data$text = sapply(a7_data$label,function(x){strsplit(x, '\\.')[[1]][1]})

a1=ggtree(a1_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.2)+geom_tiplab(aes(label=text),size=0.2)+scale_fill_npg()+scale_color_manual(values = cols)+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2label.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 60, units='mm', dpi = 600)



dist_df = data.frame(dist=c(1,8,3, 7, 2,1, 2),
                     method=c('SCTC', 'MEDALT', 'Sitka', 'MEDICC2', 'NJ', 'MP', 'ML'))
library(forcats)
a=dist_df %>%
  #mutate(method = fct_reorder(method, dist)) %>%
  mutate(method = factor(method, c('SCTC', 'MP', 'ML', 'NJ', 'Sitka', 'MEDICC2', 'MEDALT'))) %>%
  ggplot(aes(x=method, y=dist))+
  geom_point(size=2)+
  geom_text(aes(label=dist), color='white', size=1)+
  labs(y='Reciprocal cells distance')+
  scale_y_continuous(expand = expansion(0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_dist.pdf'),
       a,
       width=24, height = 26, units='mm', dpi = 600)


auc_res = c()
###
cut_tree_clone <- function(tree){
  all_path = nodepath(tree)
  max_path = max(sapply(all_path, length))
  rt_data = fortify(tree) %>% as.data.frame()
  all_layer = c()
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
    ##
    true_clone = c()
    ki=1
    for(tp in tmp_pos){
      if(tp<=length(tree$tip.label)){
        lvs = rt_data[rt_data$node==tp,'label']
      }else{
        lvs = extract.clade(tree,tp)$tip.label
      }
      true_clone = rbind(true_clone, data.frame(clone=ki, name=lvs))
      ki = ki+1
    }
    true_clone$cut_layer = i
    all_layer = rbind(all_layer, true_clone)
  }
  return(all_layer)
}
calc_right_rate <-function(ctc_res, true_label, cut_num=2){
  clone_uniq = unique(true_label$clone)
  clone_right_rate = c()
  for(cu in clone_uniq){
    tmp_bc = true_label[true_label$clone==cu, 'cell']
    cut_k_uniq = unique(ctc_res$cut_layer)
    for(k in cut_k_uniq){
      if(k<=cut_num){
        next
      }
      tmp_ctc_res1 = ctc_res[ctc_res$cut_layer==k, ]
      tmp_ctc_res1$true_label = true_label[match(tmp_ctc_res1$name,true_label$cell), 'clone']
      # clone_kind_rate = tmp_ctc_res1 %>%
      #   group_by(clone) %>%
      #   summarise(unique(true_label)) %>%
      #   group_by(clone) %>%
      #   summarise(n=1/n()) %>% as.data.frame()
      # rownames(clone_kind_rate) = clone_kind_rate$clone
      tmp_ctc_res = tmp_ctc_res1[tmp_ctc_res1$name%in%tmp_bc, ]
      right_rate = table(tmp_ctc_res$clone)
      if(max(right_rate)>=2){
        right_rate = right_rate / nrow(tmp_ctc_res)
        #right_rate = right_rate * clone_kind_rate[names(right_rate), 'n']
        max_cut_num_clone = names(right_rate)[which.max(right_rate)][1]
        nrow_clone = tmp_ctc_res1[tmp_ctc_res1$clone==max_cut_num_clone,]
        right_rate = max(right_rate)/length(unique(nrow_clone$true_label))#/nrow(tmp_ctc_res1[tmp_ctc_res1$clone==max_cut_num_clone,]) #矫正clone所在亚群大小
        clone_right_rate = rbind(clone_right_rate, c(cu, k, right_rate))
      }
    }

  }
  # clone_right_rate = c()
  # for(k in cut_k_uniq){
  #   tmp_ctc_res1 = ctc_res[ctc_res$cut_layer==k, ]
  #   tmp_ctc_res1$true_label = true_label[match(tmp_ctc_res1$name,true_label$cell), 'clone']
  #   ari = pdfCluster::adj.rand.index(tmp_ctc_res1$true_label, tmp_ctc_res1$clone)
  #   clone_right_rate = rbind(clone_right_rate, c(1, k, ari))
  # }
  clone_right_rate = as.data.frame(clone_right_rate)
  colnames(clone_right_rate) = c('true_clone', 'cut_num', 'rate')
  # 矫正线性进化对clone比例的影响
  clone_right_rate = clone_right_rate %>%
    left_join(
      clone_right_rate %>%
        group_by(true_clone) %>%
        summarise(gradien=mean(abs(diff(rate)))))

  cut_right_rate = c()
  cut_k_uniq = unique(ctc_res$cut_layer)
  for(k in cut_k_uniq){
    tmp_ctc_res = ctc_res[ctc_res$cut_layer==k, ]
    tmp_clone = unique(tmp_ctc_res$clone)
    for(tc in tmp_clone){
      tmp_bc = tmp_ctc_res[tmp_ctc_res$clone==tc, 'name']
      tmp_bc_label = true_label[true_label$cell%in%tmp_bc, 'clone']
      cut_right_rate = rbind(cut_right_rate, c(k, tc, vegan::diversity(tmp_bc_label)))
    }
  }
  cut_right_rate = as.data.frame(cut_right_rate)
  colnames(cut_right_rate) = c('cut_num', 'cut_clone', 'diversity')

  return(list(clone_right_rate=clone_right_rate,
              cut_right_rate=cut_right_rate))
}


true_label = data.frame(cell=c(c1,c2,c3,c4,other_c),
                        clone=c(rep(1, length(c1)),rep(1, length(c2)),rep(2, length(c3)),rep(3, length(c4)),rep(4, length(other_c))))

tmp_tree_res1 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sctc),true_label)
tmp_tree_res2 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$medalt),true_label)
tmp_tree_res3 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sitka),true_label)
tmp_tree = all_trees[[file_name]]$medicc2
tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
tmp_tree_res4 = calc_right_rate(cut_tree_clone(tmp_tree),true_label)
tmp_tree_res5 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$NJ),true_label)
tmp_tree_res6 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$MP),true_label)
tmp_tree_res7 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$ML),true_label)

tmp_tree_res1$clone_right_rate$method='SCTC'
tmp_tree_res2$clone_right_rate$method='MEDALT'
tmp_tree_res3$clone_right_rate$method='Sitka'
tmp_tree_res4$clone_right_rate$method='MEDICC2'
tmp_tree_res5$clone_right_rate$method='NJ'
tmp_tree_res6$clone_right_rate$method='MP'
tmp_tree_res7$clone_right_rate$method='ML'

all_clone_right_rate = rbind(tmp_tree_res1$clone_right_rate,
                             tmp_tree_res2$clone_right_rate,
                             tmp_tree_res3$clone_right_rate,
                             tmp_tree_res4$clone_right_rate,
                             tmp_tree_res5$clone_right_rate,
                             tmp_tree_res6$clone_right_rate,
                             tmp_tree_res7$clone_right_rate)
#reshape2::dcast(tmp_tree_res1$clone_right_rate, true_clone+method~cut_num, value.var='rate')
cols = c('1'="#F0E685FF", '2'="#466983FF",'3'= "#BA6338FF",'4'="#5DB1DDFF",'e'="#802268FF", 'virtual'='gray')
a=all_clone_right_rate %>%
  filter(!is.infinite(rate))%>%
  #filter(rate>0.17&rate<1) %>%
  #filter(cut_num>2) %>%
  group_by(true_clone, method) %>%
  summarise(rate=mean(rate*gradien)) %>%
  mutate(method=factor(method, c('SCTC', 'MEDICC2', 'MP', 'NJ', 'ML', 'MEDALT', 'Sitka'))) %>%
  ggplot(aes(x=method, y=rate))+
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_point(aes(fill=as.character(true_clone)), width=0.2, size=0.8, shape=21, stroke=0.1)+
  stat_compare_means(comparisons = list(c('SCTC', 'MEDICC2')), size=1)+
  scale_fill_manual(values = cols)+
  labs(y='Lineage Consistency Score')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_Lineage.pdf'),
       a,
       width=30, height = 30, units='mm', dpi = 600)

# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               border_gp = gpar(col = "gray", lwd=0.05),
               height = nrow(all_cell_cnv)*unit(1, "mm"), # 高度
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.3))),
               col = cn_col,
               column_split = as.numeric(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1])),
               row_split = factor(rownames(all_cell_cnv),rownames(all_cell_cnv)),
               column_gap = unit(0, "mm"),
               row_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               row_title = NULL,
               use_raster = F
)
width_in = 90 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

######## Fig3a  ###########
file_name = 'Fig3a_sub_dataset_P9T_20190409_1.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')

l1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup',
       'SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup',
       'SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup',
       'SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')

c1=c('SC11.P9T.09042019_dedup','SC10.P9T.09042019_dedup','SC17.P9T.09042019_dedup')
c2=c('SC13.P9T.09042019_dedup','SC23.P9T.09042019_dedup')

l1 = c2
other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c2))
other_c = grep('virtual|root', other_c, invert = T, value = T)
#all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
#tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp=all_trees[[file_name]]$sitka
tmp = drop.tip(tmp, setdiff(tmp$tip.label, all_trees[[file_name]]$sctc$tip.label))
a3_data= fortify(tmp)
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
tmp =all_trees[[file_name]]$NJ
tmp$edge.length = tmp$edge.length+1
a5_data= fortify(tmp)
a6_data= fortify(all_trees[[file_name]]$MP)
tmp =all_trees[[file_name]]$ML
tmp$edge.length = tmp$edge.length+1
a7_data= fortify(tmp)

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else if(x %in% c3){'c'}else if(x %in% c4){'d'}else if(x %in% c5){'e'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, l1), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, l1), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, l1), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, l1), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, l1), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, l1), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, l1), 'a', 'b')})

a1_data$text = sapply(a1_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a2_data$text = sapply(a2_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a3_data$text = sapply(a3_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a4_data$text = sapply(a4_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a5_data$text = sapply(a5_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a6_data$text = sapply(a6_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a7_data$text = sapply(a7_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})

cols = c('c'="#F0E685FF", 'b'="#466983FF",'a'= "#BA6338FF",'d'="#5DB1DDFF",'e'="#802268FF", 'virtual'='gray')
a1=ggtree(a1_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 30, units='mm', dpi = 600)


true_label = data.frame(cell=c(c2,other_c),
                        clone=c(rep(1, length(c2)),rep(2, length(other_c))))


tmp_tree_res1 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sctc),true_label,cut_num=1)
tmp_tree_res2 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$medalt),true_label,cut_num=0)
tmp_tree_res3 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sitka),true_label,cut_num=0)
tmp_tree = all_trees[[file_name]]$medicc2
tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
tmp_tree_res4 = calc_right_rate(cut_tree_clone(tmp_tree),true_label,cut_num=0)
tmp_tree_res5 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$NJ),true_label,cut_num=0)
tmp_tree_res6 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$MP),true_label,cut_num=0)
tmp_tree_res7 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$ML),true_label,cut_num=0)

tmp_tree_res1$clone_right_rate$method='SCTC'
tmp_tree_res2$clone_right_rate$method='MEDALT'
tmp_tree_res3$clone_right_rate$method='Sitka'
tmp_tree_res4$clone_right_rate$method='MEDICC2'
tmp_tree_res5$clone_right_rate$method='NJ'
tmp_tree_res6$clone_right_rate$method='MP'
tmp_tree_res7$clone_right_rate$method='ML'

all_clone_right_rate_fig3a = rbind(tmp_tree_res1$clone_right_rate,
                             tmp_tree_res2$clone_right_rate,
                             tmp_tree_res3$clone_right_rate,
                             tmp_tree_res4$clone_right_rate,
                             tmp_tree_res5$clone_right_rate,
                             tmp_tree_res6$clone_right_rate,
                             tmp_tree_res7$clone_right_rate)
#reshape2::dcast(tmp_tree_res1$clone_right_rate, true_clone+method~cut_num, value.var='rate')
cols = c('2'="#466983FF",'1'= "#BA6338FF")
all_clone_right_rate_fig3a$gradien[is.na(all_clone_right_rate_fig3a$gradien)] = 1
a=all_clone_right_rate_fig3a %>%
  filter(!is.infinite(rate))%>%
  #filter(rate>0.17&rate<1) %>%
  #filter(cut_num>2) %>%
  group_by(true_clone, method) %>%
  summarise(rate=mean(rate)) %>%
  mutate(method=factor(method, c('SCTC', 'MEDICC2', 'NJ', 'ML','MEDALT','MP',   'Sitka'))) %>%
  ggplot(aes(x=method, y=rate))+
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_point(aes(fill=as.character(true_clone)), width=0.2, size=0.8, shape=21, stroke=0.1)+
  stat_compare_means(comparisons = list(c('SCTC', 'MEDICC2')), size=1)+
  scale_fill_manual(values = cols)+
  labs(y='Lineage Consistency Score')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_Lineage.pdf'),
       a,
       width=30, height = 30, units='mm', dpi = 600)







dist_df = data.frame(dist=c(1,5,3, 1, 1,3, 1),
                     method=c('SCTC', 'MEDALT', 'Sitka', 'MEDICC2', 'NJ', 'MP', 'ML'))
library(forcats)
a=dist_df %>%
  #mutate(method = fct_reorder(method, dist)) %>%
  mutate(method = factor(method, c('SCTC', 'MEDICC2',  'ML', 'NJ', 'Sitka','MP', 'MEDALT'))) %>%
  ggplot(aes(x=method, y=dist))+
  geom_point(size=2)+
  geom_text(aes(label=dist), color='white', size=1)+
  labs(y='Reciprocal cells distance')+
  scale_y_continuous(expand = expansion(0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_dist.pdf'),
       a,
       width=24, height = 26, units='mm', dpi = 600)
# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               border_gp = gpar(col = "gray", lwd=0.05),
               height = nrow(all_cell_cnv)*unit(1, "mm"), # 高度
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.3))),
               col = cn_col,
               column_split = as.numeric(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1])),
               row_split = factor(rownames(all_cell_cnv),rownames(all_cell_cnv)),
               column_gap = unit(0, "mm"),
               row_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               row_title = NULL,
               use_raster = F
)
width_in = 90 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

####### Fig4a #######
file_name = 'Fig4a_sub_dataset_P9T_4_20200304.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c()#c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c()#c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c()#c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c()#c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')

l1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup',
       'SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup',
       'SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup',
       'SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')

c2=c('X4.SC18.04032020_dedup','X4.SC03.04032020_dedup','X4.SC15.04032020_dedup')
c1=c('X4.SC14.04032020_dedup','X4.SC04.04032020_dedup')

l1 = c2
other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c2))
other_c = grep('virtual|root', other_c, invert = T, value = T)
#all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
#tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp=all_trees[[file_name]]$sitka
tmp = drop.tip(tmp, setdiff(tmp$tip.label, all_trees[[file_name]]$sctc$tip.label))
a3_data= fortify(tmp)
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
tmp =all_trees[[file_name]]$NJ
tmp$edge.length = tmp$edge.length+1
a5_data= fortify(tmp)
a6_data= fortify(all_trees[[file_name]]$MP)
tmp =all_trees[[file_name]]$ML
tmp$edge.length = tmp$edge.length+1
a7_data= fortify(tmp)

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, l1), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, l1), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, l1), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, l1), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, l1), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, l1), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, l1), 'a', 'b')})

a1_data$text = sapply(a1_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a2_data$text = sapply(a2_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a3_data$text = sapply(a3_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a4_data$text = sapply(a4_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a5_data$text = sapply(a5_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a6_data$text = sapply(a6_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})
a7_data$text = sapply(a7_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][2]}else if(x %in% c1){''}else{''}})

cols = c('b'="#466983FF",'a'= "#BA6338FF")
a1=ggtree(a1_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.2, mapping=aes(color=line))+scale_color_manual(values = cols)+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 30, units='mm', dpi = 600)

true_label = data.frame(cell=c(c2,other_c),
                        clone=c(rep(1, length(c2)),rep(2, length(other_c))))


tmp_tree_res1 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sctc),true_label,cut_num=1)
tmp_tree_res2 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$medalt),true_label,cut_num=0)
tmp_tree_res3 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$sitka),true_label,cut_num=0)
tmp_tree = all_trees[[file_name]]$medicc2
tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
tmp_tree_res4 = calc_right_rate(cut_tree_clone(tmp_tree),true_label,cut_num=0)
tmp_tree_res5 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$NJ),true_label,cut_num=0)
tmp_tree_res6 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$MP),true_label,cut_num=0)
tmp_tree_res7 = calc_right_rate(cut_tree_clone(all_trees[[file_name]]$ML),true_label,cut_num=0)

tmp_tree_res1$clone_right_rate$method='SCTC'
tmp_tree_res2$clone_right_rate$method='MEDALT'
tmp_tree_res3$clone_right_rate$method='Sitka'
tmp_tree_res4$clone_right_rate$method='MEDICC2'
tmp_tree_res5$clone_right_rate$method='NJ'
tmp_tree_res6$clone_right_rate$method='MP'
tmp_tree_res7$clone_right_rate$method='ML'

all_clone_right_rate_fig4a = rbind(tmp_tree_res1$clone_right_rate,
                                   tmp_tree_res2$clone_right_rate,
                                   tmp_tree_res3$clone_right_rate,
                                   tmp_tree_res4$clone_right_rate,
                                   tmp_tree_res5$clone_right_rate,
                                   tmp_tree_res6$clone_right_rate,
                                   tmp_tree_res7$clone_right_rate)
#reshape2::dcast(tmp_tree_res1$clone_right_rate, true_clone+method~cut_num, value.var='rate')
cols = c('2'="#466983FF",'1'= "#BA6338FF")
all_clone_right_rate_fig4a$gradien[is.na(all_clone_right_rate_fig4a$gradien)] = 1
a=all_clone_right_rate_fig4a %>%
  filter(!is.infinite(rate))%>%
  #filter(rate>0.17&rate<1) %>%
  #filter(cut_num>2) %>%
  group_by(true_clone, method) %>%
  summarise(rate=mean(rate)) %>%
  mutate(method=factor(method, c('SCTC','MP', 'MEDICC2', 'NJ', 'ML','MEDALT',   'Sitka'))) %>%
  ggplot(aes(x=method, y=rate))+
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_point(aes(fill=as.character(true_clone)), width=0.2, size=0.8, shape=21, stroke=0.1)+
  stat_compare_means(comparisons = list(c('SCTC', 'MP')), size=1)+
  scale_fill_manual(values = cols)+
  labs(y='Lineage Consistency Score')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_Lineage.pdf'),
       a,
       width=30, height = 30, units='mm', dpi = 600)



dist_df = data.frame(dist=c(1,5,3, 1, 1,3, 1),
                     method=c('SCTC', 'MEDALT', 'Sitka', 'MEDICC2', 'NJ', 'MP', 'ML'))
library(forcats)
a=dist_df %>%
  #mutate(method = fct_reorder(method, dist)) %>%
  mutate(method = factor(method, c('SCTC', 'MEDICC2',  'ML', 'NJ', 'Sitka','MP', 'MEDALT'))) %>%
  ggplot(aes(x=method, y=dist))+
  geom_point(size=2)+
  geom_text(aes(label=dist), color='white', size=1)+
  labs(y='Reciprocal cells distance')+
  scale_y_continuous(expand = expansion(0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_dist.pdf'),
       a,
       width=24, height = 26, units='mm', dpi = 600)
# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               border_gp = gpar(col = "gray", lwd=0.05),
               height = nrow(all_cell_cnv)*unit(1, "mm"), # 高度
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.3))),
               col = cn_col,
               column_split = as.numeric(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1])),
               row_split = factor(rownames(all_cell_cnv),rownames(all_cell_cnv)),
               column_gap = unit(0, "mm"),
               row_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               row_title = NULL,
               use_raster = F
)
width_in = 90 / 25.4
height_in = 40 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()

###
all_clone_right_rate$sample = 'P19BT_20181009'
all_clone_right_rate$true_clone[all_clone_right_rate$true_clone==1] = 6
all_clone_right_rate$true_clone[all_clone_right_rate$true_clone==3] = 1
all_clone_right_rate$true_clone[all_clone_right_rate$true_clone==6] = 3

all_clone_right_rate_fig3a$sample = 'P9T_20190409_1'
all_clone_right_rate_fig4a$sample = 'P9T_4_20200304'
cols = c('3'="#F0E685FF", '2'="#466983FF",'1'= "#BA6338FF",'4'="#5DB1DDFF")

all_rate = rbind(all_clone_right_rate, all_clone_right_rate_fig3a, all_clone_right_rate_fig4a)
a=all_rate %>%
  filter(!is.infinite(rate))%>%
  #filter(rate>0.17&rate<1) %>%
  #filter(cut_num>2) %>%
  group_by(true_clone, sample, method) %>%
  summarise(rate=mean(rate*gradien)) %>%
  mutate(method=factor(method, c('SCTC','MP', 'MEDICC2', 'NJ', 'ML','MEDALT',   'Sitka'))) %>%
  ggplot(aes(x=method, y=rate))+
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_point(aes(color=as.character(true_clone), shape=sample), width=0.2, size=1.5,stroke=0.1)+
  stat_compare_means(comparisons = list(c('SCTC', 'MP')), size=2)+
  scale_color_manual(values = cols)+
  labs(y='Lineage Consistency Score')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/Lineageall.pdf'),
       a,
       width=80, height = 40, units='mm', dpi = 600)


######## 评估谱系划分一致性 ######
min_subtree <-function(tree, tips){
  #tree = tmp_tree
  #tips = c1
  st = subtrees(tree)
  min_tree = tree
  min_nodes = length(min_tree$tip.label)
  for(i in st){
    diff_node = setdiff(tips, i$tip.label)
    #print(diff_node)
    if(length(diff_node)==0 & length(i$tip.label)<min_nodes){
      min_tree = i
      min_nodes = length(i$tip.label)
    }
  }
  return(min_tree)
}
all_tmp_dist = c()
##
file_name = 'Fig1c_sub_dataset_P19BT_20181009.txt'
# pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1 = c('SC31.P19BT.09102018_dedup', 'SC30.P19BT.09102018_dedup')
c2 = c('SC62.P19BT.09102018_dedup','SC38.P19BT.09102018_dedup')
c3 = c('SC02.P19BT.09102018_dedup','SC36.P19BT.09102018_dedup','SC37.P19BT.09102018_dedup','SC43.P19BT.09102018_dedup','SC53.P19BT.09102018_dedup','SC65.P19BT.09102018_dedup')
c4 = c('SC18.P19BT.09102018_dedup','SC25.P19BT.09102018_dedup','SC47.P19BT.09102018_dedup','SC51.P19BT.09102018_dedup','SC63.P19BT.09102018_dedup','SC56.P19BT.09102018_dedup')
c5 = c('SC23.P19BT.09102018_dedup','SC14.P19BT.09102018_dedup')
other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c1,c2,c3,c4))
other_c = grep('virtual|root', other_c, invert = T, value = T)

for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, c(c1,c2))
  all_tmp_dist = rbind(all_tmp_dist, c(length(c(c1,c2))/length(mt$tip.label), f, file_name, '3'))
  mt = min_subtree(tmp_tree, c3)
  all_tmp_dist = rbind(all_tmp_dist, c(length(c3)/length(mt$tip.label), f,file_name, '2'))
  mt = min_subtree(tmp_tree, c4)
  all_tmp_dist = rbind(all_tmp_dist, c(length(c4)/length(mt$tip.label), f,file_name, '1'))
  mt = min_subtree(tmp_tree, other_c)
  all_tmp_dist = rbind(all_tmp_dist, c(length(other_c)/length(mt$tip.label), f,file_name, '4'))
}
##
file_name = 'Fig3a_sub_dataset_P9T_20190409_1.txt'
c1=c('SC11.P9T.09042019_dedup','SC10.P9T.09042019_dedup','SC17.P9T.09042019_dedup')
c2=c('SC13.P9T.09042019_dedup','SC23.P9T.09042019_dedup')
other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c2))
other_c = grep('virtual|root', other_c, invert = T, value = T)
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, other_c)
  all_tmp_dist = rbind(all_tmp_dist, c(length(other_c)/length(mt$tip.label), f, file_name, '2'))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(c2)/length(mt$tip.label), f,file_name, '1'))
}
##
file_name = 'Fig4a_sub_dataset_P9T_4_20200304.txt'
c2=c('X4.SC18.04032020_dedup','X4.SC03.04032020_dedup','X4.SC15.04032020_dedup')
c1=c('X4.SC14.04032020_dedup','X4.SC04.04032020_dedup')
other_c = setdiff(rownames(NG_list[[file_name]]$orig.data$all_node_data), c(c2))
other_c = grep('virtual|root', other_c, invert = T, value = T)
for(f in names(all_trees[[file_name]])){
  tmp_tree = all_trees[[file_name]][[f]]
  tmp_tree$tip.label = gsub('_NA_', '', tmp_tree$tip.label)
  mt = min_subtree(tmp_tree, other_c)
  all_tmp_dist = rbind(all_tmp_dist, c(length(other_c)/length(mt$tip.label), f, file_name, '2'))
  mt = min_subtree(tmp_tree, c2)
  all_tmp_dist = rbind(all_tmp_dist, c(length(c2)/length(mt$tip.label), f,file_name, '1'))
}
##

all_tmp_dist = as.data.frame(all_tmp_dist)
colnames(all_tmp_dist) = c('lineage_rate', 'method', 'file', 'clone')
all_tmp_dist$lineage_rate = as.numeric(all_tmp_dist$lineage_rate)
library(forcats)
all_tmp_dist %>%
  mutate(method=fct_reorder(method, min_size, median)) %>%
  mutate(group=ifelse(method=='sctc', 'sctc', 'other')) %>%
  ggplot(aes(x=method,  y=min_size, fill=group)) +
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_jitter(width = 0.2, size=0.5)+
  scale_fill_manual(values = c('sctc'='red', 'other'='gray')) +
  stat_compare_means(comparisons = list(c('sctc', 'MP')), label = 'p.format', method = 't.test')+
  theme_classic()+
  theme(legend.position = 'none')

cols = c('3'="#F0E685FF", '2'="#466983FF",'1'= "#BA6338FF",'4'="#5DB1DDFF")

a=all_tmp_dist %>%
  mutate(method=factor(method, c('sctc', 'NJ', 'ML','medicc2','MP', 'medalt', 'sitka'))) %>%
  ggplot(aes(x=method, y=lineage_rate))+
  geom_boxplot(outlier.shape = NA, size=0.1)+
  geom_jitter(aes(color=as.character(clone), shape=file),  size=1,stroke=0.1)+
  stat_compare_means(comparisons = list(c('sctc', 'NJ')), size=2)+
  scale_color_manual(values = cols)+
  labs(y='Lineage Consistency Score')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/Lineage准确性.pdf'),
       a,
       width=80, height = 40, units='mm', dpi = 600)

df.summary <- all_tmp_dist %>%
  group_by(method) %>%
  summarise(
    sd = sd(lineage_rate, na.rm = TRUE),
    lineage_rate = mean(lineage_rate)
  )

a=ggplot()+
  #stat_compare_means(data=event_coverage %>%
  #                     filter(coverage>0) ,
  #                   mapping = aes(x=chr, y=coverage, fill=type, group=type),
  #                   size=2, label = 'p.format', label.y=0.3)+
  geom_jitter(aes(x=method, y=lineage_rate,color=as.character(clone), shape=file),
              data=all_tmp_dist,  size=1,stroke=0.1, width=0.2)+
  geom_pointrange(data=df.summary, mapping = aes(x=method,y=lineage_rate,
                                                 ymin = lineage_rate-sd, ymax = lineage_rate+sd),
                  position = position_dodge(0.5), size=0.05, color='red',linewidth=0.2)+
  scale_color_manual(values = cols)+
  labs(y='lineage_rate')+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6))
a
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/Lineage准确性2.pdf'),
       a,
       width=80, height = 50, units='mm', dpi = 600)












######### Fig1b #######
file_name = 'Fig1b_dataset_P9T_20190331.txt'
pheatmap::pheatmap(NG_list[[file_name]]$orig.data$all_node_data, cluster_cols = F, cluster_rows = F)
c1=c('SC06.P9T.31032019_dedup', 'SC09.P9T.31032019_dedup')
c2=c('SC22.P9T.31032019_dedup','SC25.P9T.31032019_dedup', 'SC28.P9T.31032019_dedup')
l1 = c2

#all_tmp_dist = c()
#a1_data= fortify(all_trees$Fig1b_dataset_P9T_20190331.txt$sctc)
tmp_tree = all_trees[[file_name]]$sctc
#tmp_tree$edge.length = (tmp_tree$edge.length - min(tmp_tree$edge.length)) / (max(tmp_tree$edge.length) - min(tmp_tree$edge.length))
a1_data = fortify(tmp_tree)
#gtree_res = ggtree(tmp_tree, size=0.2)
#gtree_res$data$group = ifelse(gtree_res$data$label%in%c1,'a','b')

a2_data= fortify(all_trees[[file_name]]$medalt)
tmp=all_trees[[file_name]]$sitka
tmp = drop.tip(tmp, setdiff(tmp$tip.label, all_trees[[file_name]]$sctc$tip.label))
a3_data= fortify(tmp)
a4_data= fortify(all_trees[[file_name]]$medicc2)
a4_data$label = gsub('_NA_', '', a4_data$label)
tmp =all_trees[[file_name]]$NJ
tmp$edge.length = tmp$edge.length+1
a5_data= fortify(tmp)
a6_data= fortify(all_trees[[file_name]]$MP)
tmp =all_trees[[file_name]]$ML
tmp$edge.length = tmp$edge.length+1
a7_data= fortify(tmp)

a1_data$group = sapply(a1_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a2_data$group = sapply(a2_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a3_data$group = sapply(a3_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a4_data$group = sapply(a4_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a5_data$group = sapply(a5_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a6_data$group = sapply(a6_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})
a7_data$group = sapply(a7_data$label,function(x){if(x %in% c1){'a'}else if(x %in% c2){'b'}else{NA}})

a1_data$line = sapply(a1_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sctc,    a1_data, l1), 'a', 'b')})
a2_data$line = sapply(a2_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medalt,  a2_data, l1), 'a', 'b')})
a3_data$line = sapply(a3_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$sitka,   a3_data, l1), 'a', 'b')})
a4_data$line = sapply(a4_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$medicc2, a4_data, l1), 'a', 'b')})
a5_data$line = sapply(a5_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$NJ,      a5_data, l1), 'a', 'b')})
a6_data$line = sapply(a6_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$MP,      a6_data, l1), 'a', 'b')})
a7_data$line = sapply(a7_data$node,function(x){ifelse(x %in% tmp_func(all_trees[[file_name]]$ML,      a7_data, l1), 'a', 'b')})

a1_data$text = sapply(a1_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a2_data$text = sapply(a2_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a3_data$text = sapply(a3_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a4_data$text = sapply(a4_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a5_data$text = sapply(a5_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a6_data$text = sapply(a6_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})
a7_data$text = sapply(a7_data$label,function(x){if(x %in% c2){strsplit(x, '\\.')[[1]][1]}else if(x %in% c1){''}else{''}})

a1=ggtree(a1_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('SCTC')+   theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a2=ggtree(a2_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDALT')+ theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a3=ggtree(a3_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('Sitka')+  theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a4=ggtree(a4_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MEDICC2')+theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a5=ggtree(a5_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('NJ')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a6=ggtree(a6_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('MP')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
a7=ggtree(a7_data, size=0.1, mapping=aes(color=line))+geom_tiplab(aes(label=text),size=1)+scale_fill_npg()+scale_color_manual(values = c('a'='black', 'b'='gray'))+ggtitle('ML')+     theme(plot.title = element_text(hjust = 0.5,size=6), legend.position = 'none')
#cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=1)
#p2 = cowplot::plot_grid(a2, a3, a4, a5, a6, a7, nrow=3)
#
#pp=cowplot::plot_grid(a1, p2, nrow=1)
ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'2.pdf'),
       cowplot::plot_grid(a1, a2, a3, a4, a5, a6, a7, nrow=2),
       width=90, height = 60, units='mm', dpi = 600)


dist_df = data.frame(dist=c(1,5,3, 1, 1,3, 1),
                     method=c('SCTC', 'MEDALT', 'Sitka', 'MEDICC2', 'NJ', 'MP', 'ML'))
library(forcats)
a=dist_df %>%
  #mutate(method = fct_reorder(method, dist)) %>%
  mutate(method = factor(method, c('SCTC', 'MEDICC2',  'ML', 'NJ', 'Sitka','MP', 'MEDALT'))) %>%
  ggplot(aes(x=method, y=dist))+
  geom_point(size=2)+
  geom_text(aes(label=dist), color='white', size=1)+
  labs(y='Reciprocal cells distance')+
  scale_y_continuous(expand = expansion(0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle = 60, hjust = 1, vjust=1),
        legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        legend.position = 'none')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/',file_name,'_dist.pdf'),
       a,
       width=24, height = 26, units='mm', dpi = 600)
# heatmap
sctc = NG_list[[file_name]]
all_cell_cnv = sctc$orig.data$all_node_data
dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
all_cell_cnv = all_cell_cnv[lab, ]

lr = rep(0, nrow(all_cell_cnv))
lr[rownames(all_cell_cnv) %in% c(c1,c2)] = 1
lr_shape = lr
lr_shape[lr_shape==0] = NA
ht2 <- Heatmap(all_cell_cnv,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = F,
               border_gp = gpar(col = "gray", lwd=0.05),
               height = nrow(all_cell_cnv)*unit(1, "mm"), # 高度
               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"),
                                                                   labels = c(1:22, 'X'),
                                                                   labels_gp = gpar(cex = 0.3))),
               left_annotation = rowAnnotation(Highlight = anno_simple(rep('0', nrow(all_cell_cnv)),
                                                                       col = c('0'='white'),
                                                                       pch = lr_shape,
                                                                       pt_gp = gpar(col = "red"))),
               col = cn_col,
               column_split = as.numeric(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]][1])),
               row_split = factor(rownames(all_cell_cnv),rownames(all_cell_cnv)),
               column_gap = unit(0, "mm"),
               row_gap = unit(0, "mm"),
               heatmap_legend_param = list(
                 title = "CNV",
                 color_bar = "discrete",
                 #at = c(0,1,2),
                 #title_position = "leftcenter-rot",
                 legend_height = unit(3, "cm")
               ),
               column_title = NULL,
               row_title = NULL,
               use_raster = F
)
width_in = 90 / 25.4
height_in = 80 / 25.4
dev.off()
pdf(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图NG/不同方法比较划相同祖先划分准确性/hp_',file_name,'.pdf'), width=width_in, height=height_in)
draw(ht2)
dev.off()








#####################

