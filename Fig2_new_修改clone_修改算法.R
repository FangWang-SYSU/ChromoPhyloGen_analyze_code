source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')
# load pathway gene list
#############
DNA_list = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc_all.rds')

tmp_file = 'huh7'
sctc = DNA_list[[tmp_file]]

dna_barcodes = rownames(sctc$orig.data$all_node_data)
dna_barcodes = dna_barcodes[!grepl('virtual|root', dna_barcodes)]
'/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/DNA_F1_bc.txt'
f1_bc = c()
f2_bc = c()
for(i in dna_barcodes){
  x = strsplit(i, '_')[[1]]
  if(x[1] =='f1'){
    f1_bc = rbind(f1_bc, x[2])
  }else{
    f2_bc = rbind(f2_bc, x[2])
  }
}
write.table(f1_bc, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/DNA_F1_bc.txt',
            quote = F, row.names = F, col.names = F)
write.table(f2_bc, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/DNA_F2_bc.txt',
            quote = F, row.names = F, col.names = F)

mtx = sctc$orig.data$all_node_data[dna_barcodes, ]
write.table(mtx, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/DNA_CNA_segment_matrix.txt',
            quote = F, row.names = T, col.names = T)


# 定义模式
# gene.loc <-  read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
scTrace_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine_all_修改算法/'

a1 = read.table(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine_all//{tmp_file}mode.txt'), header = T, sep=',', row.names = "X")
a2 = read.table(glue('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/MyCellLine_all2//{tmp_file}mode.txt'), header = T, sep=',', row.names = "X")

CNA_mechnism = read.table(glue('{scTrace_dir}/{tmp_file}mode.txt'), header = T, sep=',', row.names = "X")
all_rearrange_score = read.table(glue('{scTrace_dir}/{tmp_file}re_score.txt'), header = T, sep=',', row.names = "X")
all_limit_prop = read.table(glue('{scTrace_dir}/{tmp_file}limit_score.txt'), header = T, sep=',', row.names = "X")
all_rearrange_score_pvalue = read.table(glue('{scTrace_dir}/{tmp_file}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")

all_seg_len = c()
for(cell in rownames(sctc$orig.data$all_node_data)){
  print(cell)
  tmp_cell_len = c()
  for(chr in gsub('X','', colnames(all_rearrange_score))){
    tmp_pos = grepl(paste0('^',chr,'_'), colnames(sctc$orig.data$all_node_data))
    tmp_dt = unlist(sctc$orig.data$all_node_dat[cell, tmp_pos])
    tmp_cell_len = c(tmp_cell_len,sum(diff(tmp_dt)!=0))
  }
  all_seg_len = rbind(all_seg_len, tmp_cell_len)
}
rownames(all_seg_len) = rownames(sctc$orig.data$all_node_data)
colnames(all_seg_len) = colnames(all_rearrange_score)

all_seg_len2 = all_seg_len[rownames(CNA_mechnism), ]
min_seg_num = 10
pvalue_thr = 0.001
CNA_mechnism$chromothripsis = rowMeans(all_rearrange_score)
CNA_mechnism$CPS_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_rearrange_score>0 & all_seg_len2>min_seg_num)
CNA_mechnism$limit_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5 & all_rearrange_score>0&all_seg_len2>min_seg_num)
CNA_mechnism$seismic_num = rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5 & all_rearrange_score>0&all_seg_len2>min_seg_num)
CNA_mechnism$BFB = log1p(CNA_mechnism$BFB) #/ ncol(sctc$orig.data$all_node_data)

a1=rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5 & all_rearrange_score>0&all_seg_len2>10)
a2=rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5 & all_rearrange_score>0)
boxplot(a1,a2)

a1=rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5 & all_rearrange_score>0&all_seg_len2>10)
a2=rowSums(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5 & all_rearrange_score>0)
boxplot(a1,a2)

#
all_rearrange_score2 = all_rearrange_score
all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5 & all_rearrange_score>0&all_seg_len2>min_seg_num)] = NA
limit_thr = min(all_rearrange_score2, na.rm = T)
CNA_mechnism$limit_score = rowMeans(all_rearrange_score2, na.rm = T)
# 识别gradual, seismic region #
all_region_limit = list()
for(i in rownames(all_rearrange_score2)){
  print(i)
  tmp_list = list()
  for(chr in 1:ncol(all_rearrange_score2)){
    if(!is.na(all_rearrange_score2[i, chr])){
      tmp_x = unlist(sctc$orig.data$all_node_data[i, grepl(paste0('^',chr,'_'), colnames(sctc$orig.data$all_node_data))])
      region_res = find_event_region(tmp_x, window_size=round(length(tmp_x)/10))
      tmp_list[[as.character(chr)]] = region_res
    }
  }
  all_region_limit[[i]] = tmp_list
}

all_rearrange_score2 = all_rearrange_score
all_rearrange_score2[!(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5 & all_rearrange_score>0&all_seg_len2>min_seg_num)] = NA
seismic_thr = min(all_rearrange_score2, na.rm = T)
CNA_mechnism$seismic_score = rowMeans(all_rearrange_score2, na.rm = T)
# 识别gradual, seismic region #
all_region_seismic = list()
for(i in rownames(all_rearrange_score2)){
  print(i)
  tmp_list = list()
  for(chr in 1:ncol(all_rearrange_score2)){
    if(!is.na(all_rearrange_score2[i, chr])){
      tmp_x = unlist(sctc$orig.data$all_node_data[i, grepl(paste0('^',chr,'_'), colnames(sctc$orig.data$all_node_data))])
      region_res = find_event_region(tmp_x, window_size=round(length(tmp_x)/10))
      tmp_list[[as.character(chr)]] = region_res
    }
  }
  all_region_seismic[[i]] = tmp_list
}

CNA_mechnism[is.na(CNA_mechnism)] = 0
# 修改clone信息
sctc$map_obj$cell_map_pos[sctc$map_obj$cell_map_pos$clone=='clone_3', 'clone'] = 'clone_2'
sctc$map_obj$cell_map_pos[sctc$map_obj$cell_map_pos$clone=='clone_4', 'clone'] = 'clone_3'
sctc$map_obj$cell_map_pos[sctc$map_obj$cell_map_pos$clone=='clone_5', 'clone'] = 'clone_4'


pde = predict_CEE(sctc)

CNA_mechnism$pde = pde[rownames(CNA_mechnism)]
CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

sctc$map_obj$cell_map_pos[, colnames(CNA_mechnism)] = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), ]

leaves_data = sctc$orig.data$all_node_data[grepl('f',rownames(sctc$orig.data$all_node_data)), ]
mean(as.vector(as.matrix(leaves_data)), na.rm=T)
mean(rowMeans(leaves_data))
##### 单细胞进化树 rearrange着色 ####
colorby = 'chromothripsis'
gtree_res = ggtree(sctc$orig.data$tree)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]

RE_p=gtree_res +
  aes(color=group, size=isTip)+
  scale_color_gradientn(colours = c('white', '#0060BF'))+
  scale_size_manual(values = c('TRUE'=0.2, 'FALSE'=0.1))+
  labs(title = 'chromothripsis')+
  theme_tree(legend.position='right')+
  labs(title = 'CEE')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )

#c("#045a8d","white","#a50f15")
RE_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树chromothripsis着色.pdf',
       RE_p, width=40, height = 60, units='mm', dpi = 600)



colorby = 'CPS_num'
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]
cps_p=gtree_res +
  aes(color=group, size=0.2)+
  scale_color_gradientn(colours = c('yellow', 'red'))+
  labs(title = 'CPS')+
  theme_tree(legend.position='right')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
cps_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树cps着色.pdf',
       cps_p, width=40, height = 60, units='mm', dpi = 600)


##### 单细胞进化树 BFB着色 ####
colorby = 'BFB'
# new_BFB_res = BFB_filter(sctc)
# rownames(new_BFB_res) = new_BFB_res$cellname
#cell_info = sctc$map_obj$cell_map_pos[, colorby, drop=F]
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]#log1p(new_BFB_res[gsub('\\.','-',gtree_res$data$label), colorby])

BFB_p=gtree_res +
  aes(color=group, size=isTip)+
  scale_color_gradientn(colours = c('white', '#F4BA19'))+
  theme_tree(legend.position='right')+
  labs(title = 'BFB')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
BFB_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树BFB着色.pdf',
       BFB_p, width=40, height = 60, units='mm', dpi = 600)

##### 单细胞进化树 limit_score着色 ####
colorby = 'limit_num'
#cell_info = sctc$map_obj$cell_map_pos[, colorby, drop=F]
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]#log1p(new_BFB_res[gsub('\\.','-',gtree_res$data$label), colorby])
## make sure the order of the input matrix is consistent with the tree
limit_p=gtree_res +
  aes(color=group, size=isTip)+
  scale_color_gradientn(colours = c('yellow', 'red'))+
  labs(title = 'Gradual')+
  #scale_color_gradientn(colours = c('#DCE9F5', 'red'))+
  theme_tree(legend.position='right')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
limit_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树Gradual_score着色.pdf',
       limit_p, width=40, height = 80, units='mm', dpi = 600)

##### 单细胞进化树 limit_score着色 ####
colorby = 'seismic_num'
#cell_info = sctc$map_obj$cell_map_pos[, colorby, drop=F]
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]#log1p(new_BFB_res[gsub('\\.','-',gtree_res$data$label), colorby])
## make sure the order of the input matrix is consistent with the tree
seismic_p=gtree_res +
  aes(color=group, size=isTip)+
  scale_color_gradientn(colours = c('yellow', 'red'))+
  labs(title = 'Seismic')+
  #scale_color_gradientn(colours = c('#DCE9F5', 'red'))+
  theme_tree(legend.position='right')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
seismic_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树seismic_score着色.pdf',
       seismic_p, width=40, height = 80, units='mm', dpi = 600)

### PDE着色 ####
colorby = 'pde'
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]
pde_p=gtree_res +
  aes(color=group, size=isTip)+
  #scale_color_gradientn(colours = c('#045a8d', 'red'), limits=c(0,0.7))+
  scale_colour_gradientn(colours = c("#3288BD" ,"#FEE08B","#9E0142"))+
  theme_tree(legend.position='right')+
  labs(title = 'CEE')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
pde_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树CEE着色.pdf',
       pde_p, width=40, height = 60, units='mm', dpi = 600)


### clone着色 ####
colorby = 'clone'
gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
gtree_res$data$group = CNA_mechnism[gsub('\\.','-',gtree_res$data$label), colorby]
clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')
clone_p=gtree_res +
  aes(color=group, size=0.2)+
  scale_color_manual(values=clone_col)+
  labs(title = 'Clone')+
  # scale_color_gradientn(colours = c('yellow', 'red'))+
  theme_tree(legend.position='right')+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
clone_p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树Clone着色.pdf',
       clone_p, width=40, height = 60, units='mm', dpi = 600)

a = cowplot::plot_grid(clone_p, RE_p,BFB_p,pde_p, nrow=1)
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树score着色.pdf',
       a, width=160, height = 60, units='mm', dpi = 600)


###
colorby = 'clone'
cell_info = sctc$map_obj$cell_map_pos[!grepl('root|virtual', rownames(sctc$map_obj$cell_map_pos)), colorby, drop=F]
gtree_res = ggtree(sctc$orig.data$tree)
gtree_res$data$group = cell_info[gtree_res$data$label, 'clone']
gtree_res = ggtree(gtree_res$data, aes(color=group), size=0.2)

clone_col = setNames(pal_igv()(length(unique(cell_info$clone))), unique(cell_info$clone))
clone_col = c(clone_col)
p=gtree_res +
  scale_color_manual(values=clone_col)+
  # scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  theme_tree(legend.position='right')+
  guides(color=guide_legend(title='Clone'))+
  #coord_flip()+
  #scale_x_reverse()+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

score_data = sctc$map_obj$cell_map_pos[, c('chromothripsis', 'wgd', 'CPS_num','BFB',
                                           'limit_num', 'seismic_num', 'limit_score','seismic_score', 'pde')]
score_data$nameid = rownames(score_data)
score_data$wgd = factor(score_data$wgd, c('WGD0', 'WGD1', 'WGD2'))
#tile_data = reshape2::melt(score_data[, c(2,5,6,7,8)], id='nameid')
score_data = score_data[p$data$label[p$data$isTip], ]
score_data$x1=0
#score_data$x2=score_data$x1+1
bar_data = reshape2::melt(score_data[, c('limit_num', 'seismic_num', 'nameid')], id='nameid')
bar_data$variable = factor(bar_data$variable, levels = c( 'seismic_num','limit_num'))
dd <- p$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
score_data = score_data[lab, ]

library(ggnewscale)
library(ggtreeExtra)
p2 = p+geom_tippoint(size=1)
p2 = p2 +
  new_scale_fill() +
  geom_fruit(data=score_data, geom=geom_tile,
             mapping=aes(y=nameid, x=x1,  fill=wgd),
             color = "grey",
             offset = 0.05,size = 0.05, width=100,
             axis.params=list(
               axis="x",
               text = "WGD",
               text.angle = 45,
               text.size = 2,
               line.size = 0,
               vjust = 1,
               hjust= 1,
               inherit.aes=F,
             ))+
  scale_fill_npg()+
  scale_fill_manual(values = c('WGD0'='#DEDEDE', 'WGD1'='#FFE53B', 'WGD2'='#FF2525'),drop = FALSE)+
  new_scale_fill() +
  geom_fruit(data=score_data, geom=geom_tile,
             mapping=aes(y=nameid, x=x1,  fill=BFB),
             color = "grey",
             offset = 0.05,size = 0.05, width=100,
             axis.params=list(
               axis="x",
               text = "BFB",
               text.angle = 45,
               text.size = 2,
               line.size = 0,
               vjust = 1,
               hjust= 1,
               inherit.aes=F,
             ))+
  scale_fill_gradientn(colours = c('white', '#F4BA19'))+
  new_scale_fill() +
  geom_fruit(data=score_data, geom=geom_tile,
             mapping=aes(y=nameid, x=x1,  fill=chromothripsis),
             color = "grey",
             offset = 0.05,size = 0.05,width=100,
             axis.params=list(
               axis="x",
               text = "RE",
               text.angle = 45,
               text.size = 2,
               line.size = 0,
               vjust = 1,
               hjust= 1,
               inherit.aes=F,
             ))+
  scale_fill_gradientn(colours = c('white', '#306DAA'))+
  new_scale_fill() +
  geom_fruit(mapping=aes(y=nameid, x=value, fill=variable),
             data=bar_data,
             geom=geom_bar,
             offset = 0.05,
             pwidth=0.1,
             size=0.01,
             color='gray',
             #fill='red',
             orientation="y",
             stat="identity",
             axis.params=list(
               axis="x",
               #text.angle = -45,
               text.size = 2,
               #line.size = 0,
               vjust = 1,
               inherit.aes=F,
             )
  )+
  scale_fill_manual(values = c('limit_num'='#08AEEA', 'seismic_num'='#B721FF'),drop = FALSE, name='chr_num')
p2 = p2+    theme(legend.text = element_text(size=4),
                  legend.key.size = unit(2, 'mm'),
                  legend.title = element_text(size=6)
)
p2 = p2+scale_y_continuous(expand = expansion(mult=0.1))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树多指标.pdf',
       p2, width=200, height = 160, units='mm', dpi = 600)


###


# ##### 单细胞进化树 is_chromothripsis 着色 ####
# chromothripsis_chr = CNA_mechnism$orig.rearrange_score$all_rearrange_score
# #chromothripsis_chr_pvalue = t(pvalue_rearrange_score(chromothripsis_chr))
# #chromothripsis_chr[chromothripsis_chr_pvalue>=0.01] = 0
# #CNA_mechnism$res$chromothripsis_chr_prop = rowSums(chromothripsis_chr_pvalue<0.05) #/ ncol(chromothripsis_chr)
#
# CNA_mechnism$res$chromothripsis_chr_prop = rowSums(chromothripsis_chr>quantile(chromothripsis_chr, 0.75))
#
# CNA_mechnism$res$chromothripsis_chr_prop = cut(CNA_mechnism$res$chromothripsis_chr_prop, breaks = c(-Inf, 4, 7, Inf),
#                                                labels = c('0-4', '5-7', '8-11'))
# # cell_info = sctc$map_obj$cell_map_pos[, colorby, drop=F]
# colorby = 'chromothripsis_chr_prop'
# gtree_res = ggtree(sctc$orig.data$tree, size=0.2)
# x = as.vector(CNA_mechnism$res[gsub('\\.','-',gtree_res$data$label), colorby])
# x[is.na(x)] = '0-4'
# gtree_res$data$group = x
# ## make sure the order of the input matrix is consistent with the tree
# p=gtree_res +
#   aes(color=group, size=isTip)+
#   #scale_color_gradientn(colours = c('yellow', 'red'))+
#   scale_color_manual(values = c('0-4'='#264399', '5-7'='#60993D', '8-11'='#BF6E4C'))+
#   theme_tree(legend.position='right')+
#   #coord_flip()+
#   #scale_x_reverse()+
#   theme(legend.text = element_text(size=4),
#         legend.key.size = unit(2, 'mm'),
#         legend.title = element_text(size=6)
#   )
# p
# ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树chromothripsis数量着色.pdf',
#        p, width=40, height = 60, units='mm', dpi = 600)
#

#### clone 进化树热图#######
# 1.进化树
# clone_col = setNames(pal_npg()(length(unique(cna_meta$clone))),unique(cna_meta$clone))
clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')

clone_names = c("clone_4", "clone_3", "clone_1", "clone_2")

y = sctc$clone_graph_obj$clone_graph
y=drop.tip(y, c("clone_3",'clone_2'),trim.internal = F)
y$tip.label = c("clone_3", "clone_4", "clone_1", "clone_2")
y <- ape::rotateConstr(y, rev(clone_names))

#更新clone分支长度
sctc$clone_graph_obj$cell_phylo_tree$data %>%
  filter(isTip==TRUE) %>%
  group_by(group) %>%
  summarise(max(branch.length))

tmp_meta = sctc$orig.data$cell_relation
tmp_meta$clone = sctc$map_obj$cell_map_pos[rownames(tmp_meta), 'clone']
tmp_meta%>%
  group_by(clone) %>%
  summarise(mean(Pseudotime_tree))
xx= y$edge.length
y$edge.length = c(375+1935, 408, 467+3960, 424, 276+2142, 258+2465)
plot(y)
##


clone_tree = ggtree(y)
clone_tree$data$group=clone_tree$data$label
clone_tree = clone_tree +
  geom_point2(fill='#BAD3E9',shape=21, size=3, color='black')+
  geom_tippoint(aes(fill=group),shape=21, size=5, color='black')+
  # geom_label(aes(x=branch, label=angle), fill='lightgreen')+
  scale_fill_manual(values = clone_col)+
  theme(legend.position = 'none')
clone_tree

# 添加BFB区域
###
##### 准备染色体信息 ######
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart

chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum


# 位点转换成基因
cna_data = sctc$orig.data$all_node_data
cna_loc = as.data.frame(t(sapply(colnames(cna_data), function(x)strsplit(x, '_')[[1]])))
cna_loc$V2 = as.numeric(cna_loc$V2)
cna_loc$V3 = as.numeric(cna_loc$V3)

gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)
gene_info = na.omit(gene_info)
new_gene_loc = c()
for(i in 1:nrow(gene_info)){
  print(i)
  gn = gene_info$gene_name[i]
  gs = gene_info$start[i]
  ge = gene_info$end[i]
  g_chr = gene_info$seqnames[i]
  gabs = as.numeric(chrinfo_chr_absstart[g_chr, 'absStart'])+as.numeric(gs)
  gabend = as.numeric(chrinfo_chr_absstart[g_chr, 'absStart'])+as.numeric(ge)
  tmp_loc = cna_loc[cna_loc$V1==g_chr, ]
  pos1 = nrow(tmp_loc[tmp_loc$V2<gabs, ])
  pos2 = nrow(tmp_loc[tmp_loc$V2<gabend, ])
  armname = paste0(g_chr, '_', gene_info$arm[i])
  new_gene_loc = rbind(new_gene_loc, data.frame(loc=rownames(tmp_loc)[pos1:pos2], gene=gn, arm=armname))

}
new_gene_loc = as.data.frame(new_gene_loc)
new_gene_loc = na.omit(new_gene_loc)

# cna gene level
cna_data_gene = t(cna_data)
cna_data_gene = as.data.frame(cna_data_gene[new_gene_loc$loc, ])
cna_data_gene$arm = new_gene_loc$arm
cna_data_gene = cna_data_gene %>%
  group_by(arm) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(cna_data_gene) = cna_data_gene[, 1]
cna_data_gene = as.data.frame(t(cna_data_gene[,-1]))

#
cell_mata = sctc$map_obj$cell_map_pos
clone_cnv = cna_data_gene
clone_cnv$clone = cell_mata[rownames(clone_cnv), 'clone']
clone_cnv = clone_cnv %>%
  group_by(clone) %>%
  summarise_all(list(mean)) %>% as.data.frame()
clone_cnv = na.omit(clone_cnv)
rownames(clone_cnv) = clone_cnv[,1]
clone_cnv = clone_cnv[,-1]
clone_cnv = round(clone_cnv)
#
tree_data = clone_tree$data %>% as.data.frame()

bfb_record = c()
node_queue = c('root')
while(length(node_queue)!=0){
  curr_node = node_queue[length(node_queue)]
  node_queue = node_queue[-length(node_queue)]
  curr_num = tree_data[tree_data$label==curr_node, 'node']
  curr_childs = tree_data[tree_data$parent==curr_num, 'label']
  curr_childs = setdiff(curr_childs, curr_node)
  if(length(curr_childs)==2){
    node_queue = c(node_queue, curr_childs)
    # for(cc in curr_childs){
    #   c1_data = clone_cnv[cc, ]
    #   tmp_BFB_pos1 = (c1_data!=clone_cnv[curr_node,])#((c1_data+c2_data)%%2)==0 & (c1_data!=c2_data)
    #   tmp_BFB_pos = colnames(clone_cnv)[tmp_BFB_pos1]
    #   tmp_bfb_record1 = c1_data[tmp_BFB_pos1] - clone_cnv[curr_node,][tmp_BFB_pos1]
    #   tmp_bfb_record1 = paste0(paste0(tmp_BFB_pos, tmp_bfb_record1, ','), collapse = ';')
    #   bfb_record = rbind(bfb_record, c(cc, tmp_bfb_record1))
    # }
    c1_data = clone_cnv[curr_childs[1], ]
    c2_data = clone_cnv[curr_childs[2], ]
    tmp_BFB_pos1 = ((c1_data+c2_data)%%2)==0 & (c1_data!=c2_data)#(c1_data!=c2_data)#
    tmp_BFB_pos = colnames(clone_cnv)[tmp_BFB_pos1]
    print(tmp_BFB_pos)
    tmp_bfb_record1 = c1_data[tmp_BFB_pos1] - c2_data[tmp_BFB_pos1]
    tmp_bfb_record2 = c2_data[tmp_BFB_pos1] - c1_data[tmp_BFB_pos1]
    tmp_bfb_record1 = paste0(paste0(tmp_BFB_pos, tmp_bfb_record1, ','), collapse = ';')
    tmp_bfb_record2 = paste0(paste0(tmp_BFB_pos, tmp_bfb_record2, ','), collapse = ';')
    bfb_record = rbind(bfb_record, c(curr_childs[1], tmp_bfb_record1))
    bfb_record = rbind(bfb_record, c(curr_childs[2], tmp_bfb_record2))
  }
}

#
clone_tree$data$BFB=bfb_record[match(clone_tree$data$label, bfb_record[,1]), 2]
clone_tree = clone_tree +
  geom_point2(fill='#BAD3E9',shape=21, size=3, color='black')+
  geom_tippoint(aes(fill=group),shape=21, size=5, color='black')+
  geom_label(aes(x=branch, label=BFB), fill='lightgreen')+
  scale_fill_manual(values = clone_col)+
  theme(legend.position = 'none')
clone_tree
##
# clonetree_layer = layer_data(clone_tree)
# clonetree_layer <- clonetree_layer$node[order(clonetree_layer$y, decreasing = T)]
# clone_tree_data = as.data.frame(clone_tree$data)
# clone_names = clone_tree_data[match(clonetree_layer,clone_tree_data$node), ]
# clone_names = grep('clone', clone_names$group, value = T)

# 2.热图
rescore = all_rearrange_score
rescore[all_rearrange_score_pvalue<0.001] = 0
rescore[all_seg_len<10] = 0
rescore = rescore[grepl('f', rownames(rescore)),]
#
cell_mata = sctc$map_obj$cell_map_pos
clone_cnv = sctc$orig.data$all_node_data
clone_cnv$clone = cell_mata[rownames(clone_cnv), 'clone']
clone_cnv = clone_cnv %>%
  group_by(clone) %>%
  summarise_all(list(mean)) %>% as.data.frame()
clone_cnv = clone_cnv[grepl('clone', clone_cnv$clone), ]
rownames(clone_cnv) = clone_cnv[,1]
clone_cnv = clone_cnv[,-1]
clone_cnv = round(clone_cnv)
#
clone_rescore = as.data.frame(all_rearrange_score)
# clone_rescore[all_rearrange_score_pvalue<0.001] = NA
# clone_rescore[all_seg_len<10] = NA
# clone_rescore = clone_rescore[grepl('f', rownames(clone_rescore)),]
clone_rescore$clone = cell_mata[rownames(clone_rescore), 'clone']
clone_rescore = clone_rescore %>%
  group_by(clone) %>%
  summarise_all(list(mean), na.rm=T) %>% as.data.frame()
clone_rescore = clone_rescore[grepl('clone', clone_rescore$clone), ]
rownames(clone_rescore) = clone_rescore[,1]
clone_rescore = clone_rescore[,-1]
old_clone_rescore = clone_rescore

clone_rescore = old_clone_rescore

# limit_mat = (all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)+0
# # limit_mat[limit_mat==1] = 2
# seimic_mat = (all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)+0
#
# clone_rescore = limit_mat+seimic_mat
# clone_rescore = as.data.frame(clone_rescore)
# clone_rescore$clone = cell_mata[rownames(clone_rescore), 'clone']
#
# tmp_func<-function(x){
#   aa = sort(table(x), decreasing = T)
#   as.numeric(names(aa)[1])
# }
# clone_rescore = clone_rescore %>%
#   dplyr::group_by(clone) %>%
#   dplyr::summarise_all(list(tmp_func)) %>% as.data.frame()
# clone_rescore = clone_rescore[grepl('clone', clone_rescore$clone), ]
# rownames(clone_rescore) = clone_rescore[,1]
# clone_rescore = clone_rescore[,-1]

# old_clone_rescore = clone_rescore
# # tmp_thr = quantile(as.matrix(clone_rescore), 0.75)
# # clone_rescore[clone_rescore<tmp_thr] = 0
# # clone_rescore[clone_rescore>0] = 1
# #
# # clone_rescore = as.matrix(old_clone_rescore)
# # clone_rescore_pvalue = pvalue_rearrange_score(clone_rescore)
# # clone_rescore[clone_rescore_pvalue>=0.001] = 0
# # clone_rescore[clone_rescore>0] = 1
#
# # limit rate
# limit_rate = as.data.frame(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop>0.5)
# limit_rate$clone = cell_mata[rownames(limit_rate), 'clone']
# limit_rate = limit_rate %>%
#   group_by(clone) %>%
#   summarise_all(list(mean)) %>% as.data.frame()
# limit_rate = limit_rate[grepl('clone', limit_rate$clone), ]
# rownames(limit_rate) = limit_rate[,1]
# limit_rate = limit_rate[,-1]
#
# # seismic rate
# seismic_rate = as.data.frame(all_rearrange_score_pvalue<pvalue_thr & all_limit_prop<=0.5)
# seismic_rate$clone = cell_mata[rownames(seismic_rate), 'clone']
# seismic_rate = seismic_rate %>%
#   group_by(clone) %>%
#   summarise_all(list(mean)) %>% as.data.frame()
# seismic_rate = seismic_rate[grepl('clone', seismic_rate$clone), ]
# rownames(seismic_rate) = seismic_rate[,1]
# seismic_rate = seismic_rate[,-1]
#
# # tmp_thr = quantile(as.matrix(limit_rate), 0.9)
# # limit_rate[limit_rate<limit_thr] = 0
# # limit_rate[limit_rate>0] = 1
# # limit_rate = as.matrix(limit_rate)
# # limit_rescore_pvalue = pvalue_rearrange_score(limit_rate)
# # limit_rate[limit_rescore_pvalue>=0.001] = 0
# # limit_rate[limit_rate>0] = 2
#
# clone_rescore = !(limit_rate==0 & seismic_rate== 0)+0
# clone_rescore = (limit_rate>seismic_rate)+clone_rescore
#clone_rescore[clone_rescore==3] = 1
# # 识别BFB区域
# ###
# c1_data = clone_cnv['clone_2', ]
# c2_data = clone_cnv['clone_3', ]
# tmp_BFB_pos1 = ((c1_data+c2_data)%%2)==0 & (c1_data!=c2_data)
# tmp_BFB_pos = colnames(cna_data)[tmp_BFB_pos1]
# tmp_BFB_pos = t(sapply(tmp_BFB_pos, function(x)strsplit(x, '_')[[1]]))
# tmp_BFB_pos = as.data.frame(tmp_BFB_pos)
# tmp_BFB_pos$V1 = as.numeric(tmp_BFB_pos$V1)
# tmp_BFB_pos$V2 = as.numeric(tmp_BFB_pos$V2)
# tmp_BFB_pos$V3 = as.numeric(tmp_BFB_pos$V3)
# tmp_BFB_pos$CN1 = c1_data[tmp_BFB_pos1]
# tmp_BFB_pos$CN2 = c2_data[tmp_BFB_pos1]
# tmp_BFB_pos23 = tmp_BFB_pos%>%arrange(V1, V2)
# ####
# c1_data = clone_cnv['clone_4', ]
# c2_data = clone_cnv['clone_5', ]
# tmp_BFB_pos1 = ((c1_data+c2_data)%%2)==0 & (c1_data!=c2_data)
# tmp_BFB_pos = colnames(cna_data)[tmp_BFB_pos1]
# tmp_BFB_pos = t(sapply(tmp_BFB_pos, function(x)strsplit(x, '_')[[1]]))
# tmp_BFB_pos = as.data.frame(tmp_BFB_pos)
# tmp_BFB_pos$V1 = as.numeric(tmp_BFB_pos$V1)
# tmp_BFB_pos$V2 = as.numeric(tmp_BFB_pos$V2)
# tmp_BFB_pos$V3 = as.numeric(tmp_BFB_pos$V3)
# tmp_BFB_pos$CN1 = c1_data[tmp_BFB_pos1]
# tmp_BFB_pos$CN2 = c2_data[tmp_BFB_pos1]
# tmp_BFB_pos45 = tmp_BFB_pos%>%arrange(V1, V2)
# ####
# c1_data = clone_cnv['clone_1', ]
# c2_data = round((clone_cnv['clone_5', ]+clone_cnv['clone_4', ])/2)
# tmp_BFB_pos1 = ((c1_data+c2_data)%%2)==0 & (c1_data!=c2_data)
# tmp_BFB_pos = colnames(cna_data)[tmp_BFB_pos1]
# tmp_BFB_pos = t(sapply(tmp_BFB_pos, function(x)strsplit(x, '_')[[1]]))
# tmp_BFB_pos = as.data.frame(tmp_BFB_pos)
# tmp_BFB_pos$V1 = as.numeric(tmp_BFB_pos$V1)
# tmp_BFB_pos$V2 = as.numeric(tmp_BFB_pos$V2)
# tmp_BFB_pos$V3 = as.numeric(tmp_BFB_pos$V3)
# tmp_BFB_pos$CN1 = c1_data[tmp_BFB_pos1]
# tmp_BFB_pos$CN2 = c2_data[tmp_BFB_pos1]
# tmp_BFB_pos1v = tmp_BFB_pos%>%arrange(V1, V2)
#
clone_plist = list()
#clone_names = rownames(clone_cnv)
for(i in clone_names){
  print(i)
  tmp_score_value = unlist(clone_rescore[i,])
  # pg = plot_genome(
  #   chr_bg_color = list(bk=c(0, 1,2),
  #                       col=c('white', '#C9E2F4', '#DFEDD7'),
  #                       alpha=0.2,
  #                       value=tmp_score_value),
  #   show_x_axis=F
  #   #show_x_axis = ifelse(i==clone_names[length(clone_names)], T, F)
  # )
  pg = plot_genome(
    chr_bg_color = list(# bk=c(min(clone_rescore), max(clone_rescore)),
                        bk=c(0, 0.25, 0.5),
                        col=c('white','#D7E0FF', '#0060BF'),
                        alpha=0,
                        value=tmp_score_value),
    show_x_axis=F,
    vline=F
    #show_x_axis = ifelse(i==clone_names[length(clone_names)], T, F)
  )
  tmp_cna_data = as.data.frame(t(clone_cnv[i, ,drop=F]))
  colnames(tmp_cna_data) = 'CNV'
  tmp_cna_data[, c('chr', 'start', 'end')] = t(sapply(rownames(tmp_cna_data), function(x)strsplit(x,'_')[[1]]))
  tmp_cna_data$CNV = as.numeric(tmp_cna_data$CNV)
  tmp_cna_data$start = as.numeric(tmp_cna_data$start)
  tmp_cna_data$end = as.numeric(tmp_cna_data$end)
  a = pg+
    geom_segment(aes(x=start,y=CNV, xend=end,yend=CNV), data=tmp_cna_data, color='black',linewidth=0.2)+
    labs(y=i)+
    lims(y=c(0,7))
  # if(i %in%c('clone_2', 'clone_3')){
  #   a=a+geom_segment(aes(x=V2,y=CN1, xend=V3,yend=CN2), data=tmp_BFB_pos23, color='#F4BA19',linewidth=0.6,
  #                  arrow = arrow(length = unit(0.05, "inches"), ends='both'))
  # }
  # if(i %in%c('clone_4', 'clone_5')){
  #   a=a+geom_segment(aes(x=V2,y=CN1, xend=V3,yend=CN2), data=tmp_BFB_pos45, color='#F4BA19',linewidth=0.6,
  #                  arrow = arrow(length = unit(0.05, "inches"), ends='both'))
  # }
  # if(i %in%c('clone_1')){
  #   a=a+geom_segment(aes(x=V2,y=CN1, xend=V3,yend=CN2), data=tmp_BFB_pos1v, color='#F4BA19',linewidth=0.6,
  #                  arrow = arrow(length = unit(0.05, "inches"), ends='both'))
  # }

  clone_plist[[i]] = a+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

}

clone_pg = cowplot::plot_grid(plotlist = clone_plist, ncol=1, align = 'v')

p=cowplot::plot_grid(clone_tree,clone_pg, nrow=1,align='h', rel_widths = c(0.2,0.7))
#p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树热图.pdf',
       p, width=120, height = 60, units='mm', dpi = 600)


##### clone间PDE、chromothripsis和BFB的比较 ####
library(ggpubr)
CNA_mechnism$clone = sctc$map_obj$cell_map_pos[rownames(CNA_mechnism), 'clone']

CNA_mechnism_melt = reshape2::melt(CNA_mechnism, id='clone')
CNA_mechnism_melt %>%
  filter(grepl('clone', clone)) %>%
  filter(!grepl('WGD', value)) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x=clone, y=value))+
  geom_violin(aes(fill=clone))+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  facet_wrap(~variable, scales = 'free',nrow = 1)+
  stat_compare_means(comparisons = list(c('clone_1', 'clone_2'),
                                        c('clone_1', 'clone_3'),
                                        c('clone_1', 'clone_4'),
                                        c('clone_1', 'clone_5')), label = 'p.format')+
  scale_fill_manual(values = clone_col)+
  theme_classic()+
  theme(strip.background = element_blank())

clone_col = c('C1'='#3B7035','C2'='#DEC35D','C3'='#DD805A', 'C4'='#34B6B5', 'C5'='#3152A3')

a1 = CNA_mechnism_melt %>%
  filter(grepl('clone', clone)) %>%
  filter(!grepl('WGD', value)) %>%
  filter(variable%in%c('BFB')) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C4','C3','C2'))) %>%
  ggplot(aes(x=clone, y=value))+
  #geom_violin(aes(fill=clone), size=0.001)+
  geom_violin(aes(fill=clone),scale = 'width', size=0.001, width=0.8)+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.001)+
  stat_compare_means(comparisons = list(c('C1', 'C4'),
                                        #c('C1', 'C4'),
                                        #c('C1', 'C5'),
                                        c('C1', 'C2')), label = 'p.format', size=2)+
  scale_fill_manual(values = clone_col)+
  labs(title='BFB')+
  scale_y_continuous(expand = expansion(mult=0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a1
a2 = CNA_mechnism_melt %>%
  filter(grepl('clone', clone)) %>%
  filter(!grepl('WGD', value)) %>%
  filter(variable%in%c('chromothripsis')) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C4','C3','C2'))) %>%
  ggplot(aes(x=clone, y=value))+
  #geom_violin(aes(fill=clone), size=0.001,bw=0.0025)+
  geom_violin(aes(fill=clone),scale = 'width', size=0.001, width=0.8)+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.001)+
  stat_compare_means(comparisons = list(c('C1', 'C4'),
                                        #c('C1', 'C4'),
                                        #c('C1', 'C5'),
                                        c('C1', 'C2')), label = 'p.format', size=2)+
  scale_fill_manual(values = clone_col)+
  labs(title='Chromothripsis')+
  scale_y_continuous(expand = expansion(mult=0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a2

a3 = CNA_mechnism_melt %>%
  filter(grepl('clone', clone)) %>%
  filter(!grepl('WGD', value)) %>%
  filter(variable%in%c('pde')) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C4','C3','C2'))) %>%
  ggplot(aes(x=clone, y=value))+
  #geom_violin(aes(fill=clone), size=0.001,bw=0.05)+
  geom_violin(aes(fill=clone),scale = 'width', size=0.001, width=0.8)+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.001)+
  stat_compare_means(comparisons = list(c('C1', 'C4'),
                                        #c('C1', 'C4'),
                                        #c('C1', 'C5'),
                                        c('C1', 'C2')), label = 'p.format', size=2)+
  scale_fill_manual(values = clone_col)+
  labs(title='CEE')+
  scale_y_continuous(expand = expansion(mult=0.1))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a3

a = cowplot::plot_grid(a1,a2,a3, nrow=1)
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_指标比较.pdf',
       a, width=80, height = 50, units='mm', dpi = 600)


ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_指标比较2a.pdf',
       a2, width=36, height = 44, units='mm', dpi = 600)
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_指标比较2b.pdf',
       a3, width=36, height = 44, units='mm', dpi = 600)

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=CPS_num))+
  geom_boxplot()+
  geom_jitter()

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=limit_num))+
  geom_boxplot()+
  geom_jitter()

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=seismic_num))+
  geom_boxplot()+
  geom_jitter()

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=limit_score))+
  geom_boxplot()+
  geom_jitter()

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=seismic_score))+
  geom_boxplot()+
  geom_jitter()

CNA_mechnism %>%
  filter(grepl('clone', clone)) %>%
  ggplot(aes(x=clone, y=pde))+
  geom_boxplot()+
  geom_jitter()

#### 全部chromothripsis中sesimic所占的比例 #####
# pie_1 = colSums(CNA_mechnism[rownames(sctc$map_obj$cell_map_pos), c('limit_num', 'seismic_num')])

pie_1 = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos),] %>%
  group_by(clone) %>%
  summarise(limit_num_clone = sum(limit_num), seismic_num_clone=sum(seismic_num)) %>%
  filter(grepl('clone', clone)) %>%
  as.data.frame() %>%
  reshape2::melt(id='clone')
colnames(pie_1) = c('clone', 'type', 'count')

a = pie_1 %>%
  mutate(type=factor(type, c( 'seismic_num_clone','limit_num_clone')))%>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C2','C4','C3'))) %>%
  ggplot(aes(x=clone, y=count, fill=type))+
  geom_bar(stat='identity', position = 'fill', color='black', size=0.01)+
  scale_fill_manual(values = c('limit_num_clone'='#DFEDD7','seismic_num_clone'='#C9E2F4'))+
  labs(y='chrom%')+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_比较chr比例.pdf',
       a, width=50, height = 40, units='mm', dpi = 600)
### 覆盖比例 ###
event_coverage = c()
for(i in rownames(sctc$map_obj$cell_map_pos)){
  limit_cover = 0
  for(reg in all_region_limit[[i]]){
    limit_cover = limit_cover + sum(sapply(reg, function(x)x[2]-x[1]))
  }
  event_coverage = rbind(event_coverage, c(i, 'gradual', limit_cover/ncol(sctc$orig.data$all_node_data)))
  seismic_cover = 0
  for(reg in all_region_seismic[[i]]){
    seismic_cover = seismic_cover + sum(sapply(reg, function(x)x[2]-x[1]))
  }
  event_coverage = rbind(event_coverage, c(i, 'seismic', seismic_cover/ncol(sctc$orig.data$all_node_data)))
}
event_coverage = as.data.frame(event_coverage)
colnames(event_coverage) = c('bc', 'type', 'coverage')
event_coverage$coverage = as.numeric(event_coverage$coverage)
event_coverage$clone = CNA_mechnism[event_coverage$bc, 'clone']
a = event_coverage %>%
  group_by(clone, type) %>%
  summarise(coverage = mean(coverage)) %>%
  filter(grepl('clone', clone)) %>%
  as.data.frame() %>%
  mutate(type=factor(type, c( 'gradual','seismic')))%>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C2','C4','C3'))) %>%
  ggplot(aes(x=clone, y=coverage, fill=type))+
  geom_bar(stat='identity', position = 'dodge', color='black', size=0.01)+
  scale_fill_manual(values = c('gradual'='#DFEDD7','seismic'='#C9E2F4'))+
  labs(y='coverage')+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改lone_修改算法/clone_比较event覆盖度.pdf',
       a, width=50, height = 40, units='mm', dpi = 600)
##  检验 添加显著性
df.summary <- event_coverage %>%
  group_by(clone, type) %>%
  summarise(
    sd = sd(coverage, na.rm = TRUE),
    coverage = mean(coverage)
  ) %>%
  filter(grepl('clone', clone)) %>%
  as.data.frame() %>%
  mutate(type=factor(type, c( 'gradual','seismic')))%>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C2','C4','C3')))
df2 = event_coverage %>%
  filter(grepl('clone', clone)) %>%
  as.data.frame() %>%
  #mutate(type=factor(type, c( 'gradual','seismic')))%>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C2','C4','C3')))%>%
  filter(type=='gradual')

a=ggplot(data=df2, mapping = aes(x=clone,y=coverage, fill=type))+
  #geom_boxplot()+
  stat_compare_means(comparisons = list(c("C1", "C2")), size=3)+
  geom_pointrange(data=df.summary, mapping = aes(x=clone,y=coverage, ymin = coverage-sd, ymax = coverage+sd, color=type),
                  position = position_dodge(0.5), size=0.5)+

  scale_color_manual(values = c('gradual'='#DFEDD7','seismic'='#C9E2F4'))+
  labs(y='coverage')+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_比较event覆盖度2.pdf',
       a, width=50, height = 40, units='mm', dpi = 600)


CNA_mechnism[rownames(sctc$map_obj$cell_map_pos),c('limit_num', 'seismic_num', 'clone')] %>%
  filter(grepl('clone', clone)) %>%
  reshape2::melt(id='clone') %>%
  ggplot(aes(x=clone, y=value,fill=variable))+
  geom_boxplot()



library(ggrepel)
my_pie_plot <- function(x, layer, value,xwd=1){
  p = ggplot()
  p_num = 1
  x$my_id = 1:nrow(x)
  x_pos = xwd*(1:length(layers))
  for(l in layers){
    tmp_p = ggplot()+geom_bar(data=x,
                              mapping = aes_string(x=0, y=value, fill='my_id'),
                              stat='identity', position = 'fill')
    tmp_pied = layer_data(tmp_p)[, c('ymin','ymax')]
    tmp_pied[, 'tmp_group'] = apply(x[, layers[1:p_num], drop=F], 1, function(yy)paste0(yy,collapse = ';;;'))
    tmp_pied = tmp_pied %>%
      group_by(tmp_group) %>%
      summarise(ymin=min(ymin), ymax=max(ymax))
    tmp_pied[, l] = sapply(tmp_pied$tmp_group, function(yy)strsplit(yy, ';;;')[[1]][p_num])
    tmp_pied$label = round(tmp_pied$ymax - tmp_pied$ymin,2)
    tmp_pied$label_x = 0.5+x_pos[p_num]
    tmp_pied$label_y = (tmp_pied$ymax+tmp_pied$ymin)/2

    p = p+
      geom_rect(data=tmp_pied,
                mapping = aes_string(fill=l, ymax='ymax', ymin='ymin', xmax=1+x_pos[p_num], xmin=0+x_pos[p_num]), inherit.aes = F)+
      geom_text(data=tmp_pied,mapping = aes(label=label,  y=label_y, x=label_x), inherit.aes = F)
    p_num = p_num+1
  }

  p=p+
    theme(aspect.ratio=1)+
    coord_polar(theta="y")+
    scale_fill_manual(values = c(clone_col, 'limit_num_clone'='#DFEDD7','seismic_num_clone'='#C9E2F4'))+
    theme_void()
  return(p)
}

my_pie_plot(pie_1, c('type', 'clone'), 'count')


### 对clone添加SV信息########
change_chr_name <-function(chr_name){
  chr_name = gsub('23','X', chr_name)
  chr_name = gsub('24','Y', chr_name)
  chr_name = gsub('chr0','', chr_name)
  chr_name = gsub('chr','', chr_name)

  new_name = c()
  for(i in chr_name){
    if(nchar(i)==1){
      if(i%in%c('X', 'Y')){
        new_name = c(new_name, paste0('chr',i))
      }else{
        new_name = c(new_name, paste0('chr0',i))
      }

    }else{
      new_name = c(new_name, paste0('chr',i))
    }
  }
  return(new_name)
}
# 准备染色体信息 #
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart

chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
tmp_fun = function(x,y){
  (min(x)+max(y))/2
}
chr_center_pos = chrinfo %>%
  group_by(chromNum) %>%
  summarise(center = tmp_fun(absStart, absEnd))

chrinfo_other = list(
  chr_center_pos = as.numeric(chr_center_pos$center), # 染色体中心坐标
  chr_name = as.vector(chr_center_pos$chromNum), # 染色体名字
  chr_band = chrinfo[chrinfo$chromStart==0, 'absStart'] # 染色体边界线
)
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
##
all_amplicons = c()
for(i in clone_names){
  print(i)
  amplicon_file = list.files(glue('/Volumes/WX_extend/BioSoftware/AmpliconArchitect/huh7_{i}_AA_results'), full.names = T)
  amplicon_file = grep('_edges.txt',amplicon_file, value = T)
  for(af in amplicon_file){
    amplicon_graph = tryCatch({read.table(af, sep='\t')},error = function(cond){return(NULL)})
    if(is.null(amplicon_graph)){
      next
    }
    amplicon_edges = as.data.frame(do.call(rbind,sapply(amplicon_graph$V1, function(x)strsplit(x, '->|:'))))
    colnames(amplicon_edges) = c('from_chr', 'from_pos', 'to_chr', 'to_pos')
    amplicon_edges$from_strand = sapply(amplicon_edges$from_pos, function(x)substr(x,nchar(x), nchar(x)))
    amplicon_edges$to_strand = sapply(amplicon_edges$to_pos, function(x)substr(x,nchar(x), nchar(x)))
    amplicon_edges$from_pos = as.numeric(sapply(amplicon_edges$from_pos, function(x)substr(x,1, nchar(x)-1)))
    amplicon_edges$to_pos = as.numeric(sapply(amplicon_edges$to_pos, function(x)substr(x,1, nchar(x)-1)))
    amplicon_edges$from_chr = change_chr_name(amplicon_edges$from_chr)
    amplicon_edges$to_chr = change_chr_name(amplicon_edges$to_chr)

    amplicon_edges$from_pos_abs = chrinfo_chr_absstart[amplicon_edges$from_chr, 'absStart'] + amplicon_edges$from_pos
    amplicon_edges$to_pos_abs = chrinfo_chr_absstart[amplicon_edges$to_chr, 'absStart'] + amplicon_edges$to_pos
    amplicon_edges[amplicon_edges$from_pos_abs==amplicon_edges$to_pos_abs, 'to_pos_abs'] = amplicon_edges[amplicon_edges$from_pos_abs==amplicon_edges$to_pos_abs, 'to_pos_abs']+1
    amplicon_edges$amplicon_id = grep('amplicon', strsplit(af, '_')[[1]], value = T)
    amplicon_edges$clone = i
    all_amplicons = rbind(all_amplicons, amplicon_edges)
  }
}

clone_plist2 = list()
for(tmp_clone in clone_names){
  tmp_amplicon = all_amplicons[all_amplicons$clone==tmp_clone,]
  tmp_amplicon$yend = ifelse((tmp_amplicon$to_pos-tmp_amplicon$to_pos)<1*1000*1000, 1, 0)
  tmp_amplicon$yend = ifelse(tmp_amplicon$from_chr!=tmp_amplicon$to_chr, 0, 1)

  tmp_pg = plot_genome(show_x_axis = F,
                       #special_region = region,
                       chr_bg_color = 'white')+
    geom_curve(mapping = aes(x = from_pos_abs, y = 0,
                             xend = to_pos_abs, yend = yend
    ),
    data = tmp_amplicon[tmp_amplicon$yend==0, , drop=F],
    color='red',
    curvature=-0.05,
    linewidth=0.2)+
    geom_curve(mapping = aes(x = from_pos_abs, y = 0,
                             xend = to_pos_abs, yend = yend
    ),
    data = tmp_amplicon[tmp_amplicon$yend==1, , drop=F],
    curvature=0,
    color='red',
    alpha = 0.1,
    linewidth=0.2)+
    lims(y=c(-1,1))+
    geom_hline(yintercept = 0, linetype='dashed', linewidth=0.2)+
    labs(y='SV')+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  clone_plist2[[tmp_clone]] = cowplot::plot_grid(tmp_pg,clone_plist[[tmp_clone]],ncol=1, align = 'v', rel_heights = c(0.3,0.7))
}

cowplot::plot_grid(plotlist = clone_plist2, ncol=1)



xstart = apply(all_amplicons[,c('from_pos_abs', 'to_pos_abs')], 1, min)
xend = apply(all_amplicons[,c('from_pos_abs', 'to_pos_abs')], 1, max)

tmp_chr = 'chr07'
tmp_clone = 'clone_3'
tmp_amplicon = all_amplicons[all_amplicons$clone==tmp_clone&all_amplicons$from_chr==tmp_chr,]

tmp_amplicon = all_amplicons[all_amplicons$clone==tmp_clone,]

tmp_amplicon = tmp_amplicon[1:3,]
tmp_chr_start = min(tmp_amplicon[, c('from_pos', 'to_pos')])
tmp_chr_end = max(tmp_amplicon[, c('from_pos', 'to_pos')])

region = paste0(tmp_chr, ':', tmp_chr_start, '-', tmp_chr_end)
plot_genome(show_x_axis = F,
            #special_region = region,
            chr_bg_color = 'white')+
  geom_curve(mapping = aes(x = tmp_amplicon$from_pos_abs, y = 0,
                           xend = tmp_amplicon$to_pos_abs, yend = 1),
             linewidth=0.2,
             alpha=0.5,
             curvature = -0)+
  labs(y='SC')


##### 单个细胞模式展示 ####
rescore = CNA_mechnism$orig.rearrange_score$all_rearrange_score
rescore2 = apply(rescore, 1, function(y){
  aa = sapply(y, function(x)t.test(y, mu = x, alternative='less')$p.value)
  p.adjust(aa)
})
score_vec = as.vector(as.matrix(rescore))
t_pvalue = sapply(score_vec, function(x)t.test(score_vec, mu = x, alternative='less')$p.value)
t_pvalue_fdr = p.adjust(t_pvalue)
#t_pvalue_mat = matrix(t_pvalue, nrow = nrow(limit1))
t_pvalue_fdr_mat = matrix(t_pvalue_fdr, nrow = nrow(rescore))
rescore[rescore2>0.001]=0

head(CNA_mechnism$orig.BFB_res)
tmp_x = sctc$orig.data$all_node_data[c('root', 'virtual_66', 'virtual_246'), ]
pdata=tmp_x[1,]
c1data=tmp_x[2,]
c2data=tmp_x[3,]

sum((((c1data+c2data)%%2) ==0) & c1data!=c2data)

x[x$label=='f1_TTTGGTTGTTGCTGAT-1',]
x[x$parent==637, ]
#rescore
bc = 'f2_TTGGCAAGTTCCATGA-1'
pg = plot_genome(
  chr_bg_color = list(bk=c(min(rescore), max(rescore)),
                      col=c('#52ACFF', '#FFE32C'),
                      alpha=0,
                      value=rescore[bc,])
)
pg
bc = 'f1_TTCTCCTAGCTAAACA-1'
prefix = 'huh7'
fn = strsplit(bc, '_')[[1]][1]
cl_path = '/Volumes/WX_extend/细胞轨迹推断/data/DNA_F1_f2/'
intcnv_path = glue('{cl_path}/{prefix}/{fn}/intCNV')
cellSeg_path = glue('{cl_path}/{prefix}/{fn}/cellSeg')
file_name = paste0(strsplit(bc, '_')[[1]][2],'.bam.hg19.50k.k50.varbin.data.txt')
file_data = read.table(glue('{intcnv_path}/{file_name}'), header = T, sep='\t')
file_data = na.omit(file_data)

seg_data = read.table(glue('{cellSeg_path}/{file_name}'), header = T, sep='\t')
seg_data = na.omit(seg_data)

cna_data = file_data[, c('chrom', 'chrompos', 'end', 'integerCNV', 'ratio')]
cna_data$bin_len = cna_data$end - cna_data$chrompos
cna_data$absEnd = cumsum(as.numeric(cna_data$bin_len))
cna_data$absStart = cna_data$absEnd-cna_data$bin_len

a =pg+
  geom_segment(aes(x=absStart,y=ratio, xend=absEnd,yend=ratio), data=cna_data, color='black',linewidth=2)+
  #geom_point(aes(x=absStart,y=ratio), data=cna_data, color='black',size=0.01)+
  geom_segment(aes(x=abspos,y=ratio, xend=absend,yend=ratio), data=seg_data, color='black',linewidth=1)+
  lims(y=c(0,5))
a

b = plot_genome(show_x_axis = F)+
  geom_segment(aes(x=absStart,y=integerCNV, xend=absEnd,yend=integerCNV), data=cna_data, color='red',linewidth=1)
b
limit_score = CNA_mechnism$orig.rearrange_score$all_limit_prop
c = plot_genome(
  chr_bg_color = list(bk=c(min(limit_score), max(limit_score)),
                      col=c('#52ACFF', '#FFE32C'),
                      alpha=0,
                      value=limit_score[bc,])
)

cowplot::plot_grid(b,a,c, ncol=1, align = 'v')



####
# Huh7细胞系热图
all_cell_cnv = sctc$orig.data$all_node_data
#Fn = sapply(rownames(sctc$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
#Fn = toupper(Fn)
#cna_meta = data.frame('Fn'=Fn, 'cell_name'=rownames(sctc$map_obj$cell_map_pos))
cna_meta = data.frame('clone'=sctc$map_obj$cell_map_pos$clone, 'cell_name'=rownames(sctc$map_obj$cell_map_pos))
rownames(cna_meta) = cna_meta$cell_name

col_meta = as.data.frame(t(sapply(colnames(all_cell_cnv), function(x)strsplit(x,'_')[[1]])))
col_meta$V1 = as.numeric(col_meta$V1)
col_meta$V2 = as.numeric(col_meta$V2)
col_meta = col_meta %>% arrange(V1, V2)

col_chr =col_meta$V1


all_cell_cnv[is.na(all_cell_cnv)] = 2
library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white",col="black",lwd=0.01),
                                               labels = 1:length(unique(col_chr)),
                                               labels_gp = gpar(cex = 0.6)))
cna_meta = cna_meta[, 'clone', drop=F]
cna_meta = cna_meta %>% arrange(clone)
cna_meta = cna_meta[grepl('clone', cna_meta$clone), , drop=F]

dd = sctc$clone_graph_obj$cell_phylo_tree$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
cna_meta = cna_meta[lab,,drop=F]
# Fn_col = c('F1'='#6300FF', 'F2'='#C75127FF')
clone_col = setNames(pal_npg()(length(unique(cna_meta$clone))),unique(cna_meta$clone))
clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')


left_anno = rowAnnotation(df=cna_meta,
                          col=list('clone' = clone_col),
                          width=unit(0.2, "mm"))
draw(left_anno)
cn_col = structure(c('blue','#91FF88', '#C6C6C6',  '#FFEB97',  '#FCCC00', '#ec9336', '#7d170e', 'darkred'),
                   names=c(0, 1,2,3,4,5,6, 7))
cn_col = structure(c('#2E6FAB','#6EA647', '#F5F5F5',  '#F8E596',  '#F4BA19', '#E37933', '#EC6A1A', '#B41D23'),
                   names=c(0, 1,2,3,4,5,6, 7))
width_in =  100/ 25.4
height_in = 70 / 25.4
dev.off()
pdf('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_DNA_cna热图.pdf', width=width_in, height=height_in)

ht = Heatmap(all_cell_cnv[rownames(cna_meta),rownames(col_meta)],
             cluster_rows = FALSE,cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = F,
             top_annotation = top_anno,
             left_annotation = left_anno,
             row_split = cna_meta$Fn,
             col = cn_col,
             # rect_gp = gpar(col = "white"),
             column_split = col_chr,
             column_gap = unit(0.1,  "mm"),
             heatmap_legend_param = list(
               title = "CNV",
               color_bar = "discrete",
               #at = c(0,1,2),
               direction = "horizontal",
               title_position = "leftcenter-rot"
               #legend_height = unit(3, "cm")
             ),
             row_title = NULL,
             column_title = NULL,
             use_raster=T
)
draw(ht)
dev.off()



### Huh7细胞精度进化树 （Fn着色） ####
obj =sctc
obj$map_obj$cell_map_pos$Fn = sapply(obj$map_obj$cell_map_pos$cell_name, function(x)strsplit(x, '_')[[1]][1])
ape_tree = obj$orig.data$tree
cna_data = obj$orig.data$all_node_data

colorby = 'Fn'
cell_info = obj$map_obj$cell_map_pos[, colorby, drop=F]
gtree_res = ggtree(ape_tree, size=0.1)
gtree_res$data$group = cell_info[gtree_res$data$label, colorby]

dd <- gtree_res$data
dd <- dd[dd$isTip,]
lab <- dd$label[order(dd$y, decreasing = T)]
## make sure the order of the input matrix is consistent with the tree
leaves_data = cna_data[lab,]

p = gtree_res+
  geom_tippoint(aes(color=group), size=0.2)+
  ggsci::scale_color_igv()+
  theme_tree(legend.position='right')+
  coord_flip()+
  scale_x_reverse()+
  guides(colour = guide_legend(override.aes = list(size=2))) # legend大小
p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Huh7_进化树Fn着色.pdf',
       p, width=140, height = 80, units='mm', dpi = 600)




### 细胞水平persist CNV与Denovo CNV相关性 （进化多样性velocity着色）####
meta = sctc$map_obj$cell_map_pos
meta2 = sctc$orig.data$cell_relation
denovo_cnv = meta2$Parent_gain_loc + meta2$Parent_loss_loc

plot_df = data.frame(denovo_cnv =denovo_cnv, #log1p(meta2$Parent_gain_loc),#denovo_cnv,
                     #persist_cnv = (meta2$Pseudotime_tree-min(meta2$Pseudotime_tree)) / (max(meta2$Pseudotime_tree)-min(meta2$Pseudotime_tree)),#log1p(meta2$Parent_loss_loc),#ncol(obj$orig.data$all_node_data)-denovo_cnv,
                     pde = meta[rownames(meta2), 'pde'],
                     row.names = rownames(meta2)
)
plot_df = na.omit(plot_df)

tree_data = sctc$clone_graph_obj$cell_phylo_tree$data %>% as.data.frame()

parent_num = tree_data[match(rownames(plot_df), tree_data$label), 'parent']
parent_label = tree_data[match(parent_num, tree_data$node), 'label']
persist_cnv = meta2[parent_label, 'Parent_gain_loc'] + meta2[parent_label, 'Parent_loss_loc']
plot_df$persist_cnv=persist_cnv/ncol(sctc$orig.data$all_node_data)
plot_df$denovo_cnv=plot_df$denovo_cnv/ncol(sctc$orig.data$all_node_data)

a = na.omit(plot_df) %>%
  ggplot(aes(x=persist_cnv, y=denovo_cnv, color=pde))+
  geom_point(size=0.6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='Persist_cnv', y='Denovo_cnv')+
  scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Denovo_persist相关pde着色.pdf',
       a, width=60, height = 50, units='mm', dpi = 600)


### 细胞水平velocity 与基因组变异比例相关性 （进化多样性velocity着色）####
meta = sctc$map_obj$cell_map_pos
meta2 = sctc$orig.data$cell_relation
genome_cnv = meta2$Parent_gain_cn + meta2$Parent_loss_cn
plot_df = data.frame(genome_cnv = genome_cnv/ncol(sctc$orig.data$all_node_data), #log1p(meta2$Parent_gain_loc),#denovo_cnv,
                     pde = meta[rownames(meta2), 'pde']
)
a = na.omit(plot_df) %>%
  ggplot(aes(x=pde, y=genome_cnv, color=pde))+
  geom_point(size=0.6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE', y='%Genome_altered(Parrent)')+
  scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_genome相关.pdf',
       a, width=60, height = 50, units='mm', dpi = 600)

########### PDE 与不同指标相关性 ########################
clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')
meta = sctc$map_obj$cell_map_pos[, c('chromothripsis','limit_num','seismic_num', 'BFB', 'pde', 'clone')]
plot_data = reshape2::melt(meta, id=c('pde', 'clone'))
a=plot_data %>%
  ggplot(aes(x=pde, y=value))+
  geom_point(aes(color=clone), size=1)+
  stat_cor()+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
  scale_color_manual(values = clone_col)+
  facet_wrap(~variable, scales = 'free_y')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
a

event_coverage = c()
for(i in rownames(sctc$map_obj$cell_map_pos)){
  limit_cover = 0
  for(reg in all_region_limit[[i]]){
    limit_cover = limit_cover + sum(sapply(reg, function(x)x[2]-x[1]))
  }
  event_coverage = rbind(event_coverage, c(i, 'gradual', limit_cover/ncol(sctc$orig.data$all_node_data)))
  seismic_cover = 0
  for(reg in all_region_seismic[[i]]){
    seismic_cover = seismic_cover + sum(sapply(reg, function(x)x[2]-x[1]))
  }
  event_coverage = rbind(event_coverage, c(i, 'seismic', seismic_cover/ncol(sctc$orig.data$all_node_data)))
}
event_coverage = as.data.frame(event_coverage)
colnames(event_coverage) = c('bc', 'type', 'coverage')
event_coverage$coverage = as.numeric(event_coverage$coverage)
event_coverage$clone = CNA_mechnism[event_coverage$bc, 'clone']

event_coverage$pde = sctc$map_obj$cell_map_pos[event_coverage$bc, 'pde']

a=event_coverage %>%
  ggplot(aes(x=pde, y=coverage))+
  geom_point(aes(color=clone), size=1)+
  stat_cor()+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
  scale_color_manual(values = clone_col)+
  facet_wrap(~type, scales = 'free_y')+
  scale_y_continuous(limits = c(0,0.8))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event相关.pdf',
       a, width=90, height = 50, units='mm', dpi = 600)

library(LSD)
event_coverage = na.omit(event_coverage)
x=event_coverage[event_coverage$type=='gradual',]
heatscatter(x$pde, x$coverage)

y=event_coverage[event_coverage$type=='seismic',]
heatscatter(y$pde, y$coverage)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- base::findInterval(x, dens$x)
  iy <- base::findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
x$den = get_density(x$pde, x$coverage, n = 200)
y$den = get_density(y$pde, y$coverage, n = 200)

new_event_coverage = rbind(x,y)
new_event_coverage$den = (new_event_coverage$den - min(new_event_coverage$den)) / (max(new_event_coverage$den)-min(new_event_coverage$den))
a=new_event_coverage %>%
  ggplot(aes(x=pde, y=coverage))+
  geom_point(aes(color=den), size=1)+
  stat_cor()+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
  scale_color_gradientn(colours = c('gray',"red"), limits=c(0,1))+
  #scale_color_manual(values = clone_col)+
  facet_wrap(~type, scales = 'free_y')+
  scale_y_continuous(limits = c(0,0.8))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event相关2.pdf',
       a, width=90, height = 50, units='mm', dpi = 600)

a=new_event_coverage %>%
  ggplot(aes(x=pde, y=coverage))+
  geom_point(aes(color=den), size=1)+
  stat_cor()+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
  scale_color_gradientn(colours = c('gray',"red"), limits=c(0,1))+
  #scale_color_manual(values = clone_col)+
  facet_wrap(~type, scales = 'free_y')+
  scale_y_continuous(limits = c(0,0.8))+
  theme_classic()+
  theme(strip.background = element_blank(),
        #legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event相关2_1.pdf',
       a, width=90, height = 50, units='mm', dpi = 600)

###########不同指标与基因组变异比例相关性 ########################
meta = sctc$map_obj$cell_map_pos[, c('rearrange_score','limit_score','seismic_score', 'BFB', 'pde', 'clone')]
meta2 = sctc$orig.data$cell_relation[rownames(meta),]

genome_cnv = rowSums((sctc$orig.data$all_node_data-2)!=0)[rownames(meta)] / ncol(sctc$orig.data$all_node_data)
# 2
genome_cnv = (meta2$Parent_gain_loc + meta2$Parent_loss_loc)/ ncol(sctc$orig.data$all_node_data)#(meta2$Root_gain_cn + meta2$Root_loss_loc) / (meta2$Root_gain_loc + meta2$Root_loss_loc)


meta$genome_cnas = genome_cnv
meta = meta[grepl('f', rownames(meta)), ]
library(ggpubr)
plot_data = reshape2::melt(meta, id=c('clone', 'genome_cnas'))
a=plot_data %>%
  ggplot(aes(x=value, y=genome_cnas))+
  geom_point(aes(color=clone), size=1)+
  stat_cor()+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed')+
  scale_color_manual(values = clone_col)+
  facet_wrap(~variable, scales = 'free_x')+
  labs(y='%genome_altered')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/genome与指标相关.pdf',
       a, width=180, height = 180, units='mm', dpi = 600)

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/genome与指标相关2.pdf',
       a, width=180, height = 180, units='mm', dpi = 600)


### 细胞水平细胞变异比例累计百分比，高低速率 shareed CNA 百分比 （进化多样性velocity着色）####
all_cell_cnv=sctc$orig.data$all_node_data
cell_genome_var_prop = rowSums((all_cell_cnv-2)!=0)/ncol(all_cell_cnv)

cell_genome_var_prop=cell_genome_var_prop[grepl('f', names(cell_genome_var_prop))]
plot_df = data.frame('genome_alt'=cell_genome_var_prop,
                     'pde' = pde[names(cell_genome_var_prop)])
plot_df = plot_df %>% arrange(genome_alt)
plot_df$pde_group = ifelse(plot_df$pde>median(plot_df$pde), 'High', 'low')

# plot_df %>% ggplot(aes(x=pde_group, y=genome_alt))+geom_boxplot()
# plot_df$genome_alt_cut = cut(plot_df$genome_alt, breaks = seq(0,1,length.out=10), labels = seq(0,1,length.out=10)[2:10])
#plot_df$genome_alt_cut = cut(1:nrow(plot_df), breaks = seq(0,nrow(plot_df),length.out=10),
#                             labels = seq(0,nrow(plot_df),length.out=10)[2:10])
#
#plot_df$bc = rownames(plot_df)
#
#ggtree_get_parent <- function(tree, node_labels){
#  tree_data = tree$data %>% as.data.frame()
#  p_nodes = tree_data[match(node_labels, tree_data$label), 'parent']
#  l_labels = tree_data[match(p_nodes, tree_data$node), 'label']
#  return(l_labels)
#}
#
# CNV_envents = list()
# leaves_cna_data = all_cell_cnv[plot_df$bc, ]
# leaves_parent_node = ggtree_get_parent(obj$clone_graph_obj$cell_phylo_tree, plot_df$bc)
# leaves_parent_data = all_cell_cnv[leaves_parent_node, ]
#CNV_envents = leaves_cna_data - leaves_parent_data
CNV_envents = all_cell_cnv
tmp_func_share<-function(x){
  tmp_data = CNV_envents[x,,drop=F]

  tmp_res = (rowMeans(t(tmp_data))-apply(t(tmp_data), 1, min))==0
  # tmp_res = apply(tmp_data, 2, function(y){
  #   return(ifelse(length(table(y))==1, 1, 0))
  # })
  return(sum(tmp_res)/ncol(tmp_data))
}
cps_events = (all_rearrange_score_pvalue<pvalue_thr & all_rearrange_score>0)+0
tmp_func_chromothripisis<-function(x){
  tmp_data = cps_events[x,,drop=F]
  tmp_data = sum(colSums(tmp_data)!=0) /ncol(cps_events)
  return(tmp_data)
}

BFB_pos_res = BFB_pos(sctc)
rownames(BFB_pos_res) = BFB_pos_res[,1]
BFB_pos_res = BFB_pos_res[,-1]
tmp_func_BFB<-function(x){
  tmp_data = BFB_pos_res[x,,drop=F] == 'TRUE'
  tmp_data = sum(colSums(tmp_data)!=0) /ncol(tmp_data)
  return(tmp_data)
}

plot_df_high = plot_df[plot_df$pde_group=='High',]
high_rate = c()
for(i in 1:nrow(plot_df_high)){
  print(i)
  bc = rownames(plot_df_high)[1:i]
  high_rate = rbind(high_rate, c(tmp_func_share(bc), tmp_func_chromothripisis(bc), tmp_func_BFB(bc)))
}
high_rate = as.data.frame(high_rate)
colnames(high_rate) = c('share_CNA', 'chromothripisis', 'BFB')
high_rate$group='High'
high_rate$BFB = high_rate$BFB / ncol(CNV_envents)
high_rate$cell_num = (1:nrow(high_rate))/nrow(high_rate)

plot_df_low = plot_df[plot_df$pde_group=='low',]
low_rate = c()
for(i in 1:nrow(plot_df_low)){
  print(i)
  bc = rownames(plot_df_low)[1:i]
  low_rate = rbind(low_rate, c(tmp_func_share(bc), tmp_func_chromothripisis(bc), tmp_func_BFB(bc)))
}
low_rate = as.data.frame(low_rate)
colnames(low_rate) = c('share_CNA', 'chromothripisis', 'BFB')
low_rate$BFB = low_rate$BFB / ncol(CNV_envents)
low_rate$group='low'
low_rate$cell_num = (1:nrow(low_rate))/nrow(low_rate)


plot_df2 = rbind(high_rate, low_rate)
saveRDS(plot_df2, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/high_low_rate.rds')
plot_df2 = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/high_low_rate.rds')
plot_df2 = reshape2::melt(plot_df2, id=c('group', 'cell_num'))

a=plot_df2 %>%
  ggplot(aes(x=cell_num, y=value, color=group))+
  geom_line(linewidth=1.5)+
  # geom_point(size=0.5)+
  # geom_smooth(linewidth=1, se=F)+
  labs(x='Cumulative percentage of cells')+
  facet_wrap(~variable, scales = 'free_y')+
  labs(y='prop')+
  scale_color_npg()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/CellProp_shareCNA.pdf', a,
       width=160, height=60, units='mm', dpi = 600, bg = 'transparent')

a=plot_df2 %>%
  filter(variable=='share_CNA') %>%
  ggplot(aes(x=cell_num, y=value, color=group))+
  geom_line(linewidth=0.5)+
  # geom_point(size=0.5)+
  # geom_smooth(linewidth=1, se=F)+
  labs(x='Cumulative percentage of cells', y='Share CNA(%)')+
  scale_color_npg()+
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

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/CellProp_shareCNA1.pdf', a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

#
plot_df_high = plot_df[plot_df$pde_group=='High',]
plot_df_high$rearrange_score = CNA_mechnism[rownames(plot_df_high), 'rearrange_score']
plot_df_high = plot_df_high %>% arrange(rearrange_score)
plot_df_high$cell_num = 1:nrow(plot_df_high)
plot_df_low = plot_df[plot_df$pde_group=='low',]
plot_df_low$rearrange_score = CNA_mechnism[rownames(plot_df_low), 'rearrange_score']
plot_df_low = plot_df_low %>% arrange(rearrange_score)
plot_df_low$cell_num = 1:nrow(plot_df_low)

rbind(plot_df_low,plot_df_high) %>%
  ggplot(aes(x=cell_num, y=rearrange_score, color=pde_group))+
  geom_line()


#### seismic limit 随着细胞比例变化 #####
# all_cell_cnv=sctc$orig.data$all_node_data
# cell_genome_var_prop = rowSums((all_cell_cnv-2)!=0)/ncol(all_cell_cnv)
#
# cell_genome_var_prop=cell_genome_var_prop[grepl('f', names(cell_genome_var_prop))]
# plot_df = data.frame('genome_alt'=cell_genome_var_prop,
#                      'pde' = pde[names(cell_genome_var_prop)])
# plot_df = plot_df %>% arrange(genome_alt)
# plot_df$pde_group = ifelse(plot_df$pde>median(plot_df$pde), 'High', 'low')
#
#
# tmp_func_seismic<-function(x){
#   return(sum(sctc$map_obj$cell_map_pos[x, 'seismic_score'], na.rm=T))
# }
# tmp_func_limit<-function(x){
#   return(sum(sctc$map_obj$cell_map_pos[x, 'limit_score'], na.rm=T))
# }
# tmp_func_BFB<-function(x){
#   return(sum(sctc$map_obj$cell_map_pos[x, 'BFB'], na.rm=T))
# }
# tmp_func_RE<-function(x){
#   return(sum(sctc$map_obj$cell_map_pos[x, 'rearrange_score'], na.rm=T))
# }
# plot_df_high = plot_df[plot_df$pde_group=='High',]
# high_rate = c()
# for(i in 1:nrow(plot_df_high)){
#   print(i)
#   bc = rownames(plot_df_high)[1:i]
#   high_rate = rbind(high_rate, c(tmp_func_seismic(bc), tmp_func_limit(bc),tmp_func_BFB(bc), tmp_func_RE(bc)))
# }
# high_rate = as.data.frame(high_rate)
# colnames(high_rate) = c('seismic', 'limit', 'BFB', 'rearrange_score')
# high_rate$group='High'
# high_rate$cell_num = (1:nrow(high_rate))/nrow(high_rate)
#
# plot_df_low = plot_df[plot_df$pde_group=='low',]
# low_rate = c()
# for(i in 1:nrow(plot_df_low)){
#   print(i)
#   bc = rownames(plot_df_low)[1:i]
#   low_rate = rbind(low_rate, c(tmp_func_seismic(bc), tmp_func_limit(bc),tmp_func_BFB(bc), tmp_func_RE(bc)))
# }
# low_rate = as.data.frame(low_rate)
# colnames(low_rate) = c('seismic', 'limit', 'BFB', 'rearrange_score')
# low_rate$group='low'
# low_rate$cell_num = (1:nrow(low_rate))/nrow(low_rate)
#
#
# plot_df2 = rbind(high_rate, low_rate)
# plot_df2 = reshape2::melt(plot_df2, id=c('group', 'cell_num'))
#
# a=plot_df2 %>%
#   #filter(variable=='limit')%>%
#   ggplot(aes(x=cell_num, y=value, color=group))+
#   geom_line()+
#   #geom_point(size=0.5)+
#   #geom_smooth(linewidth=1, se=F)+
#   facet_wrap(~variable, scales = 'free_y')+
#   labs(x='Cumulative percentage of cells')+
#   labs(y='Limit_score')+
#   scale_color_npg()+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank())
# a
#
# ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/CellProp_Limit_score.pdf', a,
#        width=80, height=50, units='mm', dpi = 600, bg = 'transparent')
#
#### 展示seismic 和limit #####
rescore = all_rearrange_score#CNA_mechnism$orig.rearrange_score$all_rearrange_score
chr=8
rescore = rescore[all_rearrange_score_pvalue[,chr]<0.001,]
rescore = rescore[rownames(all_seg_len[all_seg_len[,chr]>10,]),]
rescore = rescore[grepl('f', rownames(rescore)),]
#rescore = rescore[grepl('f', rownames(rescore)),]
pheatmap::pheatmap(rescore)
# 选择8号染色体
rescore = rescore[order(-rescore[,chr]),]
pheatmap::pheatmap(rescore, cluster_rows = F, cluster_cols = F)

#
bc = tail(rownames(rescore),1)[10]
bc = rownames(rescore)[1]
#'#C9E2F4', '#DFEDD7'
prefix = 'huh7'
fn = strsplit(bc,'_')[[1]][1]
cl_path = '/Volumes/WX_extend/细胞轨迹推断/data/DNA_F1_f2/'
intcnv_path = glue('{cl_path}/{prefix}/{fn}/intCNV')
cellSeg_path = glue('{cl_path}/{prefix}/{fn}/cellSeg')
file_name = paste0(strsplit(bc, '_')[[1]][2],'.bam.hg19.50k.k50.varbin.data.txt')
file_data = read.table(glue('{intcnv_path}/{file_name}'), header = T, sep='\t')
file_data = na.omit(file_data)

seg_data = read.table(glue('{cellSeg_path}/{file_name}'), header = T, sep='\t')
seg_data = na.omit(seg_data)

cna_data = file_data[, c('chrom', 'chrompos', 'end', 'integerCNV', 'ratio')]
cna_data$bin_len = cna_data$end - cna_data$chrompos
cna_data$absEnd = cumsum(as.numeric(cna_data$bin_len))
cna_data$absStart = cna_data$absEnd-cna_data$bin_len
# pg = plot_genome_chr(chr_num = paste0('chr0', chr),bg_color = '#C9E2F4')#, start = '59900000', end='127100000')
pg=ggplot()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill='#C9E2F4'),
  )
pg
tmp_cna_data = cna_data[cna_data$chrom==as.character(chr),]
tmp_seg_data = seg_data[seg_data$chr==as.character(chr),]
for(i in 1:nrow(tmp_seg_data)){
  tmp_start = tmp_seg_data[i, 'start']
  tmp_end = tmp_seg_data[i, 'end']

  tmp_ratio = subset(tmp_cna_data, chrompos>=tmp_start&end<tmp_end)
  new_ratio = mean(tmp_ratio$ratio)
  tmp_seg_data[i, 'ratio'] = new_ratio
}
a1 =pg+
  geom_segment(aes(x=chrompos/1000/1000,y=ratio+0.5, xend=end/1000/1000,yend=ratio+0.5), data=tmp_cna_data, color='black',linewidth=.3)+
  #geom_point(aes(x=absStart,y=ratio), data=cna_data, color='black',size=0.01)+
  geom_segment(aes(x=start/1000/1000,y=ratio+0.5, xend=end/1000/1000,yend=ratio+0.5), data=tmp_seg_data, color='black',linewidth=0.5)+
  scale_x_continuous(limits=c(0, max(tmp_seg_data$end)/1000/1000),
                     expand = expansion(mult = c(0,0)))+
  labs(y='ratio')+
  lims(y=c(0,5))+
  theme(#axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.15),
    axis.ticks.length=unit(0.025, "cm"),
    axis.text = element_text(size=6))

a1

b1 = plot_genome_chr(chr_num = paste0('chr0', chr),bg_color = 'white')+#, start = '59900000', end='127100000')+
  labs(y='CN')+
  lims(y=c(0,7))+
  geom_segment(aes(x=chrompos,y=integerCNV, xend=end,yend=integerCNV),
               data=cna_data[cna_data$chrom==as.character(chr),], color='red',linewidth=0.5)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6))


p_seismic = cowplot::plot_grid(b1,a1, ncol=1, align = 'v', rel_heights = c(0.3,0.6))
p_seismic


ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Seismic展示',bc,'chr',chr,'.pdf'),
       p_seismic,
       width=80, height=60, units='mm', dpi = 600, bg = 'transparent')


# limit
rescore = all_limit_prop#CNA_mechnism$orig.rearrange_score$all_limit_prop
chr=4
rescore = rescore[all_rearrange_score_pvalue[,chr]<0.001,]
rescore = rescore[rownames(all_seg_len[all_seg_len[,chr]>10,]),]
rescore = rescore[grepl('f', rownames(rescore)),]
#rescore = rescore[grepl('f', rownames(rescore)),]
pheatmap::pheatmap(rescore)
# 选择8号染色体
#rescore = rescore[order(-rescore[,12]),]
#pheatmap::pheatmap(rescore, cluster_rows = F, cluster_cols = F)

rescore = rescore[order(-rescore[,chr]),]
pheatmap::pheatmap(rescore, cluster_rows = F, cluster_cols = F)

#
bc = rownames(rescore)[28]
bc
#'#C9E2F4', '#DFEDD7'
prefix = 'huh7'
fn = strsplit(bc,'_')[[1]][1]
cl_path = '/Volumes/WX_extend/细胞轨迹推断/data/DNA_F1_f2/'
intcnv_path = glue('{cl_path}/{prefix}/{fn}/intCNV')
cellSeg_path = glue('{cl_path}/{prefix}/{fn}/cellSeg')
file_name = paste0(strsplit(bc, '_')[[1]][2],'.bam.hg19.50k.k50.varbin.data.txt')
file_data = read.table(glue('{intcnv_path}/{file_name}'), header = T, sep='\t')
file_data = na.omit(file_data)

seg_data = read.table(glue('{cellSeg_path}/{file_name}'), header = T, sep='\t')
seg_data = na.omit(seg_data)

cna_data = file_data[, c('chrom', 'chrompos', 'end', 'integerCNV', 'ratio')]
cna_data$bin_len = cna_data$end - cna_data$chrompos
cna_data$absEnd = cumsum(as.numeric(cna_data$bin_len))
cna_data$absStart = cna_data$absEnd-cna_data$bin_len
#pg = plot_genome_chr2(chr_num = 'chr04',bg_color = '#DFEDD7')#, start = '123800000', end='170100000')

pg=ggplot()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill='#DFEDD7'))
# pg = plot_genome_chr2(chr_num = paste0('chr0', chr),bg_color = '#DFEDD7', start = 110*1000*1000)

# plot_genome_chr(chr_num = paste0('chr0', chr),bg_color = 'white', start = 110*1000*1000)
tmp_seg_data = seg_data[seg_data$chr==as.character(chr),]
tmp_cna_data = cna_data[cna_data$chrom==as.character(chr),]
#tmp_seg_data = rbind(tmp_seg_data, c(4, 123800000, 143019057, 0.8021812))
#tmp_seg_data = rbind(tmp_seg_data, c(4, 164026166, 170100000, 0.7909728))
for(i in 1:nrow(tmp_seg_data)){
  tmp_start = tmp_seg_data[i, 'start']
  tmp_end = tmp_seg_data[i, 'end']

  tmp_ratio = subset(tmp_cna_data, chrompos>=tmp_start&end<tmp_end)
  new_ratio = mean(tmp_ratio$ratio)
  tmp_seg_data[i, 'ratio'] = new_ratio
}
s=57558840
e=99254497
a =pg+
  geom_segment(aes(x=chrompos/1000/1000,y=ratio+0.5, xend=end/1000/1000,yend=ratio+0.5),
               data=cna_data[cna_data$chrom==as.character(chr),], color='black',linewidth=0.3)+
  #geom_point(aes(x=absStart,y=ratio), data=cna_data, color='black',size=0.01)+
  geom_segment(aes(x=start/1000/1000,y=ratio+0.5, xend=end/1000/1000,yend=ratio+0.5), data=tmp_seg_data, color='black',linewidth=0.5)+
  scale_x_continuous(#limits=c(0, max(tmp_seg_data$end)/1000/1000),
                     limits=c(s/1000/1000, e/1000/1000),
                     expand = expansion(mult = c(0,0)))+
  labs(y='ratio')+
  lims(y=c(0,5))+
  theme(#axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.15),
    axis.ticks.length=unit(0.025, "cm"),
    axis.text = element_text(size=6))

a

b = plot_genome_chr(chr_num = paste0('chr0', chr),bg_color = 'white', start = s, end=e)+
  #labs(y='CN')+
  lims(y=c(0,7))+
  geom_segment(aes(x=chrompos,y=integerCNV, xend=end,yend=integerCNV),
               data=cna_data[cna_data$chrom==as.character(chr),], color='red',linewidth=0.5)+
  scale_x_continuous(#limits=c(0, max(tmp_seg_data$end)/1000/1000),
    limits=c(s, e),
    expand = expansion(mult = c(0,0)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6))
#bb = cowplot::plot_grid(b1,b, nrow=1)
#aa = cowplot::plot_grid(a1,a, nrow=1)
#bb+aa
#cowplot::plot_grid(bb,aa,ncol=1, align = 'v')
p_gradual = cowplot::plot_grid(b,a, ncol=1, align = 'v', rel_heights = c(0.3,0.6))
p_gradual

# b+geom_vline(xintercept = 59853470, color='blue')

ggsave(paste0('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Seismic展示',bc,'chr',chr,'.pdf'),
       p_gradual,
       width=80, height=60, units='mm', dpi = 600, bg = 'transparent')



ab = cowplot::plot_grid(p_gradual,p_seismic, nrow=1, align = 'h')
ab
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/模式展示.pdf', ab,
       width=100, height=40, units='mm', dpi = 600, bg = 'transparent')

### cell-ITH 与 PDE相关性 ######
all_cell_cnv = obj$orig.data$all_node_data
cell_pos = tmp_class$map_obj$cell_map_pos
cell_ITH_mat = c()
for(clone in unique(cell_pos$clone)){
  if(grepl('clone', clone)){
    print(clone)
    tmp_clone_cells = rownames(cell_pos[cell_pos$clone==clone, ])

    for(ii in 1:(length(tmp_clone_cells)-1)){
      print(ii)
      i = tmp_clone_cells[ii]
      for(jj in (ii+1):length(tmp_clone_cells)){
        print(jj)
        j = tmp_clone_cells[jj]
        x = unlist(all_cell_cnv[i,,drop=T])
        y = unlist(all_cell_cnv[j,,drop=T])

        up_value = sum((x-mean(x))*(y-mean(y)))
        down_value = sqrt(sum((x-mean(x))^2)) * sqrt(sum((y-mean(y))^2))
        cell_ITH_mat = rbind(cell_ITH_mat, c(clone, i, j, 1-up_value/down_value))
      }
    }
  }
}

cell_ITH_mat = as.data.frame(cell_ITH_mat)
colnames(cell_ITH_mat) = c('clone','c1', 'c2', 'cell_ITH')
cell_ITH_mat$cell_ITH = as.numeric(cell_ITH_mat$cell_ITH)

saveRDS(cell_ITH_mat, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Huh7_clone_ITH.rds')
cell_ITH_mat = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/Huh7_clone_ITH.rds')
library(forcats)
cell_ITH_mat[cell_ITH_mat$clone=='clone_3', 'clone'] = 'clone_2'
cell_ITH_mat[cell_ITH_mat$clone=='clone_4', 'clone'] = 'clone_3'
cell_ITH_mat[cell_ITH_mat$clone=='clone_5', 'clone'] = 'clone_4'


clone_col = c('C1'='#3B7035','C2'='#DEC35D','C3'='#DD805A', 'C4'='#34B6B5', 'C5'='#3152A3')
a = cell_ITH_mat %>%
  mutate(clone=fct_reorder(clone, -cell_ITH, median)) %>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  mutate(clone=factor(clone, levels = c('C1', 'C4','C3','C2'))) %>%
  ggplot(aes(x=clone, y=cell_ITH))+
  geom_violin(aes(fill=clone), size=0.001, width=0.8)+
  geom_boxplot(width=0.1, outlier.shape = NA, size=0.001)+
  stat_compare_means(comparisons = list(c('C1', 'C4'),
                                        c('C1', 'C2')), label = 'p.format', size=2)+
  scale_fill_manual(values = clone_col)+
  labs(title='ITH')+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/Clone_ITH.pdf', a,
       width=36, height=44, units='mm', dpi = 600, bg = 'transparent')



clone_ITH = cell_ITH_mat %>%
  group_by(clone, c1) %>%
  summarise(cell_ITH=median(cell_ITH)) %>% as.data.frame()
clone_ITH$PDE = sctc$map_obj$cell_map_pos[clone_ITH$c1, 'pde']

a = clone_ITH %>%
  mutate(clone=gsub('lone_','', clone))%>%
  mutate(clone=toupper(clone))%>%
  #mutate(clone=factor(clone, levels = c('C1', 'C3','C5','C4','C2'))) %>%
  ggplot(aes(x=PDE, y=cell_ITH))+
  geom_point(aes(color=clone), size=1)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE', y='cell_ITH')+
  stat_cor(size=2)+
  scale_color_manual(values = clone_col)+
  # scale_fill_igv()+
  theme_classic()+
  #theme(panel.grid = element_blank())+
  #theme(legend.text = element_text(size=4),
  #      legend.key.size = unit(1, 'mm'),
  #      legend.title = element_text(size=6)
  #)+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))

a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_ITH_相关clone.pdf',
       a, width=40, height = 40, units='mm', dpi = 600)

### 高低PDE两组细胞之间 ######
library(ggpubr)
rownames(clone_ITH) = clone_ITH$c1
plot_df = CNA_mechnism[rownames(sctc$map_obj$cell_map_pos),]
plot_df = plot_df[grepl('^f', rownames(plot_df)),]
plot_df$pde_group = ifelse(plot_df$pde>median(plot_df$pde,na.rm = T), 'High', 'low')
plot_df$cell_ITH = clone_ITH[rownames(plot_df), 'cell_ITH']

a=plot_df %>%
  ggplot(aes(x=pde_group, y=cell_ITH))+
  geom_violin(aes(fill=pde_group), size=0.1)+
  geom_boxplot(width=0.1, size=0.1, outlier.shape = NA)+
  #geom_jitter()+
  stat_compare_means(comparisons = list(c('High', 'low')), label = 'p.signif')+
  labs(x='PDE_F2', y='Cell_ITH')+
  scale_fill_npg()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        axis.title.x = element_blank()
  )
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//PDE_Group_ITH.pdf',a,
       width=60, height=60, units='mm', dpi = 600, bg = 'transparent')

# 0
a = plot_df[,c('pde_group', 'limit_num', 'seismic_num')] %>%
  reshape2::melt(id='pde_group') %>%
  mutate(variable=factor(variable, levels = c('seismic_num','limit_num'))) %>%
  ggplot(aes(x=pde_group,  y=value, fill=variable))+
  geom_bar(stat='identity', position = 'fill', width=0.6)+
  scale_fill_manual(values = c('limit_num'='#DFEDD7','seismic_num'='#C9E2F4'))+

  # coord_polar(theta="y")+
  #facet_wrap(~variable, scales = 'free')+
  #scale_fill_igv()+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6),
        plot.title = element_text(hjust=0.5, size=8),
        axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text = element_text(size=6))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event_prop.pdf',a,
       width=40, height=40, units='mm', dpi = 600, bg = 'transparent')

# 1
a=plot_df[,c('pde_group', 'limit_num', 'seismic_num')] %>%
  reshape2::melt(id='pde_group') %>%
  mutate(value=cut(value, breaks = c(-Inf,5,10,15,20,Inf), labels=c('<5', '<10', '<15', '<20', '<24'))) %>%
  ggplot(aes(x=pde_group,  fill=value))+
  geom_bar(stat='count', position = 'fill', width=0.6)+
  scale_fill_manual(values = c('<5'='#D8E4F0','<10'='#B5C3DA','<15'='#8B9FC1', '<20'='#607FAA', '<24'='#2D6699'))+
  # coord_polar(theta="y")+
  facet_wrap(~variable, scales = 'free')+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event_prop2.pdf',a,
       width=80, height=46, units='mm', dpi = 600, bg = 'transparent')

# 1.1
a = plot_df[,c('pde_group', 'limit_num', 'seismic_num')] %>%
  reshape2::melt(id='pde_group') %>%
  mutate(variable=factor(variable, levels = c('seismic_num','limit_num'))) %>%
  mutate(value2=cut(value, breaks = c(-Inf,5,10,Inf), labels=c('>0', '>5', '>10'))) %>%
  group_by(pde_group, value2, variable) %>%
  summarise(value=sum(value)) %>%
  mutate(variable=factor(variable, levels = c('seismic_num','limit_num'))) %>%
  ggplot(aes(x=pde_group,  y=value, color=variable, fill=value2))+
  geom_bar(stat='identity', position = 'fill', width=0.6, size=1)+
  scale_fill_manual(values = c('>0'='#FFECB3FF','>5'='#FFD54FFF','>10'='#FFA000FF', '>15'='#FF6F00FF', '<24'='#2D6699'))+
  scale_color_manual(values = c('limit_num'='#DFEDD7','seismic_num'='#C9E2F4'))+
  # coord_polar(theta="y")+
  #facet_wrap(~variable, scales = 'free')+
  #scale_fill_igv()+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_event_prop3.pdf',a,
       width=70, height=40, units='mm', dpi = 600, bg = 'transparent')

# 2
plot_df[,c('pde_group', 'CPS_num', 'limit_num', 'seismic_num')] %>%
  reshape2::melt(id='pde_group') %>%
  mutate(value=cut(value, breaks = c(-Inf,5,10,15,20,Inf))) %>%
  ggplot(aes(x=value,  fill=pde_group))+
  geom_bar(stat='count', position = 'fill')+
  # coord_polar(theta="y")+
  facet_wrap(~variable, scales = 'free')+
  scale_fill_igv()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

# 3
plot_df[,c('pde_group', 'CPS_num', 'limit_num', 'seismic_num', 'clone')] %>%
  reshape2::melt(id=c('pde_group','clone')) %>%
  mutate(value=cut(value, breaks = c(-Inf,5,10,15,20,Inf))) %>%
  ggplot(aes(x=pde_group,  fill=clone))+
  geom_bar(stat='count', position = 'fill')+
  # coord_polar(theta="y")+
  facet_wrap(~variable, scales = 'free')+
  scale_fill_manual(values = clone_col)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# 4
plot_df[,c('pde_group', 'CPS_num', 'limit_num', 'seismic_num', 'clone')] %>%
  reshape2::melt(id=c('pde_group','clone')) %>%
  mutate(value=cut(value, breaks = c(-Inf,5,10,15,20,Inf))) %>%
  ggplot(aes(x=clone,  fill=pde_group))+
  geom_bar(stat='count', position = 'fill')+
  # coord_polar(theta="y")+
  facet_wrap(~variable, scales = 'free')+
  scale_fill_npg()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# 5
a=plot_df[,c('pde_group', 'CPS_num', 'limit_num', 'seismic_num', 'clone')] %>%
  reshape2::melt(id=c('pde_group','clone')) %>%
  mutate(value=cut(value, breaks = c(-Inf,5,10,15,20,Inf), label=c('<5', '<10', '<15', '<20', '<24'))) %>%
  filter(variable=='CPS_num') %>%
  ggplot(aes(x=clone,  fill=value))+
  geom_bar(stat='count', position = 'fill')+
  # coord_polar(theta="y")+
  facet_wrap(~variable, scales = 'free')+
  # scale_fill_igv()+
  scale_fill_manual(values = c('<5'='#D8E4F0','<10'='#B5C3DA','<15'='#8B9FC1', '<20'='#607FAA', '<24'='#2D6699'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
layer_data(a)
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_CPS_构成.pdf',
       a, width=110, height = 80, units='mm', dpi = 600)

# 6
plot_df[,c('pde_group', 'limit_score', 'seismic_score', 'BFB', 'rearrange_score')] %>%
  reshape2::melt(id='pde_group') %>%
  filter(value>0) %>%
  ggplot(aes(x=pde_group,y=value))+
  geom_violin(aes(fill=pde_group))+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  facet_wrap(~variable, scales = 'free', nrow=1)+
  stat_compare_means(comparisons = list(c('High', 'low')), label='p.signif')+
  scale_fill_npg()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

# 7:clone中高低pde构成
a=plot_df[,c('pde_group', 'clone')] %>%
  mutate(clone=factor(clone, levels = c('clone_1','clone_3','clone_5','clone_4','clone_2'))) %>%
  ggplot(aes(x=clone,  fill=pde_group))+
  geom_bar(stat='count', position = 'fill')+
  scale_fill_npg()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
layer_data(a)
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/clone_PDE_构成.pdf',
       a, width=110, height = 80, units='mm', dpi = 600)

#### RNA-seq 映射clonealign #####

## clone align
# run clonealign_huh7.R
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

srt = readRDS('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/cell_line_srt_with_SC3_cytotrace_clonealign_huh7.RDS')
#srt = RunUMAP(srt, dims = 1:30, min.dist = 1)
# DimPlot(srt, group.by = 'cell_line')
srt = subset(srt, cell_line=='Huh7'&clonealign!='unassigned')
srt = process_srt(srt)

srt@meta.data[srt$clonealign=='clone_4', 'clonealign'] = 'clone_3'
srt@meta.data[srt$clonealign=='clone_5', 'clonealign'] = 'clone_4'
clone_col = c('clone_1'='#3B7035','clone_2'='#DEC35D','clone_3'='#DD805A', 'clone_4'='#34B6B5', 'clone_5'='#3152A3')
a = DimPlot(srt, group.by = 'clonealign', pt.size=2, cols=clone_col)
xx = as.data.frame(srt@reductions$umap@cell.embeddings)
xx$clonealign = gsub('huh7.txt_', '', srt$clonealign)

a=xx %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, fill=clonealign))+
  geom_point(shape=21, stroke=0.1)+
  scale_fill_manual(values =clone_col)+
  theme_void()
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/huh7_clonealign_UMAP.pdf',
       a, width=70, height = 50, units='mm', dpi = 600)

RNA_bc = as.data.frame(colnames(srt))
write.table(RNA_bc, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/RNA_bc.txt',
            quote = F, row.names = F, col.names = F)

mtx = srt@assays$RNA@counts
write.table(mtx, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/fitPhylo_wx/提交GEO数据/RNA_count_matrix.txt',
            quote = F, row.names = T, col.names = T)


clone_names = c("clone_4", "clone_3", "clone_1", "clone_2")
y = sctc$clone_graph_obj$clone_graph
y=drop.tip(y, c("clone_3",'clone_2'),trim.internal = F)
y$tip.label = c("clone_3", "clone_4", "clone_1", "clone_2")
y <- ape::rotateConstr(y, rev(clone_names))

#更新clone分支长度
sctc$clone_graph_obj$cell_phylo_tree$data %>%
  filter(isTip==TRUE) %>%
  group_by(group) %>%
  summarise(max(branch.length))

tmp_meta = sctc$orig.data$cell_relation
tmp_meta$clone = sctc$map_obj$cell_map_pos[rownames(tmp_meta), 'clone']
tmp_meta%>%
  group_by(clone) %>%
  summarise(mean(Pseudotime_tree))
xx= y$edge.length
y$edge.length = c(375+1935, 408, 467+3960, 424, 276+2142, 258+2465)
plot(y)

findPointD <- function(x1, y1, x2, y2, M, K) {
  # 计算AB向量的分量
  dx <- x2 - x1
  dy <- y2 - y1

  # 计算AB的长度L
  L1 <- sqrt(dx^2 + dy^2)
  L2 = K*L1/M
  # 计算单位向量UAB
  UAB <- c(dx / L1, dy / L1)

  # 计算点D的坐标
  Dx <- x2 + L2 * (dx / L1)
  Dy <- y2 + L2 * (dy / L1)

  # 返回点D的坐标
  return(c(Dx, Dy))
}

prefix = 'huh7'
tmp_class = DNA_list[[prefix]]
tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
cell_pos = tmp_class$map_obj$cell_map_pos
clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label
# c(375, 310+2288, 358+1684, 408, 467+3960, 424, 276+2142, 258+2465)

clone_pos['clone_2', c('x', 'y')] = findPointD(clone_pos['clone_2', 'xend'],clone_pos['clone_2', 'yend'], clone_pos['clone_2', 'x'],clone_pos['clone_2', 'y'], 375,1935)
clone_pos['clone_3', c('x', 'y')] = findPointD(clone_pos['clone_3', 'xend'],clone_pos['clone_3', 'yend'], clone_pos['clone_3', 'x'],clone_pos['clone_3', 'y'], 276,2142)#358,1684)
clone_pos['clone_1', c('x', 'y')] = findPointD(clone_pos['clone_1', 'xend'],clone_pos['clone_1', 'yend'], clone_pos['clone_1', 'x'],clone_pos['clone_1', 'y'], 467,3960)
clone_pos['clone_4', c('x', 'y')] = findPointD(clone_pos['clone_4', 'xend'],clone_pos['clone_4', 'yend'], clone_pos['clone_4', 'x'],clone_pos['clone_4', 'y'], 258,2465)#276,2142)
#clone_pos['clone_5', c('x', 'y')] = findPointD(clone_pos['clone_5', 'xend'],clone_pos['clone_5', 'yend'], clone_pos['clone_5', 'x'],clone_pos['clone_5', 'y'], 258,2465)

clone_tree = ggtree(y)
clone_tree$data$group=clone_tree$data$label
clone_tree = clone_tree +
  geom_point2(fill='#BAD3E9',shape=21, size=3, color='black')+
  geom_tippoint(aes(fill=group),shape=21, size=5, color='black')+
  # geom_label(aes(x=branch, label=angle), fill='lightgreen')+
  scale_fill_manual(values = clone_col)+
  theme(legend.position = 'none')

clone_pos = clone_tree$data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label

cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]
#
tmp_srt = subset(srt, low_CL==prefix&clonealign!='unassigned')
srt_align_clone = sapply(unique(tmp_srt$clonealign), function(x)gsub(paste0(prefix,'.txt_'), '', x))
tmp_new_xy = c()
for(i in names(srt_align_clone)){
  tmp_xy = clone_pos[srt_align_clone[i],c('x', 'y')]
  cell_xy = subset(tmp_srt, clonealign==i)@reductions$umap@cell.embeddings
  tmp_center = colMeans(cell_xy)
  # 移动
  delt = tmp_center - tmp_xy
  new_xy = t(t(cell_xy)-unlist(delt))
  # 收缩
  delt =  t(t(new_xy)-unlist(tmp_xy))
  new_xy = as.data.frame(new_xy)# +0.975*delt)
  new_xy[,1] = new_xy[,1]+200*delt[,1]
  new_xy[,2] = new_xy[,2]-0.85*delt[,2]
  new_xy$clone = srt_align_clone[i]
  tmp_new_xy = rbind(tmp_new_xy, new_xy)
}

#
# p = ggplot()
# p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
#                      color='#dddddd',
#                      alpha=1,
#                      linewidth=1,
#                      #lwd=pmax(10/nchar(g$branches), 1),
#                      linetype='solid',
#                      lineend = "round",
#                      linejoin='round',
# )
p = ggtree(y)
p$data$group=p $data$label
clone_pos2 = clone_pos[srt_align_clone, ]
p=p+geom_text_repel(aes(x=x,y=y, label=label, color=label),
                    nudge_x=0.2,
                    size=2,
                    data=clone_pos2)
p = p+geom_point(aes(x=UMAP_1, y=UMAP_2, fill=clone),shape=21, stroke=0.1, size=1,data=tmp_new_xy)
p=p+
  #scale_color_igv()+
  scale_fill_manual(values = clone_col)+
  scale_color_manual(values = clone_col)+
  theme_void() +
  labs(title=prefix)+
  guides(alpha=FALSE)+
  #theme(plot.title = element_text(hjust=0.5), legend.position = 'none')+
  theme(plot.title = element_text(hjust=0.5, size=6), legend.position = 'none')+
  scale_y_continuous(expand = expansion(mult = 0.2))+
  scale_x_continuous(expand = expansion(mult = 0.2))#+
#coord_flip()#+scale_x_reverse()
p

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_clonealign_Tree2.pdf',p,
       width=60, height=60, units='mm', dpi = 600, bg = 'transparent')



## DNA RNA clone在Fn比例 ####
meta = sctc$map_obj$cell_map_pos

DNA_clone_num = meta[grepl('clone', meta$clone), ] %>%
  group_by(clone) %>%
  summarise(cell_num=n()) %>%
  mutate(group = 'DNA') %>% as.data.frame()

RNA_clone_num = srt@meta.data %>%
  mutate(clonealign = gsub('huh7.txt_', '', clonealign)) %>%
  group_by(clonealign) %>%
  summarise(cell_num = n()) %>%
  mutate(group='RNA')%>% as.data.frame()

clone_num = as.data.frame(rbind(as.matrix(DNA_clone_num), as.matrix(RNA_clone_num)))
clone_num$cell_num = as.numeric(clone_num$cell_num)

p = clone_num %>%
  ggplot(aes(x=group, y=cell_num, fill=clone))+
  geom_bar(stat='identity', position='fill', width=0.8)+
  scale_fill_manual(values = clone_col)+
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
p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_clone_比例柱图.pdf',p,
       width=40, height=40, units='mm', dpi = 600, bg = 'transparent')


### 高低PDE 比较转录组diversity ########
meta = sctc$map_obj$cell_map_pos
cell = meta[grepl('clone', meta$clone), ]
# 计算转录组diversity
srt = subset(srt, clonealign%in%c('clone_1', 'clone_2')) # 去掉细胞小于50的clone
# gene_expression = FetchData(srt, vars = VariableFeatures(srt), slot = 'counts')
gene_expression = FetchData(srt, vars = VariableFeatures(srt))
gene_expression = GetAssayData(srt)
res = apply(gene_expression, 2, function(y){
  tmp_up = quantile(y, 3/4)
  tmp_dw = quantile(y, 1/4)
  # pos1 = y>tmp_up
  # pos2 = y<tmp_dw
  # x = rep(NA, length(y))
  # x[pos1] = 2
  # x[pos2] = 1
  # names(x) = names(y)
  # if(tmp_dw==tmp_up){
  #   x = cut(y, breaks = 3, labels = F)
  # }else{
  #   x = cut(y, breaks = c(-Inf,tmp_dw, tmp_up,  Inf), labels = F)
  # }
  x = cut(y, breaks = 5, labels = F)
  x = as.numeric(as.vector(x))
  names(x) = names(y)
  return(x)
})
cell_ITH = apply(res, 1, function(a)vegan::diversity(a[a!=1]))
# cell_ITH = apply(res, 1, function(a)var(table(a[a!=1])))
#
# # 计算转录组dist
# srt_pca_center = as.data.frame(srt@reductions$pca@cell.embeddings)
# # srt_pca_center = FetchData(srt, vars = rownames(srt))
#
# srt_pca_center$clone = srt$clonealign
# srt_pca_center = srt_pca_center %>%
#   group_by(clone) %>%
#   summarise_all(list(mean)) %>% as.data.frame()
# rownames(srt_pca_center) = srt_pca_center[,1]
# srt_pca_center = srt_pca_center[,-1]
#
c1_name = rownames(srt@meta.data[srt$clonealign=='clone_1',])
c2_name = rownames(srt@meta.data[srt$clonealign=='clone_2',])

srt_pca = as.data.frame(srt@reductions$pca@cell.embeddings)
srt_pca = as.data.frame(srt@reductions$umap@cell.embeddings)
srt_pca = t(as.matrix(gene_expression))
srt_pca = gene_expression
#srt_pca = FetchData(srt, vars = rownames(srt))
#
# srt_pca = rbind(srt_pca, srt_pca_center)#'center'=colMeans(srt_pca))
# euc_dist = as.matrix(dist(srt_pca,method = "euclidean"))
euc_dist = 1-cor(t(srt_pca))
#euc_dist = as.matrix(dist(srt_pca,method = "euclidean"))
cell_ITH_dist = euc_dist[colnames(srt),colnames(srt)]
cell_ITH_1 = apply(cell_ITH_dist[c1_name, c1_name], 1, mean)
cell_ITH_2 = apply(cell_ITH_dist[c2_name, c2_name], 1, mean)
cell_ITH = c(cell_ITH_1, cell_ITH_2)

####
# filter_cell = rownames(srt_pca[rowSums(srt_pca!=0)>(ncol(srt_pca)*0.10), ])
# srt_pca_filter = srt_pca[, colSums(srt_pca!=0)>(269*0.10)]
# srt_pca_center = as.data.frame(srt_pca_filter)
# srt_pca_center$clone = srt$clonealign
# srt_pca_center = srt_pca_center %>%
#   group_by(clone) %>%
#   summarise_all(list(mean)) %>% as.data.frame()
# rownames(srt_pca_center) = srt_pca_center[,1]
# srt_pca_center = srt_pca_center[,-1]
# srt_pca2 = rbind(srt_pca_filter, srt_pca_center)
#
# #euc_dist = 1-cor(t(srt_pca2))
# euc_dist = as.matrix(dist(srt_pca2,method = "euclidean"))
# cell_ITH_1 = euc_dist['clone_1', c1_name]
# cell_ITH_2 = euc_dist['clone_2', c2_name]
# cell_ITH = c(cell_ITH_1, cell_ITH_2)
#
# x1=cell_ITH_1[cell_ITH_1>quantile(cell_ITH_1, 0.8)]
# x2=cell_ITH_2[cell_ITH_2>quantile(cell_ITH_2, 0.8)]
#
# length(x2)/length(cell_ITH_2)
# length(x1)/length(cell_ITH_1)

###
#filter_cell = rownames(srt_pca[rowSums(srt_pca!=0)>(ncol(srt_pca)*0.10), ])
#srt_pca_filter = srt_pca[, colSums(srt_pca!=0)>(269*0.01)]
#euc_dist = 1-cor(t(srt_pca))
#gene_expression
tmp_func <- function(x){
  #sd(x)/mean(x)
  sd(x)
}
gene_expression = FetchData(srt, vars = VariableFeatures(srt))
srt_pca = gene_expression
euc_dist = as.matrix(dist(srt_pca,method = "euclidean"))
cell_ITH_dist = euc_dist[colnames(srt),colnames(srt)]
cell_ITH_1 = apply(cell_ITH_dist[c1_name, c1_name], 1, tmp_func)
cell_ITH_2 = apply(cell_ITH_dist[c2_name, c2_name], 1, tmp_func)
cell_ITH = c(cell_ITH_1, cell_ITH_2)
###
cell_ITH = as.data.frame(cell_ITH)
cell_ITH$clone = srt$clonealign[rownames(cell_ITH)]
cell_ITH$clone_new = cell_ITH$clone#'other'
#cell_ITH[cell_ITH$clone%in%c('clone_1'), 'clone_new'] = 'clone_1'

a = cell_ITH %>%
  ggplot(aes(x=clone_new,y=cell_ITH))+
  geom_violin(aes(fill=clone_new), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_jitter()+
  stat_compare_means(comparisons = list(c('clone_1', 'clone_2')), size=2)+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  scale_fill_manual(values = clone_col)+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_ITH.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')


# 2.2000 genes expression
gene_expression = FetchData(srt, vars = VariableFeatures(srt))
srt_pca = gene_expression
euc_dist = as.matrix(dist(srt_pca,method = "euclidean"))
cell_ITH_dist = euc_dist[colnames(srt),colnames(srt)]
cell_ITH_1 = apply(cell_ITH_dist[c1_name, c1_name], 1, tmp_func)
cell_ITH_2 = apply(cell_ITH_dist[c2_name, c2_name], 1, tmp_func)
cell_ITH = c(cell_ITH_1, cell_ITH_2)

cell_ITH = as.data.frame(cell_ITH)
cell_ITH$clone = srt$clonealign[rownames(cell_ITH)]
cell_ITH$clone_new = cell_ITH$clone#'other'
a = cell_ITH %>%
  ggplot(aes(x=clone_new,y=cell_ITH))+
  geom_violin(aes(fill=clone_new), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_jitter()+
  stat_compare_means(comparisons = list(c('clone_1', 'clone_2')), size=2)+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  scale_fill_manual(values = clone_col)+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_ITH.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')


cell_ITH %>%
  ggplot(aes(x=cell_ITH, color=clone_new))+
  geom_density()+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  #scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  #scale_fill_manual(values = clone_col)+
  theme_bw()


# cytotrace in clone1 vs other
library(ggridges)
tmp_meta = srt@meta.data
tmp_meta = subset(tmp_meta, clonealign%in%c('clone_1', 'clone_2')) # 去掉细胞小于50的clone
tmp_meta$clone_new = tmp_meta$clonealign#'other'
# tmp_meta[tmp_meta$clonealign%in%c('clone_1'), 'clone_new'] = 'clone_1'
tmp_meta %>%
  ggplot(aes(x=cytotrace_score, y=clone_new, fill=clone_new))+
  #geom_density(bw=0.2)+
  geom_density_ridges(bandwidth=0.05)+
  labs(x='', y='cytotrace_score')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = clone_col)+
  # scale_fill_igv()+
  #scale_fill_manual(values = clone_col)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6)
  )
library(ggbeeswarm)
a =tmp_meta %>%
  ggplot(aes(x=clone_new,y=cytotrace_score))+
  geom_violin(aes( fill=clone_new), size=0.01)+
  geom_boxplot(width=0.1, size=0.01)+
  #geom_quasirandom(color='white', size=0.1)+
  #geom_jitter()+
  labs(x='', y='cytotrace_score')+
  stat_compare_means(comparisons = list(c('clone_1', 'clone_2')))+
  scale_fill_manual(values = clone_col)+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_cytotrace.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

# 增殖评分
srt = AddModuleScore(srt, features = list(c(cc.genes$s.genes, cc.genes$g2m.genes)), name='CC')
srt[['clone_new']] = srt$clonealign#'clone_2'
srt@meta.data[srt@meta.data$clonealign%in%c('clone_1'), 'clone_new'] = 'clone_1'
a =srt@meta.data %>%
  ggplot(aes(x=clone_new,y=CC1))+

  geom_violin(aes( fill=clone_new), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_quasirandom(color='white', size=0.1)+
  #geom_jitter()+
  labs(x='', y='CellCycle')+
  stat_compare_means(comparisons = list(c('clone_1', 'clone_2')))+
  scale_fill_manual(values = clone_col)+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_CC.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

# 识别差异基因clone1 vs other
#srt[['clone_new']] = 'other'
#srt@meta.data[srt@meta.data$clonealign%in%c('clone_2'), 'clone_new'] = 'clone_1'
Idents(srt) = srt$clonealign
DimPlot(srt)
clone_marker = FindMarkers(srt, ident.1 = 'clone_1', ident.2 = 'clone_2')
clone_marker = clone_marker %>% arrange(-avg_log2FC, -pct.1)
srt2 = subset(srt, clonealign%in%c('clone_1', 'clone_2'))
a=DoHeatmap(srt2, features = c(rownames(clone_marker)[1:10], tail(rownames(clone_marker), 10)),
            group.colors=clone_col, raster=F,label=F)+
  scale_fill_gradient2(low='#3471a9',  high='#dc2221')+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_text(size=6),
        axis.text = element_text(size=6),)
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_alignHP.pdf',a,
       width=100, height=50, units='mm', dpi = 600, bg = 'transparent')

# GO
library(clusterProfiler)
library(org.Hs.eg.db)
tmp_marker = clone_marker[clone_marker$avg_log2FC>0&clone_marker$p_val_adj<0.05,]
test1 = bitr(rownames(tmp_marker)[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
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

a = dotplot(tmp_go)
a=a +
  scale_size_continuous(range = c(0.5,2))+
  theme_classic()+
  scale_x_continuous(expand = expansion(0.1))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone1_GO.pdf',a,
       width=120, height=40, units='mm', dpi = 600, bg = 'transparent')
###
tmp_marker = tail(clone_marker, 50)
tmp_marker = tmp_marker[tmp_marker$avg_log2FC<0&tmp_marker$p_val_adj<0.05,]
test1 = bitr(rownames(tmp_marker)[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
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

a = dotplot(tmp_go)
a=a +
  scale_size_continuous(range = c(0.5,2))+
  theme_classic()+
  scale_x_continuous(expand = expansion(0.1))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone2_GO.pdf',a,
       width=120, height=40, units='mm', dpi = 600, bg = 'transparent')


### monocle ######
srt2 = readRDS('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/cell_line_srt_with_SC3_cytotrace_clonealign_huh7.RDS')
srt2 = subset(srt2, cell_line=='Huh7'&clonealign!='unassigned')
srt2 = process_srt(srt2)
srt2 = AddModuleScore(srt2, features = list(c(cc.genes$s.genes, cc.genes$g2m.genes)), name='CC')


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
FC_cds = run_monocle2(data=srt2@assays$RNA@counts,
                      metadata= srt2@meta.data,
                      ordering_genes = VariableFeatures(srt2),
                      FormulaStr = NULL)

a = plot_cell_trajectory(FC_cds, color_by = "clonealign", cell_size = 1)+theme(legend.position = 'right') + scale_color_manual(values = clone_col)
b = plot_cell_trajectory(FC_cds, color_by = "cytotrace_score", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
c = plot_cell_trajectory(FC_cds, color_by = "CC1", cell_size = 1)+theme(legend.position = 'right')+ scale_color_gradientn(colours = c('gray', 'red'))
d = plot_cell_trajectory(FC_cds, color_by = "State", cell_size = 1)+theme(legend.position = 'right')+ scale_color_d3()
a+b+c+d


cell_ITH$State = pData(FC_cds)[rownames(cell_ITH), 'State']
a = cell_ITH %>%
  ggplot(aes(x=State,y=cell_ITH))+
  geom_violin(aes(fill=State), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_jitter()+
  stat_compare_means(comparisons = list(c('1', '2'),c('1', '3'),c('3', '2')), label = 'p.signif')+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  scale_fill_igv()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_ITH_Monocle.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

a =pData(FC_cds) %>%
  ggplot(aes(x=State,y=CC1))+
  geom_violin(aes( fill=State), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_quasirandom(color='white', size=0.1)+
  #geom_jitter()+
  labs(x='', y='CellCycle')+
  stat_compare_means(comparisons = list(c('1', '2'),c('1', '3'),c('3', '2')), label = 'p.signif')+
  scale_fill_igv()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_CC_monocle.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')

a =pData(FC_cds) %>%
  ggplot(aes(x=State,y=cytotrace_score))+
  geom_violin(aes( fill=State), size=0.01)+
  geom_boxplot(width=0.1, size=0.01, outlier.shape = NA)+
  #geom_quasirandom(color='white', size=0.1)+
  #geom_jitter()+
  labs(x='', y='cytotrace_score')+
  stat_compare_means(comparisons = list(c('1', '2'),c('1', '3'),c('3', '2')), label = 'p.signif')+
  scale_fill_igv()+
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0))+
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
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//huh7_Clone_cytotrace_moncle.pdf',a,
       width=50, height=40, units='mm', dpi = 600, bg = 'transparent')




srt2[['State']] = pData(FC_cds)[colnames(srt2), 'State']
srt2$State[srt2$State=='3'] = '2'
Idents(srt2) = srt2$State
DimPlot(srt2)
mono_marker = FindMarkers(srt2, ident.1 = '1', ident.2 = '2')
mono_marker = mono_marker %>% arrange(-avg_log2FC, -pct.1)

tmp_marker = mono_marker[mono_marker$avg_log2FC>0&mono_marker$p_val_adj<0.05,]
test1 = bitr(rownames(tmp_marker)[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
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

dotplot(tmp_go)



#######
clone_pde$group = sapply(clone_pde$PDE > median(clone_pde$PDE), function(x)ifelse(x, 'high','low'))

library(ggpubr)
a = clone_pde %>%
  ggplot(aes(x=group, y=ITH))+
  geom_boxplot(width=0.5)+
  geom_point(aes(fill=clone),shape=21, size=2, width=0.2)+
  #geom_point(shape=21, size=2)+
  # geom_smooth(method = 'loess', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  scale_fill_manual(values = setNames(pal_igv()(8)[c(2:5,7)], sort(unique(xx$clonealign))))+
  stat_compare_means(comparisons = list(c('high', 'low')), label='p.format')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_ITH_相关clone_RNA.pdf',
       a, width=70, height = 50, units='mm', dpi = 600)


cell_ITH_mat = c()
srt_pca = as.data.frame(srt@reductions$pca@cell.embeddings)
srt$clone = gsub('huh7.txt_', '', srt$clonealign)
for(clone in unique(cell_pos$clone)){
  if(grepl('clone', clone)){
    print(clone)
    tmp_clone_cells = rownames(srt@meta.data[srt$clone==clone, ])

    for(ii in 1:(length(tmp_clone_cells)-1)){
      print(ii)
      i = tmp_clone_cells[ii]
      for(jj in (ii+1):length(tmp_clone_cells)){
        print(jj)
        j = tmp_clone_cells[jj]
        x = unlist(srt_pca[i,,drop=T])
        y = unlist(srt_pca[j,,drop=T])

        up_value = sum((x-mean(x))*(y-mean(y)))
        down_value = sqrt(sum((x-mean(x))^2)) * sqrt(sum((y-mean(y))^2))
        cell_ITH_mat = rbind(cell_ITH_mat, c(clone, i, j, 1-up_value/down_value))
      }
    }
  }
}
cell_ITH_mat = as.data.frame(cell_ITH_mat)
colnames(cell_ITH_mat) = c('clone','c1', 'c2', 'cell_ITH')
cell_ITH_mat$cell_ITH = as.numeric(cell_ITH_mat$cell_ITH)

clone_ITH = cell_ITH_mat %>%
  group_by(clone) %>%
  summarise(cell_ITH=mean(cell_ITH)) %>% as.data.frame()
rownames(clone_ITH) = clone_ITH$clone

clone_pde$ITH = clone_ITH[clone_pde$clone, 'cell_ITH']

library(ggpubr)
a = clone_pde %>%
  ggplot(aes(x=group, y=ITH))+
  geom_boxplot(width=0.5)+
  geom_point(aes(fill=clone),shape=21, size=2)+
  #geom_point(shape=21, size=2)+
  # geom_smooth(method = 'loess', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE_F2', y='Transcriptome_ITH')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  scale_fill_manual(values = setNames(pal_igv()(8)[c(2:5,7)], sort(unique(xx$clonealign))))+
  stat_compare_means(comparisons = list(c('high', 'low')), label='p.format', method = 't.test',method.args=list(alternative='greater'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )+
  scale_y_continuous(expand = expansion(mult = 0.2))

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_ITH_相关clone_RNA.pdf',
       a, width=60, height = 50, units='mm', dpi = 600)

#### cytotrace 评分 #######
cyto = srt@meta.data %>%
  group_by(clone) %>%
  summarise(cytotrace=median(cytotrace_score))%>% as.data.frame()
rownames(cyto) = cyto$clone

clone_pde$cytotrace = cyto[clone_pde$clone, 'cytotrace']
a = clone_pde %>%
  ggplot(aes(x=group, y=cytotrace))+
  geom_boxplot(width=0.5)+
  geom_point(aes(fill=clone),shape=21, size=2)+
  #geom_point(shape=21, size=2)+
  # geom_smooth(method = 'loess', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE_F2', y='Cytotrace')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  scale_fill_manual(values = setNames(pal_igv()(8)[c(2:5,7)], sort(unique(xx$clonealign))))+
  stat_compare_means(comparisons = list(c('high', 'low')), label='p.format', method = 't.test',method.args=list(alternative='greater'))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )+
  scale_y_continuous(expand = expansion(mult = 0.2))

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_Cytotrace_相关clone_RNA.pdf',
       a, width=60, height = 50, units='mm', dpi = 600)

srt@meta.data %>%
  ggplot()+
  geom_histogram(aes(x=cytotrace_score, fill=clone))+
  scale_fill_manual(values = setNames(pal_igv()(8)[c(2:5,7)], sort(unique(xx$clonealign))))















### clone tree 着色PDE #####
prefix = 'huh7'
print(prefix)
tmp_class = DNA_list[[prefix]]
tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
cell_pos = tmp_class$map_obj$cell_map_pos
clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label
cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

#
cell_pos$predict_rate = predict_PDE(tmp_class)[rownames(cell_pos)]
x=boxplot(cell_pos$predict_rate)
cell_pos = cell_pos[!cell_pos$predict_rate%in%x$out, ]
p = ggplot()
p = p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=clone_pos, inherit.aes = F,
                     color='#dddddd',
                     alpha=1,
                     linewidth=2,
                     #lwd=pmax(10/nchar(g$branches), 1),
                     linetype='solid',
                     lineend = "round",
                     linejoin='round',
)
p = p+ggshadow::geom_glowline(mapping = aes(x=x, y=y, group=clone, color=predict_rate),
                              data=cell_pos[cell_pos$clone!='clone_0', ],
                              alpha=1,
                              lwd=0.4,
                              #color='gray',
                              lineend='round',linejoin='round',
                              #size=10,
)

p = p + geom_point(aes(x=x, y=y, fill=!!sym('predict_rate')),data=cell_pos,
                   shape=21, size=1,color='#00000000',inherit.aes = F)
p=p +
  scale_color_gradientn(colors = c('blue', 'white', 'red'), values = c(0,0.5,1))+#, limits=c(0,0.5))+
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), values = c(0,0.5,1))#, limits=c(0,0.5))
p=p+  theme_void() +
  labs(title=prefix)+
  guides(alpha=FALSE)+
  #theme(plot.title = element_text(hjust=0.5), legend.position = 'none')+
  theme(plot.title = element_text(hjust=0.5, size=6), legend.position = 'none')+
  scale_y_continuous(expand = expansion(mult = 0.2))+
  scale_x_continuous(expand = expansion(mult = 0.2))+
  coord_flip()+scale_x_reverse()
p
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法//clone_tree_PDE着色.pdf',p,
       width=30, height=30, units='mm', dpi = 600, bg = 'transparent')



### clone水平velocity 与基因组变异比例相关性 （进化多样性PDE着色）####
tmp_class = DNA_list[[prefix]]
tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
cell_pos = tmp_class$map_obj$cell_map_pos
clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label
cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

meta = obj$map_obj$cell_map_pos
meta2 = tmp_class$orig.data$cell_relation[rownames(cell_pos), ]
genome_cnv = meta2$Parent_gain_cn + meta2$Parent_loss_cn
plot_df = data.frame(genome_cnv = genome_cnv, #log1p(meta2$Parent_gain_loc),#denovo_cnv,
                     pde = meta[rownames(cell_pos), 'pde'],
                     clone=cell_pos$clone
)

a = na.omit(plot_df) %>%
  group_by(clone) %>%
  summarise(genome_cnv=mean(genome_cnv), pde=mean(pde)) %>%
  ggplot(aes(x=pde, y=genome_cnv/ncol(tmp_class$orig.data$all_node_data), fill=clone))+
  geom_point(shape=21, size=4)+
  #geom_smooth(method = 'loess', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE', y='%Genome_altered')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_genome_相关clone.pdf',
       a, width=70, height = 50, units='mm', dpi = 600)


### clone水平F2比例与基因组变异比例相关性 （进化多样性PDE着色）####
tmp_class = DNA_list[[prefix]]
tmp_class$map_obj$cell_map_pos$Fn = sapply(rownames(tmp_class$map_obj$cell_map_pos), function(x)strsplit(x, '_')[[1]][1])
tmp_class$map_obj$cell_map_pos$Fn = toupper(tmp_class$map_obj$cell_map_pos$Fn)
cell_pos = tmp_class$map_obj$cell_map_pos
clone_pos = tmp_class$map_obj$clone_seg_data %>% as.data.frame()
rownames(clone_pos) = clone_pos$label
cell_pos = cell_pos[cell_pos$is_virtual=="FALSE", ]

meta = obj$map_obj$cell_map_pos
plot_df = data.frame(Fn = cell_pos$Fn, #log1p(meta2$Parent_gain_loc),#denovo_cnv,
                     pde = meta[rownames(cell_pos), 'pde'],
                     clone=cell_pos$clone
)
tmp_fun <- function(x){
  y = table(x)
  return(ifelse(is.na(y['F2']/sum(y)), 0, y['F2']/sum(y)))
}
a = na.omit(plot_df) %>%
  group_by(clone) %>%
  summarise(F2_prop=tmp_fun(Fn), pde=mean(pde)) %>%
  ggplot(aes(x=pde, y=F2_prop, fill=clone))+
  geom_point(shape=21, size=2)+
  geom_smooth(method = 'loess', formula = 'y ~ x', se=F, linetype='dashed', color='black', linewidth=0.5)+
  labs(x='PDE', y='%F2')+
  #scale_color_gradientn(colours = c('blue','yellow', 'red'))+
  scale_fill_manual(values = setNames(pal_igv()(7)[2:7], sort(unique(plot_df$clone))))+
  # scale_fill_igv()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.text = element_text(size=4),
        legend.key.size = unit(1, 'mm'),
        legend.title = element_text(size=6)
  )

a

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_huh7_修改clone_修改算法/PDE_F2_相关clone.pdf',
       a, width=70, height = 50, units='mm', dpi = 600)





