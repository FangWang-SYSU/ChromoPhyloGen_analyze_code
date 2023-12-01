#### 展示seismic 和limit #####
rescore = all_rearrange_score#CNA_mechnism$orig.rearrange_score$all_rearrange_score
rescore = rescore[grepl('f', rownames(rescore)),]
pheatmap::pheatmap(rescore)
# 选择8号染色体
rescore = rescore[order(-rescore[,5]),]
pheatmap::pheatmap(rescore, cluster_rows = F, cluster_cols = F)

pvalue=all_rearrange_score_pvalue

bc = rownames(pvalue[pvalue[, 'X5']>0.01, ])[1]
#
bc = tail(rownames(rescore),1)[1]
bc = rownames(rescore)[20]
bc
#'#C9E2F4', '#DFEDD7'
prefix = 'huh7'
fn = 'f2'
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
pg = plot_genome_chr(chr_num = 'chr05',bg_color = '#C9E2F4')#, start = '59900000', end='127100000')
pg=ggplot()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill='#C9E2F4'),
  )
pg
tmp_cna_data = cna_data[cna_data$chrom=='5',]
tmp_seg_data = seg_data[seg_data$chr=='5',]
for(i in 1:nrow(tmp_seg_data)){
  tmp_start = tmp_seg_data[i, 'start']
  tmp_end = tmp_seg_data[i, 'end']

  tmp_ratio = subset(tmp_cna_data, chrompos>=tmp_start&end<tmp_end)
  new_ratio = mean(tmp_ratio$ratio)
  tmp_seg_data[i, 'ratio'] = new_ratio
}
a1 =pg+
  geom_segment(aes(x=chrompos/1000/1000,y=ratio, xend=end/1000/1000,yend=ratio), data=tmp_cna_data, color='black',linewidth=.3)+
  #geom_point(aes(x=absStart,y=ratio), data=cna_data, color='black',size=0.01)+
  geom_segment(aes(x=start/1000/1000,y=ratio, xend=end/1000/1000,yend=ratio), data=tmp_seg_data, color='black',linewidth=0.5)+
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

b1 = plot_genome_chr(chr_num = 'chr05',bg_color = 'white')+#, start = '59900000', end='127100000')+
  labs(y='CN')+
  lims(y=c(0,7))+
  geom_segment(aes(x=chrompos,y=integerCNV, xend=end,yend=integerCNV),
               data=cna_data[cna_data$chrom=='5',], color='red',linewidth=0.5)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"),
        axis.text.y = element_text(size=6))

p_seismic = cowplot::plot_grid(b1,a1, ncol=1, align = 'v', rel_heights = c(0.3,0.6))
p_seismic






# limit
rescore = all_limit_prop#CNA_mechnism$orig.rearrange_score$all_limit_prop
chr=5
rescore = rescore[all_rearrange_score_pvalue[,chr]<0.001,]
rescore = rescore[rownames(all_seg_len[all_seg_len[,chr]>5,]),]
rescore = rescore[grepl('f', rownames(rescore)),]
pheatmap::pheatmap(rescore)

rescore = rescore[order(-rescore[,chr]),]
pheatmap::pheatmap(rescore, cluster_rows = F, cluster_cols = F)

#
"f2_ATTCTACTCCTACCCA-1"
rescore["f2_ATTCTACTCCTACCCA-1",]
bc = rownames(tail(rescore[rescore[,chr]>0.7,]))[1]
bc = rownames(rescore[rescore[,chr]<=0.4,])[1]
bc = rownames(rescore)[1]
bc

all_seg_len[bc,]
#'#C9E2F4', '#DFEDD7'
prefix = 'huh7'
fn = 'f2'
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
a =pg+
  geom_segment(aes(x=chrompos/1000/1000,y=ratio, xend=end/1000/1000,yend=ratio), data=cna_data[cna_data$chrom==as.character(chr),], color='black',linewidth=0.3)+
  #geom_point(aes(x=absStart,y=ratio), data=cna_data, color='black',size=0.01)+
  geom_segment(aes(x=start/1000/1000,y=ratio, xend=end/1000/1000,yend=ratio), data=tmp_seg_data, color='black',linewidth=0.5)+
  scale_x_continuous(limits=c(0, max(tmp_seg_data$end)/1000/1000),
                     expand = expansion(mult = c(0,0)))+
  labs(y='ratio')+
  lims(y=c(0,5))+
  theme(#axis.text.x = element_text(angle = 60, hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y=element_blank(),
    #axis.line = element_line(linewidth = 0.15),
    axis.ticks.x = element_line(linewidth = 0.15),
    axis.ticks.length=unit(0.025, "cm"),
    axis.text = element_text(size=6))

a

b = plot_genome_chr(chr_num = paste0('chr0', chr),bg_color = 'white')+#, start = '123800000', end='170100000')+
  #labs(y='CN')+
  lims(y=c(0,7))+
  geom_segment(aes(x=chrompos,y=integerCNV, xend=end,yend=integerCNV),
               data=cna_data[cna_data$chrom==as.character(chr),], color='red',linewidth=0.5)+
  theme(axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),#element_line(linewidth = 0.15),
        axis.ticks.length=unit(0.025, "cm"))

p_gradual = cowplot::plot_grid(b,a, ncol=1, align = 'v', rel_heights = c(0.3,0.6))
p_gradual



