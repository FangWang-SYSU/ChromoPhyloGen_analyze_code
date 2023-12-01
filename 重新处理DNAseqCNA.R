library(glue)
library(ggplot2)

cl_path = '/Volumes/WX_extend/细胞轨迹推断/data/DNA_F1_f2/'


prefix = 'huh7'
fn = 'f1'
intcnv_path = glue('{cl_path}/{prefix}/{fn}/intCNV')
cellSeg_path = glue('{cl_path}/{prefix}/{fn}/cellSeg')

for(file_name in list.files(intcnv_path)){
  bc = strsplit(file_name, '\\.')[[1]][1]
  bc = glue("{fn}_{bc}")

  file_data = read.table(glue('{intcnv_path}/{file_name}'), header = T, sep='\t')
  file_data = na.omit(file_data)

  seg_data = read.table(glue('{cellSeg_path}/{file_name}'), header = T, sep='\t')
  seg_data = na.omit(seg_data)

}

x = data.frame(start=c(0,4,0,8,9),
               end = c(4,5,8,9,11))


# 绘制染色体
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

chrstart_min = 0 # 其实坐标
chrstart_max = max(chrinfo$absEnd) # 终止坐标
p = ggplot()+
  labs(x='chromosome')+
  theme_bw()+
  theme(panel.grid = element_blank())

p=p+
  geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
  scale_x_continuous(breaks=chr_center_pos$center, labels = chrinfo_other$chr_name,
                     limits=c(chrstart_min, chrstart_max),
                     expand = expansion(mult = c(0,0)))
p
# 添加背景颜色
library(RColorBrewer)
bg_color <- rasterGrob(rep(c('red', 'blue'), 12),
                       width=unit(1,"npc"), height = unit(1,"npc"),
                interpolate = TRUE)

p+annotation_custom(bg_color,xmin=0,xmax=1100000000,ymin=-Inf,ymax=Inf)
bg_color = rep(c('gray', 'white'), 12)

col_fun = circlize::colorRamp2(c(1, 24), c("blue", "red"),transparency = 0.8)
col_num = 1:24

x = p
for(i in 1:nrow(chrinfo_chr_absstart)){
  tmp_grop = rasterGrob(col_fun(col_num[i]),
                                    width=unit(1,"npc"), height = unit(1,"npc"),
                                    interpolate = TRUE)
  x = x + annotation_custom(tmp_grop,
                           xmin=chrinfo_chr_absstart$absStart[i],
                           xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
}
x
x+
  geom_segment(aes(x=absStart,y=ratio, xend=absEnd,yend=ratio), data=cna_data, color='black',linewidth=2)+
  lims(y=c(0,3))
# 绘制拷贝数
cna_data = file_data[, c('chrom', 'chrompos', 'end', 'integerCNV', 'ratio')]
cna_data$bin_len = cna_data$end - cna_data$chrompos
cna_data$absEnd = cumsum(as.numeric(cna_data$bin_len))
cna_data$absStart = cna_data$absEnd-cna_data$bin_len

a =p+
  geom_segment(aes(x=absStart,y=ratio, xend=absEnd,yend=ratio), data=cna_data, color='red',linewidth=2)+
  lims(y=c(0,3))


a = a+
  geom_segment(aes(x=abspos,y=ratio, xend=absend,yend=ratio), data=seg_data, color='blue',linewidth=1)+
  lims(y=c(0,3))


b = p+
  geom_segment(aes(x=absStart,y=integerCNV, xend=absEnd,yend=integerCNV), data=cna_data, color='red',linewidth=1)

cowplot::plot_grid(a,b,ncol=1, align='v')

edge_start = 234795382
edge_end = 1233657027+116509016



c = p+geom_curve(aes(x = edge_start, y = 0,
                 xend = edge_end, yend = 0),    # 设置曲率
             curvature = -0.2)+
  lims(y=c(0,0.5))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
cowplot::plot_grid(c,a,b,ncol=1, align='v')
library(ggplot2)

# 创建数据
amplicon_graph = read.table('/Volumes/WX_extend/BioSoftware/AmpliconArchitect/huh7_AA_results/huh7_amplicon3_edges.txt', sep='\t')
amplicon_edges = as.data.frame(do.call(rbind,sapply(amplicon_graph$V1, function(x)strsplit(x, '->|:'))))
colnames(amplicon_edges) = c('from_chr', 'from_pos', 'to_chr', 'to_pos')
amplicon_edges$from_strand = sapply(amplicon_edges$from_pos, function(x)substr(x,nchar(x), nchar(x)))
amplicon_edges$to_strand = sapply(amplicon_edges$to_pos, function(x)substr(x,nchar(x), nchar(x)))
amplicon_edges$from_pos = as.numeric(sapply(amplicon_edges$from_pos, function(x)substr(x,1, nchar(x)-1)))
amplicon_edges$to_pos = as.numeric(sapply(amplicon_edges$to_pos, function(x)substr(x,1, nchar(x)-1)))

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
amplicon_edges$from_chr = change_chr_name(amplicon_edges$from_chr)
amplicon_edges$to_chr = change_chr_name(amplicon_edges$to_chr)


amplicon_edges$from_pos_abs = chrinfo_chr_absstart[amplicon_edges$from_chr, 'absStart'] + amplicon_edges$from_pos
amplicon_edges$to_pos_abs = chrinfo_chr_absstart[amplicon_edges$to_chr, 'absStart'] + amplicon_edges$to_pos
amplicon_edges=amplicon_edges[amplicon_edges$from_pos_abs!=amplicon_edges$to_pos_abs, ]

xstart = apply(amplicon_edges[,c('from_pos_abs', 'to_pos_abs')], 1, min)
xend = apply(amplicon_edges[,c('from_pos_abs', 'to_pos_abs')], 1, max)
c = p+geom_curve(mapping = aes(x = xstart, y = 0,
                               xend = xend, yend = 0),
                 #data=amplicon_edges[3,,drop=F],
                 curvature = -0.3)+
  #lims(y=c(0,0.5))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
c
cowplot::plot_grid(c,a,b,ncol=1, align='v')


# chr7
tmp_chr='chr07'
chr_start = min(chrinfo[chrinfo$chromNum==tmp_chr, 'absStart'])
chr_end = max(chrinfo[chrinfo$chromNum==tmp_chr, 'absEnd'])

chr_start = min(amplicon_edges[amplicon_edges$to_chr==tmp_chr&amplicon_edges$from_chr==tmp_chr, 'from_pos_abs'])
chr_end = max(amplicon_edges[amplicon_edges$to_chr==tmp_chr&amplicon_edges$from_chr==tmp_chr, 'to_pos_abs'])

chr_start=113290429+1233657027
chr_end=116610225+1233657027

a2 = a+scale_x_continuous(limits = c(chr_start, chr_end))
b2 = b+scale_x_continuous(limits = c(chr_start, chr_end))
c2 = c+scale_x_continuous(limits = c(chr_start, chr_end))
cowplot::plot_grid(c2,a2,b2,ncol=1, align='v')


library(DNAcopy)
library(circlize)
library(ggplot2)
plot_genome <- function(
    chr_bg_color=NULL,
    show_x_axis=TRUE,
                        special_region=NULL){
  ##### 准备染色体信息 ######
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

  chrstart_min = 0 # 起始坐标
  chrstart_max = max(chrinfo$absEnd) # 终止坐标

  ##### 绘制画布 ######
  plot_panel = ggplot()+
    theme_bw()+
    theme(panel.grid = element_blank())

  ##### 绘制染色体 ######
  plot_chr=plot_panel+
    geom_vline(xintercept = chrinfo_other$chr_band, color='gray', linewidth=0.5, linetype='dashed')+
    scale_x_continuous(breaks=chr_center_pos$center, labels = chrinfo_other$chr_name,
                       limits=c(chrstart_min, chrstart_max),
                       expand = expansion(mult = c(0,0)))+
    labs(x='chromosome')

  if(!show_x_axis){
    ##### 绘制染色体，没有text，trick ######
    plot_chr = plot_chr+
      theme(axis.text.x = element_blank(),
            axis.ticks.x  = element_blank(),
            axis.title.x = element_blank())


  }
  ##### 添加背景颜色 ######
  if(is.null(chr_bg_color)){
    bg_color = rep(c('gray', 'white'), 12)
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                xmin=chrinfo_chr_absstart$absStart[i],
                                xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }else if(is.character(chr_bg_color)){
    bg_color = rep(chr_bg_color, 24)
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }else if(is.vector(chr_bg_color)){
    bg_color = chr_bg_color
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
    }else{
    bk = chr_bg_color$bk
    col = chr_bg_color$col
    alpha = chr_bg_color$alpha
    col_num = chr_bg_color$value
    col_fun = circlize::colorRamp2(bk, col,transparency = alpha)

    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(col_fun(col_num[i]),
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                xmin=chrinfo_chr_absstart$absStart[i],
                                xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }
  ##### 绘制特定区域(chr:start-end) ######
  if(!is.null(special_region)){
    tmp_chr = strsplit(special_region, ':')[[1]][1]
    tmp_start = strsplit(special_region,':|-')[[1]][2]
    tmp_end = strsplit(special_region,':|-')[[1]][3]
    if(tmp_start=='none'){
      tmp_start=0
    }else{
      tmp_start = as.numeric(tmp_start)
    }
    if(tmp_end=='none'){
      tmp_end = chrinfo_chr_absstart[tmp_chr, 'absEnd'] - chrinfo_chr_absstart[tmp_chr, 'absStart']
    }else{
      tmp_end = as.numeric(tmp_end)
    }
    chr_start = chrinfo_chr_absstart[tmp_chr, 'absStart'] + tmp_start
    chr_end = chrinfo_chr_absstart[tmp_chr, 'absStart'] + tmp_end
    #print(c(chr_start,chr_end))
    plot_chr = plot_chr+scale_x_continuous(limits = c(chr_start, chr_end))
  }

  return(plot_chr)
}




##### 添加信息 ######
# 绘制拷贝数
cna_data = file_data[, c('chrom', 'chrompos', 'end', 'integerCNV', 'ratio')]
cna_data$bin_len = cna_data$end - cna_data$chrompos
cna_data$absEnd = cumsum(as.numeric(cna_data$bin_len))
cna_data$absStart = cna_data$absEnd-cna_data$bin_len

plot_genome(chr_bg_color=NULL,
            show_x_axis=TRUE,
            special_region=NULL)+
  geom_segment(aes(x=absStart,y=ratio, xend=absEnd,yend=ratio), data=cna_data, color='red',linewidth=2)+
  lims(y=c(0,3))





all_amplicons = c()
amplicon_file = list.files(glue('/Volumes/WX_extend/BioSoftware/AmpliconArchitect/huh7_AA_results'), full.names = T)
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


tmp_amplicon = all_amplicons
tmp_amplicon$yend = ifelse((tmp_amplicon$to_pos-tmp_amplicon$to_pos)<1*1000*1000, 1, 0)
tmp_amplicon$yend = ifelse(tmp_amplicon$from_chr!=tmp_amplicon$to_chr, 0, 1)

tmp_pg = plot_genome(show_x_axis = T,
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

tmp_pg


