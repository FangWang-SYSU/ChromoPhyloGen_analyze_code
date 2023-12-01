setwd('/Volumes/WX_extend/细胞轨迹推断/RNAseq/')
library(reticulate)
use_python("/Users/lab/.virtualenvs/clonealign_env/bin/python")
use_virtualenv("clonealign_env")
library(tensorflow)
#install_tensorflow()
#library(tensorflow)
# tf$constant("Hello Tensorflow!")
tf$print(
  tf$constant('Hello, TensorFlow!')
)
# devtools::install_github("kieranrcampbell/clonealign")
library(clonealign)
################
chrinfo = DNAcopy::cytoBand
chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
  group_by(chromNum) %>%
  summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
chrinfo_chr_absstart$chromNum = gsub('chr0', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('chr', '', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('X', '23', chrinfo_chr_absstart$chromNum)
chrinfo_chr_absstart$chromNum = gsub('Y', '24', chrinfo_chr_absstart$chromNum)
rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum

library(Seurat)
srt = readRDS('./output/cell_line_srt_with_SC3_cytotrace_new.RDS')

expdata = FetchData(srt, vars = VariableFeatures(srt), slot = 'counts')
# cna
cna_dir = '/Users/lab/wangxin/MEDALT/SCTMS/example/Mydata_tree/'
file_names = c(
  'huh7.txt'
  #'plc.txt',
  #'hep3b.txt',
  #'smc7721.txt',
  #'hepg2.txt',
  #'mhcc97h.txt',
  #'mhcc97l.txt'
  #'lm3.txt'
)
DNA_sctc = readRDS('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_sctc_all.rds')

cna_data = DNA_sctc$huh7$orig.data$all_node_data
cna_data = cna_data[grepl('f2', rownames(cna_data)), ]

srt[['low_CL']] = tolower(srt$cell_line)

srt = FindVariableFeatures(srt, nfeatures = 1000)

tmp_srt = subset(srt, low_CL=='huh7')
tmp_srt = seurat_process(tmp_srt)
expdata = FetchData(tmp_srt, vars = rownames(tmp_srt), slot = 'counts')


expdata = AverageExpression(tmp_srt, features = VariableFeatures(tmp_srt), group.by = 'seurat_clusters')$RNA
expdata = t(expdata)
# 位点转换成基因
cna_loc = as.data.frame(t(sapply(colnames(cna_data), function(x)strsplit(x, '_')[[1]])))
cna_loc$V2 = as.numeric(cna_loc$V2)
cna_loc$V3 = as.numeric(cna_loc$V3)


gene_info = read.table('/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', header = T)
gene_info$seqnames = gsub('chr','', gene_info$seqnames)
gene_info$seqnames = gsub('X','23', gene_info$seqnames)
gene_info$seqnames = gsub('Y','24', gene_info$seqnames)

gene_info = gene_info[match(colnames(expdata),gene_info$gene_name), ]
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
  new_gene_loc = rbind(new_gene_loc, data.frame(loc=rownames(tmp_loc)[pos1:pos2], gene=gn))

}
new_gene_loc = as.data.frame(new_gene_loc)
new_gene_loc = na.omit(new_gene_loc)

# cna gene level
cna_data_gene = t(cna_data)
cna_data_gene = as.data.frame(cna_data_gene[new_gene_loc$loc, ])
cna_data_gene$gene = new_gene_loc$gene
cna_data_gene = cna_data_gene %>%
  group_by(gene) %>%
  summarise_all(list(mean)) %>% as.data.frame()

rownames(cna_data_gene) = cna_data_gene[, 1]
cna_data_gene = as.data.frame(t(cna_data_gene[,-1]))
# cna_data_gene = round(cna_data_gene)

cna_data_gene$clone = DNA_list$huh7$map_obj$cell_map_pos[gsub('-','\\.',rownames(cna_data_gene)), 'clone']
cna_data_gene = na.omit(cna_data_gene)
cna_data_gene = cna_data_gene[grepl('clone', cna_data_gene$clone), ]

#cna_data_gene = cna_data_gene[grepl('f2', rownames(cna_data_gene)), ]
cna_data_gene_clone = cna_data_gene %>%
  group_by(clone) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(cna_data_gene_clone) = cna_data_gene_clone[, 1]
cna_data_gene_clone = as.data.frame(t(cna_data_gene_clone[,-1]))
#cna_data_gene_clone = round(cna_data_gene_clone)

cnadata_filter = cna_data_gene_clone[rowSums(cna_data_gene_clone==0)==0, ]

same_gene = intersect(colnames(expdata), rownames(cnadata_filter))
expdata = as.matrix(expdata[,same_gene])
cnadata_filter = cnadata_filter[same_gene,]

cal <- run_clonealign(expdata, cnadata_filter, n_repeats = 2)
clones <- cal$clone
names(clones) = rownames(expdata)
table(clones)


tmp_srt[['clonealign']] = clones[colnames(tmp_srt)]
VlnPlot(tmp_srt, features = 'cytotrace_score', group.by = 'clonealign', pt.size = 1)
y=FeaturePlot(tmp_srt, features = 'cytotrace_score', pt.size = 1)
x=DimPlot(tmp_srt, group.by = 'clonealign')+scale_color_igv()
x+y
DimPlot(tmp_srt, group.by = 'sc3_cluster')
##
DimPlot(tmp_srt, group.by = 'clonealign')#+scale_color_igv()
saveRDS(tmp_srt,'./output/cell_line_srt_with_SC3_cytotrace_clonealign_huh7.RDS')




seurat_process <-function(obj, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:20)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:20,reduction = "pca", n.components=n.components)
  return(obj)
}

file_names = c(
  'huh7.txt',
  'plc.txt',
  #'hep3b.txt',
  'smc7721.txt',
  'hepg2.txt',
  'mhcc97h.txt',
  'mhcc97l.txt'
  #'lm3.txt'
)
cl_list = list()
cl_plot_list = list()
for(i in file_names){
  CL = strsplit(i, '\\.')[[1]][1]
  tmp_srt = subset(srt, low_CL==CL&clonealign!='unassigned')
  tmp_srt  = seurat_process(tmp_srt)
  cl_list[[i]] = tmp_srt
  cl_plot_list[[i]] = DimPlot(tmp_srt, group.by = 'clonealign')+scale_color_igv()
}
cowplot::plot_grid(plotlist = cl_plot_list, nrow=2)

saveRDS(srt,'./output/cell_line_srt_with_SC3_cytotrace_clonealign.RDS')




