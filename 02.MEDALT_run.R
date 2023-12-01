library(igraph)
library(ape)
library(pROC)
library(phangorn)
library(ggtree)
library(phangorn)
scripts_dir='/Users/lab/wangxin/MEDALT/'
python_path='/usr/bin/python2'

sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num//'

output_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_medalt/'

fl = sapply(list.files(sim_100_dir), function(x)substring(x,1,16))
fl = unique(fl)

for(f in fl){
    tmp_data = paste0(sim_100_dir, f, '.txt')
    cmd = paste0(
      python_path, ' ',
      scripts_dir, '/scTree.py',  ' ',
      '-P ', scripts_dir,  ' ',
      '-I ', tmp_data,  ' ',
      '-D D ',
      '-G hg19 ',
      '-O ', output_dir, f,' '
    )
    system(cmd)
}

for(f in fl){
  file_path = paste0(output_dir, f, '/CNV.tree.txt')
  tcn = read.table(file_path, header = T)
  if('root' %in% tcn$from){
    root_d = tcn[tcn$from=='root',]
    new_root = root_d[order(root_d$dist),][1,2]
    tcn[tcn$from=='root','from'] = new_root
    tcn = tcn[tcn$to!=new_root,]
  }
  net <- graph_from_data_frame(tcn, directed=T)
  wc <- cluster_walktrap(net)
  phylo = as.phylo.hclust(as.hclust(wc))
  write.tree(phylo, paste0(sim_100_dir, f, '_medalt_tree.txt'))
}


# run sitka
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num/'
output_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_sitka/'
fl = sapply(list.files(sim_100_dir), function(x)substring(x,1,16))
fl = unique(fl)

sitka_sh = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/Rscripts/sitka_pipline.sh'
tmp_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sitka_tmp_dir/'

options(scipen=200)
for(f in fl){
  tmp_data = paste0(sim_100_dir, f, '.txt')
  tmp_data = read.table(tmp_data, header = T)
  end = as.character(1:nrow(tmp_data) * 500000)
  start = as.character(0:(nrow(tmp_data)-1) * 500000 + 1)
  rownames(tmp_data) = paste0(tmp_data$chr, '_', start, '_', end)
  tmp_data = tmp_data[, -c(1,2)]


  write.table(tmp_data, paste0(tmp_dir, f, '.csv'), quote = F, sep=',')
  tmp_input_dir = paste0(tmp_dir, f, '.csv')
  # x = read.delim(tmp_input_dir, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  cmd = paste0(
    'sh ', sitka_sh, ' ',
    tmp_input_dir,  ' ',
    output_dir,  ' ',
    tmp_dir
  )
  system(cmd)
  system(paste0('mv ', output_dir, 'tree.newick ', output_dir, f, '_sitka.newick'))
}


# run medicc2
sim_100_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_RFdist_cell_num/'
output_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_tree_RFdist_medicc2/'
fl = sapply(list.files(sim_100_dir), function(x)substring(x,1,16))
fl = unique(fl)
options(scipen=200)

tmp_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/medicc2_tmp_dir/'

for(f in fl){
  print(f)
  tryCatch({
    tmp_data = paste0(sim_100_dir, f, '.txt')
    tmp_data = read.table(tmp_data, header = T)
    tmp_data[,3:ncol(tmp_data)] = tmp_data[, sample(3:ncol(tmp_data), ncol(tmp_data)-2)]
    end = as.character(1:nrow(tmp_data) * 1)
    start = as.character(0:(nrow(tmp_data)-1) * 1)

    tmp_data[,2] = start
    tmp_cn_data = tmp_data[,3:ncol(tmp_data)]#tmp_data[,sample(3:ncol(tmp_data), 5)]#
    colnames(tmp_cn_data) = paste0(colnames(tmp_cn_data), 1:ncol(tmp_cn_data))

    dup = rownames(t(tmp_cn_data)[duplicated(t(tmp_cn_data)), ])
    x = matrix(0, nrow=nrow(tmp_cn_data),ncol=length(dup))
    for(i in 1:length(dup)){
      spos = sample(1:nrow(tmp_cn_data), 5)
      x[spos, i] = 1
    }
    tmp_cn_data[, dup] = tmp_cn_data[, dup] + x

    tmp_data = cbind(tmp_data[,1:2], end, tmp_cn_data)
    cl = as.integer(nrow(tmp_data)/3)
    c1 = 1:cl
    c2 = (cl+1):(2*cl)
    c3 = (2*cl+1):(nrow(tmp_data))
    tmp_data[c2, 'chr'] = 2
    tmp_data[c3, 'chr'] = 3

    tmp_data = reshape2::melt(tmp_data, id=c('chr', 'absposstart', 'end'))
    tmp_data = tmp_data[, c(4,1,2,3,5)]
    colnames(tmp_data) = c('sample_id','chrom', 'start', 'end', 'cn')
    tmp_data$chrom = paste0('chrom', tmp_data$chrom)

    write.table(tmp_data, paste0(tmp_dir, f, '.tsv'), quote = F, sep='\t', row.names = F)
    tmp_input_dir = paste0(tmp_dir, f, '.tsv')
    # x = read.delim(tmp_input_dir, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
    cmd = paste0(
      '/Users/lab/anaconda3/bin/medicc2 ',
      tmp_input_dir,  ' ',
      output_dir, f, ' ',
      '--input-type tsv ',
      '--total-copy-numbers ',
      '-a cn ',
      '-j 4 '
    )
    system(cmd)
    },error = function(cond){return(NA)})

}


####################### crispr #################
# run sitka
output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sitka/Crispr/'
sitka_sh = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/Rscripts/sitka_pipline.sh'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sitka//tmp/'

options(scipen=200)

for(f in c('A549','786O','U251','Huh7')){
  tmp_data = paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/',f, '_infercnv_int_gene_filter.txt')
  tmp_data = t(read.table(tmp_data, header = T))
  set.seed(1994)
  tmp_data = tmp_data[sample(1:nrow(tmp_data), 500),]
  end = as.character(1:nrow(tmp_data) * 500000)
  start = as.character(0:(nrow(tmp_data)-1) * 500000 + 1)
  rownames(tmp_data) = paste0(1, '_', start, '_', end)
  end = as.character(1:nrow(tmp_data) * 500000)

  write.table(tmp_data, paste0(tmp_dir,'infercnv_int_gene.csv'), quote = F, sep=',')
  tmp_input_dir = paste0(tmp_dir,'infercnv_int_gene.csv')
  # x = read.delim(tmp_input_dir, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
  cmd = paste0(
    'sh ', sitka_sh, ' ',
    tmp_input_dir,  ' ',
    output_dir,  ' ',
    tmp_dir
  )
  system(cmd)
  system(paste0('mv ', output_dir, 'tree.newick ', output_dir, f, '_sitka.newick'))

}

###
# run medicc2
output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medicc2/Crispr'
options(scipen=200)

tmp_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/medicc2_tmp_dir/'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medicc2/tmp/'
for(f in c('Huh7','A549','786O','U251')){
  print(f)
  tryCatch({
    tmp_data = paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/',f, '_infercnv_int_gene_filter.txt')
    tmp_data = t(read.table(tmp_data, header = T))
    tmp_data = cbind(rep(1, nrow(tmp_data)),rep(0, nrow(tmp_data)), tmp_data)
    # tmp_data[,3:ncol(tmp_data)] = tmp_data[, sample(3:ncol(tmp_data), ncol(tmp_data)-2)]
    end = as.character(1:nrow(tmp_data) * 1)
    start = as.character(0:(nrow(tmp_data)-1) * 1)

    tmp_data[,2] = start
    tmp_cn_data = tmp_data[,3:ncol(tmp_data)]
    tmp_data = cbind(tmp_data[,1:2], end, tmp_cn_data)
    cl = as.integer(nrow(tmp_data)/3)
    c1 = 1:cl
    c2 = (cl+1):(2*cl)
    c3 = (2*cl+1):(nrow(tmp_data))
    tmp_data[c2, 1] = 2
    tmp_data[c3, 1] = 3
    colnames(tmp_data)[1:3] = c('chr', 'absposstart', 'end')
    tmp_data = reshape2::melt(as.data.frame(tmp_data), id=c('chr', 'absposstart', 'end'))
    tmp_data = tmp_data[, c(4,1,2,3,5)]
    colnames(tmp_data) = c('sample_id','chrom', 'start', 'end', 'cn')
    tmp_data$chrom = paste0('chrom', tmp_data$chrom)

    write.table(tmp_data, paste0(tmp_dir, f, '.tsv'), quote = F, sep='\t', row.names = F)
    tmp_input_dir = paste0(tmp_dir, f, '.tsv')
    # x = read.delim(tmp_input_dir, check.names = FALSE, stringsAsFactors = FALSE, sep = ",")
    cmd = paste0(
      '/Users/lab/anaconda3/bin/medicc2 ',
      tmp_input_dir,  ' ',
      output_dir, f, ' ',
      '--input-type tsv ',
      '--total-copy-numbers ',
      '-a cn ',
      '-j 4 '
    )
    system(cmd)
  },error = function(cond){return(NA)})
}


###
# run medalt
library(igraph)
library(ape)
library(pROC)
library(phangorn)
library(ggtree)
library(phangorn)


output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medalt/Crispr/'
output_dir =  '/Volumes/WX_extend/BioSoftware/MEDALT/Crispr/'
tmp_dir = '/Volumes/WX_extend/BioSoftware/MEDALT/tmp/'
scripts_dir='/Users/lab/wangxin/MEDALT/'
scripts_dir='/Volumes/WX_extend/BioSoftware/MEDALT/'
python_path='/usr/bin/python2'

for(f in c('Huh7','A549','786O','U251')){
  print(f)
  tmp_data = paste0('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/',f, '_infercnv_int_gene_filter.txt')
  tmp_data = t(read.table(tmp_data, header = T))
  gene.loc <-  read.table('/Volumes/WX_extend/细胞轨迹推断/RNAseq/output/infercnv_crispr2/geneFile.txt')
  colnames(gene.loc) = c('gene', 'chr', 'start', 'end')
  rownames(gene.loc) = gene.loc$gene
  same_gene = intersect(rownames(tmp_data), gene.loc$gene)
  gene.loc = gene.loc[match(same_gene, rownames(gene.loc)), c('chr', 'start')]

  tmp_data = cbind(gene.loc, tmp_data[rownames(gene.loc), ])
  colnames(tmp_data)[1:2] = c('chr', 'absposstart')
  tmp_data[, 1] = gsub('chr', '', tmp_data[, 1])

  tmp_data = tmp_data[!duplicated(paste0(tmp_data[,1], tmp_data[,2])),]

  write.table(tmp_data, paste0(tmp_dir, f, '.tsv'), quote = F, sep='\t', row.names = F)
  tmp_input_dir = paste0(tmp_dir, f, '.tsv')

  cmd = paste0(
    python_path, ' ',
    scripts_dir, '/scTree.py',  ' ',
    '-P ', scripts_dir,  ' ',
    '-I ', tmp_input_dir,  ' ',
    '-D D ',
    '-G hg19 ',
    '-O ', output_dir, f,' '
  )
  system(cmd)
}


medalt_output_dir='/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/medalt/Crispr/'
for(f in c('Huh7','A549','786O','U251')){
  file_path = paste0(output_dir, f, '/CNV.tree.txt')
  tcn = read.table(file_path, header = T)
  if('root' %in% tcn$from){
    root_d = tcn[tcn$from=='root',]
    new_root = root_d[order(root_d$dist),][1,2]
    tcn[tcn$from=='root','from'] = new_root
    tcn = tcn[tcn$to!=new_root,]
  }
  net <- graph_from_data_frame(tcn, directed=T)
  wc <- cluster_walktrap(net)
  phylo = as.phylo.hclust(as.hclust(wc))
  write.tree(phylo, paste0(medalt_output_dir, f, '_medalt_tree.txt'))
}

################# 运行huh7 ###############
# run medalt
library(igraph)
library(ape)
library(pROC)
library(phangorn)
library(ggtree)
library(phangorn)


output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_huh7/'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_huh7/tmp/'
scripts_dir='/Users/lab/wangxin/MEDALT/'
scripts_dir='/Volumes/WX_extend/BioSoftware/MEDALT/'
python_path='/usr/bin/python2'

tmp_input_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/huh7_medalt.tsv'
cmd = paste0(
  python_path, ' ',
  scripts_dir, '/scTree.py',  ' ',
  '-P ', scripts_dir,  ' ',
  '-I ', tmp_input_dir,  ' ',
  '-D D ',
  '-G hg19 ',
  '-O ', output_dir
)
system(cmd)

file_path = paste0(output_dir, '/CNV.tree.txt')
tcn = read.table(file_path, header = T)
if('root' %in% tcn$from){
  root_d = tcn[tcn$from=='root',]
  new_root = root_d[order(root_d$dist),][1,2]
  tcn[tcn$from=='root','from'] = new_root
  tcn = tcn[tcn$to!=new_root,]
}
net <- graph_from_data_frame(tcn, directed=T)
wc <- cluster_walktrap(net)
phylo = as.phylo.hclust(as.hclust(wc))
write.tree(phylo, paste0(output_dir, 'medalt_tree.txt'))


### medicc2
tmp_input_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/huh7_medicc2.tsv'
output_dir='/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_huh7/medicc2/'
cmd = paste0(
  '/Users/lab/anaconda3/bin/medicc2 ',
  tmp_input_dir,  ' ',
  output_dir, ' ',
  '--input-type tsv ',
  '--total-copy-numbers ',
  '-a cn ',
  '-j 4 '
)
system(cmd)



# run sitka
output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/DNA_huh7/sitka/'
sitka_sh = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/Rscripts/sitka_pipline.sh'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sitka//tmp/'
tmp_input_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/huh7_sitka.tsv'

cmd = paste0(
  'sh ', sitka_sh, ' ',
  tmp_input_dir,  ' ',
  output_dir,  ' ',
  tmp_dir
)
system(cmd)
system(paste0('mv ', output_dir, 'tree.newick ', output_dir, 'sitka.newick'))



########### 运行NG数据 #############
# run medalt
library(igraph)
library(ape)
library(pROC)
library(phangorn)
library(ggtree)
library(phangorn)


output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/tmp/'
scripts_dir='/Users/lab/wangxin/MEDALT/'
scripts_dir='/Volumes/WX_extend/BioSoftware/MEDALT/'
python_path='/usr/bin/python2'

dna_file = c('Fig1b_dataset_P9T_20190331.txt',
             'Fig1c_dataset_P19BT_20181009.txt',
             'Fig2c_dataset_P9T_20190409_2.txt',
             'Fig3a_dataset_P9T_20190409_1.txt',
             'Fig4a_dataset_P9T_4_20200304.txt',
             'Fig4b_dataset_photoconverted_M3_M4.txt',
             'Fig3a_sub_dataset_P9T_20190409_1.txt',
             'Fig1c_sub_dataset_P19BT_20181009.txt',
             'Fig4a_sub_dataset_P9T_4_20200304.txt'
)
data_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/Colorectal/处理/'

for(prefix in dna_file){
  print(prefix)
  tmp_data = read.table(paste0(data_dir, prefix), sep=' ', header = T)
  tmp_data = tmp_data[, -3]
  colnames(tmp_data)[1:2] = c('chr', 'absposstart')
  tmp_data[,1] = gsub('X', 23, tmp_data[,1])

  duplicated(tmp_data[,1:2])
  write.table(tmp_data, paste0(data_dir, 'medalt_',prefix), sep='\t', quote = F, row.names = F)
  cmd = paste0(
    python_path, ' ',
    scripts_dir, '/scTree.py',  ' ',
    '-P ', scripts_dir,  ' ',
    '-I ', paste0(data_dir, 'medalt_',prefix),  ' ',
    '-D D ',
    '-G hg19 ',
    '-O ', output_dir
  )
  system(cmd)

  file_path = paste0(output_dir, '/CNV.tree.txt')
  tcn = read.table(file_path, header = T)
  if('root' %in% tcn$from){
    root_d = tcn[tcn$from=='root',]
    new_root = root_d[order(root_d$dist),][1,2]
    tcn[tcn$from=='root','from'] = new_root
    tcn = tcn[tcn$to!=new_root,]
  }
  net <- graph_from_data_frame(tcn, directed=T)
  wc <- cluster_walktrap(net)
  phylo = as.phylo.hclust(as.hclust(wc))
  write.tree(phylo, paste0(output_dir, 'medalt_', prefix, '_tree.txt'))
}



### medicc2
dna_file = c(#'Fig1b_dataset_P9T_20190331.txt',
             #'Fig1c_dataset_P19BT_20181009.txt',
             #'Fig2c_dataset_P9T_20190409_2.txt',
             #'Fig3a_dataset_P9T_20190409_1.txt',
             #'Fig4a_dataset_P9T_4_20200304.txt',
             'Fig4b_dataset_photoconverted_M3_M4.txt'
             #'Fig3a_sub_dataset_P9T_20190409_1.txt',
             #'Fig1c_sub_dataset_P19BT_20181009.txt',
             #'Fig4a_sub_dataset_P9T_4_20200304.txt'
)
data_dir = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/Colorectal/处理/'

output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/medicc2/'
for(prefix in dna_file){
  print(prefix)
  tmp_data = read.table(paste0(data_dir, prefix), sep=' ', header = T)
  tmp_data[,1] = gsub('X', 23, tmp_data[,1])
  tmp_data[,1] = paste0('chrom', tmp_data[,1])
  tmp_data = reshape2::melt(tmp_data, id=c("seqnames","start","end"))
  tmp_data = tmp_data[, c(4,1,2,3,5)]
  colnames(tmp_data) = c('sample_id', 'chrom', 'start', 'end', 'cn')
  write.table(tmp_data, paste0(data_dir, 'medicc2_',prefix), sep='\t', quote = F, row.names = F)


  cmd = paste0(
    '/Users/lab/anaconda3/bin/medicc2 ',
    paste0(data_dir, 'medicc2_',prefix),  ' ',
    output_dir, ' ',
    '--input-type tsv ',
    '--total-copy-numbers ',
    '-a cn ',
    '-j 4 '
  )
  system(cmd)

}

# run sitka
dna_file = c('Fig1b_dataset_P9T_20190331.txt',
             'Fig1c_dataset_P19BT_20181009.txt',
             'Fig2c_dataset_P9T_20190409_2.txt',
             'Fig3a_dataset_P9T_20190409_1.txt',
             'Fig4a_dataset_P9T_4_20200304.txt',
             'Fig4b_dataset_photoconverted_M3_M4.txt',
             'Fig3a_sub_dataset_P9T_20190409_1.txt',
             'Fig1c_sub_dataset_P19BT_20181009.txt',
             'Fig4a_sub_dataset_P9T_4_20200304.txt'
)
sitka_sh = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/Rscripts/sitka_pipline.sh'
tmp_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/sitka/tmp/'
output_dir = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/NG_other/sitka/'
options(scipen=200)
for(prefix in dna_file){
  print(prefix)
  tmp_data = read.table(paste0(data_dir, prefix), sep=' ', header = T)
  tmp_data[,1] = gsub('X', '23', tmp_data[,1])
  #500000
  for(chr in unique(tmp_data[,1])){
    tmp_x = tmp_data[tmp_data[,1]==chr, 2:3]
    tmp_x = data.frame(((1:nrow(tmp_x))-1)*500000+1, (1:nrow(tmp_x))*500000)
    tmp_data[tmp_data[,1]==chr, 2:3]  = tmp_x
  }

  rownames(tmp_data) = paste0(tmp_data[,1],'_', as.character(tmp_data[,2]), '_', as.character(tmp_data[,3]))


  tmp_data = tmp_data[, -c(1:3)]
  write.table(tmp_data, paste0(data_dir, 'sitka_',prefix), sep=',', quote = F, row.names = T)


  cmd = paste0(
    'sh ', sitka_sh, ' ',
    paste0(data_dir, 'sitka_',prefix),  ' ',
    output_dir,  ' ',
    tmp_dir
  )

  system(cmd)
  system(paste0('mv ', output_dir, 'tree.newick ', output_dir, 'sitka_',prefix,'.newick'))
}












