library(glue)
# NG 文献数据
# Reconstructing single-cell karyotype alterations in colorectal cancer identifies punctuated and gradual diversification patterns
# https://zenodo.org/record/4732372

# /Volumes/WX_extend/细胞轨迹推断/data/公开数据/Colorectal
data_path = '/Volumes/WX_extend/细胞轨迹推断/data/公开数据/Colorectal/'

# NG:Fig1b,62 cells, ID=31032019, dataset_P9T_20190331.csv
tmp_data = read.csv(glue('{data_path}/dataset_P9T_20190331.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
tmp_data2 = tmp_data[,4:ncol(tmp_data)]
tmp_data2[tmp_data2>6] = 6
tmp_data[,4:ncol(tmp_data)] = tmp_data2
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F, border_color = 'black')
write.table(tmp_data, glue('{data_path}/处理/Fig1b_dataset_P9T_20190331.txt'), quote=F, row.names = F)

select_cells = c(1,2,8,9,10,16:25, 31, 32)
pheatmap::pheatmap(t(tmp_data[,select_cells+3]), cluster_cols = F, cluster_rows = F, border_color = 'black')
tmp_data2 = tmp_data[,select_cells+3]
pseudo = tmp_data[, -c(select_cells+3)]
pseudo = pseudo[,4:ncol(pseudo)]
#tmp_data2 = cbind(tmp_data2, other=round(rowMeans(pseudo), 0))
tmp_data = cbind(tmp_data[,1:3], tmp_data2)
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F, border_color = 'black')
write.table(tmp_data, glue('{data_path}/处理/Fig1b_sub_dataset_P9T_20190331.txt'), quote=F, row.names = F)

# NG:Fig1c,58 cells, ID=09102018, dataset_P19BT_20181009.csv
tmp_data = read.csv(glue('{data_path}/dataset_P19BT_20181009.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
tmp_data2 = tmp_data[,4:ncol(tmp_data)]
tmp_data2[tmp_data2>7] = 7
tmp_data[,4:ncol(tmp_data)] = tmp_data2
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F)
write.table(tmp_data, glue('{data_path}/处理/Fig1c_dataset_P19BT_20181009.txt'), quote=F, row.names = F)

select_cells = c(1:4,10:15,21:26,49:58)
pheatmap::pheatmap(t(tmp_data[,select_cells+3]), cluster_cols = F, cluster_rows = F, border_color = 'black')
tmp_data2 = tmp_data[,select_cells+3]
pseudo = tmp_data[, -c(select_cells+3)]
pseudo = pseudo[,4:ncol(pseudo)]
#tmp_data2 = cbind(tmp_data2, other=round(rowMeans(pseudo), 0))
tmp_data = cbind(tmp_data[,1:3], tmp_data2)
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = T, border_color = 'black')
write.table(tmp_data, glue('{data_path}/处理/Fig1c_sub_dataset_P19BT_20181009.txt'), quote=F, row.names = F)


# NG:Fig2c,12 cells, ID=09102018, dataset_P9T_20190409_2.csv
tmp_data = read.csv(glue('{data_path}/dataset_P9T_20190409_2.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
tmp_data2 = tmp_data[,4:ncol(tmp_data)]
tmp_data2[tmp_data2>6] = 6
tmp_data[,4:ncol(tmp_data)] = tmp_data2
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F)
write.table(tmp_data, glue('{data_path}/处理/Fig2c_dataset_P9T_20190409_2.txt'), quote=F, row.names = F)
write.table(tmp_data, glue('{data_path}/处理/Fig2c_sub_dataset_P9T_20190409_2.txt'), quote=F, row.names = F)


# NG:Fig3a, 25 cells, ID=09042019_1, dataset_P9T_20190409_1.csv
tmp_data = read.csv(glue('{data_path}/dataset_P9T_20190409_1.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
tmp_data2 = tmp_data[,4:ncol(tmp_data)]
tmp_data2[tmp_data2>6] = 6
tmp_data[,4:ncol(tmp_data)] = tmp_data2
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F)
write.table(tmp_data, glue('{data_path}/处理/Fig3a_dataset_P9T_20190409_1.txt'), quote=F, row.names = F)

select_cells = c(1:7)
pheatmap::pheatmap(t(tmp_data[,select_cells+3]), cluster_cols = F, cluster_rows = F, border_color = 'black')
tmp_data2 = tmp_data[,select_cells+3]
pseudo = tmp_data[, -c(select_cells+3)]
pseudo = pseudo[,4:ncol(pseudo)]
#tmp_data2 = cbind(tmp_data2, other=round(rowMeans(pseudo), 0))
tmp_data = cbind(tmp_data[,1:3], tmp_data2)
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = T, border_color = 'black')
write.table(tmp_data, glue('{data_path}/处理/Fig3a_sub_dataset_P9T_20190409_1.txt'), quote=F, row.names = F)


# NG:Fig4a, 18 cells, ID=04032020, dataset_P9T_4_20200304.csv
tmp_data = read.csv(glue('{data_path}/dataset_P9T_4_20200304.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = F)
write.table(tmp_data, glue('{data_path}/处理/Fig4a_dataset_P9T_4_20200304.txt'), quote=F, row.names = F)

select_cells = c(2:8)
pheatmap::pheatmap(t(tmp_data[,select_cells+3]), cluster_cols = F, cluster_rows = F, border_color = 'black')
tmp_data2 = tmp_data[,select_cells+3]
pseudo = tmp_data[, -c(select_cells+3)]
pseudo = pseudo[,4:ncol(pseudo)]
#tmp_data2 = cbind(tmp_data2, other=round(rowMeans(pseudo), 0))
tmp_data = cbind(tmp_data[,1:3], tmp_data2)
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F, cluster_rows = T, border_color = 'black')
write.table(tmp_data, glue('{data_path}/处理/Fig4a_sub_dataset_P9T_4_20200304.txt'), quote=F, row.names = F)


# NG:Fig4b, 18 cells, ID=04032020, dataset_photoconverted_M3_M4.csv
tmp_data = read.csv(glue('{data_path}/dataset_photoconverted_M3_M4.csv'), sep = ';')
colnames(tmp_data)[1:3] = c('seqnames', 'start', 'end')
tmp_data = tmp_data[, -4]
tmp_data2 = tmp_data[,4:ncol(tmp_data)]
tmp_data2[tmp_data2>6] = 6
tmp_data[,4:ncol(tmp_data)] = tmp_data2
pheatmap::pheatmap(t(tmp_data[,4:ncol(tmp_data)]), cluster_cols = F)
write.table(tmp_data, glue('{data_path}/处理/Fig4b_dataset_photoconverted_M3_M4.txt'), quote=F, row.names = F)


