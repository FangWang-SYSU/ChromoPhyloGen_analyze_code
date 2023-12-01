library(Seurat)
library(ComplexHeatmap)
library(ggsci)
# ccRCC

cnv_obj = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNV_ccRCC/tumorCNVstate_mtx_noTherapy.rds')

cnv_obj_gse159115 = readRDS('/Volumes/WX_extend/细胞轨迹推断/data/公开数据/wangying/CNV_ccRCC/infercnv_outs/GSE159115/run.final.infercnv_obj')

cnv_obj_gse159115@expr.data[1:3,1:3]
