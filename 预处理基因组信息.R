library(GenomicFeatures)
library(rtracklayer)

# hg38
gtfFile <- "/Volumes/WX_extend/BioSoftware/gene_info/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- import(gtfFile)
gtf_df <- as.data.frame(gtf)

select_col = c('seqnames', 'start', 'end', 'strand', 'type', 'gene_id', 'gene_name')
gtf_df = gtf_df[,select_col]
gtf_df = gtf_df[gtf_df$type=='gene', ]
gtf_df$seqnames = as.vector(gtf_df$seqnames)
gtf_df = gtf_df[grepl('chr', gtf_df$seqnames), ]
gtf_df = gtf_df[!grepl('chrM', gtf_df$seqnames), ]
rownames(gtf_df) = 1:nrow(gtf_df)
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
cytoband = read.table('/Volumes/WX_extend/BioSoftware/gene_info/cytoBand_hg38.txt', sep='\t')

for(i in 1:nrow(gtf_df)){
  print(i)
  s = gtf_df[i, 'start']
  e = gtf_df[i, 'end']
  c = as.vector(gtf_df[i, "seqnames"])
  tmp_band = cytoband[cytoband$V1==c&s>=cytoband$V2&e<=cytoband$V3, ]
  if(nrow(tmp_band)==0){
    band1 = cytoband[cytoband$V1==c&s>=cytoband$V2, ]
    #band2 = cytoband[cytoband$V1==c&e<=cytoband$V3, ]
    gtf_df[i, 'arm'] = substr(band1[nrow(band1), 'V4'] ,1,1)
  }else{
    gtf_df[i, c('band_start', 'band_end', 'band_name', 'band_type')] = tmp_band[1, 2:5]
    gtf_df[i, 'arm'] = substr(tmp_band[1,'V4'] ,1,1)
  }
}
write.table(gtf_df, '/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg38.txt', sep='\t', row.names = F, quote = F)


# hg19
gtfFile <- "/Volumes/WX_extend/BioSoftware/gene_info/gencode.v43lift37.annotation.gtf.gz"
gtf <- import(gtfFile)
gtf_df <- as.data.frame(gtf)

select_col = c('seqnames', 'start', 'end', 'strand', 'type', 'gene_id', 'gene_name')
gtf_df = gtf_df[,select_col]
gtf_df = gtf_df[gtf_df$type=='gene', ]
gtf_df$seqnames = as.vector(gtf_df$seqnames)
gtf_df = gtf_df[grepl('chr', gtf_df$seqnames), ]
gtf_df = gtf_df[!grepl('chrM', gtf_df$seqnames), ]
rownames(gtf_df) = 1:nrow(gtf_df)
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
cytoband = read.table('/Volumes/WX_extend/BioSoftware/gene_info/cytoBand_hg19.txt', sep='\t')

for(i in 1:nrow(gtf_df)){
  print(i)
  s = gtf_df[i, 'start']
  e = gtf_df[i, 'end']
  c = as.vector(gtf_df[i, "seqnames"])
  tmp_band = cytoband[cytoband$V1==c&s>=cytoband$V2&e<=cytoband$V3, ]
  if(nrow(tmp_band)==0){
    band1 = cytoband[cytoband$V1==c&s>=cytoband$V2, ]
    #band2 = cytoband[cytoband$V1==c&e<=cytoband$V3, ]
    gtf_df[i, 'arm'] = substr(band1[nrow(band1), 'V4'] ,1,1)
  }else{
    gtf_df[i, c('band_start', 'band_end', 'band_name', 'band_type')] = tmp_band[1, 2:5]
    gtf_df[i, 'arm'] = substr(tmp_band[1,'V4'] ,1,1)
  }
}
write.table(gtf_df, '/Volumes/WX_extend/BioSoftware/gene_info/gene_info_hg19.txt', sep='\t', row.names = F, quote = F)


