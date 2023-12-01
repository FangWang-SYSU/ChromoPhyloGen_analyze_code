###### BFB_WGD ######
sim_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_bfb_wgd2/'
sim_dir_output = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sim_for_bfb_wgd2_out/'
sim_dir_output2 = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sim_for_bfb_wgd2_out/'
fl = sapply(list.files(sim_dir), function(x)substring(x,1,16))
fl = unique(fl)

res_wgd = c()
for(f in fl){
  print(f)
  sctc = load_scTrace(sim_dir_output, f)
  CNA_mechnism = read.table(glue('{sim_dir_output2}/{f}mode.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score = read.table(glue('{sim_dir_output2}/{f}re_score.txt'), header = T, sep=',', row.names = "X")
  all_limit_prop = read.table(glue('{sim_dir_output2}/{f}limit_score.txt'), header = T, sep=',', row.names = "X")
  all_rearrange_score_pvalue = read.table(glue('{sim_dir_output2}/{f}re_score_pvalue.txt'), header = T, sep=',', row.names = "X")
  wgd_prop = read.table(glue('{sim_dir_output2}/{f}wgd_density.txt'), header = T, sep=',', row.names = "X")
  colnames(wgd_prop) = c('wgd0_prop', 'wgd1_prop', 'wgd2_prop')
  rownames(wgd_prop) = rownames(CNA_mechnism)
  # 加载真实的拷贝数变异信息
  true_sv = c()
  true_sv_names = c()
  handle = readLines(paste0(sim_dir, f, '.param'))
  for(line in handle){
    tmp_sv = c()
    if(grepl('aneuploidy_rate', line)){
      aneuploidy_rate = strsplit(line, '=')[[1]][2]
    }
    if(grepl('root_comment', line)){
      root_comment = strsplit(line, '=')[[1]][2]
    }

    line_split = strsplit(line, '\t')[[1]]
    if(grepl('virtual',line_split[[1]])){
      if(length(line_split)==1){
        true_sv = rbind(true_sv, c(0,0,0,0,0, 0))
        true_sv_names = c(true_sv_names, line_split[[1]])
        next
      }
      tmp_info = strsplit(line_split[[2]], ';')[[1]]

      cell_num = length(tmp_info)
      evt_num = 0
      evt_length = 0
      evt_aneu_num = 0
      evt_re_num = 0
      evt_coverage = c()
      for(events in tmp_info){
        tmp_evt = strsplit(events, '\\[|\\]|\\]\\[')[[1]]
        evt_name = tmp_evt[1]
        if(evt_name=='aneu'){
          evt_aneu_num = evt_aneu_num +1
        }else{
          evt_re_num = evt_re_num+1
        }
        evt_locs = as.numeric(strsplit(tmp_evt[2], ',| |, ')[[1]])
        evt_len = as.numeric(strsplit(tmp_evt[3], ',| |, ')[[1]])
        evt_num = evt_num+length(evt_locs)
        evt_length = evt_length+sum(evt_len)
        for(el in 1:length(evt_locs)){
          if(is.na(evt_len[el])){
            next
          }
          evt_coverage = c(evt_coverage, evt_locs[el]:(evt_locs[el]+evt_len[el]))
        }
      }
      true_sv = rbind(true_sv, c(cell_num, evt_num, evt_length, evt_aneu_num, evt_re_num,
                                 length(unique(evt_coverage))/100))
      true_sv_names = c(true_sv_names, line_split[[1]])
    }
  }
  true_sv = as.data.frame(true_sv)
  true_sv$aneuploidy_rate = aneuploidy_rate
  true_sv$root_comment = root_comment
  rownames(true_sv) = true_sv_names
  colnames(true_sv) = c('cell_num', 'event_num', 'event_len', 'aneu_num', 're_num', 'coverage', 'real_aneu', 'real_wgd')

  true_sv[, 'rearrange_score'] = rowSums(all_rearrange_score[rownames(true_sv), ])
  true_sv[, 'WGD_label'] = CNA_mechnism[rownames(true_sv), 'wgd']
  true_sv[, 'ploidy'] = rowMeans(sctc$all_node_data)[rownames(true_sv)]
  true_sv[, c('wgd0_prop', 'wgd1_prop', 'wgd2_prop')] = wgd_prop[rownames(true_sv), ]
  true_sv$file = f
  true_sv$sim_cell_num = nrow(true_sv)
  res_wgd = rbind(res_wgd, true_sv)
}


res_wgd2 = res_wgd
res_wgd2$real_wgd[res_wgd2$real_wgd>3]='W2'
res_wgd2$real_wgd[res_wgd2$real_wgd==2]='W0'
res_wgd2$real_wgd[res_wgd2$real_wgd==3]='W1'


real_0 = res_wgd2$real_wgd=='W0'
pROC::plot.roc(real_0, res_wgd2$wgd0_prop)
real_1 = res_wgd2$real_wgd=='W1'
pROC::plot.roc(real_1, res_wgd2$wgd1_prop)
real_2 = res_wgd2$real_wgd=='W2'
pROC::plot.roc(real_2, res_wgd2$wgd2_prop)

aa=ggroc(list('WGD0'=roc(real_0, res_wgd2$wgd0_prop),
           'WGD1' = roc(real_1, res_wgd2$wgd1_prop),
           'WGD2' = roc(real_2, res_wgd2$wgd2_prop)
           ),
      legacy.axes = TRUE, size = 0.5)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw()+ggtitle('WGD ROC')+
  #scale_color_manual(values=c('WGD0'='#F9E4D4', 'WGD1'='#EAB78B', 'WGD2'='#E37931'))+
  scale_color_igv()+
  #scale_color_manual(values=c('WGD0'='#F9E4D4', 'WGD1'='#EAB78B', 'WGD2'='#E37931'))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.6))+
  annotate("text",x=0.75,y=0.25,label=paste0('WGD0 AUC:', round(auc(real_0, res_wgd2$wgd0_prop),4)), color='black')+
  annotate("text",x=0.75,y=0.15,label=paste0('WGD1 AUC:', round(auc(real_1, res_wgd2$wgd1_prop),4)), color='black')+
  annotate("text",x=0.75,y=0.05,label=paste0('WGD2 AUC:', round(auc(real_2, res_wgd2$wgd2_prop),4)), color='black')

ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图_仿真/WGD_AUC.pdf',
       aa, width=110, height=80, units='mm', dpi = 600, bg = 'transparent')

