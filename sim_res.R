source('/Users/lab/wangxin/PyProjects/scTracePlot/R/force_map.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/get_clone_graph.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/load_scTrace.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/scTraceClass.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/plot_cells.R')
source('/Users/lab/wangxin/PyProjects/scTracePlot/R/define_mode_func.R')



sim_dir = '/Users/lab/wangxin/MEDALT/SCTMS/analysis/simdata/sim_for_rearrange/'
sim_dir_output = '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/sim_for_rearrange_out/'
fl = sapply(list.files(sim_dir), function(x)substring(x,1,16))
fl = unique(fl)
sim_chr = data.frame('chr'= c(rep(1,50), rep(2,50)),
                     'start'= c((0:49)*1*1000*1000+2,(0:49)*1*1000*1000+2),
                     'end' = c((1:50)*1*1000*1000+1,(1:50)*1*1000*1000+1)
)

all_res = c()
for(f in fl[82:140]){
  print(f)
  sctc = load_scTrace(sim_dir_output, f)
  sctc = list(orig.data=sctc)
  colnames(sctc$orig.data$all_node_data) = paste0(sim_chr$chr, '_', sim_chr$start, '_', sim_chr$end)
  sim_mechnism = define_CNA_mechnism(sctc, 'ALL', 1,2,0.75, 0.5*1000*1000)


  # 加载真实的拷贝数变异信息
  true_sv = c()
  true_sv_names = c()
  handle = readLines(paste0(sim_dir, f, '.param'))
  for(line in handle){
    tmp_sv = c()
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
      true_sv = rbind(true_sv, c(cell_num,evt_num,evt_length,evt_aneu_num,evt_re_num, length(unique(evt_coverage))/100))
      true_sv_names = c(true_sv_names, line_split[[1]])
    }
  }
  true_sv = as.data.frame(true_sv)
  rownames(true_sv) = true_sv_names
  colnames(true_sv) = c('cell_num', 'event_num', 'event_len', 'aneu_num', 're_num', 'coverage')

  true_sv[, c('rearrange_score', 'limit_score', 'BFB_score', 'wgd')] = sim_mechnism$res[rownames(true_sv), c('rearrange_score', 'limit_score', 'BFB_score', 'wgd')]
  true_sv$file = f
  true_sv$sim_cell_num = nrow(true_sv)
  all_res = rbind(all_res, true_sv)
}

write.table(all_res, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/data/all_sim_events_record.txt', sep='\t', quote = F)
head(all_res)
library(ggpubr)

dens <- MASS::kde2d(all_res$rearrange_score^3,all_res$event_num)
gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
names(gr) <- c("xgr", "ygr", "zgr")
mod <- loess(zgr~xgr*ygr, data=gr)
all_res$pointdens_rearrange <- predict(mod, newdata=data.frame(xgr=all_res$rearrange_score^3, ygr=all_res$event_num))

all_res %>%
  #mutate(event_len = (event_len-min(event_len))/(max(event_len)-min(event_len))) %>%
  #mutate(rearrange_score = (rearrange_score-min(rearrange_score))/(max(rearrange_score)-min(rearrange_score))) %>%
  ggplot(aes(x=rearrange_score^3, y=event_num)) +
  geom_point(aes(color=pointdens_rearrange),size=0.5)+
  scale_colour_gradientn(colours = rainbow(2))
  stat_cor()+
  theme_bw()
library(LSD)
x = heatscatter(all_res$rearrange_score^3,all_res$event_num,cor=T)

all_res %>%
  ggplot(aes(x=wgd, y=coverage)) +
  geom_boxplot()+
  lims(y=c(0.75,1))+
  theme_bw()

newick_to_phylo <- function(newick_dir){
  newick_file = file(newick_dir,open='r')
  newick_str = readLines(newick_file, n = 1)
  close(newick_file)
  newick_str = paste0('(', newick_str, ');')
  newick_str <- read.tree(text = newick_str)
  return(newick_str)
}
all_rate_res = c()
for(f in fl){
  print(f)
  sctc = load_scTrace(sim_dir_output, f)
  # sctc$cell_relation$aneu_rate
  infer_aneu_rate = c()
  tmp_data = sctc$all_node_data
  tree_data = fortify(sctc$tree) %>% as.data.frame()
  leaves_node = tree_data[tree_data$isTip==TRUE, ]
  for(p in unique(leaves_node$parent)){
    tmp_child = leaves_node[leaves_node$parent==p,,drop=F]
    if(nrow(tmp_child)==2){
      c1 = unlist(tmp_child[1, 'label'])
      c2 = unlist(tmp_child[2, 'label'])
      c1_data = tmp_data[c1, ]
      c2_data = tmp_data[c2, ]
      p_data = tmp_data[p, ]
      # aneu rate
      tmp_aneu_rate1 = sum(c1_data!=p_data) / 100
      tmp_aneu_rate2 = sum(c2_data!=p_data) / 100

      tmp_BFB_rate = sum(((c1_data+c1_data)%%2)==0 & (c1_data!=c2_data)) / 100
      infer_aneu_rate = rbind(infer_aneu_rate, c(c1, tmp_aneu_rate1, tmp_BFB_rate))
      infer_aneu_rate = rbind(infer_aneu_rate, c(c2, tmp_aneu_rate2, tmp_BFB_rate))
    }else{
      c1 = unlist(tmp_child[1, 'label'])
      c1_data = tmp_data[c1, ]
      p_data = tmp_data[p, ]
      tmp_aneu_rate1 = sum(c1_data!=p_data) / 100
      infer_aneu_rate = rbind(infer_aneu_rate, c(c1, tmp_aneu_rate1, NA))
    }
  }
  infer_aneu_rate = as.data.frame(infer_aneu_rate)
  colnames(infer_aneu_rate) = c('cellname', 'aneu_rate', 'bfb_rate')
  infer_aneu_rate$aneu_rate = as.numeric(infer_aneu_rate$aneu_rate)
  infer_aneu_rate$bfb_rate = as.numeric(infer_aneu_rate$bfb_rate)
  ######
  true_aneu_rate = c()
  tmp_data = read.table(paste0(sim_dir, f, '.txt'), header = T)
  tmp_data = t(tmp_data[,2:ncol(tmp_data)])
  real_tree = newick_to_phylo(paste0(sim_dir,f,'.newick'))

  tree_data = fortify(real_tree) %>% as.data.frame()
  leaves_node = tree_data[tree_data$isTip==TRUE, ]
  for(p in unique(leaves_node$parent)){
    tmp_child = leaves_node[leaves_node$parent==p,,drop=F]
    if(nrow(tmp_child)==2){
      c1 = unlist(tmp_child[1, 'label'])
      c2 = unlist(tmp_child[2, 'label'])
      c1_data = tmp_data[c1, ]
      c2_data = tmp_data[c2, ]
      p_data = round((c1_data+c2_data) / 2)
      # aneu rate
      tmp_aneu_rate1 = sum(c1_data!=p_data) / 100
      tmp_aneu_rate2 = sum(c2_data!=p_data) / 100

      tmp_BFB_rate = sum(((c1_data+c1_data)%%2)==0 & (c1_data!=c2_data)) / 100
      true_aneu_rate = rbind(true_aneu_rate, c(c1, tmp_aneu_rate1, tmp_BFB_rate))
      true_aneu_rate = rbind(true_aneu_rate, c(c2, tmp_aneu_rate2, tmp_BFB_rate))
    }else{
      c1 = unlist(tmp_child[1, 'label'])
      c1_data = tmp_data[c1, ]
      p_data = c1_data
      tmp_aneu_rate1 = sum(c1_data!=p_data) / 100
      true_aneu_rate = rbind(true_aneu_rate, c(c1, tmp_aneu_rate1, NA))
    }
  }
  true_aneu_rate = as.data.frame(true_aneu_rate)
  colnames(true_aneu_rate) = c('cellname', 'real_aneu_rate', 'real_bfb_rate')
  true_aneu_rate$real_aneu_rate = as.numeric(true_aneu_rate$real_aneu_rate)
  true_aneu_rate$real_bfb_rate = as.numeric(true_aneu_rate$real_bfb_rate)

  aneu_rate = merge(true_aneu_rate, infer_aneu_rate, by='cellname')
  aneu_rate$file = f
  aneu_rate$cell_num = nrow(tmp_data)
  all_rate_res = rbind(all_rate_res, aneu_rate)
}

all_rate_res %>%
  ggplot(aes(x=real_aneu_rate, y=aneu_rate)) +
  geom_point(aes(color='red'),size=0.5)+
  stat_cor()+
  theme_bw()

library(pROC)
library(plotROC)
all_rate_res_bfb = all_rate_res#[!is.na(all_rate_res$bfb_rate), ]
all_rate_res_bfb = all_rate_res#[!is.na(all_rate_res$real_bfb_rate), ]

all_rate_res_bfb$is_bfb = all_rate_res_bfb$real_bfb_rate > 0.4
p=all_rate_res_bfb %>%
  mutate(cell_num = as.factor(cell_num-1)) %>%
  ggplot(aes(d = is_bfb, m = bfb_rate, color=cell_num)) +
  geom_roc(n.cuts = 0, size=0.5) +
  style_roc()+
  ggsci::scale_color_lancet()+
  theme(panel.grid = element_blank())
auc_value = calc_auc(p)
for(i in 1:nrow(auc_value)){
  p = p+annotate("text",x = .5, y = 0.5-0.05*i,
                 label = paste("AUC of ",auc_value[i, 'cell_num']," = ", round(auc_value[i, 'AUC'], 2)))

}
p

all_rate_res %>%
  ggplot(aes(x=real_bfb_rate, y=bfb_rate)) +
  geom_point(aes(color='red'),size=0.5)+
  stat_cor()+
  theme_bw()




