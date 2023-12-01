
x = log1p(abs(rnorm(100, 1, 5)))
y= x^2+rnorm(100, 1, 1)

pd = data.frame('x'=x, 'y'=y)

pd$obs = 'Virtual'
pd[sample(1:100, 40),'obs'] = 'Real'
write.table(pd, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图.txt')
pd = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图.txt')
a=pd %>%
  ggplot(aes(x=abs(x),y=abs(y)))+
  geom_point(aes(color=obs), size=0.2)+
  geom_smooth(se=F, color='black', linetype='dashed', linewidth=0.25)+
  scale_color_manual(values = c('Real'='#EA5514', 'Virtual'='#329393'))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size=4),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25))
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图.pdf',
       a, width=30, height=30, units='mm', dpi = 600, bg = 'transparent')


y = runif(100, -1,1)
x= y^2+runif(100, -0.2,0.2)

pd = data.frame('x'=x, 'y'=y)

pd$rate = pd$x
aa = sample(1:100, 40)
pd[aa,'rate'] = pd[aa,'y']

write.table(pd, '/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图clone.txt')
pd = read.table('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图clone.txt')
a=pd %>%
  ggplot(aes(x=x,y=y))+
  geom_point(aes(color=rate), size=0.2)+
  scale_color_gradientn(colours = c('blue', 'red'))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size=4),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25))
a
ggsave('/Users/lab/wangxin/单细胞进化轨迹推断/手稿文件/第一次/原图/Fig1_散点示意图clone.pdf',
       a, width=30, height=30, units='mm', dpi = 600, bg = 'transparent')
