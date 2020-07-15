#zlab
#zhoujian-2020
#linux转biom到TXT
#biom convert -i otu_table.biom -o otu.txt --to-tsv
#准备好环境因子表格
#devtools安装ggcor,建议通过原生R安装，R-studio安装会报错。
if(!require(devtools))
  install.packages("devtools")
if(!require(ggcor))
  devtools::install_github("houyunhuang/ggcor")
#读取otu表和环境因子表
s<- read.table("otu.txt", sep = "\t", header = T, row.names = 1, check.names = F)
env <- read.table("env.txt", sep = "\t", header = T, row.names = 1, check.names = F)
spec<-as.data.frame(t(s))
#加载依赖关系包
library(vegan)
library(ggcor)
library(dplyr)
library(ggplot2)
#计算环境因子相关系数矩阵
corr <- fortify_cor(env, type = "upper", show.diag = TRUE,
                    cor.test = TRUE, cluster.type = "all")
head(corr)
#不同群落与单个环境因子的偏mantel test
mantel<-mantel_test(spec, env, mantel.fun = "mantel.randtest",
                       spec.select = list(OTUs = 1:30,
                                          mOTUs = 31:40,
                                          genes = 41:50))
#对分析结果按R和P值重新赋值
mantel<-mantel_test(spec, env,
                        spec.select = list(OTUs =1:30,
                                           mOTUs =31:40,
                                           genes =41:50)) %>%
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                 labels = c("<0.25", "0.25-0.5", ">=0.5"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))
#绘制组合图
p1<-quickcor(env, type = "upper") + geom_square(inherit.aes = TRUE) +
  anno_link(mantel, mapping = aes(colour = p.value, size = r)) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_diag_label(angle = 0,geom = "text", remove.axis = TRUE)
p1

#自定义热图颜色
p2<-p1+scale_fill_gradient2(midpoint = 0, low = "#56B4E9", mid = "white",
                            high = "#77C034", space = "Lab" )
p2
#自定义连线
p3<-p2+scale_color_manual(values=c("#77C034", "#E69F00", "#56B4E9"))
p3
#定义图例标题
p3+guides(size=guide_legend(title="Mantel's R",override.aes=list(colour="grey35"),order=2),
          colour=guide_legend(title="Mantel's P",override.aes = list(size=3),order=1),
          fill=guide_colorbar(title="Pearson's R",order=3))
