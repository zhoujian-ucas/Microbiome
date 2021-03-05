####简单线性回归观测数据
#读取鱼类物种丰度和水体环境数据
dat <- read.delim('fish_data.txt', sep = '\t', row.names = 1)

#绘制二维散点图观测各环境变量与鱼类物种丰度的关系
library(ggplot2)

dat_plot <- reshape2::melt(dat, id = 'fish')

p <- ggplot(dat_plot, aes(value, fish)) +
geom_point(color="red",shape = 16,size=2) +
facet_wrap(~variable, ncol = 3, scale = 'free') +
geom_smooth(method = 'lm',color="blue",size=3)
p

#拟合各环境变量与鱼类物种丰度的一元回归，并提取各自 R2 和 p 值
env <- c('acre', 'do2', 'depth', 'no3', 'so4', 'temp')
R2_adj <- c()
p_value <- c()

for (i in env) {
    fit_stat <- summary(lm(dat[['fish']]~dat[[i]]))  #一元线性回归
    R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #提取校正后 R2
    p_value <- c(p_value, fit_stat$coefficients[2,4])  #提取显著性 p 值
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  #数据框中存储了各环境变量与鱼类物种丰度的一元回归的 R2 和 p 值

#在散点图中添加各环境变量与鱼类物种丰度的一元回归的 R2 和 p 值作为标识
#注：该 R2 和 p 值仅为单个变量一元回归的 R2 和 p 值，和下文即将提到的多元回归的 R2 和 p 值存在区别
env_stat$fish <- max(dat$fish) * 0.8
for (i in env) env_stat[i,'value'] <- max(dat[[i]]) * 0.8  #这句和上一句，定义文字在图中的展示坐标，这里仅粗略设置下
env_stat$variable <- rownames(env_stat)
env_stat$label <- paste('R2.adj =', round(env_stat$R2_adj, 3), '\np =', round(env_stat$p_value, 3))  #文字标签

env_stat <- env_stat[c('fish', 'variable', 'value', 'label')]
dat_plot$label <- NA
dat_plot <- rbind(dat_plot, env_stat)  #和先前的数据框合并，便于作图

p + geom_text(data = dat_plot, aes(label = label), size = 3)

####多元线性回归
#首先不妨使用全部环境变量拟合与鱼类物种丰度的多元线性回归
#自变量之间用“+”相连
fit1 <- lm(fish~acre+do2+depth+no3+so4+temp, data = dat)
summary(fit1)  #展示拟合方程的简单统计

#如果数据集中所有预测变量都使用到，且不考虑交互作用，则上式可以简化如下
#使用“.”代替所有预测变量
fit1 <- lm(fish~., data = dat)
summary(fit1)  #展示拟合方程的简单统计

##提取或统计重要数值，例如
coefficients(fit1)  #各变量的斜率和截距项
 
#也可以 names(summary(fit)) 后查看主要的内容项，然后从中提取，例如
summary(fit1)$adj.r.squared  #校正后 R2

#作图，三维空间中最多展示 1 个响应变量和 2-3 个预测变量间的线性关系，例如
library(car)
scatter3d(fish~acre+no3, data = dat)

#要么就考虑使用多个一元回归散点图替代，但是它和多元线性回归是存在一些区别的

####模型比较和选择
#只考虑使用 3 个和鱼类物种丰度线性关系较为明显的环境变量拟合多元线性回归
#这里将处于线性关系临界值的 depth 也算在内
fit2 <- lm(fish~acre+depth+no3, data = dat)
summary(fit2)  #展示拟合方程的简单统计

#anova() 比较两个嵌套模型的拟合优度
anova(fit2, fit1)

#AIC 比较两个回归方程
AIC(fit2, fit1)

