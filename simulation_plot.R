##############################################################################
#加载包
##批量加载包
{library(ggplot2)#绘图工具
  library(reshape2)#融合数据
  library(cowplot)
  library(openxlsx)
}
##############################################################################
##############################################################################
##############################################################################
#折线图 true与imputed的RMSE/PCC
#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/折线图/PCC_1.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","DropoutRate","value")#更改列名
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)

p1<- ggplot(data = mydata,aes(x=Method,y=value,group = DropoutRate,color=DropoutRate,shape=DropoutRate))+
  geom_point()+
  geom_line()+
  xlab("800×1000")+#横坐标名称
  ylab("PCC")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    #text = element_text(family = "STXihei"),#设置中文字体的显示
    #legend.position = c(.075,.075),#更改图例的位置，放至图内部的左上角
    #legend.box.background = element_rect(color="black")
  )+#为图例田间边框线
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061"))

#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/折线图/PCC_2.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","DropoutRate","value")#更改列名
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)

p2<- ggplot(data = mydata,aes(x=Method,y=value,group = DropoutRate,color=DropoutRate,shape=DropoutRate))+
  geom_point()+
  geom_line()+
  xlab("1000×1000")+#横坐标名称
  ylab("PCC")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    #text = element_text(family = "STXihei"),#设置中文字体的显示
    #legend.position = c(.075,.075),#更改图例的位置，放至图内部的左上角
    #legend.box.background = element_rect(color="black")
  )+#为图例田间边框线
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061"))


#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/折线图/PCC_3.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","DropoutRate","value")#更改列名
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)

p3<- ggplot(data = mydata,aes(x=Method,y=value,group = DropoutRate,color=DropoutRate,shape=DropoutRate))+
  geom_point()+
  geom_line()+
  xlab("1000×5000")+#横坐标名称
  ylab("PCC")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    #text = element_text(family = "STXihei"),#设置中文字体的显示
    #legend.position = c(.075,.075),#更改图例的位置，放至图内部的左上角
    #legend.box.background = element_rect(color="black")
  )+#为图例田间边框线
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061"))


#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/折线图/PCC_4.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","DropoutRate","value")#更改列名
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)

p4<- ggplot(data = mydata,aes(x=Method,y=value,group = DropoutRate,color=DropoutRate,shape=DropoutRate))+
  geom_point()+
  geom_line()+
  xlab("1000×800")+#横坐标名称
  ylab("PCC")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    #text = element_text(family = "STXihei"),#设置中文字体的显示
    #legend.position = c(.075,.075),#更改图例的位置，放至图内部的左上角
    #legend.box.background = element_rect(color="black")
  )+#为图例田间边框线
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061"))


#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/折线图/PCC_5.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","DropoutRate","value")#更改列名
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)

p5<- ggplot(data = mydata,aes(x=Method,y=value,group = DropoutRate,color=DropoutRate,shape=DropoutRate))+
  geom_point()+
  geom_line()+
  xlab("5000×1000")+#横坐标名称
  ylab("PCC")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    #text = element_text(family = "STXihei"),#设置中文字体的显示
    #legend.position = c(.075,.075),#更改图例的位置，放至图内部的左上角
    #legend.box.background = element_rect(color="black")
  )+#为图例田间边框线
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061"))


library(lemon)
pp<-grid_arrange_shared_legend(p1,p2,p3,p4,p5,
                               ncol = 2, nrow = 3,position='top')
pp
ggsave(pp,filename = "/Users/kangkangyue/Desktop/plot/PCC_new1.pdf",width = 12,height = 10,dpi=300)



##############################################################################
##############################################################################
##############################################################################
#cell-cell/gene-gene 相似性  分组柱状图

#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/分组柱状图/gene_PCC_1.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
#mydata <- mydata[order(mydata[,1]),] 
colnames(mydata) <- c("Method","Dropout","value")#更改列名
mydata$Dropout<- factor(mydata$Dropout,
                        levels = c('20%','30%','40%','50%','60%'),ordered = TRUE)
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)
#绘图
a <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),#使柱子并排放置
           color="white", width=.9) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+#消除x轴与绘图区的间隙
  labs(x = '800×1000', y = 'PCC')+
  theme_bw()+
  theme(
    # panel.grid.major =element_blank(),  # 去除边框
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    # text=element_text(family="Songti SC",size=12,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    # axis.text.x = element_text(size = 10, color = "black"),##设置x轴字体大小
    # axis.text.y = element_text(size = 10, color = "black")##设置y轴字体大小
  )  +   
  guides(fill = guide_legend(title = NULL))+ # 删掉图例名称
  scale_fill_manual(values = c("#C10025","#E75646","#FFA27B","#FFDAC4","#CDE5F1","#84C7E1","#1A95C8","#0064B2","#053061"))
#scale_fill_brewer(palette="RdYlBu")


#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/分组柱状图/gene_PCC_2.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
#mydata <- mydata[order(mydata[,1]),] 
colnames(mydata) <- c("Method","Dropout","value")#更改列名
mydata$Dropout<- factor(mydata$Dropout,
                        levels = c('20%','30%','40%','50%','60%'),ordered = TRUE)
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)
#绘图
b <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),#使柱子并排放置
           color="white", width=.9) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+#消除x轴与绘图区的间隙
  labs(x = '1000×1000', y = 'PCC')+
  theme_bw()+
  theme(
    # panel.grid.major =element_blank(),  # 去除边框
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    # text=element_text(family="Songti SC",size=12,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    # axis.text.x = element_text(size = 10, color = "black"),##设置x轴字体大小
    # axis.text.y = element_text(size = 10, color = "black")##设置y轴字体大小
  )  +   
  guides(fill = guide_legend(title = NULL))+ # 删掉图例名称
  scale_fill_manual(values = c("#C10025","#E75646","#FFA27B","#FFDAC4","#CDE5F1","#84C7E1","#1A95C8","#0064B2","#053061"))
#scale_fill_brewer(palette="RdBu")

#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/分组柱状图/gene_PCC_3.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
#mydata <- mydata[order(mydata[,1]),] 
colnames(mydata) <- c("Method","Dropout","value")#更改列名
mydata$Dropout<- factor(mydata$Dropout,
                        levels = c('20%','30%','40%','50%','60%'),ordered = TRUE)
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)
#绘图
c <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),#使柱子并排放置
           color="white", width=.9) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+#消除x轴与绘图区的间隙
  labs(x = '1000×5000', y = 'PCC')+
  theme_bw()+
  theme(
    # panel.grid.major =element_blank(),  # 去除边框
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    # text=element_text(family="Songti SC",size=12,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    # axis.text.x = element_text(size = 10, color = "black"),##设置x轴字体大小
    # axis.text.y = element_text(size = 10, color = "black")##设置y轴字体大小
  )  +   
  guides(fill = guide_legend(title = NULL))+ # 删掉图例名称
  scale_fill_manual(values = c("#C10025","#E75646","#FFA27B","#FFDAC4","#CDE5F1","#84C7E1","#1A95C8","#0064B2","#053061"))
#scale_fill_brewer(palette="RdBu")

#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/分组柱状图/gene_PCC_4.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
#mydata <- mydata[order(mydata[,1]),] 
colnames(mydata) <- c("Method","Dropout","value")#更改列名
mydata$Dropout<- factor(mydata$Dropout,
                        levels = c('20%','30%','40%','50%','60%'),ordered = TRUE)
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)
#绘图
d <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),#使柱子并排放置
           color="white", width=.9) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+#消除x轴与绘图区的间隙
  labs(x = '1000×800', y = 'PCC')+
  theme_bw()+
  theme(
    # panel.grid.major =element_blank(),  # 去除边框
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    # text=element_text(family="Songti SC",size=12,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    # axis.text.x = element_text(size = 10, color = "black"),##设置x轴字体大小
    # axis.text.y = element_text(size = 10, color = "black")##设置y轴字体大小
  )  +   
  guides(fill = guide_legend(title = NULL))+ # 删掉图例名称
  scale_fill_manual(values = c("#C10025","#E75646","#FFA27B","#FFDAC4","#CDE5F1","#84C7E1","#1A95C8","#0064B2","#053061"))
#scale_fill_brewer(palette="RdBu")

#读入数据
mydata <- read.csv("/Users/kangkangyue/Desktop/plot/分组柱状图/gene_PCC_5.csv")
names(mydata)<-c("Method","60%","50%","40%","30%","20%")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
#mydata <- mydata[order(mydata[,1]),] 
colnames(mydata) <- c("Method","Dropout","value")#更改列名
mydata$Dropout<- factor(mydata$Dropout,
                        levels = c('20%','30%','40%','50%','60%'),ordered = TRUE)
mydata$Method<- factor(mydata$Method,
                       levels = c('dropout','scINRB','drimpute','scimpute','scrabble','viper','SAVER','MAGIC','CMF'),ordered = TRUE)
#绘图
e <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),#使柱子并排放置
           color="white", width=.9) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+#消除x轴与绘图区的间隙
  labs(x = '5000×1000', y = 'PCC')+
  theme_bw()+
  theme(
    # panel.grid.major =element_blank(),  # 去除边框
    # panel.grid.minor = element_blank(),
    # panel.background = element_blank(),
    # text=element_text(family="Songti SC",size=12,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
    axis.title =  element_text(size=10,face = "bold"),#设置标题字体大小
    # axis.text.x = element_text(size = 10, color = "black"),##设置x轴字体大小
    # axis.text.y = element_text(size = 10, color = "black")##设置y轴字体大小
  )  +   
  guides(fill = guide_legend(title = NULL))+ # 删掉图例名称
  scale_fill_manual(values = c("#C10025","#E75646","#FFA27B","#FFDAC4","#CDE5F1","#84C7E1","#1A95C8","#0064B2","#053061"))
#scale_fill_brewer(palette="RdBu")


#######合并多张图#########
# library(cowplot)
# pp<-plot_grid(a,b,c,d,e, labels = LETTERS[1:5], ncol = 3)

library(lemon)
pp<-grid_arrange_shared_legend(a, b, 
                               c, d,
                               e,
                               ncol = 2, nrow = 3,position='top')

#保存图片
ggsave(pp,filename = "/Users/kangkangyue/Desktop/plot/gene_PCC.pdf",width = 12,height = 10,dpi=300)

#library(eoffice)
# topptx(pp,filename = "/Users/kangkangyue/Desktop/plot/gene_PCC.pptx")

# 将这个ggplot对象转化成可编辑的对象保存在ppt
library(officer)
editable_graph <- dml(ggobj = pp)
doc <- read_pptx() %>%
  
  add_slide() %>%
  
  ph_with(value = editable_graph,location = ph_location_fullsize()) %>%
  
  print(target = "/Users/kangkangyue/Desktop/plot/pp.pptx") # 这里对导出步骤做了简化，本质上和方法1一样





##############################################################################
#imputed聚类热图

#rm(list=ls())#环境设置，清理环境中的变量
#install.packages('pheatmap')#安装热图所需的R包
require(pheatmap)#加载热图所需R包
require(tidyverse)
install.packages("RColorBrewer")
library(RColorBrewer)
#display.brewer.all()

#df <- read.csv("/Users/kangkangyue/Desktop/plot/1_20.csv",row.names = 1)

mydata <- read.csv("/Users/kangkangyue/Desktop/plot/1_ARI.csv")
#对这一数据进行融合操作
mydata <- melt(mydata,id="Method")
colnames(mydata) <- c("Method","Dropout","ARI")#更改列名

ggplot(mydata,aes(Method,Dropout,fill=ARI))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90 ))
#scale_fill_gradient(low="#E64B35B2", high="#4DBBD5B2")
# scale_y_discrete(position="right") +
# theme_minimal()+
# theme(panel.grid.major=element_blank())

mydata <- as.matrix(read.csv("/Users/kangkangyue/Desktop/plot/1_ARI.csv",row.names = 1))
pheatmap(mydata,cluster_row = FALSE,cluster_cols = FALSE
         ,cellwidth = 42, cellheight = 12,gaps_row = c(10,13)
         ,breaks = c(seq(0.4,1,by = 0.02))
         ,display_numbers = TRUE, number_format = "%s"
         ,border_color = "black",color = colorRampPalette(c( "white", "firebrick3"))(30)
         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
         ,main = "Pearson to P")






##############################################################################
#矩阵热图
##############################################################################

path="/Users/kangkangyue/Desktop/论文/kangyue/LPLS9.28/simulation_data_changerate" #声明test2.R所在位置
setwd(path)  #把工作路径设置到path
data <- readRDS("1_60%.rds")
data_dropout <- as.matrix(data$data_dropout)
data_true <- as.matrix(data$data_true)
data_bulk <- as.matrix(data$data_bulk)
impute_data<- read.csv("/Users/kangkangyue/Desktop/论文/kangyue/LPLS9.28/imputation_FPimpute_changerate_data/FPimpute_1_60%.csv",header = F)
impute_data<-as.matrix(impute_data)



path="/Users/kangkangyue/Desktop/论文/kangyue/analysis" #声明test2.R所在位置
setwd(path)
source('analysis_plot_methods.R')
library("RColorBrewer")
library("reshape2")
library("ggplot2")

pl <- list()
pl[[1]] <- plot_data(log(data_true+1),"True Data")
pl[[2]] <- plot_data(log(data_dropout+1),"Drop-out Data")
pl[[3]] <- plot_data(log(impute_data+1),"impute_data")
pl[[4]] <- plot_data(log(impute_data+1),"impute_data")
pl[[5]] <- plot_data(log(impute_data+1),"impute_data")
pl[[6]] <- plot_data(log(impute_data+1),"impute_data")
pl[[7]] <- plot_data(log(impute_data+1),"impute_data")

main <- gridExtra::grid.arrange(grobs = pl, ncol = 3, top = "")



heatmap(sg, Colv = NA, Rowv = NA, scale="column")

heatmap(sc)


##############################################################################
#cell/gene correlation pheatmap
##############################################################################
#install.packages("corrplot")
#install.packages("pheatmap")
rm(list=ls())
library("pheatmap")
library("corrplot")
library("RColorBrewer")
library("reshape2")
library("ggplot2")
library(gridExtra)
library(lattice)
library(aplot)
library(tidyr)
path="/Users/kangkangyue/Desktop/scINRB" #声明test2.R所在位置
setwd(path)  #把工作路径设置到path
source('functions.R')#“预装“函数
source('analysis_methods.R')

# change_rate=20
# seed_value=1

for(change_rate in c(20,30,40,50)){
  
  for(seed_value in c(1:3,5,6)){
    
    data_cor = get_cor_data(change_rate, seed_value)
    
    data_cell = data_cor[[1]]
    
    data_gene = data_cor[[2]]
    
    cell_min<-min(min(data_cell[[1]]),min(data_cell[[2]]),min(data_cell[[3]]),min(data_cell[[4]]),min(data_cell[[5]]),min(data_cell[[6]]),min(data_cell[[7]]))
    
    cell_min<-round(cell_min,2)
    
    cell_min<-floor(cell_min*10)/10
    
    gene_min<-min(min(data_gene[[1]]),min(data_gene[[2]]),min(data_gene[[3]]),min(data_gene[[4]]),min(data_gene[[5]]),min(data_gene[[6]]),min(data_gene[[7]]))
    
    gene_min<-round(gene_min,2)
    
    gene_min<-floor(gene_min*10)/10
    
    cell_bk= c(seq(cell_min,(1+cell_min)/2-0.01,by=0.01),seq((1+cell_min)/2,1,by=0.01))
    
    gene_bk= c(seq(gene_min,(1+gene_min)/2-0.01,by=0.01),seq((1+gene_min)/2,1,by=0.01))    
    
    name<-list('true','dropout','FPimpute','drimpute','scimpute','scrabble','viper','SAVER','magic','CMF')
    
    gene_pl <- list()
    
    cell_pl <- list()
    
    for(j in c(1:10)){
      
      cell_pl[[j]] <- pheatmap(data_cell[[j]], 
                               scale = "none",
                               legend_breaks = c(seq(cell_min,1,by=0.1)), 
                               legend_labels = c(seq(cell_min,1,by=0.1)),
                               breaks = cell_bk,
                               color = c(colorRampPalette(colors = c("navy","white"))(length(cell_bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(cell_bk)/2)),
                               show_rownames = F,show_colnames = F,border=F,cluster_cols = F,cluster_rows = F,main=name[[j]])
      
      filename= paste0("/cell_correlation/cell_correlation_",name[[j]],"_",seed_value,'_',change_rate,".png")
      
      ggsave(filename,cell_pl[[j]])
      
      
      gene_pl[[j]] <- pheatmap(data_gene[[j]], 
                               scale = "none",
                               legend_breaks = c(seq(gene_min,1,by=0.1)), 
                               legend_labels = c(seq(gene_min,1,by=0.1)),
                               breaks = gene_bk,
                               color = c(colorRampPalette(colors = c("navy","white"))(length(gene_bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(gene_bk)/2)),
                               show_rownames = F,show_colnames = F,border=F,cluster_cols = F,cluster_rows = F,main=name[[j]])
      
      filename= paste0("/gene_correlation/gene_correlation_",name[[j]],"_",seed_value,'_',change_rate,".png")
      
      ggsave(filename,gene_pl[[j]])
      
    }
    
  }
}


#####################################################################
########批量模拟数据tsne降维聚类 画图
for(seed_value in c(1)){
  data <- list()
  p <- list()
  asn <- list()
  for(change_rate in c(60)){
    
    result <- get_simulation_data(change_rate, seed_value)
    data <- result[[1]]
    realcluster<- result[[2]]
    res <- tsen_plot(data,realcluster)
    p <- res[[1]]
    asn <- res[[2]]
    
  }
  
}

# for(k in 1:10){
#   print(k)
#   print(asn[[k]])
# }

library(Rmisc) # 加载包
library(ggpubr) # 加载包
library(ggplot2) # 加载包
#图形组合
# multiplot(p[[1]], p[[2]], p[[3]], p[[4]],p[[5]], p[[6]], p[[7]], p[[8]],p[[9]], p[[10]],   # 要布局的图形
#           cols= 2)
library(lemon)
pp<-grid_arrange_shared_legend(p[[1]], p[[2]], p[[3]], p[[4]],p[[5]], p[[6]], p[[7]], p[[8]],p[[9]], p[[10]],
                               ncol = 10, nrow = 1,position='top')




ggsave(pp,filename = "/Users/kangkangyue/Desktop/plot/tsne/1_60.pdf",dpi=600)

#######################################################################
##############真实数据tsne降维聚类 画图
data <- 
  realcluster<- 
  res <- tsen_plot(data,realcluster)
p <- res[[1]]
asn <- res[[2]]


