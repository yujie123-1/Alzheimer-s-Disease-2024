#有3个血液的数据
#下载文章给出的原始的结果，然后检查下自己算出的count和原来文章作者给出的结果的异同
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
setwd("C:/Users/siat/Downloads")
a <- read.csv("DemoData.csv",header = T,row.names = 1)
a <- as.matrix(a)
a <- t(a)

pheatmap(a)
pdf("heatmap_output.pdf", width = 8, height = 8)
pheatmap(a, scale = "none",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cellwidth = 15,cellheight = 15,
         cluster_rows = FALSE,
         cluster_cols = FALSE)  
dev.off()









#下载原文的血液的数据，批量加载文件，然后确认这3个基因在这3个数据集中的表达量#
#GSE168813 不用处理直接下载
#GSE153881 需要处理
#GSE161199 需要处理
rm(list=ls())
setwd("C:/Users/siat/Downloads/临时文件夹GSE161199_RAW")
# 设置文件夹路径
folder_path <- "C:/Users/siat/Downloads/临时文件夹GSE161199_RAW"

# 获取文件夹中所有 .txt 文件的文件名
txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# 使用 lapply 读取所有表格文件并保存到列表中，同时为每个数据框的列名加上文件名前缀
data_list <- lapply(seq_along(txt_files), function(i) {
  data <- read.table(txt_files[i], header = TRUE,row.names = 1)[,c(1,8:10)]
  # 获取文件名（不带路径和扩展名）
  file_name <- tools::file_path_sans_ext(basename(txt_files[i]))
  # 修改列名，在每个列名前加上文件名
  colnames(data) <- paste(file_name, colnames(data), sep = "_")
  return(data)
})

# 假设每个文件都有一个共同的列（例如“ID”）用于合并
# 将该共同列恢复为统一名称
data_list <- lapply(data_list, function(df) {
  colnames(df)[1] <- "annotation.gene_id"
  return(df)
})

# 使用 Reduce 和 merge 函数将所有数据框合并
merged_data <- Reduce(function(x, y) merge(x, y, by = "annotation.gene_id", all = TRUE), data_list)

# 查看合并后的数据
head(merged_data)[1:5,1:5]
#将生成的文件进行保存
write.csv(merged_data,file = "GSE161199_merged_data.csv",quote = F,row.names = F)
#只提取出来这3个基因的数据
#Gene: SNAP25 (ENSG00000132639) - Summary
#Gene: ENC1 (ENSG00000171617) - Summary
#Gene: VSNL1 (ENSG00000163032) - Summary
SUB_data <- subset(merged_data,merged_data$annotation.gene_id %in% c("ENSG00000132639","ENSG00000171617","ENSG00000163032"))
write.csv(SUB_data,file = "GSE161199_sub_data.csv",quote = F,row.names = F)




#找出第二个数据集的数据
setwd("C:/Users/siat/Downloads/临时文件夹GSE153881_RAW")
# 设置文件夹路径
folder_path <- "C:/Users/siat/Downloads/临时文件夹GSE153881_RAW"

# 获取文件夹中所有 .txt 文件的文件名
txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# 使用 lapply 读取所有表格文件并保存到列表中
data_list <- lapply(txt_files, function(x) read.table(x, header = TRUE))

# 通过 Reduce 和 merge 函数将所有数据框合并
# 假设文件中都有一个共同的列名 "ID" 作为合并依据
merged_data <- Reduce(function(x, y) merge(x, y, by = "Tracking_ID", all = TRUE), data_list)

# 查看合并后的数据
head(merged_data)[1:5,1:5]
#将生成的文件进行保存
write.csv(merged_data,file = "GSE153881_merged_data.csv",quote = F,row.names = F)
#只提取出来这3个基因的数据
#Gene: SNAP25 (ENSG00000132639) - Summary
#Gene: ENC1 (ENSG00000171617) - Summary
#Gene: VSNL1 (ENSG00000163032) - Summary
SUB_data <- subset(merged_data,merged_data$Tracking_ID %in% c("ENSG00000132639","ENSG00000171617","ENSG00000163032"))
write.csv(SUB_data,file = "GSE153881_sub_data.csv",quote = F,row.names = F)




































#ggplot2中设置边框尺寸的方法
#设置外面边框的粗细
# theme(panel.border = element_rect(fill=NA,color="black", size=5, linetype="solid"))
#设置外边框的长和宽


# to compile RawData statistics
#统计AD和Healthy中两个组的差别，大家都是基于统一的标准
#所以就绘制堆积的柱状图
# rm(list = ls())
# getwd()
# setwd("C:/Users/yujie/Desktop/datacollect")
# #加载原始数据的信息
# a <- read.csv("sum.csv",header = T,fill = TRUE)
# colnames(a)
# #加载帮助配色的R包
# library(ggsci)
# #对脑区的数据进行统计
# b <- a[,c("Brain_Region","Severity")]
# library(dplyr)
# d <- b |> group_by(Severity,Brain_Region) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# #对脑区的名字进行更改
# colnames(d)
# d <- d%>%mutate(Brain_Region=case_when(d$Brain_Region==0~"Hippo",
#                                        d$Brain_Region==1~"PC",
#                                        d$Brain_Region==2~"Blood",
#                                        d$Brain_Region==3~"TL",
#                                        d$Brain_Region==4~"TLPC",
#                                        d$Brain_Region==5~"PL",
#                                        d$Brain_Region==6~"OB",
#                                        d$Brain_Region==7~"CB"))
# 
# #将数据转化为因子类型
# d$Severity <- as.factor(d$Severity)
# d$Brain_Region <- as.factor(d$Brain_Region)
# #绘制堆积的柱状图
# library(ggplot2)
# p1 <- ggplot(d,aes(x=Severity,y=freq,fill=Brain_Region))+
#   geom_bar(stat="identity")+
#   #挑选配色vignette( "ggsci")+scale_fill_npg()
#   #scale_fill_aaas()+scale_fill_nejm()+scale_fill_lancet()+scale_fill_jama()
#   #scale_fill_jco()+scale_fill_ucscgb()+scale_fill_d3()+scale_fill_locuszoom()
#   #scale_fill_igv()+scale_fill_uchicago()+scale_fill_startrek()+scale_fill_tron()
#   #scale_fill_futurama()+scale_fill_rickandmorty()+scale_fill_simpsons()
#   scale_fill_npg()+
#   xlab("") +
#   ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of Region") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# p1
# library(egg)
# ggsave("p1.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# colnames(a)
# #对脑区的数据进行统计
# b <- a[,c("Bank_Location","Severity")]
# d <- b |> group_by(Severity,Bank_Location) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# d <- d%>%mutate(Bank_Location=case_when(d$Bank_Location==0~"America",
#                                         d$Bank_Location==1~"Britain",
#                                         d$Bank_Location==2~"China",
#                                         d$Bank_Location==3~"Australia",
#                                         d$Bank_Location==4~"Italy",
#                                         d$Bank_Location==5~"Japan",
#                                         d$Bank_Location==6~"Colombia"))
# 
# 
# #将数据转化为因子类型
# d$Severity <- as.factor(d$Severity)
# d$Bank_Location <- as.factor(d$Bank_Location)
# #绘制堆积的柱状图
# library(ggplot2)
# p2 <- ggplot(d,aes(x=Severity,y=freq,fill=Bank_Location))+
#   geom_bar(stat="identity")+
#   scale_fill_npg()+
#   xlab("") +
#   ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of Bank Location") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# # library(ggpubr)
# # p <- ggarrange(p1,p2,nrow = 2)
# # ggsave("catagory1.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
# #        width = 8, height = 7, units = 'in', dpi = 600)
# ggsave("p2.pdf", egg::set_panel_size(p2, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 
# #对AD的病理特征进行描述
# #对BraakStage进行统计
# rm(list = ls())
# setwd("C:/Users/yujie/Desktop/datacollect")
# #加载原始数据的信息
# a <- read.csv("sum.csv",header = T,fill = TRUE)
# colnames(a)
# #对脑区的数据进行统计
# b <- a[,c("Braak_Stage","Severity")]
# d <- b |> group_by(Severity,Braak_Stage) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# d <- d%>%mutate(Braak_Stage=case_when(d$Braak_Stage==0~"0",
#                                       d$Braak_Stage==1~"1~2",
#                                       d$Braak_Stage==2~"3~4",
#                                       d$Braak_Stage==3~"5~6"))
# 
# #将数据转化为因子类型
# d$Severity <- as.factor(d$Severity)
# d$Braak_Stage <- as.factor(d$Braak_Stage)
# #绘制堆积的柱状图
# library(ggplot2)
# p3 <- ggplot(d,aes(x=Severity,y=freq,fill=Braak_Stage))+
#   geom_bar(stat="identity")+
#   # scale_fill_npg()+
#   scale_fill_manual(values =  rev(c('#D73227', '#47ADCB', '#209073', "#EA866B"))) +
#   # geom_text(aes(label=freq),position=position_stack(vjust=0.5), size=4)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
#   xlab("") +ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of Braak Stage") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# p3
# ggsave("p3.pdf", egg::set_panel_size(p3, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# #对CERAD进行统计，最新的修改版的文章直接去掉CERAD的数据，不使用这个数据
# 
# #对APOE的数据进行统计
# setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
# #加载原始数据的信息
# a <- read.csv("BackgroundInformationAPOE.csv",header = T,fill = TRUE)
# colnames(a)
# #对脑区的数据进行统计
# b <- a[,c("Phenotypedit","Severity")]
# d <- b |> group_by(Severity,Phenotypedit) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# #将数据转化为因子类型
# d$Severity <- as.factor(d$Severity)
# d$Phenotypedit <- as.factor(d$Phenotypedit)
# colnames(d)[2] <- "APOE"
# #绘制堆积的柱状图
# library(ggplot2)
# p5 <- ggplot(d,aes(x=Severity,y=freq,fill=APOE))+
#   geom_bar(stat="identity")+
#   scale_fill_manual(values =  rev(c('#D73227', '#47ADCB', '#209073', "#EA866B"))) +
#   xlab("") +ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of APOE") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# p5
# library(ggpubr)
# ggsave("p7.pdf", egg::set_panel_size(p5, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 
# 
# #对数据的年龄进行统计
# rm(list = ls())
# setwd("C:/Users/yujie/Desktop/datacollect")
# #加载原始数据的信息
# a <- read.csv("sum.csv",header = T,fill = TRUE)
# colnames(a)
# #对Age进行划分
# b <- a[,c("Age","Severity")]
# d <- b |> group_by(Severity,Age) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# 
# d <-d %>% mutate(Age=case_when(40 < d$Age & d$Age <= 50 ~"40~50",
#                                50 < d$Age & d$Age <= 60 ~"50~60",
#                                60 < d$Age & d$Age <= 70 ~"60~70",
#                                70 < d$Age & d$Age <= 80 ~"70~80",
#                                80 < d$Age & d$Age <= 90 ~"80~90",
#                                90 < d$Age & d$Age <= 100 ~"90~100",
#                                d$Age >= 100 ~">100"))
# #统计每一个类有多少个样本
# m <- d |> group_by(Severity,Age) |> summarise(sum = sum(freq))
# m <- as.data.frame(m)
# colnames(m)
# #将数据转化为因子类型
# m$Severity <- as.factor(m$Severity)
# m$Age <- as.factor(m$Age)
# m <- arrange(m,Severity,Age)
# #绘制堆积的柱状图
# library(ggplot2)
# p1 <- ggplot(m,aes(x=Severity,y=sum,fill=Age))+
#   geom_bar(stat="identity")+
#   scale_fill_manual(values =  rev(c('#D73227', '#47ADCB', '#209073', "#EA866B",
#                                     '#2E4276','#717EA5','#d7a227'))) +
#   xlab("") +ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of Age") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# p1
# ggsave("p5.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 
# colnames(a)
# #对年龄的数据进行统计
# b <- a[,c("Sex","Severity")]
# d <- b |> group_by(Severity,Sex) |> summarise(freq=n()) 
# #对变量的名字进行更改
# d <- as.data.frame(d)
# d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
#                                    d$Severity==1~"AD",
#                                    d$Severity==2~"MCI"))
# 
# d <- d%>%mutate(Sex=case_when(d$Sex==1~"F",
#                               d$Sex==0~"M"))
# 
# #将数据转化为因子类型
# d$Severity <- as.factor(d$Severity)
# d$Sex <- as.factor(d$Sex)
# #绘制堆积的柱状图
# library(ggplot2)
# p2 <- ggplot(d,aes(x=Severity,y=freq,fill=Sex))+
#   geom_bar(stat="identity")+
#   scale_fill_manual(values =  rev(c('#D73227', '#47ADCB', '#209073', "#EA866B"))) +
#   xlab("") +ylab("Number")+
#   theme_bw()+#theme_classic()+
#   coord_flip()+
#   ggtitle("Category of Gender") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title = element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# p2
# # library(ggpubr)
# # p <- ggarrange(p1,p2,nrow = 2)
# 
# ggsave("p6.pdf", egg::set_panel_size(p2, width=unit(5, "in"), height=unit(3, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 


####表型数据的分析#####
#统计一下sum表格中每一个表型数据的缺失率
setwd("C:/Users/yujie/Desktop/datacollect")
sum <- read.csv("sum.csv",header = T,fill = T)
missing_rates <- sapply(sum, function(x) mean(is.na(x)))
# 将缺失率结果转化为数据框
missing_rates_df <- data.frame(Column = names(missing_rates), MissingRate = missing_rates)
print(missing_rates_df)
# Column MissingRate
# Brain_Region   Brain_Region  0.00000000
# GEO                     GEO  0.00000000
# GSM_number       GSM_number  0.00000000
# Length               Length  0.00000000
# LibraryLayout LibraryLayout  0.00000000
# Severity           Severity  0.00000000
# CERAD                 CERAD  0.77464789
# Braak_Stage     Braak_Stage  0.19718310
# Age                     Age  0.11907810
# Sex                     Sex  0.08706786
# Bank_Location Bank_Location  0.00000000






# KNN填充缺失值(真实的数据)
#将所有的数据放在一块进行填充，不对数据进行分类
library(DMwR2)
setwd("C:/Users/yujie/Desktop/datacollect")
sum <- read.csv("sum.csv",header = T,fill = T)
sum_compelete=sum[complete.cases(sum),]
nrow(sum_compelete)
sum <- sum[!duplicated(sum$GSM_number),]
colnames(sum)
#这里在对数据进行筛选的时候没有考虑结果变量
#没有考虑是AD还是是健康人这一列，去掉CERAD这一列的数据
#第一版的时候使用的情况，结果缺失率太高，会影响准确率，所以采用第二种情况
# sumKnn <- knnImputation(sum[,c(1,6:10)])
#第二种情况
sumKnn <- knnImputation(sum[,c(1,6,8:10)])
sumKnn <- round(sumKnn)
#检查一下样本行的顺序有没有变化
which(sumKnn$Brain_Region != sum$Brain_Region)
sumKnn <- cbind(sumKnn,sum[,c(2:5,11)])
summary(sumKnn)
anyNA(sumKnn)
write.csv(sumKnn,file = "sumknnImputation20240820.csv",row.names = F,quote = F)



###对数据重新进行填充
rm(list=ls())
library(VIM)
setwd("C:/Users/yujie/Desktop/datacollect")
sum <- read.csv("sum.csv",header = T,fill = T)
sum_compelete=sum[complete.cases(sum),]
nrow(sum_compelete)
sum_compelete <- sum_compelete[order(sum_compelete$GSM_number),]
#给sum_compelete添加APOE以及batch的信息
version <- read.csv("BackgroundInformationAPOE.csv",header = T)
colnames(version)
merged_file <- merge(sum_compelete, version[, c("GSM_number", "Sample","Batch","Phenotypedit")], by = "GSM_number")
merged_file <- merged_file[order(merged_file$GSM_number),]
merged_file <- na.omit(merged_file)
#检查一下完整的数据中，哪个因素对数据的影响最大
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
FPKM <- as.data.frame(t(FPKM))
FPKM <- FPKM[rownames(FPKM) %in% merged_file$GSM_number,]
FPKM <- FPKM[order(rownames(FPKM)),]
identical(rownames(FPKM),merged_file$GSM_number)
merged_file <- merged_file[,-c(3,1)]
colnames(merged_file)
# 将数据转化为因子
merged_file$Brain_Region <- as.factor(merged_file$Brain_Region)
merged_file$Batch <- factor(merged_file$Batch)
merged_file$Braak_Stage <- factor(merged_file$Braak_Stage)
merged_file$CERAD <- factor(merged_file$CERAD)
merged_file$Sex <- factor(merged_file$Sex)
merged_file$Bank_Location <- factor(merged_file$Bank_Location)
merged_file$LibraryLayout <- factor(merged_file$LibraryLayout)
merged_file$Severity <- factor(merged_file$Severity)
#因为APOE存在缺失值，所以要去除这几个样本
merged_file <- na.omit(merged_file)


# 将样本信息转化为因子
#使用多因素方差分析从batchFactor里面计算出来那个因素的影响比较大

# adonis2分析(多因素方差分析)
library(vegan)
# ?adonis()
# adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information)
# 如果你希望变量的顺序不影响结果，那么需要使用adonis2，并且设置参数by="margin"
colnames(merged_file)
#因为是随机置换，在未指定随机数种子时，每次执行的结果都会略有不同，但通常对结论没有影响。也可以如下设置随机数种子，则结果稳定
set.seed(1)
a <- adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Batch+Phenotypedit, 
             data = merged_file,
             permutations=999, by="margin",method = "euclidean")
#将多因素方差分析的结果绘制柱状图
a
b <- as.data.frame(a)
b <- na.omit(b)
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Brain_Region"~"Region",
                                   b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age",
                                   b$category=="Sex"~"Gender"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=0.2,angle = 0, size=4)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("MANOVA") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("MANOVA_Complete.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)



#通过上面的分析，我知道对数据影响最大的是脑区，因此在进行数据填充的时候k的个数设置为脑区的个数
#k指的是使用中心点附近的多少个点进行填充而不是说中心点的个数
rm(list=ls())
library(VIM)
setwd("C:/Users/yujie/Desktop/datacollect")
sum <- read.csv("sum.csv",header = T,fill = T)
sum_compelete=sum[complete.cases(sum),]
nrow(sum_compelete)
sum_compelete <- sum_compelete[order(sum_compelete$GSM_number),]
sum <- sum[!duplicated(sum$GSM_number),]
colnames(sum)
sum <- sum[,-c(2:5,7,10)]
#检查缺失值的个数
sapply(sum, function(x) sum(is.na(x)))
# Brain_Region      Severity   Braak_Stage           Age Bank_Location 
# 0             0           154            93             0 
#检查数据的类型，连续性数据和分类数据使用的是不同的方法
str(sum)
sum$Braak_Stage <- as.factor(sum$Braak_Stage)
# sum$Sex <- as.factor(sum$Sex)
sum$Brain_Region <- as.factor(sum$Brain_Region)
sum$Severity <- as.factor(sum$Severity)
sum$Bank_Location <- as.factor(sum$Bank_Location)

table(sum$Brain_Region)
#k的个数设置为8
set.seed(1)
sumKnn <- VIM::kNN(sum,trace = T,k = 15)
sumKnn$Age <- round(sumKnn$Age)
sumKnn$Age <- as.integer(sumKnn$Age)
#尝试了这种填充方法，得到的数据的分布很奇怪，因此最后选择不更改任何的参数
# sumKnn <- VIM::kNN(sum,trace = T, weights = "auto",addRF = TRUE)
# which(sumKnn$CERAD != sumKnn2$CERAD)
which(sumKnn$Brain_Region != sum$Brain_Region)
sum <- read.csv("sum.csv",header = T,fill = T)

#检查一下样本名是否对应
data <- cbind(sum[,c(2:5)],sumKnn[,1:5])
summary(data)
anyNA(data)
setwd("C:/Users/yujie/Desktop/datacollect")
write.csv(data,file = "sumkNN20240827_removeSex.csv",row.names = F,quote = F)




#检查填充的数据是否准确，手动挖掉一些点，然后看填充之后的结果和本来的结果之间的差距
#使用sum_compelete的数据
#
sum_compelete <- sum_compelete[order(sum_compelete$GSM_number),]
sum_compelete <- sum_compelete[!duplicated(sum_compelete$GSM_number),]
#remove sex
sum_compelete <- sum_compelete[,-c(2:5,7,10)]
# Braak_Stage     Braak_Stage  0.19718310
# Age                     Age  0.11907810
sum_compelete$Braak_Stage <- as.factor(sum_compelete$Braak_Stage)
# sum_compelete$Sex <- as.factor(sum_compelete$Sex)
sum_compelete$Brain_Region <- as.factor(sum_compelete$Brain_Region)
sum_compelete$Severity <- as.factor(sum_compelete$Severity)
sum_compelete$Bank_Location <- as.factor(sum_compelete$Bank_Location)
sum_compelete$Age <- as.numeric(sum_compelete$Age)
colnames(sum_compelete)
#随机对数据设置缺失值
set.seed(1)  # 设置种子以保证结果可重复
# 设置要添加缺失值的列及其对应的缺失值比例
missing_ratios <- list(Braak_Stage = 0.2, Age = 0.12 ) 
# k BA_accuracy      mse      mae
# TRUE 15   0.8857143 102.5238 7.380952
missing_ratios <- list(Braak_Stage = 0.05, Age = 0.05) 
# k BA_accuracy      mse      mae
# TRUE 15   0.8888889 61.44444 6.333333
missing_ratios <- list(Braak_Stage = 0.09, Age = 0.09 ) 
# k BA_accuracy     mse    mae
# TRUE 15      0.9375 62.6875 6.4375
missing_ratios <- list(Braak_Stage = 0.1, Age = 0.1) 
# k BA_accuracy  mse      mae
# TRUE 15   0.9444444 63.5 6.611111
missing_ratios <- list(Braak_Stage = 0.12, Age = 0.12) 
# k BA_accuracy      mse      mae
# TRUE 15   0.9047619 105.6667 7.666667
missing_ratios <- list(Braak_Stage = 0.15, Age = 0.15) 
# k BA_accuracy      mse      mae
# TRUE 15   0.9230769 105.9615 8.038462
missing_ratios <- list(Braak_Stage = 0.2, Age = 0.2 ) 
# k BA_accuracy      mse      mae
# TRUE 15   0.8857143 84.42857 7.285714
missing_ratios <- list(Braak_Stage = 0.25, Age = 0.25)
# k BA_accuracy      mse      mae
# TRUE 15   0.8636364 96.34091 7.931818
missing_ratios <- list(Braak_Stage = 0.30, Age = 0.30) 
# k BA_accuracy      mse mae
# TRUE 15   0.8679245 108.6792   8

df <- sum_compelete
# 对每列单独进行缺失值处理
for (col in names(missing_ratios)) {
  num_values <- nrow(df)  # 该列中的值数量
  num_missing <- round(missing_ratios[[col]] * num_values)  # 计算需要替换为 NA 的数量
  # 随机选择要设置为 NA 的行索引
  set.seed(1)
  missing_indices <- sample(1:num_values, num_missing)
  # 将选中的位置替换为 NA（针对指定列）
  df[missing_indices, col] <- NA
}
# 查看结果
# print(df)
#查看每一列数据的类型
# str(df)
#

# 定义不同的 k 值
k_values <- c(5,8,10,15, 20,25, 30,35,40,45,50)
k_values <- 15
results <- data.frame(k = integer(),
                      BA_accuracy = numeric(),
                      # Sex_accuracy = numeric(),
                      mse = numeric(),
                      mae = numeric(),
                      stringsAsFactors = FALSE)
# 找出 df 数据框中所有 NA 值的索引
na_indices <- which(is.na(df), arr.ind = TRUE)
# 遍历每个 k 值
for (k in k_values) {
  # 再次填充缺失值
  set.seed(1)
  df_Knn <- VIM::kNN(df, trace = F, k = k)
  df_Knn$Age <- round(df_Knn$Age, 0)
  
  # 提取 df_Knn 中对应的 NA 值的数据
  na_values_from_Knn <- df_Knn[na_indices]
  
  # 根据 col 的值分成两列
  Braak_Stage_Filled <- na_values_from_Knn[na_indices[, "col"] == 3]
  Age_Filled <- na_values_from_Knn[na_indices[, "col"] == 4]
  # Sex_Filled <- na_values_from_Knn[na_indices[, "col"] == 5]
  # 获取原始数据
  na_values_from_sum <- sum_compelete[na_indices]
  Braak_Stage_original <- na_values_from_sum[na_indices[, "col"] == 3]
  Age_original <- na_values_from_sum[na_indices[, "col"] == 4]
  # Sex_original <- na_values_from_sum[na_indices[, "col"] == 5]
  
  # 计算准确率
  BA_accuracy <- table(Braak_Stage_original == Braak_Stage_Filled)[2] / length(Braak_Stage_original)
  # Sex_accuracy <- table(Sex_original == Sex_Filled)[2] / length(Sex_original)
  
  # 计算Age的 MSE 和 MAE
  Age_Filled <- as.numeric(Age_Filled)
  Age_original <- as.numeric(Age_original)
  
  mse_age_filled <- mean((Age_Filled - Age_original)^2, na.rm = TRUE)
  mae_age_filled <- mean(abs(Age_Filled - Age_original), na.rm = TRUE)
  
  # 将结果存储到数据框中
  results <- rbind(results, data.frame(k = k,
                                       BA_accuracy = BA_accuracy,
                                       # Sex_accuracy = Sex_accuracy,
                                       mse = mse_age_filled,
                                       mae = mae_age_filled))
  
}

# 查看结果
print(results)

#remove sex
# k BA_accuracy       mse      mae
# TRUE    5   0.8571429  73.38095 7.095238
# TRUE1   8   0.8000000  86.52381 6.904762
# TRUE2  10   0.8000000  90.80952 7.095238
# TRUE3  15   0.8857143 102.52381 7.380952
# TRUE4  20   0.8857143  87.38095 7.380952
# TRUE5  25   0.8285714  86.95238 7.428571
# TRUE6  30   0.8571429  84.23810 7.380952
# TRUE7  35   0.8857143  91.00000 7.857143
# TRUE8  40   0.8571429  91.42857 7.523810
# TRUE9  45   0.7714286  87.80952 7.714286
# TRUE10 50   0.7142857  85.76190 7.285714



#检查填充的数据是否正确
# 连续性数据的检查
# data_combined <- data.frame(
#   value = c(sum$Age, sumKnn$Age),
#   dataset = rep(c("Original", "Imputed"), each = length(sum$Age)))
# 
# # 绘制密度曲线
# library(ggplot2)
# p1 <- ggplot(data_combined, aes(x = value, color = dataset)) +
#       geom_density(size = 1) +
#       labs(title = "Distribution of Age: Original vs. Imputed Data",
#           x = "Age",
#           y = "Density") +
#      theme_minimal() +
#      scale_color_manual(values = c("blue", "red"))
# p1
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# ggsave("density_Age.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 
# #分类数据的检查
# #方法1：
# # 对比原始数据与填充数据的类别频率
# table_original_braakStage <- table(sum$Braak_Stage)
# table_imputed_braakStage <- table(sumKnn$Braak_Stage)
# prop_original_braakStage <- prop.table(table_original_braakStage)
# prop_imputed_braakStag <- prop.table(table_imputed_braakStage)
# #age
# table_original_sex <- table(sum$Sex)
# table_imputed_sex <- table(sumKnn$Sex)
# prop_original_sex <- prop.table(table_original_sex)
# prop_imputed_sex <- prop.table(table_imputed_sex)
# # 可视化对比
# pdf("Two_Barplots.pdf", width = 8, height = 4)
# # 设置绘图区域为1行2列
# layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE), widths = c(2, 1),heights = c(1, 1))
# # 绘制第一张条形图
# barplot(rbind(prop_original_braakStage, prop_imputed_braakStag), beside = TRUE,
#         col = c("blue", "red"), legend = c("Original", "Imputed"),
#         main = "Braak Stage")
# 
# # 绘制第二张条形图（假设有另一组数据）
# barplot(rbind(prop_original_sex, prop_imputed_sex), beside = TRUE,
#         col = c("blue", "red"), legend = c("Original", "Imputed"),
#         main = "Gender")
# 
# # 关闭图形设备并保存文件
# dev.off()



# #方法2
# # Cohen's Kappa 一致性检验
# #Kappa 统计量：计算填充前后数据的Cohen's Kappa值。如果Kappa值较低，表示填充后的类别数据与原始数据之间一致性不高。
# library(irr)
# kappa_result <- kappa2(cbind(sum$Braak_Stage, sumKnn$Braak_Stage))
# print(kappa_result)
# #Kappa 值的范围从 -1 到 1，数值越接近 1 表示一致性越高
# # 显示 Cohen's Kappa 值
# library(ggplot2)
# ggplot() +
#   annotate("text", x = 1, y = 1, label = sprintf("Kappa = %.2f", kappa_result$value), size = 6) +
#   labs(title = "Cohen's Kappa Value") +
#   theme_void() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 





#方法3
#交叉验证，
#混淆矩阵：生成填充前后的混淆矩阵，观察哪些类别的变化最大。如果某些类别频繁被误填充为其他类别，这可能是一个警示信号。
# 混淆矩阵
library(caret)
conf_matrix <- confusionMatrix(sumKnn$Braak_Stage[!is.na(sum$Braak_Stage)], sum$Braak_Stage[!is.na(sum$Braak_Stage)])
print(conf_matrix)
conf_matrix_table <- as.table(conf_matrix)
# 将混淆矩阵转化为数据框
conf_matrix_df <- as.data.frame(as.table(conf_matrix_table))

# 绘制混淆矩阵图
ggplot(conf_matrix_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Reference Category", y = "Predicted Category", title = "Confusion Matrix") +
  theme_minimal()


####
# 提取混淆矩阵的各项数据
conf_matrix_df <- as.data.frame(conf_matrix$table)
conf_matrix_df$Type <- factor(
  with(conf_matrix_df, ifelse(Prediction == Reference, "True Positive", "False Positive")),
  levels = c("True Positive", "False Positive")
)

# 绘制条形图
ggplot(conf_matrix_df, aes(x = Prediction, fill = Type)) +
  geom_bar(position = "dodge") +
  labs(x = "Predicted Category", y = "Count", title = "Confusion Matrix by Type") +
  theme_minimal()


# 绘制堆叠条形图
ggplot(conf_matrix_df, aes(x = Reference, fill = Prediction)) +
  geom_bar(position = "fill") +
  labs(x = "True Category", y = "Proportion", title = "Confusion Matrix (Stacked Bar)") +
  theme_minimal()

# 计算精确度和召回率
library(PRROC)
roc_curve <- roc.curve(scores.class0 = conf_matrix$byClass[, "Sensitivity"], 
                       scores.class1 = conf_matrix$byClass[, "Specificity"], 
                       curve = TRUE)

# 绘制精确度与召回率图
plot(roc_curve)


#####进行基因表达layer的分析#####
# 将count转换为FPKM 
#计算FPKM的时候要求使用的是一个基因所有的外显子的长度
#而不是一个基因在基因组中所有的长度
#因为不清楚自己之前使用的是什么，所以要重新做一次分析
#经过检查之后发现自己之前使用的就是一个基因所有的外显子的长度
##将count值转换为FPKM的值#
#首先将featurecount产生的每一个文件通过循环加载进来
setwd("E:/AD_Patient/STARFeatureCount/04.STARcount")
filenames = list.files("E:/AD_Patient/STARFeatureCount/04.STARcount")
dir <- paste0("E:/AD_Patient/STARFeatureCount/04.STARcount/",filenames)
head(dir)

n = length(dir) 
merge.data <- read.table(file = dir[1],header = T)
for (i in 2:n){
  new.data = read.table(file = dir[i],header = T)
  merge.data = merge(merge.data,new.data,by = "Geneid", all = FALSE)
}

a1 <- colnames(merge.data)[2:30]
a2 <- colnames(merge.data)[31:143]
a3 <- colnames(merge.data)[144:812]
b1 <- substr(a1,46,55)
b2 <- substr(a2,63,72)
b3 <- substr(a3,46,55)
col_name <- c("Geneid",b1,b2,b3)
colnames(merge.data) <- col_name
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(merge.data,file = "ADPatientCount.csv",row.names=F,quote=F)




#导入所有reads的count值,开始计算FPKM值
gene_Length <- read.table(file = "HumanGeneLength.txt",header = T)
count <- read.csv(file="ADPatientCount.csv",header = T)
merge <- merge(count,gene_Length,by = 'Geneid') #匹配gene的count数据与gene的length数据进行合并
dim(merge)
count <- merge[,1:(dim(merge)[2]-1)]#此时已经去掉了基因长度这一列
gene_num <- dim(merge)[1]
sample_num <- dim(merge)[2]-2 #减去gene_name和gene_length列
#从第二列开始，将每个样本的count循环转化成FPKM
i <- 2
repeat{
  mapped_reads <- sum(merge[1:gene_num,i])#计算每个样本的mapped reads数
  FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])#计算FPKM值
  FPKM <- pmax(FPKM,0)#去掉矫正带来的负值,但是不可能会出现负值
  count = data.frame(count[1:gene_num,],FPKM)#添加到count表格后面i
  i <- i + 1
  
  if(i > sample_num+1){
    break
  }
}

#生成表格列名称
head(count)
head(colnames(count))
count_colname <- read.csv("ADPatientCount.csv",header = F,nrow = 1,as.is=TRUE)
FPKM_colname <- paste(count_colname[1,],"_FPKM",sep="")
head(FPKM_colname)
colname <- c(count_colname[1,],FPKM_colname[2:length(FPKM_colname)])
names(count) <- colname
head(count)
FPKM <- count[,c(1,(sample_num+2):(sample_num*2+1))]
FPKM <- cbind(FPKM[,1],round(FPKM[,2:812],3))
colnames(FPKM)[1] <- "Geneid"
#生成FPKM表
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(FPKM,"ADPatientFPKM.csv",row.names = FALSE, quote = FALSE)





# PCA图查看样本的分布情况 
rm(list = ls())
#使用sumknnImputation的信息对每个样本根据AD和Healthy进行编号
setwd("C:/Users/yujie/Desktop/datacollect")
SampleInformtion <- read.csv("sumkNN20240827_removeSex.csv",header = T)
#将APOE，Sample，Batch的信息添加到文件的后面,因为这个是之前经常使用到的编号，所以需要保持一致
version <- read.csv("BackgroundInformationAPOE.csv",header = T)
colnames(version)
merged_file <- merge(SampleInformtion, version[, c("GSM_number", "Sample","Batch","Phenotypedit")], by = "GSM_number")
write.csv(merged_file,file = "BackgroundInformationkNN20240827_k_15.csv",row.names = F,quote = F)


#检查一下自己新填充的数据和之前填充好的数据之间存在多大的差别
setwd("C:/Users/yujie/Desktop/datacollect")
old_version <- read.csv("BackgroundInformationAPOE.csv",header = T)
old_version <- old_version[order(old_version$GSM_number),]
colnames(old_version)
new_version <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
new_version <- new_version[order(new_version$GSM_number),]
colnames(new_version)
identical(old_version$GSM_number,new_version$GSM_number)
#检查不一致的的填充的数据所在数据框中的位置,以及个数
length(which(old_version$Braak_Stage != new_version$Braak_Stage))
# 111
length(which(old_version$Age != new_version$Age))
# 93
# length(which(old_version$Sex != new_version$Sex))
# 42
length(which(old_version$Bank_Location != new_version$Bank_Location))
# 0


#使用FPKM计算PCA查看各个样本之间的距离
# 用FPKM计算样本之间的距离,加载描述样本的信息
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
# colnames(FPKM) <- gsub("_FPKM","",colnames(FPKM))
FPKM[1:5,1:5]
dim(FPKM)
FPKM <- as.matrix(FPKM)
FPKM[FPKM<=1] <- 1
FPKM <- as.data.frame(FPKM)
FPKM <- na.omit(FPKM)

# 在R中删除矩阵中含有0的行
# #方法一：
# A = sapply(1:nrow(A),function(x) if(all(A[x,])!=0) A[x,])
# #方法二：
# x[!as.logical(rowSums(dat==0)), ]
# #方法三：
# dat[dat==0] <- NA
# na.omit(dat)


FPKM <- log10(FPKM)
head(colnames(FPKM))
#加载背景信息
setwd("C:/Users/yujie/Desktop/datacollect")
SampleInformtion <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")

FPKM <- as.data.frame(t(FPKM))
FPKM$GSM_number <- rownames(FPKM)
sum <- merge(FPKM,SampleInformtion,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
library(factoextra)
library(ggpubr)
library(gmodels)
pca.info <- prcomp(sum)
head(pca.info)
a <- summary(pca.info) #能看到每个PC具体的解释比例
b <- as.data.frame(a[["importance"]])
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
# 构建数据框
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
colnames(pca.data)[1] <- "Sample"
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="confidence")
b <- merge(SampleInformtion,pca.data,by="Sample")
b$Brain_Region <- as.factor(b$Brain_Region)
#将数字0转换为Healthy,1转换为AD,2用MCI表示
library(tidyverse)
pca.data <- pca.data%>%mutate(Type=case_when(pca.data$Type==0~"Healthy",
                                             pca.data$Type==1~"AD",
                                             pca.data$Type==2~"MCI"))
str(pca.data$Type)
pca.data$Type <- as.factor(pca.data$Type)
library(ggforce)
# shape =Brain_Region
p1 <- ggplot(data=b,aes(x=PC1,y=PC2,color=Type))+
  geom_point(size=3)+
  #scale_x_continuous(limits = c(-800000, 6000000))+
  #scale_y_continuous(limits = c(-1700000,900000))+
  scale_color_manual(values = c("#3366CC","#911F27","#009933"))+
  ggtitle("All Sample PCA Plot")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 6,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",25,"%"),
       y=paste0("PCA2 ",23,"%"))
# geom_mark_ellipse(aes(color = Type),
#                 expand = unit(3, "mm"),
#                 stat = "identity",
#                 position = "identity")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AllSamplePCAFPKM.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)


# 做所有样本的多因素方差分析
# 在进行多因素方差分析之前，需要将sample的FPKM和sample的样本信息对应
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
FPKM <- as.data.frame(t(FPKM))
FPKM <- FPKM[order(rownames(FPKM)),]
# 将样本信息转化为因子
#使用多因素方差分析从batchFactor里面计算出来那个因素的影响比较大
#使用用KNN填充之后的数据
setwd("C:/Users/yujie/Desktop/datacollect")
sumknnImputation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sumknnImputation <- sumknnImputation[order(sumknnImputation$GSM_number),]
identical(rownames(FPKM),sumknnImputation$GSM_number)
sumknnImputation <- sumknnImputation[,-c(1,2)]
colnames(sumknnImputation)
# 将数据转化为因子
sumknnImputation$Brain_Region <- as.factor(sumknnImputation$Brain_Region)
sumknnImputation$Batch <- factor(sumknnImputation$Batch)
sumknnImputation$Braak_Stage <- factor(sumknnImputation$Braak_Stage)
# sumknnImputation$Sex <- factor(sumknnImputation$Sex)
sumknnImputation$Bank_Location <- factor(sumknnImputation$Bank_Location)
sumknnImputation$LibraryLayout <- factor(sumknnImputation$LibraryLayout)
sumknnImputation$Severity <- factor(sumknnImputation$Severity)


# adonis2分析(多因素方差分析)
library(vegan)
# ?adonis()
# adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information)
# 如果你希望变量的顺序不影响结果，那么需要使用adonis2，并且设置参数by="margin"
colnames(sumknnImputation)
#因为是随机置换，在未指定随机数种子时，每次执行的结果都会略有不同，但通常对结论没有影响。也可以如下设置随机数种子，则结果稳定
set.seed(1)
a <- adonis2(FPKM ~ Severity+Braak_Stage+Brain_Region+Age+Batch+Bank_Location, 
             data = sumknnImputation,
             permutations=999, by="margin",method = "gower")
a
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = FPKM ~ Severity + Braak_Stage + Brain_Region + Age + Batch + Bank_Location, data = sumknnImputation, permutations = 999, method = "gower", by = "margin")
# Df SumOfSqs      R2        F Pr(>F)    
# Severity        2  0.00327 0.00131   1.3147  0.185    
# Braak_Stage     3  0.00550 0.00220   1.4757  0.075 .  
# Brain_Region    1  0.14823 0.05920 119.2054  0.001 ***
#   Age             1  0.00154 0.00061   1.2369  0.227    
# Batch           3  0.12295 0.04911  32.9594  0.001 ***
#   Bank_Location   0  0.00000 0.00000     -Inf           
# Residual      758  0.94253 0.37644                    
# Total         780  2.50380 1.00000                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = FPKM ~ Severity + Braak_Stage + Brain_Region + Age + Batch + Bank_Location, data = sumknnImputation, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2      F Pr(>F)  
# Severity        2  1.3970e+11 0.00048 1.5197  0.228  
# Braak_Stage     3  7.5419e+10 0.00026 0.5470  0.674  
# Brain_Region    1  2.1374e+11 0.00074 4.6502  0.023 *
#   Age             1  7.7408e+09 0.00003 0.1684  0.770  
# Batch           3  3.5053e+11 0.00121 2.5422  0.057 .
# Bank_Location   0 -3.3387e+07 0.00000   -Inf         
# Residual      758  3.4840e+13 0.12013                
# Total         780  2.9003e+14 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 


#将多因素方差分析的结果绘制柱状图
b <- as.data.frame(a)
b <- na.omit(b)
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Brain_Region"~"Region",
                                   b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=0.2,angle = 0, size=4)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("All sample MANOVA") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("all_samples_MANOVA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)

#对不同脑区的数据量进行统计
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sample <- read.csv("SampleDistribution.csv",header = T)
colnames(sample)
p <- ggplot(sample,aes(x=reorder(Region,-Number),y=Number,fill=Region))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Number),vjust=-0.2,angle = 0, size=3)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("Region") +ylab("Number")+
  theme_bw()+#theme_classic()+
  ggtitle("All Sample") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("SampleDistribution.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
#将上面的PCA图和多因素方差分析的图和样本脑区分布的图，共3张图放在一个pdf文件中


#####分脑区对AD病人的数据进行分析，首先分析颞叶的数据#####
# 首先分析颞叶的数据，使用FPKM计算一下样本之间的距离
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
#加载背景信息
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
#取出来颞叶的数据
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
#根据BackgroundInformation提取出来想要看的样本
FPKM <- FPKM[,colnames(FPKM) %in% BackgroundInformation$GSM_number]
#使用AD和healthy作为分组的标准
FPKM <- as.data.frame(t(FPKM))
FPKM$GSM_number <- rownames(FPKM)
sum <- merge(FPKM,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
# 使用层次聚类，看一下样本之间的距离
# dist_mat <- dist(sum)
# clustering <- hclust(dist_mat, method = "complete")
# plot(clustering)

# 颞叶的数据做PCA检查一下样本之间的差异性
library(gmodels)
pca.info <- fast.prcomp(sum)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
# 构建数据框
# 构建数据框
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="confidence")

#将数字0转换为Healthy,1转换为AD,2用MCI表示
library(dplyr)
pca.data <- pca.data%>%mutate(Type=case_when(pca.data$Type==0~"Healthy",
                                             pca.data$Type==1~"AD",
                                             pca.data$Type==2~"MCI"))
str(pca.data$Type)
pca.data$Type <- as.factor(pca.data$Type)
library(ggforce)
p1 <- ggplot(data=pca.data,aes(x=PC1,y=PC2,color=Type))+
  geom_point(size=3)+
  scale_x_continuous(limits = c(-90000, 3000000))+
  scale_y_continuous(limits = c(-300000,300000))+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
  ggtitle("Temporal Lobe PCA Plot")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  geom_mark_ellipse(aes(color = Type),
                    expand = unit(3, "mm"),
                    stat = "identity",
                    position = "identity")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240820")
ggsave("TemporalLobePCAFPKM.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
#并不能看出来有很明显的差异
# 所以决定使用一下DESeq2里面的count值计算一下样本之间的差异性
# 首选选出来所需要的数据
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")


#加载FPKM的数据
myFPKM <- read.csv("ADPatientFPKM.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- subset(myFPKM,rownames(myFPKM) %in% BackgroundInformation$GSM_number)
myFPKM <- myFPKM[order(rownames(myFPKM)),]

#对颞叶的数据要进行多因素方差分析
setwd("C:/Users/yujie/Desktop/datacollect")
sumKnn <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sumKnn <- subset(sumKnn,sumKnn$GSM_number %in% BackgroundInformation$GSM_number)
sumKnn <- sumKnn[order(sumKnn$GSM_number),]

#检查样本表型数据的顺序和样本FPKM值的顺序是否一致
identical(rownames(myFPKM),sumKnn$GSM_number)
sample_information <- sumKnn
colnames(sample_information)
# 将表型类型转化为因子
sample_information$Brain_Region <- as.factor(sample_information$Brain_Region)
# sample_information$CERAD <- factor(sample_information$CERAD)
sample_information$Braak_Stage <- factor(sample_information$Braak_Stage)
# sample_information$Sex <- factor(sample_information$Sex)
sample_information$Bank_Location <- factor(sample_information$Bank_Location)
sample_information$LibraryLayout <- factor(sample_information$LibraryLayout)
sample_information$Severity <- factor(sample_information$Severity)
sample_information$Batch <- factor(sample_information$Batch)

library(vegan)
# ?adonis()
# adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information)
# 如果你希望变量的顺序不影响结果，那么需要使用adonis2，并且设置参数by="margin"
#使用count数据
# a1 <- adonis2(mycount ~ Severity+CERAD+Braak_Stage+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information,
#         permutations=999, by="margin",method = "euclidean")

#使用FPKM数据
set.seed(1)
a <- adonis2(myFPKM ~ Severity+Braak_Stage+Age+Batch+Bank_Location, data = sample_information,
             permutations=999, by="margin",method = "euclidean")
a
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = myFPKM ~ Severity + Braak_Stage + Age + Batch + Bank_Location, data = sample_information, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2       F Pr(>F)    
# Severity        1  7.9662e+09 0.00078  0.9813  0.336    
# Braak_Stage     3  7.4155e+10 0.00722  3.0449  0.033 *  
#   Age             1  1.4129e+10 0.00138  1.7405  0.158    
# Batch           2  2.0771e+11 0.02022 12.7934  0.001 ***
#   Bank_Location   0 -1.5791e+06 0.00000    -Inf           
# Residual      345  2.8007e+12 0.27259                   
# Total         353  1.0274e+13 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#将多因素方差分析的结果绘制柱状图
b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age",
                                   # b$category=="CERAD"~"CERAD",
                                   b$category=="Sex"~"Gender"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=0.3,angle = 0, size=4)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("GE Raw MANOVA") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_Raw_MANOVA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)






#使用DEseq2做PCA和差异表达分析
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount$GSM_number <- rownames(mycount)
#加载背景信息
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
#只取出来颞叶的数据
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) #构建condition因子，作为基因差异表达的自变量（即实验组和对照组的对比）。i和j分别为control和exp的个数
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  #将database中各个实验名称加上condition标签，即说明是control还是exp
coldata
library(DESeq2)  #加载DESeq2包
#构建dds及利用results函数得到最终结果
max(mycount)
# which(is.na(mycount))

dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  #countData用于说明数据来源，colData用于说明不同组数据的实验操作类型，design用于声明自变量，即谁和谁进行对比
vst <- vst(dds, blind=T)#在vst处理大样本的时候表现的很好，rlog处理小样本的时候表现的很好
head(assay(vst), 3)
library(ggforce)
plotPCA(vst)

#下面将自己的PCA图画上标签,所以首先需要做的就是将PCA画图的data导出来
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
#将Severity转换为AD和Healthy
colnames(p1data)
p1data <- p1data%>%mutate(Severity=case_when(p1data$Severity==0~"Healthy",
                                             p1data$Severity==1~"AD",
                                             p1data$Severity==2~"MCI"))
###将部分的表型数据转换为因子，方便用颜色显示出来
p1data$Severity <- as.factor(p1data$Severity)
p1data$Batch <- as.factor(p1data$Batch)
# p1data$Sex <- as.factor(p1data$Sex)
p1data$Bank_Location <- as.factor(p1data$Bank_Location)
library(ggrepel)
library(ggplot2)
library(ggforce)
# ,shape =Bank_Location
p1 <- ggplot(data=p1data,aes(x=PC1,y=PC2,color=Severity,shape =Batch))+
  geom_point(size=3)+
  scale_x_continuous(limits = c(-80, 20))+
  scale_y_continuous(limits = c(-20,40))+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
  ggtitle("GE Raw PCA Plot")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",32,"%"),
       y=paste0("PCA2 ",17,"%"))
# geom_mark_ellipse(aes(color = group),
#               expand = unit(3, "mm"),
#               stat = "identity",
#               position = "identity")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_Raw_PCACount.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)








#将所有的样本合并在一起去除批次效应#
###使用ComBat函数来去除批次效应，使用的data是FPKM
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
BackgroundInformation$Severity <- as.factor(BackgroundInformation$Severity)
BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
# BackgroundInformation$Sex <- as.factor(BackgroundInformation$Sex)

#加载count的数据
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount <- subset(mycount,rownames(mycount) %in% BackgroundInformation$GSM_number)
mycount <- mycount[order(rownames(mycount)),]
identical(rownames(mycount),BackgroundInformation$GSM_number)
mycount <- as.matrix(t(mycount))


# #加载FPKM的数据
# myFPKM <- read.csv("ADPatientFPKM.csv",header = T,row.names = 1)
# myFPKM <- as.data.frame(t(myFPKM))
# myFPKM <- subset(myFPKM,rownames(myFPKM) %in% BackgroundInformation$GSM_number)
# myFPKM <- myFPKM[order(rownames(myFPKM)),]
# #将FPKM转化为combat去除批次需要的格式
# identical(rownames(myFPKM),BackgroundInformation$GSM_number)
# myFPKM <- as.data.frame(t(myFPKM))
# #对FPKM进行过滤，因为当使用所有的数据的时候，FPKM会出现负值
# #所以对FPKM进行ID转换
# a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myFPKM))
# head(a)
# rownames(myFPKM) <- a
# myFPKM$ENSEMBL <- rownames(myFPKM)
# ####将所有的差异基因做GO和KEGG
# library(clusterProfiler)
# library(org.Hs.eg.db)
# df1 <- bitr(a, fromType = "ENSEMBL",
#             toType = c("SYMBOL"),
#             OrgDb = org.Hs.eg.db)
# 
# colnames(df1)
# 
# filterFPKM <- merge(df1,myFPKM,by="ENSEMBL")
# #去重
# filterFPKM <- filterFPKM[!duplicated(filterFPKM$SYMBOL),]
# rownames(filterFPKM) <- filterFPKM$SYMBOL
# filterFPKM <- filterFPKM[,-c(1,2)]

#使用combat函数校正批次效应
library(sva)
# batch <- BackgroundInformation$Batch
# #full mode包含的是校正的变量和感兴趣的变量，full mode用mod表示
# #null mode包含的是校正的变量，但不包括感兴趣的变量，null mode用mod0表示
# #这里采用null mode,只包含校正的变量
# modcombat = model.matrix(~1, data=BackgroundInformation)
# combat_edata <- ComBat(dat=filterFPKM, batch=batch, mod=modcombat)
# #将去了批次之后的数据转换为矩阵
# combat_edata <- as.data.frame(combat_edata)
#去除了批次之后的数据，大概有2000个基因的FPKM出现了负值

#所以在去除批次效应的时候，采用count值
# 输入的mycount的matrix的数据的格式应该是行名是基因，列名是样品的名字
# covar_mod里面存放的是想要看到的生物学变异，但是性别和年龄并不是我想要看到的生物学变异
# 使用ComBat_seq去除批次效应
#当把bank_location放入到去除批次的模型里面的时候，就会出现报错
# # At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq
# dim(mycount)
# BackgroundInformation$Severity <- as.factor(BackgroundInformation$Severity)
# BackgroundInformation$Sex <- as.factor(BackgroundInformation$Sex)
# BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
# #计算这些协变量之间的相关性，防止出现共线性问题
# str(BackgroundInformation[, c( "Age", "Sex")])
# # 将因子变量转换为虚拟变量
# # 使用 model.matrix 函数来处理因子型变量，并排除掉常数列（Intercept）
design_matrix <- model.matrix(~ Age -1, data = BackgroundInformation)
cor_matrix <- cor(design_matrix)
print(cor_matrix)
# #分别计算了age与性别为0和性别为1的相关性
# # Age       Sex0       Sex1
# # Age   1.0000000 -0.2232475  0.2232475
# # Sex0 -0.2232475  1.0000000 -1.0000000
# # Sex1  0.2232475 -1.0000000  1.0000000
# # 计算分类变量之间的相关性，使用的是卡方检验
# # 创建 Sex 和 Batch 的交叉表
contingency_table <- table(BackgroundInformation$Sex, BackgroundInformation$Batch)
# 进行卡方检验
chi_test_result <- chisq.test(contingency_table)
# 输出检验结果
chi_test_result
# # Pearson's Chi-squared test
# # 
# # data:  contingency_table
# X-squared = 16.959, df = 3, p-value = 0.0007206
# 结果的解释
# 如果卡方检验的 p 值非常小（例如小于 0.05），则表示 Sex 和 Batch 之间存在显著关联，说明可能存在混淆关系。
# 如果你要检查分类变量（如 Sex 和 Batch）与连续变量（如 Age）之间的关系，可以使用方差分析（ANOVA）或 Kruskal-Wallis 检验（非参数检验）。

# 检查 Batch 和 Age 之间的关系（ANOVA）
anova_result <- aov(Age ~ Batch, data = BackgroundInformation)
# 输出 ANOVA 结果
summary(anova_result)
# # 如果 p 值较小，说明 Age 和 Batch 之间存在显著关联。
# # 应用 ComBat_seq 函数
# corrected_data <- ComBat_seq(
#   counts = mycount,
#   batch = BackgroundInformation$Batch,  # 批次信息
#   group = BackgroundInformation$Severity,  # 生物学条件
#   covar_mod = design_matrix,  # 设计矩阵
#   full_mod = TRUE,  # 包括生物学条件
#   shrink = FALSE  # 根据需要决定是否应用收缩
# )


# full_mod=TRUE：包括生物学条件（如实验组别），同时进行批次效应校正。适用于需要同时校正批次效应和生物学差异的情况。
# full_mod=FALSE：不包括生物学条件，仅进行批次效应校正。适用于只需要去除批次效应的情况。
adjust_counts <- ComBat_seq(mycount, batch=BackgroundInformation$Batch,group=BackgroundInformation$Severity,covar_mod=NULL,full_mod=FALSE)
# design_matrix <- model.matrix(~ Age + Sex, data = BackgroundInformation)
# fit <- lmFit(log(adjust_counts+1), design_matrix)
head(adjust_counts)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(adjust_counts,file = "GE_Count_RemoveBatch.csv",quote = F)
mycount <- as.data.frame(t(adjust_counts))
mycount$GSM_number <- rownames(mycount)


#加载尝试之前去除了批次的数据
# setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
# test <- read.csv("Temporal_LobeCount_filterAdjustCounts.csv",header = T,row.names = 1)
# mycount <- as.data.frame(t(test))
# mycount$GSM_number <- rownames(mycount)
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) #构建condition因子，作为基因差异表达的自变量（即实验组和对照组的对比）。i和j分别为control和exp的个数
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  #将database中各个实验名称加上condition标签，即说明是control还是exp
coldata
library(DESeq2)  #加载DESeq2包
#构建dds及利用results函数得到最终结果
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  #countData用于说明数据来源，colData用于说明不同组数据的实验操作类型，design用于声明自变量，即谁和谁进行对比
vst <- vst(dds, blind=T)
library(ggforce)
plotPCA(vst)
#下面将自己的PCA图画上标签,所以首先需要做的就是将PCA画图的data导出来
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"

p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
p1data$Severity <- as.factor(p1data$Severity)

p1data <- p1data%>%mutate(Severity=case_when(p1data$Severity==0~"Healthy",
                                             p1data$Severity==1~"AD",
                                             p1data$Severity==2~"MCI"))
p1data$Batch <- as.factor(p1data$Batch)
p1data$Sex <- as.factor(p1data$Sex)
# p1data <- p1data%>%mutate(Sex=case_when(p1data$Sex==0~"M",
#                                         p1data$Sex==1~"F"))
p1data$Bank_Location <- as.factor(p1data$Bank_Location)
p1data <- p1data%>%mutate(Bank_Location=case_when(p1data$Bank_Location==0~"America",
                                                  p1data$Bank_Location==3~"Australia"))

p1data$Braak_Stage <- as.factor(p1data$Braak_Stage)
p1data <- p1data%>%mutate(Braak_Stage=case_when(p1data$Braak_Stage==0~"0",
                                                p1data$Braak_Stage==1~"transentorhinal_stage",
                                                p1data$Braak_Stage==2~"limbic_stage",
                                                p1data$Braak_Stage==3~"isocortical_stage"))

# p1data$CERAD <- as.factor(p1data$CERAD)
# p1data <- p1data%>%mutate(CERAD=case_when(p1data$CERAD==0~"C0",
#                                           p1data$CERAD==1~"C1",
#                                           p1data$CERAD==2~"C2",
#                                           p1data$CERAD==3~"C3"))
p1data$Age <- as.numeric(p1data$Age)
library(ggrepel)
library(ggplot2)
library(ggforce)
# display.brewer. all#展示所有的配色

# col2<-brewer.pal(5, "YlOrRd") #选色盘的几种
# pal<-colorRampPalette(col2)
# colors2<-pal(100)  #分成多少种颜色
# # mycolor3<-brewer.pal(9, "YlOrRd")


# ,shape=Bank_Location,
p1 <- ggplot(data=p1data,aes(x=PC1,y=PC2,color=Severity,shape=Batch))+
  geom_point(size=3)+
  scale_x_continuous(limits = c(-30, 30))+
  scale_y_continuous(limits = c(-25,30))+
  ggtitle("GE RemoveBatch PCA")+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position = "right")+
  labs(x=paste0("PCA1 ",16,"%"),
       y=paste0("PCA2 ",10,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_RemoveBatch_PCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)





#将校正之后的数据计算FPKM，然后计算多因素方差分析
dim(adjust_counts)
count <- as.data.frame(adjust_counts)
count$Geneid <- rownames(count)
#导入所有reads的count值,开始计算FPKM值
gene_Length <- read.table(file = "HumanGeneLength.txt",header = T)
merge <- merge(count,gene_Length,by = 'Geneid') #匹配gene的count数据与gene的length数据进行合并
dim(merge)
count <- merge[,1:(dim(merge)[2]-1)]#此时已经去掉了基因长度这一列
gene_num <- dim(merge)[1]
sample_num <- dim(merge)[2]-2 #减去gene_name和gene_length列
#从第二列开始，将每个样本的count循环转化成FPKM
i <- 2
repeat{
  mapped_reads <- sum(merge[1:gene_num,i])#计算每个样本的mapped reads数
  FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])#计算FPKM值
  FPKM <- pmax(FPKM,0)#去掉矫正带来的负值,但是不可能会出现负值
  count = data.frame(count[1:gene_num,],FPKM)#添加到count表格后面i
  i <- i + 1
  
  if(i > sample_num+1){
    break
  }
}

#生成表格列名称
head(count)
head(colnames(count))
FPKM_colname <- colnames(merge[,1:(dim(merge)[2]-1)])
head(FPKM_colname)
colname <- c(colnames(merge[,1:(dim(merge)[2]-1)]),FPKM_colname[2:length(FPKM_colname)])
names(count) <- colname
head(count)
FPKM <- count[,c(1,(sample_num+2):(sample_num*2+1))]
FPKM <- cbind(FPKM[,1],round(FPKM[,2:355],3))
colnames(FPKM)[1] <- "Geneid"
#生成FPKM表
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(FPKM,"GE_RemoveBatch_FPKM.csv",row.names = FALSE, quote = FALSE)




#加载FPKM的数据
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myFPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- myFPKM[order(rownames(myFPKM)),]



#对颞叶的数据要进行多因素方差分析
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]



# 检查一下是否顺序就是一致的
identical(rownames(myFPKM),BackgroundInformation$GSM_number)
BackgroundInformation$Brain_Region <- as.factor(BackgroundInformation$Brain_Region)
# BackgroundInformation$CERAD <- factor(BackgroundInformation$CERAD)
BackgroundInformation$Braak_Stage <- factor(BackgroundInformation$Braak_Stage)
# BackgroundInformation$Sex <- factor(BackgroundInformation$Sex)
BackgroundInformation$Bank_Location <- factor(BackgroundInformation$Bank_Location)
BackgroundInformation$LibraryLayout <- factor(BackgroundInformation$LibraryLayout)
BackgroundInformation$Severity <- factor(BackgroundInformation$Severity)
BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
str(BackgroundInformation)
colnames(BackgroundInformation)
library(vegan)
# ?adonis()
# adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information)
# 如果你希望变量的顺序不影响结果，那么需要使用adonis2，并且设置参数by="margin"
set.seed(1)
a <- adonis2(myFPKM ~ Severity+Braak_Stage+Age+Bank_Location+Batch, data = BackgroundInformation,
             permutations=999, by="margin",method = "euclidean")
a
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = myFPKM ~ Severity + Braak_Stage + Age + Bank_Location + Batch, data = BackgroundInformation, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2      F Pr(>F)  
# Severity        1  6.8525e+09 0.00875 3.0984  0.069 .
# Braak_Stage     3  4.2270e+09 0.00539 0.6371  0.633  
# Age             1  5.5967e+09 0.00714 2.5306  0.105  
# Bank_Location   0 -1.5581e+05 0.00000   -Inf         
# Batch           2  3.3566e+09 0.00428 0.7589  0.470  
# Residual      345  7.6301e+11 0.97376                
# Total         353  7.8357e+11 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="Sex"~"Gender",
                                   # b$category=="CERAD"~"CERAD",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("GE Remove Batch MANOVA") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_RemoveBatch_MANOVA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 9, height = 7, units = 'in', dpi = 600)


#使用去除了批次效应的count数据使用线性回归模型校正性别和年龄的影响
# rm(list=ls())
# setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
# mycount <- read.csv("Temporal_LobeRawCount_RemoveBatch20220831.csv",header = T,row.names = 1)
# mycount <- as.data.frame(t(mycount))
# mycount <- mycount[order(rownames(mycount)),]
# 
# #将背景信息加载进来进行校正
# BackgroundInformation <- read.csv("BackgroundInformationAPOE.csv",header = T)
# BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
# BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
# 
# #判断两个数据对应的行是否是一样的样本呢
# identical(rownames(mycount),BackgroundInformation$GSM_number)
# adjustFactor <- BackgroundInformation[,-c(1,2,8,9,10)]
# #将APOE的数据转换为数字
# adjustFactor <- adjustFactor%>%mutate(Phenotypedit=case_when(adjustFactor$Phenotypedit=="noE4"~"0",
#                                                              adjustFactor$Phenotypedit=="E4carrier"~"1",
#                                                              adjustFactor$Phenotypedit=="E4/4"~"2"))
# 
# colnames(adjustFactor)
# #使用线性回归进行校正
# m <- lm(mycount~CERAD+Braak_Stage+Age+Sex+Bank_Location+Phenotypedit,data =adjustFactor )
# 




#使用去除了批次之后的数据进行基因的表达差异分析
rm(list=ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
mycount <- read.csv("GE_Count_RemoveBatch.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount$GSM_number <- rownames(mycount)

setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]

sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
mycount <- as.data.frame(t(sum))
#对基因进行ID转化，去除那些没有办法转换的基因
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(mycount))
rownames(mycount) <- a
####将所有的差异基因做GO和KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
mycount$ENSEMBL <- rownames(mycount)
df <- merge(df1,mycount,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
mycount <- df[,-c(1:2)]

condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) #构建condition因子，作为基因差异表达的自变量（即实验组和对照组的对比）。i和j分别为control和exp的个数
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  #将database中各个实验名称加上condition标签，即说明是control还是exp
coldata
library(DESeq2)  #加载DESeq2包
#构建dds及利用results函数得到最终结果
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  #countData用于说明数据来源，colData用于说明不同组数据的实验操作类型，design用于声明自变量，即谁和谁进行对比
#过滤低表达的基因
library(edgeR)
# keep <- rowSums(cpm(mycount)>0.5) >= 10#海马
keep <- rowSums(cpm(mycount)>0.5) >= 20#颞叶
#rowSums函数计算矩阵每一行的加和并生成一个新的向量
#cpm函数是edgeR中，提供了一种名为CPM的定量方式，全称为count-per-millon。
#cpm <- apply(count ,2, function(x) { x/sum(x)*1000000 })
#原始的表达量除以该样本表达量的总和，在乘以一百万就得到了CPM值 。从公式可以看出， CPM其实就是相对丰度，只不过考虑到测序的reads总量很多，所以总的reads数目以百万为单位。
# data=data[which(apply(data,1,function(x){return(sum(x>10))})>ncol(data)*0.25),]
table(keep)
# keep
# FALSE  TRUE 
# 18777 18219  
dds <- dds[keep,]
dds_norm <- DESeq(dds)#dds标准化,不剔除outliers; 与cookscutoff结果相同,设置参数：,minReplicatesForReplace = Inf
dds_norm

resultsNames(dds_norm)
# res1 <- results(dds_norm,contrast = c("condition","KO","WT"),cooksCutoff = FALSE)#指的是将KO和WT进行比较lgfc<0，KO的量少，,cooksCutoff = FALSE
res1 <- results(dds_norm,contrast = c("condition","1","0"))#指的是将KO和WT进行比较lgfc<0，KO的量少，,cooksCutoff = FALSE
head(res1)
res1 <- res1[order(res1$padj,decreasing = F),]
head(res1)
sum(res1$padj<0.01, na.rm = TRUE)
summary(res1,alpha = 0.01)

res1 <-na.omit(res1)
res1 <- res1[!duplicated(rownames(res1)),]#去重，会保留第一个，而去掉后面出现的重复
dim(res1)
res1 <- as.data.frame(res1)
table(res1$padj<0.05 )
# FALSE  TRUE 
# 6893 11326 
table(res1$padj<0.01 )
# FALSE  TRUE 
# 8756  9463
AD_Healthyres <- res1
head(AD_Healthyres)
dim(AD_Healthyres)
# [1] 18219     6
AD_Healthyres$change <- ifelse(AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >= 0.585, 
                               ifelse(AD_Healthyres$log2FoldChange > 0.585 ,'Up','Down'),'Stable')

table(AD_Healthyres$change)
# Down Stable     Up 
# 471  17221    527 
#差异基因的筛选标准是padj < 0.05
table(AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange)>0.585)    
# FALSE  TRUE 
# 17221   998
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres,file = "GE_RemoveBatch_DEseq2.csv")
# AD_Healthyres$label <- ifelse(AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >=1,
#                          as.character(rownames(AD_Healthyres)), "")

AD_Healthyres_def <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585 )

dim(AD_Healthyres_def)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_def,file = "GE_RemoveBatch_DEseq2_Significant.csv")
####将所有的差异基因做GO和KEGG
df1 <- bitr(rownames(AD_Healthyres_def), fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
GOgene=as.character(df1$ENTREZID)#得出所有的ENTREZID
GOgene <- na.omit(GOgene)#得出所有的ENTREZID，就可以进行GO和kegg的分析
AD_Healthyres_GOplot <- enrichGO(GOgene, 
                                 OrgDb = org.Hs.eg.db, 
                                 ont='All',
                                 pAdjustMethod = 'BH',
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.1,
                                 keyType = 'ENTREZID')
##得到的结果是一个大的enrichresult
AD_Healthyres_GOplot_genelist<-setReadable(AD_Healthyres_GOplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_GOplot_genelist <- as.data.frame(AD_Healthyres_GOplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_GOplot_genelist, file='GE_GO_genelist.csv')

#将GO分开进行分析
AD_Healthyres_GOBP <- enrichGO(GOgene, 
                               OrgDb = org.Hs.eg.db, 
                               ont='BP',
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.1, 
                               qvalueCutoff = 0.1,
                               keyType = 'ENTREZID')


##得到的结果是一个大的enrichresult
library(stringr)
p1 <- barplot(AD_Healthyres_GOBP,
              x = "Count",
              color = "p.adjust",
              showCategory = 15,
              size = NULL,
              font.size = 12,
              title = "AD_Healthy padj<0.05 |lgfc|>0.585 GO BP",
              label_format = 30)+
  scale_size(range=c(2, 12))+
  theme_bw()+
  # theme_classic()+
  labs(x="Count",y="GO BP") +
  theme(#panel.border = element_blank(),
    #axis.text.y = element_blank()+
    panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
    plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
    axis.title.x = element_text(color = "black",size = 10,face = "plain"),
    axis.title.y = element_text(color = "black",size = 10,face = "plain"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black",size = 10,face = "plain"),
    axis.text.y = element_text(color = "black",size = 10,face = "plain"),
    legend.position = "right",
    legend.title = element_text(color = "black",size = 6,face = "plain"),
    legend.text = element_text(color = "black",size = 6,face = "plain"),
    legend.background = element_rect(fill = "transparent",size = 0.4,linetype = "solid",colour = "transparent"))

AD_Healthyres_GOBP_genelist<-setReadable(AD_Healthyres_GOBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_GOBP_genelist <- as.data.frame(AD_Healthyres_GOBP_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_GOBP_genelist, file='GOBP_genelist.csv')



AD_Healthyres_GOCC <- enrichGO(GOgene, 
                               OrgDb = org.Hs.eg.db, 
                               ont='CC',
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.1, 
                               qvalueCutoff = 0.1,#0.3
                               keyType = 'ENTREZID')

##得到的结果是一个大的enrichresult
p2 <- barplot(AD_Healthyres_GOCC,
              x = "Count",
              color = "p.adjust",
              showCategory = 15,
              size = NULL,
              font.size = 12,
              title = "AD_Healthy padj<0.05 |lgfc|>0.585 GO CC",
              label_format = 30)+
  scale_size(range=c(2, 12))+
  theme_bw()+
  # theme_classic()+
  labs(x="Count",y="GO CC") +
  theme(#panel.border = element_blank(),
    #axis.text.y = element_blank()+
    panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
    plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
    axis.title.x = element_text(color = "black",size = 10,face = "plain"),
    axis.title.y = element_text(color = "black",size = 10,face = "plain"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black",size = 10,face = "plain"),
    axis.text.y = element_text(color = "black",size = 10,face = "plain"),
    legend.position = "right",
    legend.title = element_text(color = "black",size = 6,face = "plain"),
    legend.text = element_text(color = "black",size = 6,face = "plain"),
    legend.background = element_rect(fill = "transparent",size = 0.4,linetype = "solid",colour = "transparent"))

AD_Healthyres_GOCC_genelist<-setReadable(AD_Healthyres_GOCC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_GOCC_genelist <- as.data.frame(AD_Healthyres_GOCC_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_GOCC_genelist, file='GOCC_genelist.csv')


AD_Healthyres_GOMF <- enrichGO(GOgene, 
                               OrgDb = org.Hs.eg.db, 
                               ont='MF',
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.1, 
                               qvalueCutoff = 0.1,
                               keyType = 'ENTREZID')

#得到的结果是一个大的enrichresult
p3 <- barplot(AD_Healthyres_GOMF,
              x = "Count",
              color = "p.adjust",
              showCategory = 15,
              size = NULL,
              font.size = 12,
              title = "AD_Healthy padj<0.05 |lgfc|>0.585 GO MF",
              label_format = 30)+
  scale_size(range=c(2, 12))+
  theme_bw()+
  # theme_classic()+
  labs(x="Count",y="GO MF") +
  theme(#panel.border = element_blank(),
    #axis.text.y = element_blank()+
    panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
    plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
    axis.title.x = element_text(color = "black",size = 10,face = "plain"),
    axis.title.y = element_text(color = "black",size = 10,face = "plain"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black",size = 10,face = "plain"),
    axis.text.y = element_text(color = "black",size = 10,face = "plain"),
    legend.position = "right",
    legend.title = element_text(color = "black",size = 6,face = "plain"),
    legend.text = element_text(color = "black",size = 6,face = "plain"),
    legend.background = element_rect(fill = "transparent",size = 0.4,linetype = "solid",colour = "transparent"))

AD_Healthyres_GOMF_genelist<-setReadable(AD_Healthyres_GOMF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_GOMF_genelist <- as.data.frame(AD_Healthyres_GOMF_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_GOMF_genelist, file='GOMF_genelist.csv')
#将GO的3张图组合在一起
library(ggpubr)
p <- ggarrange(p1, p2,p3, ncol = 3, nrow = 1)+
  theme(plot.margin = margin(t = 20,r = 20,b = 20,l = 20))
# ggsave("GO15_0.05_0.585.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(5, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)


#绘制合在一块的GO图
go_enrich_df <- data.frame(
  ID=c(AD_Healthyres_GOBP$ID[1:15], AD_Healthyres_GOCC$ID[1:15], AD_Healthyres_GOMF$ID[1:15]),
  Description=c(AD_Healthyres_GOBP$Description[1:15],AD_Healthyres_GOCC$Description[1:15],AD_Healthyres_GOMF$Description[1:15]),
  GeneNumber=c(AD_Healthyres_GOBP$Count[1:15], AD_Healthyres_GOCC$Count[1:15], AD_Healthyres_GOMF$Count[1:15]),
  type=factor(c(rep("biological process", 15), 
                rep("cellular component", 15),
                rep("molecular function", 15)), 
              levels=c("biological process", "cellular component","molecular function" )))

go_enrich_df <- arrange(go_enrich_df,type,GeneNumber)
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# write.csv(go_enrich_df, file='go_enrich_dfGE.csv')
##开始绘制GO柱状图
###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#0D6CA6","#099963", "#911F27")#设定颜色

# p4 <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
#   geom_segment(aes(x=type_order,
#                    xend=type_order,
#                    y=0,
#                    yend=GeneNumber))+
#   geom_point(aes(size=GeneNumber,color=type)) +
#   scale_color_manual(values=COLS)+
#   # geom_bar(stat="identity", width=0.8)+
#   # scale_fill_manual(values = COLS) +
#   theme_bw() + 
#   xlab("GO Term") + 
#   ylab("Num of Genes") + 
#   labs(title = "GO")+ 
#   scale_size(range=c(2, 8))+
# theme(
#   plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#   plot.margin = margin(t = 40,r = 80,b = 40,l = 80),
#   axis.title.x = element_text(color = "black",size = 10,face = "plain"),
#   axis.title.y = element_text(color = "black",size = 10,face = "plain"),
#   panel.grid = element_blank(),
#   axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 70,vjust = 1, hjust = 1 ),
#   axis.text.y = element_text(color = "black",size = 8,face = "plain"),
#   legend.position = "right",
#   legend.title = element_text(color = "black",size = 6,face = "plain"),
#   legend.text = element_text(color = "black",size = 6,face = "plain"),
#   legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"))

library(ggpubr)
p4 <- ggdotchart(go_enrich_df, x = "type_order", y = "GeneNumber",
                 color = "type",                                # 按照cyl填充颜色
                 palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                 sorting = "descending",                      
                 add = "segments",                             # 添加棒子
                 add.params = list(color = "type", size = 1.3),#改变棒子参数
                 # rotate = TRUE,                                # 方向转为垂直
                 group = "type",                                
                 dot.size = "GeneNumber",                                 # 改变点的大小
                 #label = round(go_enrich_df$GeneNumber),                       # 添加label
                 font.label = list(color = "black", size = 7,
                                   vjust = 0.5),               # 设置label参数
                 ggtheme = theme_pubr(),                        # 改变主题
                 xlab="GO Term",
                 ylab="GeneNumber",
                 title = "GE GO")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 14,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 15,r = 10,b = 10,l = 20),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 70,vjust = 1, hjust = 1 ),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.position = c(0.9,0.68),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"))
p4
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_GO.pdf", egg::set_panel_size(p4, width=unit(16, "in"), height=unit(5, "in")), 
       width = 20, height = 12, units = 'in', dpi = 600)

#GO图的另外一种展现形式
# a <- AD_Healthyres_GO_genelist$geneID[15]
# a <- unlist(strsplit(a,split='/'))
# b <- subset(AD_Healthyres,rownames(AD_Healthyres) %in% a)
# table(b$change)
# library(reshape2)
# library(knitr)
# GO_term <- read.csv("GO_trem.csv",header = T,check.names = F)
# kable(GO_term,format="markdown") 
# mydata<-melt(GO_term,id.vars=c("ï»¿ID","Description"),variable.name="Change",value.name="Number")
# kable(mydata,format='markdown')
# p1 <- ggplot(mydata,aes(Description,Number)) + 
#   geom_bar(aes(fill=factor(Change)),stat="identity", colour="white",width=0.99)+ 
#   scale_fill_manual(values = c("#008080","firebrick3"))+
#   geom_text(aes(label=abs(Number)),color="black", size=3.5,
#             position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
#   coord_flip()+
#   theme_bw()+
#   ggtitle("GO BP Plot")+
#   # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
#   xlab("GO BP")+ylab("")+
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 35))+
#   theme(
#     legend.key = element_blank(),
#     plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#     axis.title.x = element_text(color = "black",size = 8,face = "plain"),
#     axis.title.y = element_text(color = "black",size = 10,face = "plain"),
#     panel.grid= element_blank(),panel.border = element_blank(),
#     axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(color = "black",size = 10,face = "plain"),
#     legend.title = element_text(color = "black",size = 10,face = "plain"),
#     legend.text = element_text(color = "black",size = 10,face = "plain"),
#     legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "white"),
#     legend.direction = 'vertical', legend.position ="right")
# ggsave(p1, filename = "AD_Healthy_BPbar0.585.png", width = 7, height = 5)







#KEGG富集分析
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
AD_Healthyres_kegg <- enrichKEGG(gene = GOgene,
                                 keyType = "kegg",
                                 organism  = 'hsa',
                                 pvalueCutoff  = 0.1,
                                 pAdjustMethod  = "BH",
                                 qvalueCutoff  = 0.1)

#用ggplot2来展示

p1 <- barplot(AD_Healthyres_kegg,
              x = "Count",
              color = "p.adjust",
              showCategory = 20,
              size = NULL,
              split = NULL,
              font.size = 12,
              title = "GE KEGG",
              label_format = 30)+
  scale_size(range=c(2, 12))+
  theme_bw()+
  # theme_classic()+
  labs(x="Count",y="KEGG") +
  theme(#panel.border = element_blank(),
    #axis.text.y = element_blank()+
    plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
    axis.title.x = element_text(color = "black",size = 10,face = "plain"),
    axis.title.y = element_text(color = "black",size = 10,face = "plain"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black",size = 10,face = "plain"),
    axis.text.y = element_text(color = "black",size = 8,face = "plain"),
    legend.position = "right",
    legend.title = element_text(color = "black",size = 6,face = "plain"),
    legend.text = element_text(color = "black",size = 6,face = "plain"),
    legend.background = element_rect(fill = "transparent",size = 0.4,linetype = "solid",colour = "transparent"))
# ggsave(p1, filename = "AD_Healthy_KEGG小于0.01下调.png", width = 4.5, height = 5.5)
AD_Healthyres_kegg<-setReadable(AD_Healthyres_kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_kegg <- as.data.frame(AD_Healthyres_kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_kegg,file = "AD_Healthyres_kegg_genelist0.05and0.585.csv",row.names = F)



#KEGG图的另外一种展现形式
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
AD_Healthyres_kegg <- enrichKEGG(gene = GOgene,
                                 keyType = "kegg",
                                 organism  = 'hsa',
                                 pvalueCutoff  = 0.1,
                                 pAdjustMethod  = "BH",
                                 qvalueCutoff  = 0.1)
kegg<-setReadable(AD_Healthyres_kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
kegg <- kegg[order(kegg$Count,decreasing = TRUE),]

#KEGG图的另外一种展现形式,展现的是前15条信号通路
#将富集到的信号通路按照富集到的基因的个数进行排序
# 定义一个函数来计算下调和上调基因的数量
table1 <- AD_Healthyres_def[,c(1,2,7)]#做GE分析的时候注意修改这里的数字
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  # 统计下调和上调基因数量
  down_genes <- sum(gene_list %in% rownames(AD_Healthyres_def)[AD_Healthyres_def$change == "Down"])
  up_genes <- sum(gene_list %in% rownames(AD_Healthyres_def)[AD_Healthyres_def$change == "Up"])
  return(c(down_genes, up_genes))
}

# 对表格2中的每一行计算down和up的数量
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))

# 将结果加到表格2中
table2$down <- -abs(result[, 1])
table2$up <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "GE_KEGG.csv",row.names = F)
#对数据进行转换画图
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:15,c(3,4,12,13)] #做GE分析的时候注意这里的数字
colnames(KEGGTerm)
#对数据格式进行转换
mydata<-melt(KEGGTerm,id.vars=c("ID","Description"),variable.name="Change",value.name="Number")
mydata$type_order=factor(rev(as.integer(rownames(mydata))),labels=rev(mydata$Description))

library(stringr)
p1 <- ggplot(mydata,aes(type_order,Number)) + 
  geom_bar(aes(fill=factor(Change)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#008080","firebrick3"))+
  geom_text(aes(label=abs(Number)),color="black", size=4,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("GE KEGG")+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Gene Number")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill = guide_legend(title = 'Change'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 15,face = "plain"),
        axis.title.y = element_text(color = "black",size = 15,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 15,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "white"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)


#绘制Volcano plot
library(dplyr)
up <- subset(AD_Healthyres, AD_Healthyres$change == 'Up')
up <- arrange(up,padj,log2FoldChange)
up <- up[order(up$padj), ][1:5, ]
down <- subset(AD_Healthyres, AD_Healthyres$change == 'Down')
down <- down[order(down$padj), ][1:5, ]
a <- rbind(up, down)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(a,file = "GE_MostSignificantGene.csv",quote = F)
library(ggplot2)
library(ggrepel)
p <- ggplot(
  # 数据、映射、颜色
  AD_Healthyres, aes(x = log2FoldChange, y = -log10(padj),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080", "gray", "firebrick3")) +
  geom_text_repel(data = rbind(up, down), aes(x = log2FoldChange, y = -log10(padj), label = rownames(a)),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
  # 坐标轴
  labs(x="log2FoldChange",y="-log10(padj)") +
  ggtitle("AD_Healthy volcano")+
  theme_bw()+
  # theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.position = "right",
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.background = element_rect(fill = "transparent",linewidth = 0.4,linetype = "solid"))+
  scale_x_continuous(breaks = seq(-10, 10, 5))+
  scale_y_continuous(breaks = seq(0, 50, 10))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_Volcano.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)




# #在文章中查到的和AD有关的基因
# ADgene <- c("CEBPA","REST","TREM2","APOE","ABCA7","GPR141","PTK2B","SPI1",
#             "MEF2C","CD2AP","SORL1","CR1","IL15","MS4A6A","MS4A4A","NME8",
#             "GPR141","CECR2","CLU","PLCG2","GFAP","SERPING1","KRT5","NRN1",
#             "ZNF568","ZBTB16","HIF3A")




#绘制热图
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AD_Healthyres <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
table(AD_Healthyres$padj<0.05 & abs(AD_Healthyres$log2FoldChange)>0.585 )    ##看下p小于0.05的有多少
# FALSE  TRUE 
# 17221   998
AD_Healthyres <- AD_Healthyres[order(AD_Healthyres$padj,decreasing = F),]
head(AD_Healthyres)
choose_gene <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585)
dim(choose_gene)
# [1] 998   7


#要标注的最显著的基因是a
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#MostSigGene指的是上调和下调的基因中padj最小的前10个基因
MostSigGene <- read.csv("GE_MostSignificantGene.csv",header = T,row.names = 1)


#从FPKM表格中提取数据，首先将fpkm中的ENSEMBL转换为SYMBOL
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myfpkm <- read.csv(file = "GE_RemoveBatch_FPKM.csv",header = T ,row.names = 1)
#将基因的名字进行转换，从ensembl转换为Symbol
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
####将所有的差异基因做GO和KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]


AD_Healthyres_heatmap <- subset(myfpkm, rownames(myfpkm) %in% rownames(choose_gene))#这个非常棒
#将GSM号转换为AD还是Healthy
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AD_Healthyres_heatmap <- as.data.frame(t(AD_Healthyres_heatmap))
AD_Healthyres_heatmap$GSM_number <- rownames(AD_Healthyres_heatmap)
sum <- merge(AD_Healthyres_heatmap,BackgroundInformation,by = "GSM_number",all= FALSE)
sum[5,999:1010]
sum[5,1:5]
rownames(sum) <- sum$Sample
#对样本进行排序，使得AD和Healthy的样本分为上下两个不同的组
sum <- arrange(sum,rownames(sum))
data <- sum[-c(1,1000:1010)]
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){
  sd(x)/mean(x)*100
})
data_df <- data[order(data$cv, decreasing = T),1:354]
dim(data_df)
a <- apply(data_df,1,scale)
# 进行scale，但是不设定最大最小值，这个时候的热图一片绿
data_scale <- as.data.frame(t(apply(data_df,1,scale))) ##Z-score标准化
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.01
#将scale之后的数据设置最大值最小值
data_scale <- as.matrix(data_scale)
#设置最大值与最小值
table((data_scale)>2)
table((data_scale)<(-2))
data_scale[data_scale>=2]=2
data_scale[data_scale<=-2]=-2
table(is.na(data_scale))
# library(pheatmap)
# pheatmap(data_scale)
library(ComplexHeatmap)
library(circlize)

#读取待展示的基因名称，并添加到热图中
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "GE_Heatmap.pdf",width =3,height = 4)
# png("p.png",res = 300)

p <- Heatmap(data_scale,name = "Expression", 
             na_col = "grey",
             cluster_rows = TRUE,
             clustering_distance_rows = "pearson",
             clustering_method_rows = "complete",
             col = colorRampPalette(c("#008080", "white", "firebrick3"))(50),
             cluster_columns = TRUE,
             clustering_distance_columns = "pearson",
             clustering_method_columns = "complete",
             row_dend_side = "left",
             show_row_names = FALSE,
             show_column_names = FALSE,
             use_raster = T,
             heatmap_width = unit(1, "npc"),
             width = NULL,
             heatmap_height = unit(1, "npc"),
             height = NULL,
             show_column_dend = TRUE,
             column_dend_height = unit(5, "mm"),
             show_row_dend = FALSE,
             column_labels = colnames(data_scale),
             column_names_side = "bottom",
             column_names_centered = TRUE,
             column_title = "",
             column_names_rot = 0,
             column_title_gp = gpar(fontsize = 0.1),
             column_names_gp = gpar(fontsize = 0.1),
             heatmap_legend_param = list( #设置标签
               color_bar = 'continuous',
               legend_direction = 'vertical',
               legend_width = unit(2, 'cm'),
               legend_height = unit(2, 'cm'),
               title_position = 'topcenter',
               title_gp = gpar(fontsize = 5, fontface = 'plain'),
               labels_gp = gpar(fontsize = 5, fontface = 'plain'))
)+ 
  rowAnnotation(link = anno_mark(at = which(rownames(data_scale) %in% rownames(MostSigGene)), 
                                 labels = rownames(data_scale[which(rownames(data_scale) %in% rownames(MostSigGene)),]), labels_gp = gpar(fontsize = 7)))

print(p)
dev.off()








#对差异基因进行WGCNA分析之后寻找基因和表型之间的关系
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sigGene <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",row.names = 1)

library(stringr)
#加载背景信息，其中包括APOE的类型
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
BackgroundInformationAPOE <- BackgroundInformationAPOE[BackgroundInformationAPOE$Brain_Region==3,]

#找出这些基因的FPKM值，然后使得行作为样本，列作为基因的名字
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
#对FPKM的基因名进行ID转换
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
####
library(clusterProfiler)
# packageVersion("clusterProfiler")
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]

##将差异基因的FPKM值提取出来
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% rownames(sigGene))
#将数据转换为行是样品名列是基因名的格式
sig_FPKM <- as.data.frame(t(sig_FPKM))
# 对样本进行检查
library(WGCNA)
gsg <- goodSamplesGenes(sig_FPKM, verbose=3)
gsg

#从表型数据中提取出来颞叶的数据
#表型的数据也要是行是样本，列是表型的数据
BackgroundInformation <- subset(BackgroundInformationAPOE,BackgroundInformationAPOE$GSM_number %in% rownames(sig_FPKM))
rownames(BackgroundInformation) <- BackgroundInformation$GSM_number
colnames(BackgroundInformation)
BackgroundInformation <- BackgroundInformation[,-c(1:6,10:11)]



#

#将APOE的类型转换为数字1，2，3
library(dplyr)
colnames(BackgroundInformation)

BackgroundInformation <- BackgroundInformation%>%mutate(Phenotypedit=case_when(BackgroundInformation$Phenotypedit=="noE4"~"0",
                                                                               BackgroundInformation$Phenotypedit=="E4carrier"~"1",
                                                                               BackgroundInformation$Phenotypedit=="E4/4"~"2"))

#将Phenotypedit转换为int类型，而不是chr类型
BackgroundInformation$Phenotypedit <- as.numeric(BackgroundInformation$Phenotypedit)
BackgroundInformation$Braak_Stage <- as.numeric(BackgroundInformation$Braak_Stage)
# BackgroundInformation$Sex <- as.numeric(BackgroundInformation$Sex)
BackgroundInformation$Bank_Location <- as.numeric(BackgroundInformation$Bank_Location)
str(BackgroundInformation)
#正式开始WGCNA分析
#计算样本之间的距离，对样本进行聚类
sampleTree <- hclust(dist(sig_FPKM), method="average")
sizeGrWindow(2, 2)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="",
     xlab="", cex.lab=1.5, 
     cex.axis=1.5, cex.main=2)



## 样本-表达-临床关联分析
sampleTree2 <- hclust(dist(sig_FPKM), method="average")
# 对每一个样本和临床指标都赋予一个颜色，所以就是134*26=3484
traitColors <- numbers2colors(BackgroundInformation, signed=F)
# signed=F设置了之后表示的就是白色表示最小值，红色表示最大值
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=colnames(BackgroundInformation), 
                    main="Sample dendrogram and trait heatmap")


# 2.2.1 确定软阈值β
# 这里我们使用pickSoftThreshold函数，来获取合适的软阈值β：

# 软阈值的预设范围
powers <- c(c(1:10), seq(from=12, to=50, by=2))

# 自动计算推荐的软阈值
library("doParallel")
sft <- pickSoftThreshold(sig_FPKM, powerVector=powers, verbose=5, networkType="unsigned")
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate
# [1] 6

# 作图
sizeGrWindow(2, 2)
par(mfrow=c(1, 2))#画出来一个1行2列的图
cex1 <- 0.9

# Scale-free topology fit index 作为评价软阈值指标
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
# h值作为筛选指标
abline(h=0.85, col="red")

# Mean connectivity 作为评价软阈值的指标
# 选择软阈值β，使其满足无尺度网络和较好的网络连接性
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")


# 如果是NA，此时就需要根据图来自己指定数值了

sft$powerEstimate <- 6



# 2.2.2 构建网络并识别模块

# 直接调用blockwiseModules函数构建网络并确定模块

cor <- WGCNA::cor
#如果对这个不设置的话会出现报错
# Error in (new("standardGeneric", .Data = function (x, y = NULL, use = "everything", : unused arguments (weights.x = NULL, weights.y = NULL, cosine = FALSE)
# I am pretty sure the problem is that the stats::cor function is being used under the hood, rather than WGCNA::cor. I have not found a fix to this issue and it seems many people are experiencing it. Can anyone recommend a workaround.
#解决办法如下
# WGCNA与其他软件包之间存在冲突。WGCNA中的cor函数与R中自带的cor在命名空间上有冲突。
# 解决方法为：在使用该函数之前暂时重新分配功能
net <- blockwiseModules(sig_FPKM, power = sft$powerEstimate,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,randomSeed=47,
                        saveTOMFileBase = "ADPatientTOM20240910", 
                        verbose = 3)
table(net$colors)
# 0   1   2 
# 155 570 273 
# 在 WGCNA 中，编号为 0 或标记为 grey 的模块通常表示没有被分配到其它任何模块的基因。
# 这些基因被认为是不共表达或与其它基因的共表达关系不明显，
# 因此它们被归类为一个“未分配”或“噪声”模块。
cor<-stats::cor

# 查看每个模块的基因数，其中0模块下为没有计算进入模块的基因数
table(net$colors)
# 0   1   2 
# 155 570 273
# 打开新的绘图窗口
sizeGrWindow(2, 2)
# 把模块编号转成颜色
mergedColors <- labels2colors(net$colors)
# 绘制树状图和下面的模块颜色
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "GE_WGCNA_Module_colors.pdf",width =5,height = 3,bg = "white")
p1 <- plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          cex.rowText = 0.9,cex.colorLabels = 0.9, cex.dendroLabels = 0.9)
print(p1)
dev.off()


# 计算模块特征向量MEs
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
# 保存数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
save(MEs, moduleLabels, moduleColors, geneTree,
     file="GE_networkConstruction-auto.RData")


# 3. 模块-临床信息关联并识别重要基因
# 3.1. 导入处理的数据：
# 导入包
library(WGCNA)
# 设置将strings不要转成factors
options(stringsAsFactors=F)

# 允许使用最大线程
allowWGCNAThreads()

# 3.2. 模块-临床信息关联定量

# 这步发分析主要是为例确定与我们感兴趣的临床特征有关的模块
# 在这项分析中，我们想要确定与定量与临床特征显著相关的模块。由于我们已经有了每个模块的信息(即，特征基因 Eigengene )，我们只是将特征基因与临床特征相关联，并寻找最重要的关联：

# 确定基因与样本个数
nGenes = ncol(sig_FPKM)
nSamples = nrow(sig_FPKM)
# 用颜色标签重新估计MEs（模块特征基因）
a <- moduleEigengenes(sig_FPKM, moduleColors)
table(moduleColors)
MEs0 = moduleEigengenes(sig_FPKM, moduleColors)$eigengenes
#求出来不同的模块喝临床特征之间的相关性
#模块的dataframe的横坐标是样本，纵坐标表示的是模块的标签
#这个相关性的计算和一般的cor函数计算得到的结果是一样的
identical(rownames(MEs0),rownames(BackgroundInformation))
#第二次做的时候这里可能会发生变化，因为少了一个表型
moduleTraitCor = cor(MEs0, BackgroundInformation,use = 'p',method = "spearman")
#计算相关性的p值
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# 由于我们有相当多的模块和临床特征，适当的图形表示将有助于理解关联表。我们用相关值对每个关联进行颜色编码:

sizeGrWindow(2,2)
# 显示相关性及其p值
textMatrix <-  paste(signif(moduleTraitCor, 4), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
colnames(BackgroundInformation)
colnames(BackgroundInformation)[4] <- "APOE"
colnames(BackgroundInformation)[1] <- "Braak Stage"
colnames(BackgroundInformation)[3] <- "Bank Location"
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "GE_WGCNA_Module_trait.pdf",width =5,height = 2.5,bg = "white")
par(mar = c(5, 7, 3, 3))
# 在热图图中显示相关值
p1 <- labeledHeatmap(Matrix = moduleTraitCor,
                     xLabels = colnames(BackgroundInformation),
                     yLabels = colnames(MEs0),
                     ySymbols = colnames(MEs0),
                     colorLabels = FALSE,
                     colors = blueWhiteRed(50),
                     textMatrix = textMatrix,
                     setStdMargins = FALSE,
                     cex.text = 0.5,
                     zlim = c(-1,1),
                     cex.main = 1, cex.lab = 0.75, cex.axis = 0.75,
                     main = paste("GE Module Trait Relationships"))
print(p1)
dev.off()


# 从临床信息中单独提取体重APOE的数据
APOE <- as.data.frame(BackgroundInformation$APOE)
colnames(APOE) = "APOE"
# 获得模块的名称（颜色）
modNames <- substring(names(MEs0), 3)#从第3个字符开始截取字符

# 计算Module Membership MM指标
#计算出来基因的表达量和模块之间的相关性
identical(rownames(MEs0),rownames(sig_FPKM))
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) = paste("p.MM", modNames, sep="")

## 计算Gene Significance GS指标
#计算出来基因的表达量和APOE之间的相关性
identical(rownames(sig_FPKM),rownames(BackgroundInformation))
#保证我提取出来的APOE的样本的顺序和样本的FPKM值对应的样本的顺序是一致的
geneTraitSignificance = as.data.frame(cor(sig_FPKM, APOE, use = "p",method = "spearman"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

colnames(geneTraitSignificance) = paste("GS.", names(APOE), sep="");
colnames(GSPvalue) = paste("p.GS.", names(APOE), sep="")
# 3.3. 模块内分析：确定高MM和GS的基因
# 
# 在这一步主要是为了展示模块中基因的GS和MM
# 
# 基于GS和MM指标，我们可以筛选出一些重要的基因，这些基因即与APOE高度相关，也在模块内很重要。比如，我们看一下与重量关联最大的棕色模块（MEbrown）。我们在MEbrown中绘制了基因重要性与模块中基因的关系的散点图：


##将相关系数和p值进行合并
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
#取出来blue模块中的基因,保证上面展示的基因的相关系数都是有意义的，
#所有同时要保证pvalue<0.05
allLLIDs <- colnames(sig_FPKM)

###绘制点状图
# 0   1   2 
# 155 570 273 
module = "blue"#273
# 获得模块的探针ID
modGenes = (moduleColors==module)
table(modGenes)
# 获得探针ID对应的Entrez ID
modLLIDs = allLLIDs[modGenes]#

b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]#绿色
#相关分为正相关和负相关，所以需要对相关系数取绝对值
m$MMblue <- abs(m$MMblue)
m$GS.APOE <- abs(m$GS.APOE)
m1 <- m
m <- subset(m,m$p.MMblue<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$GS.APOE) > 0.3)
#最终筛选完之后没有剩余任何的基因
blue_APOE <- n$Row.names
#0个基因
write.csv(n,file = "GE_nblueAPOE.csv")
colnames(n)
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMblue) , y = abs(GS.APOE))) +
  geom_point(size=4,colour="blue") +
  # geom_text_repel(data =n, aes(x = abs(MMblue), y = abs(GS.APOE), label = Row.names),
  #                 size = 2,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
  #                 color = "black",max.overlaps=100)+
  labs(x="MMblue",y="GS.APOE") +
  ggtitle("MM vs GS")+
  geom_smooth(method = 'lm',se = F, color = '#800000')+
  stat_cor(data=m1, method = "spearman",label.x = 0.1, label.y = 0.38,digits = 1,size = 7)+
  theme_bw()+
  geom_vline(xintercept=0.8,lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = 0.3,lty=4,col="#666666",lwd=0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))+
  scale_y_continuous(breaks = seq(0, 1, 0.2))+
  scale_x_continuous(breaks = seq(0, 1, 0.2))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_GS_MMblueAPOE.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)

#绘制另外一个模块的饿点状图
# 0   1   2 
# 155 570 273 
module = "turquoise"#570
# 获得模块的探针ID
modGenes = (moduleColors==module)
# 获得探针ID对应的Entrez ID
modLLIDs = allLLIDs[modGenes]#

b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]#青色
m$MMturquoise <- abs(m$MMturquoise)
m$GS.APOE <- abs(m$GS.APOE)
m1 <- m
m <- subset(m,m$p.MMturquoise<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$GS.APOE) > 0.3)
colnames(n)
turquoise_APOE <- n$Row.names
# [1] "ADCYAP1"   "LINC01182" "RPH3A" 
write.csv(n,file = "GE_nturquoiseAPOE.csv")
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMturquoise) , y = abs(GS.APOE))) +
  geom_point(size=4,colour="turquoise") +
  # geom_text_repel(data =n, aes(x = abs(MMturquoise), y = abs(GS.APOE), label = Row.names),
  #                 size = 2,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
  #                 color = "black",max.overlaps=100)+
  labs(x="MMturquoise",y="GS.APOE") +
  ggtitle("MM vs GS")+
  geom_smooth(method = 'lm',se = F, color = '#800000')+
  stat_cor(data=m1, method = "spearman",label.x = 0.1, label.y = 0.35,digits = 1,size = 7)+
  theme_bw()+
  geom_vline(xintercept=0.8,lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = 0.3,lty=4,col="#666666",lwd=0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_GS_MMturquoiseAPOE.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)



#绘制第二个部分的临床信息，选择的是BraakStage
# 从临床信息中单独提取体重BraakStage的数据
Braak_Stage <- as.data.frame(BackgroundInformation$`Braak Stage`)
colnames(Braak_Stage) = "Braak Stage"
# 获得模块的名称（颜色）
modNames <- substring(names(MEs0), 3)#从第3个字符开始截取字符

## 计算Module Membership MM指标
#计算出来基因的表达量和模块之间的相关性
identical(rownames(MEs0),rownames(sig_FPKM))
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) = paste("p.MM", modNames, sep="")

## 计算Gene Significance GS指标
#计算出来基因的表达量和APOE之间的相关性
identical(rownames(sig_FPKM),rownames(BackgroundInformation))
geneTraitSignificance = as.data.frame(cor(sig_FPKM, Braak_Stage, use = "p",method = "spearman"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

colnames(geneTraitSignificance) = paste("GS.", names(Braak_Stage), sep="");
colnames(GSPvalue) = paste("p.GS.", names(Braak_Stage), sep="")
# 3.3. 模块内分析：确定高MM和GS的基因
# 
# 在这一步主要是为了展示模块中基因的GS和MM
# 
# 基于GS和MM指标，我们可以筛选出一些重要的基因，这些基因即与APOE高度相关，也在模块内很重要。比如，我们看一下与重量关联最大的棕色模块（MEbrown）。我们在MEbrown中绘制了基因重要性与模块中基因的关系的散点图：


##将相关系数和p值进行合并
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
#取出来blue模块中的基因,保证上面展示的基因的相关系数都是有意义的，

###绘制点状图
# 0   1   2 
# 155 570 273 
module = "blue"#273
# 获得模块的探针ID
modGenes = (moduleColors==module)
# 获得探针ID对应的Entrez ID
modLLIDs = allLLIDs[modGenes]#

b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]#绿色
m$MMblue <- abs(m$MMblue)
m$`GS.Braak Stage` <- abs(m$`GS.Braak Stage`)
m1 <- m
m <- subset(m,m$p.MMblue<0.05 & m$`p.GS.Braak Stage` <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$`GS.Braak Stage`) > 0.3)
write.csv(n,file = "GE_nblueBraak_Stage.csv")
colnames(n)
blue_BraakStage <- n$Row.names
# [1] "ABCC2"   "C8orf88" "CLDN15"  "CLMN"    "CPM"     "PPFIBP2" "PRKX"    "SGK1"    "TNS1" 
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMblue) , y = abs(`GS.Braak Stage`))) +
  geom_point(size=4,colour="blue") +
  # geom_text_repel(data =n, aes(x = abs(MMblue), y = abs(`GS.Braak Stage`), label = Row.names),
  #                 size = 2,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
  #                 color = "black",max.overlaps=100)+
  labs(x="MMblue",y="GS.Braak Stage") +
  ggtitle("MM vs GS")+
  geom_smooth(method = 'lm',se = F, color = '#800000')+
  stat_cor(data=m1, method = "spearman",label.x = 0.1, label.y = 0.38,digits = 1,size = 7)+
  theme_bw()+
  geom_vline(xintercept=0.8,lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = 0.3,lty=4,col="#666666",lwd=0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_GS_MMblueBraakStage.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)

#绘制另外一个模块的饿点状图
# 0   1   2 
# 155 570 273 
module = "turquoise"#570
# 获得模块的探针ID
modGenes = (moduleColors==module)
# 获得探针ID对应的Entrez ID
modLLIDs = allLLIDs[modGenes]#

b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]#青色
m$MMturquoise <- abs(m$MMturquoise)
m$`GS.Braak Stage` <- abs(m$`GS.Braak Stage`)
m1 <- m
m <- subset(m,m$p.MMturquoise<0.05 & m$`p.GS.Braak Stage` <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$`GS.Braak Stage`) > 0.3)
colnames(n)
turquoise_BraakStage <- n$Row.names
# [1] "ABCC12"         "ADCYAP1"        "AGBL1"          "AGBL4"          "ANO3"          
# [6] "ARL4D"          "ARMCX5-GPRASP2" "ATOH7"          "ATP2B1-AS1"     "BDNF"          
# [11] "BEX1"           "BEX5"           "C2orf80"        "C3orf80"        "CALY"          
# [16] "CAP2"           "CARTPT"         "CCDC184"        "CDK5R2"         "CHN1"          
# [21] "CPLX1"          "CRYM"           "DGKB"           "ENC1"           "FAM81A"        
# [26] "FLRT2-AS1"      "GABRA1"         "GAD1"           "GAD2"           "GAP43"         
# [31] "GCNT4"          "GFRA2"          "GNG13"          "GPRASP2"        "HPRT1"         
# [36] "HSPB3"          "HTR3B"          "IPCEF1"         "KCNB2"          "KCNV1"         
# [41] "KRT5"           "LINC00460"      "LINC00507"      "LINC01007"      "LINC01182"     
# [46] "LINC01331"      "LINC02347"      "LINC02458"      "LNCBRM"         "LOC101927189"  
# [51] "LOC101927531"   "LOC105376121"   "LOC124902439"   "LOC124902537"   "LY86-AS1"      
# [56] "LYRM9"          "MAL2"           "MAS1"           "MCHR1"          "MDH1"          
# [61] "MLIP"           "MLIP-IT1"       "MSC-AS1"        "NAP1L2"         "NAP1L5"        
# [66] "NCALD"          "NEFL"           "NEUROD6"        "NPTXR"          "NRN1"          
# [71] "NXPH2"          "OLFM3"          "OTOGL"          "PART1"          "PAX7"          
# [76] "PCP4"           "PCSK1"          "PPEF1"          "PRMT8"          "PTER"          
# [81] "RAB3C"          "RGS4"           "RPH3A"          "SLC30A3"        "SNAP25"        
# [86] "SNX10"          "SOWAHB"         "SPTSSB"         "ST6GALNAC5"     "STAT4"         
# [91] "STMN2"          "SV2B"           "SVOP"           "SYP"            "SYT1"          
# [96] "TAC1"           "TESPA1"         "TMEM132E-DT"    "TUBB2A"         "VGF"           
# [101] "VSNL1" 
write.csv(n,file = "GE_nturquoiseBraak_Stage.csv")
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMturquoise) , y = abs(`GS.Braak Stage`))) +
  geom_point(size=4,colour="turquoise") +
  # geom_text_repel(data =n, aes(x = abs(MMturquoise), y = abs(`GS.Braak Stage`), label = Row.names),
  #                 size = 2,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
  #                 color = "black",max.overlaps=100)+
  labs(x="MMturquoise",y="GS.Braak Stage") +
  ggtitle("MM vs GS")+
  geom_smooth(method = 'lm',se = F, color = '#800000')+
  stat_cor(data=m1, method = "spearman",label.x = 0.1, label.y = 0.4,digits = 1,size = 7)+
  theme_bw()+
  geom_vline(xintercept=0.8,lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = 0.3,lty=4,col="#666666",lwd=0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_GS_MMturquoiseBraakStage.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(3, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)


#绘制这些基因的小提琴图
#绘制青色模块中AD和Healthy中选择到的基因的表达量
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#第二次分析得出的结果
turquoise_mark_gene <- append(turquoise_APOE,turquoise_BraakStage)
gene <- turquoise_mark_gene
gene <- gene[!duplicated(gene)]
# #第一次分析得出的结果
# gene <- read.table("turquoise_mark_gene.txt",header = T)
# table(gene$GeneName)
# gene <- gene[!duplicated(gene),]
ADvsHealthy_DEseq2 <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
#按照基因的平均表达水平对数据排序
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
#选择表达量最高的前5个基因
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]

#提取出这几个基因的FPKM值
gene <- ADvsHealthy_DEseq2$X
turquoise_gene <- ADvsHealthy_DEseq2$X
#找出这些基因的FPKM值，然后使得行作为样本，列作为基因的名字
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
#对FPKM的基因名进行ID转换
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
####将所有的差异基因做GO和KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]

##将差异基因的FPKM值提取出来
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)

#将样本的名字转换为AD和Healthy
#ADCYAP1,NRN1.RPH3A,SCG3,VGF--Braak_Stage
#ADCYAP1,RPH3A,VGF,SCG3--CERAD
#ADCYAP1,RPH3A--Phenotypedit
#blue模块的数据HIP1
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,11)]
library(tidyverse)
gene_group <- gene_group%>%mutate(Severity=case_when(gene_group$Severity==0~"Healthy",
                                                     gene_group$Severity==1~"AD",
                                                     gene_group$Severity==2~"MCI"))
gene_group$Severity <- as.factor(gene_group$Severity)
str(gene_group)
# gene_group <- gene_group[,order(colnames(gene_group))]
#绘制小提琴图
colnames(gene_group)
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]
table(gene_group$Severity)
# AD Healthy 
# 257      97 
library(ggsignif)
library(ggpubr)
library(glue)
plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_violin(aes_string(fill = .x),trim = FALSE)+
                   geom_signif(comparisons = list(c("AD", "Healthy")), 
                               tip_length = 0.02,
                               margin_top = 0.15,size = 0.8,textsize = 8,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#008080","#CC6600"))+
                   stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                geom="pointrange", color = "red")+
                   theme_bw() + 
                   xlab("Group") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))

setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'GE_turquoise_GeneViolin.pdf')





#画出不同的boxplot绘制出不同的基因，在不同的病理情况下基因表达的变化
# 绘制Braak Stage模块的图
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,12)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
str(gene_group)
colnames(gene_group)
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]
table(gene_group$Braak_Stage)
# 0   1   2   3 
# 17  24 164 149 
#批量绘制boxplot
plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+
                   geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),
                                                  c("1","2"),c("1","3"),c("2","3")), 
                               tip_length = 0.02,step_increase = 0.1,
                               margin_top = 0.001,size = 0.8,textsize = 4,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#3A5489","#F29B7F","#8491B4","#93CCC0"))+
                   theme_bw() + 
                   xlab("Braak Stage") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain" ),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'GE_BraakStage_turquoise.pdf')




#
# 绘制APOE模块的图
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,17)]
#去除NA值的那一行
gene_group <- na.omit(gene_group)
#对APOE那一列的变量名进行更改
colnames(gene_group)
gene_group <- gene_group%>%mutate(Phenotypedit=case_when(gene_group$Phenotypedit=="noE4"~"0",
                                                         gene_group$Phenotypedit=="E4carrier"~"1",
                                                         gene_group$Phenotypedit=="E4/4"~"2"))



gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
str(gene_group)
colnames(gene_group)
table(gene_group$Phenotypedit)
# 0   1   2 
# 194 122  33 
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]

#批量绘制boxplot
plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+
                   geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("1","2")), 
                               tip_length = 0.02,step_increase = 0.1,
                               margin_top = 0.001,size = 0.8,textsize = 4,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#3A5489","#F29B7F","#8491B4","#93CCC0"))+
                   theme_bw() + 
                   xlab("APOE4 Number") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain" ),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'GE_APOE4_turquoise.pdf')






#对绿色模块的数据进行绘制
#第二次分析得出的结果
blue_mark_gene <- append(blue_APOE,blue_BraakStage)
gene <- blue_mark_gene 
gene <- gene[!duplicated(gene)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ADvsHealthy_DEseq2 <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
#按照基因的平均表达水平对数据排序
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
#选择表达量最高的前5个基因
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]

#提取出这几个基因的FPKM值
gene <- ADvsHealthy_DEseq2$X
blue_gene <- ADvsHealthy_DEseq2$X
#找出这些基因的FPKM值，然后使得行作为样本，列作为基因的名字
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
#对FPKM的基因名进行ID转换
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
####将所有的差异基因做GO和KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]

##将差异基因的FPKM值提取出来
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)

#将样本的名字转换为AD和Healthy
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,11)]

gene_group <- gene_group%>%mutate(Severity=case_when(gene_group$Severity==0~"Healthy",
                                                     gene_group$Severity==1~"AD",
                                                     gene_group$Severity==2~"MCI"))
gene_group$Severity <- as.factor(gene_group$Severity)
str(gene_group)
table(gene_group$Severity)
# AD Healthy 
# 257      97
# gene_group <- gene_group[,order(colnames(gene_group))]
#绘制小提琴图
colnames(gene_group)
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]

plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_violin(aes_string(fill = .x),trim = FALSE)+
                   geom_signif(comparisons = list(c("AD", "Healthy")), 
                               tip_length = 0.02,
                               margin_top = 0.05,size = 0.8,textsize = 8,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#008080","#CC6600"))+
                   stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                geom="pointrange", color = "red")+
                   theme_bw() + 
                   xlab("Group") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain" ),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",filename = 'GE_Blue_GeneViolin.pdf')




#画出不同的boxplot绘制出不同的基因，在不同的病理情况下基因表达的变化
# 绘制Braak Stage模块的图
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,12)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
str(gene_group)
table(gene_group$Braak_Stage)
# 0   1   2   3 
# 17  24 164 149
colnames(gene_group)
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]
table(gene_group$Braak_Stage)
#批量绘制boxplot
plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+
                   geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),
                                                  c("1","2"),c("1","3"),c("2","3")), 
                               tip_length = 0.02,step_increase = 0.1,
                               margin_top = 0.0001,size = 0.8,textsize = 4,
                               map_signif_level=TRUE)+  
                   scale_fill_manual(values = c("#3A5489","#F29B7F","#8491B4","#93CCC0"))+
                   theme_bw() + 
                   xlab("Braak Stage") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain" ),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'GE_Blue_BraakStage.pdf')


# 绘制APOE模块的图
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,17)]
#去除NA值的那一行
gene_group <- na.omit(gene_group)
#对APOE那一列的变量名进行更改
colnames(gene_group)
gene_group <- gene_group%>%mutate(Phenotypedit=case_when(gene_group$Phenotypedit=="noE4"~"0",
                                                         gene_group$Phenotypedit=="E4carrier"~"1",
                                                         gene_group$Phenotypedit=="E4/4"~"2"))



gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
str(gene_group)
colnames(gene_group)
table(gene_group$Phenotypedit)
# 0   1   2 
# 194 122  33
# 同时对多个文件进行输出
x = names(gene_group)[6]
y = names(gene_group)[-6]

#批量绘制boxplot
plot_list = map2(x, y, 
                 ~ gene_group %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+
                   geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("1","2")), 
                               tip_length = 0.02,step_increase = 0.1,
                               margin_top = 0.0001,size = 0.8,textsize = 4,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#3A5489","#F29B7F","#8491B4","#93CCC0"))+
                   theme_bw() + 
                   xlab("APOE4 Number") + 
                   ylab("FPKM") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain" ),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",size = 0,linetype = "solid",colour = "black")))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'GE_Blue_APOE.pdf')

#对发现的基因进行ROC分析，数据下面的代码还没有进行更改
#注意使用去完批次之后的数据

#绘制ROC曲线的时候要使用的是去除完批次之后的数据
#对所有有差异的intron建立ROC曲线分析
#取出有差异的APA在不同的样本中的表达情况
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#加载去除完批次之后的数据
gene <- c("SYT1","CHN1","SNAP25","VSNL1","ENC1","TNS1","SGK1","CPM","CLMN","PPFIBP2")
# gene <- c(turquoise_gene,blue_gene)
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
#对FPKM的基因名进行ID转换
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
####将所有的差异基因做GO和KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

# 'select()' returned 1:many mapping between keys and columns
# Warning message:
#   In bitr(a, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db) :
#   41.77% of input gene IDs are fail to map...
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]

##将差异基因的FPKM值提取出来
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)

#将样本的名字转换为AD和Healthy
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:11,16)]

gene_group <- gene_group%>%mutate(Severity=case_when(gene_group$Severity==0~"Healthy",
                                                     gene_group$Severity==1~"AD",
                                                     gene_group$Severity==2~"MCI"))
gene_group$Severity <- as.factor(gene_group$Severity)
str(gene_group)
table(gene_group$Severity)
# AD Healthy 
# 257      97
b <- gene_group
colnames(b)[11] <- "Group"
#只取出来第一个数字,在制作分组信息的时候不能使用字母，要使用数字
library(dplyr)
b <- b%>%mutate(Group=case_when(b$Group=="Healthy"~"0",
                                b$Group=="AD"~"1"))

#进行ROC分析的基因的个数
n <- dim(b)[1]
y <- b$Group
all <- b

###开始随机抽样
set.seed(11)
require(caret)
folds <- createFolds(y,k=6)
library(ROCR)
library(magrittr) # pipe operator
library(plyr)
auc_value<-as.numeric()
#将所有的fitvalue取出来
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- all[-folds[[i]],] # 剩下的数据作为训练集
  model <- glm(as.numeric(fold_train$Group)~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- predict(model,type='response',newdata=fold_test)
  pred <- prediction(predictions = fold_predict, labels = fold_test$Group)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
auc_value
# [1] 0.7703488 0.8221289 0.6497093 0.7296512 0.7674419 0.6962209
mean(auc_value)
# [1] 0.7392502

#去掉数据框rocData的第一行
rocData <- rocData[-1,]
#更改数据框rocData的行名
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
#计算均值
# write.csv(rocData,file = "ASroc.csv",quote = F,row.names = F)
#将数据框中出现的NA用1进行补全
# rocData[is.na(rocData)] <- 1
#添加Xmean和Ymean的值
colnames(rocData)
rocData$Xmean <- rowMeans(rocData[,1:6])
rocData$Ymean <- rowMeans(rocData[,7:12])
rocData$Gene <- c(rep("GE", 60))
rocData <- rocData[,13:15]
#计算95%的置信区间
#计算不出来
# library(magrittr) # pipe operator
# rocit_emp <- rocit(score = rocData$Ymean, class = all$BackgroundInformation[test], method = "emp")
# summary(rocit_emp)
# #计算CI值
# ciAUC(rocit_emp)

# Xmean指的是所有X的平均值
# Ymean指的是所有Y的平均值
#构造样本ROC的终点数据
#构造中间虚线的数据
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("Control", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
#合并数据
ROC <- rbind(rocData,inn)
#画图
library(ggrepel)
library(ggplot2)
a<- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene)) + 
  geom_line(aes(linetype=Gene,color=Gene),size=1)+
  scale_linetype_manual(values = c(2,7))+
  scale_color_manual(values=c("black","#00872D"))+
  xlab("False Positive Rate") +
  ylab("True Positive Rate")+
  theme_bw()+
  ggtitle("GE ROC AUC=0.73925") +
  #scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=22),
        legend.title = element_text(size=8),
        legend.position = c(0.8,0.25),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain"),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"))
a
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)








#####在可变剪切的分析方面，表型数据的变化会影响多因素方差分析的结果#####
#将leafcutter产生的结果文件加载进来
rm(list = ls())
getwd()
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
a <- load("E:/AD_Patient/ASanalysiaLeafcutter/testtemporal_lobeADpatient.Rdata")
# ASgene <- as.vector(clusters$gene)
# library(stringr)
# m <- gsub("<i>","",ASgene)
# clusters$gene <- gsub("</i>","",m)

#将基因的count去除批次
#首先将count文件加载进来，并对文件的格式进行修改
#count文件每一行表示一个内含子，每一列表示一个样本
head(counts)
#对count文件的列名进行修改
a <- colnames(counts)
b <- substr(a,1,10)
colnames(counts) <- b
head(colnames(counts))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(counts,file = "RawIntroRead.csv",quote = F)
#保存exon的信息
exons_table <- as.data.frame(exons_table)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(exons_table,file = "exons_table.csv",quote = F)
intron <- introns[,1:7]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(intron,file = "Intron_To_ClusterNames.csv",quote = F)


#加载背景文件，将GSM号转换为AD还是Healthy
# APOE的基因型缺失的就认为是缺失值，没有进行填充
rm(list = ls())
getwd()
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
#从rowsum的结果中可以看出，leafcutter不会在cluster中包含0个read的行，leafcutter有自己的过滤标准，最少是30个reads
a <- as.data.frame(rowSums(counts))
head(counts[1:20,1:5],10)
#这里没有更改因为表型信息只是起到一个标签的作用
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
colnames(Temporal_Lobe)
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
#转置产生的sum文件太大，没法在R中展示出来，所以展示一小部分
head(colnames(sum),15)
head(rownames(sum),5)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
head(colnames(sum),15)
head(rownames(sum),5)
#再将转换之后的文件进行转置
counts <- sum
head(colnames(counts))
head(rownames(counts))
#将处理之后的文件进行PCA分析
#使用deseq中的pca函数进行分析
counts <- as.data.frame(t(counts))
condition <- factor(substring(colnames(counts), 1,1))
condition
coldata <- data.frame(row.names =colnames(counts),condition) 
coldata
library(DESeq2)  
dds <- DESeqDataSetFromMatrix(countData = counts,colData = coldata,design = ~ condition) 
vst <- vst(dds, blind=T)
head(assay(vst), 3)
library(ggforce)
library(ggrepel)
plotPCA(vst,intgroup = "condition")
p1data <- plotPCA(vst,returnData = T,intgroup = "condition")
#将p1的data与表型数据进行结合，在图中表示不同的表型数据
colnames(p1data)[5] <- "Sample"
pca_data <- merge(p1data,Temporal_Lobe,by="Sample")
pca_data$condition <- as.factor(pca_data$condition)
pca_data$Batch <- as.factor(pca_data$Batch)
p1 <- ggplot(data=pca_data,aes(x=PC1,y=PC2,color=condition,shape =Batch))+
  geom_point(size=3)+
  # stat_ellipse(level = 0.95) +
  # ggforce::geom_mark_ellipse(aes(color = group), 
  #                            fill = NA,  # 去除填充色
  #                            alpha = 0.2) +  # 设置边框透明度
  # geom_text_repel(data = p1data, aes(x = PC1, y = PC2, label = rownames(p1data)),
  #                 size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
  #                 color = "black",max.overlaps=100)+
  scale_x_continuous(limits = c(-30, 40))+
  scale_y_continuous(limits = c(-70,30))+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#FF9900","#FF3333","#66FF66","#3399CC","#663300"))+
  ggtitle("AS PCA")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",30,"%"),
       y=paste0("PCA2 ",19,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_RawCount_PCA.pdf", egg::set_panel_size(p1, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)



##采用多因素方差分析的方法对数据进行分析
df <- as.data.frame(t(counts))#保证df数据框的行是样本列是基因
head(df[1:5,1:5])
dim(df)
min(df)
max(df)
#将数据格式转换为因子
colnames(Temporal_Lobe)
Temporal_Lobe$Brain_Region <- as.factor(Temporal_Lobe$Brain_Region)
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
#年龄的分布区间很大，不能将其转换为因子
# Temporal_Lobe$Sex <- as.factor(Temporal_Lobe$Sex)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
##APOE是遗传上的风险因素
#进行多因素方差分析
str(Temporal_Lobe)
library(vegan)
set.seed(1)
#注意这里df表示的是count值，是没有经过标准化的数据
a <- adonis2(df ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
#将多因素方差的结果绘图
a#结果
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df ~ Braak_Stage + Age + Bank_Location + Batch + Severity, data = Temporal_Lobe, permutations = 999, method = "euclidean", by = "margin")
# Df   SumOfSqs      R2       F Pr(>F)    
# Braak_Stage     3 4.5945e+10 0.01624  3.4758  0.004 ** 
#   Age             1 4.2672e+09 0.00151  0.9685  0.408    
# Bank_Location   0 0.0000e+00 0.00000    -Inf           
# Batch           2 3.5721e+11 0.12624 40.5347  0.001 ***
#   Severity        1 2.0385e+09 0.00072  0.4627  0.830    
# Residual      345 1.5201e+12 0.53725                   
# Total         353 2.8295e+12 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="Age"~"Age",
                                   # b$category=="CERAD"~"CERAD",
                                   b$category=="Sex"~"Gender"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("AS Raw Count") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_RawCount.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)

#使用combat_seq对数据去除批次效应
# rm(list = ls())
# # 示例数据
# set.seed(123)
# data1 <- rnbinom(100, size = 10, mu = 20) # 生成负二项分布的数据
# data2 <- rnorm(100, mean = 50, sd = 10)   # 生成正态分布的数据
# # 负二项分布拟合
# nb_fit <- fitdistr(data1, "negative binomial")
# print(nb_fit)
# # 正态分布拟合
# norm_fit <- fitdistr(data2, "normal")
# print(norm_fit)
# # 绘制直方图
# hist(data1, breaks = 20, probability = TRUE, main = "Histogram of Two Distributions", 
#      xlab = "Value", col = rgb(0, 0, 1, 0.5), xlim = c(min(data1, data2), max(data1, data2)),
#      ylim = c(0, 0.05))
# # 绘制正态分布直方图在同一张图中
# hist(data2, breaks = 20, probability = TRUE, col = rgb(1, 0, 0, 0.5), add = TRUE)
# # 负二项分布拟合曲线
# curve(dnbinom(x, size = nb_fit$estimate["size"], mu = nb_fit$estimate["mu"]), 
#       col = "blue", lwd = 2, add = TRUE)
# # 正态分布拟合曲线
# curve(dnorm(x, mean = norm_fit$estimate["mean"], sd = norm_fit$estimate["sd"]), 
#       col = "red", lwd = 2, add = TRUE)
# # 添加图例
# legend("topright", legend = c("Negative Binomial", "Normal"), 
#        col = c("blue", "red"), lwd = 2)
# 
# 
# #第二种方法可以使用Shapiro-Wilk测试和Kolmogorov-Smirnov测试来检验数据是否符合正态分布
# #
# #使用模拟的负二项分布的数据检验正态性
# # Shapiro-Wilk 正态性检验
# shapiro_test <- shapiro.test(data1)
# print(shapiro_test)
# # Kolmogorov-Smirnov 正态性检验
# ks_test <- ks.test(data1, "pnorm", mean=norm_fit$estimate["mean"], sd=norm_fit$estimate["sd"])
# print(ks_test)
# 
# #使用正态分布的数据检验数据的正态性
# # Shapiro-Wilk 正态性检验
# shapiro_test <- shapiro.test(data2)
# print(shapiro_test)
# # Kolmogorov-Smirnov 正态性检验
# ks_test <- ks.test(data2, "pnorm", mean=norm_fit$estimate["mean"], sd=norm_fit$estimate["sd"])
# print(ks_test)
# 
# ###检验我的数据是否属于负二项分布#
# # 检验数据是否符合负二项分布可以通过以下几种方法
# # 拟合负二项分布并与数据比较： 使用最大似然估计（MLE）拟合负二项分布，并将拟合的分布与实际数据进行比较。这可以通过 fitdistr 函数来完成。
# # 绘制拟合曲线与数据对比： 比较数据的直方图与负二项分布的拟合曲线，以直观地检查数据是否符合负二项分布。
# # 卡方检验： 使用卡方检验比较数据的频率分布与负二项分布的理论分布。
# # 1. 拟合负二项分布并与数据比较
# # 生成示例数据
# rm(list = ls())
# set.seed(123)
# data <- rnbinom(100, size=10, mu=20) # 生成负二项分布的数据
# library(MASS)
# nb_fit <- fitdistr(data, "negative binomial")
# # 绘制数据的直方图
# hist(data, breaks=20, probability=TRUE, main="Histogram with Negative Binomial Fit", xlab="Value", col=rgb(0, 0, 1, 0.5))
# # 绘制负二项分布的拟合曲线
# curve(dnbinom(x, size=nb_fit$estimate["size"], mu=nb_fit$estimate["mu"]), add=TRUE, col="red", lwd=2)
# # 2. 卡方检验
# # 将数据分组
# breaks <- seq(min(data), max(data), length.out=10)
# data_cut <- cut(data, breaks=breaks, include.lowest=TRUE)
# observed_freq <- table(data_cut)
# # 计算理论频率
# # 生成每个区间的中点值
# interval_midpoints <- (breaks[-length(breaks)] + breaks[-1]) / 2
# theoretical_prob <- dnbinom(interval_midpoints, size=nb_fit$estimate["size"], mu=nb_fit$estimate["mu"])
# expected_freq <- theoretical_prob * length(data)
# # 执行卡方检验
# # 使用理论概率而不是理论频率
# chisq_test <- chisq.test(observed_freq, p=theoretical_prob / sum(theoretical_prob))
# # 打印检验结果
# print(chisq_test)
# 
# 
# 
# ####使用其中的一列数据进行检验#
# # 过滤掉负数或0
# counts_filtered <- counts[,1][counts[,1] > 0]
# nb_fit <- fitdistr(counts_filtered, "negative binomial")
# #但是我的数据出现了报错，
# # Error in stats::optim(x = c(1L, 3L, 3L, 5L, 4L, 2L, 3L, 1L, 1L, 3L, 4L,  : 
# #                               non-finite finite-difference value [1]
# #                             此外: Warning message:
# #                               In densfun(x, parm[1], parm[2], ...) : 产生了NaNs
# #然后发现出现这个报错的原因是这一列的数据中出现过多的0,因此准备换一个模型进行拟合
# #零膨胀负二项分布
# library(pscl)
# # 使用 zeroinfl 函数拟合零膨胀负二项分布
# zinb_fit <- zeroinfl(counts[,1] ~ 1 | 1, dist = "negbin")
# # 查看拟合结果
# summary(zinb_fit)
# #结果显示数据的负二项分布拟合的非常好
# # Count model coefficients (negbin with log link):
# #   Estimate Std. Error z value Pr(>|z|)    
# # (Intercept)  2.338521   0.008564   273.1   <2e-16 ***
# #   Log(theta)  -2.710692   0.004556  -595.0   <2e-16 ***
# 
# 
# #检查一下count数据的分布，决定使用combat还是combat_seq
# #负二项分布选择combat_seq
# #正态分布选择combat
# rm(list = ls())
# setwd("E:/AD_Patient/ASanalysiaLeafcutter")
# counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
# #绘制一下数据的分布——首先使用一列数据进行分析
# # 绘制直方图
# hist(counts[,1], breaks = 50, probability = TRUE, main = "Histogram of Two Distributions", 
#      xlab = "Value", col = rgb(0, 0, 1, 0.5), xlim = c(min(data1, data2), max(data1, data2)),
#      ylim = c(0, 0.05))
# curve(dnbinom(x, size = nb_fit$estimate["size"], mu = nb_fit$estimate["mu"]), 
#       col = "blue", lwd = 2, add = TRUE)
# #使用统计学方法检验数据的正态性和负二项性
# # Shapiro-Wilk 正态性检验，适用于少量的样本，样本数在3-5000之间
# shapiro_test <- shapiro.test(counts[,1])
# print(shapiro_test)
# # Kolmogorov-Smirnov 正态性检验——适用于大样本
# # 计算样本数据的均值和标准差
# mean_data <- mean(counts[,1])
# sd_data <- sd(counts[,1])
# # 执行 Kolmogorov-Smirnov 检验
# ks_test <- ks.test(counts[,1], "pnorm", mean=mean_data, sd=sd_data)
# # 打印检验结果
# print(ks_test)
# 
# #上面是进行一列数据分析的结果，下面对每一列的数据都进行分析。
# # Shapiro-Wilk 正态性检验
# shapiro_results <- list()
# # Kolmogorov-Smirnov 正态性检验
# ks_results <- list()
# # 遍历所有列进行正态性检验
# for (i in 1:ncol(counts)) {
#   data <- counts[,i]
#   # Shapiro-Wilk 检验
#   if (length(data) <= 5000 && !any(is.na(data))) {
#     shapiro_test <- shapiro.test(data)
#     shapiro_results[[paste("Column", i)]] <- shapiro_test
#   } else {
#     shapiro_results[[paste("Column", i)]] <- "Not applicable (sample size too large or contains NA values)"
#   }
#   # Kolmogorov-Smirnov 检验
#   mean_data <- mean(data, na.rm = TRUE)
#   sd_data <- sd(data, na.rm = TRUE)
#   ks_test <- ks.test(data, "pnorm", mean = mean_data, sd = sd_data)
#   ks_results[[paste("Column", i)]] <- ks_test
# }
# # 打印 Shapiro-Wilk 检验结果
# print("Shapiro-Wilk Test Results:")
# print(shapiro_results)
# # 打印 Kolmogorov-Smirnov 检验结果
# print("Kolmogorov-Smirnov Test Results:")
# print(ks_results)
# ###由于结果不容易看出来，所以以数据框的形式展现出来
# # 创建一个数据框来存储每列数据的 p 值
# ks_pvalues_df <- data.frame(
#   Column = names(ks_results),
#   P_Value = sapply(ks_results, function(result) {
#     if (inherits(result, "htest")) {
#       return(result$p.value)
#     } else {
#       return(NA)  # 如果结果不是检验对象，则返回 NA
#     }
#   })
# )
# 
# # 打印数据框
# print("Kolmogorov-Smirnov Test p-values for each column:")
# print(ks_pvalues_df)
# 
# #因此说明数据不是正态分布，检查数据是否符合负二项分布
# library(pscl)
# # 创建一个空的数据框，用于存储每列的检验结果
# results <- data.frame(Column = character(), Intercept_Estimate = numeric(),
#                       Intercept_Std_Error = numeric(), Intercept_P_value = numeric(),
#                       Theta = numeric(), Log_theta_p_value = numeric(),
#                       LogLik = numeric(), stringsAsFactors = FALSE)
# # 遍历 counts 数据框的每一列
# # for (i in 1:3) { #如果是只是想测试的话，先使用前3列的数据进行计算
# for (i in 1:ncol(counts[,1:3])) { #正规情况下使用所有的数据进行计算
#   # 对当前列数据进行零膨胀负二项分布拟合
#   zinb_fit <- tryCatch({
#     zeroinfl(counts[,i] ~ 1 | 1, dist = "negbin")
#   }, error = function(e) {
#     return(NULL)  # 如果某一列无法拟合，返回 NULL
#   })
#   # 如果拟合成功，提取拟合结果
#   if (!is.null(zinb_fit)) {
#     # 提取计数模型 (负二项分布) 的截距系数、标准误差、p 值
#     intercept_estimate <- coef(zinb_fit)[1]
#     intercept_se <- summary(zinb_fit)$coefficients$count[1, "Std. Error"]
#     intercept_p_value <- summary(zinb_fit)$coefficients$count[1, "Pr(>|z|)"]
#     Log_theta_p_value <- summary(zinb_fit)$coefficients$count[2, "Pr(>|z|)"]
#     # 提取 Theta 参数
#     theta <- zinb_fit$theta
#     # 提取对数似然值
#     logLik_value <- logLik(zinb_fit)
#     # 将结果添加到数据框中
#     results <- rbind(results, data.frame(Column = colnames(counts)[i], 
#                                          Intercept_Estimate = intercept_estimate,
#                                          Intercept_Std_Error = intercept_se, 
#                                          Intercept_P_value = intercept_p_value,
#                                          Theta = theta, 
#                                          Log_theta_p_value = Log_theta_p_value,
#                                          LogLik = as.numeric(logLik_value)))
#   } else {
#     # 如果拟合失败，将 NA 添加到结果中
#     results <- rbind(results, data.frame(Column = colnames(counts)[i], 
#                                          Intercept_Estimate = NA, 
#                                          Intercept_Std_Error = NA, 
#                                          Intercept_P_value = NA, 
#                                          Theta = NA, 
#                                          Log_theta_p_value= NA,
#                                          LogLik = NA))
#   }
# }
# # 查看所有列的拟合结果
# print(results)
# write.csv(results, file = "results_output.csv", row.names = FALSE)
# 

###综上可得，上面的数据不是正态分布的数据，属于负二项分布,因此采用combat_seq对数据进行去批次的处理
# rm(list = ls())
# setwd("E:/AD_Patient/ASanalysiaLeafcutter")
# counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
# dim(counts)
# counts <- counts[rowSums(is.na(counts)) == 0, ]  # 去除含有 NA 值的行
# counts <- counts[rowSums(counts)!=0, ]  # 去除全零行
# dim(counts)
# head(counts)
# #去除低表达量的基因
# # 设置一个最小表达阈值
# min_count_threshold <- 5  # 或者你可以设置为更高，比如 10
# sample_proportion_threshold <- 0.2  # 在30%的样本中基因的表达值超过阈值
# # 计算每个基因的样本中超过阈值的样本数
# genes_to_keep <- rowSums(counts > min_count_threshold) >= (sample_proportion_threshold * ncol(counts))
# table(genes_to_keep)
# head(genes_to_keep)
# # 过滤低表达基因
# filtered_counts <- counts[genes_to_keep, ]
# # 查看过滤后的数据
# dim(filtered_counts)  # 检查过滤后基因的数量
# filtered_counts <- filtered_counts[,order(colnames(filtered_counts))]
# #加载表型数据
# setwd("C:/Users/yujie/Desktop/datacollect")
# BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
# setwd("E:/AD_Patient/ASanalysiaLeafcutter")
# Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
# colnames(Temporal_Lobe)
# Temporal_Lobe <- Temporal_Lobe[order(Temporal_Lobe$GSM_number),]
# identical(colnames(filtered_counts),Temporal_Lobe$GSM_number)
# str(Temporal_Lobe)
# #将自己的数据格式转换为因子
# Temporal_Lobe$Brain_Region <- as.factor(Temporal_Lobe$Brain_Region)
# Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
# Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
# Temporal_Lobe$Sex <- as.factor(Temporal_Lobe$Sex)
# Temporal_Lobe$Bank_Location <-as.factor(Temporal_Lobe$Bank_Location)
# Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
# str(Temporal_Lobe)
# # 检查 batch 和 group 变量是否有NA
# length(Temporal_Lobe$Batch)
# length(Temporal_Lobe$Severity)
# library(sva)
# set.seed(1)
# 
# adjust_counts <- ComBat_seq(as.matrix(filtered_counts), batch=Temporal_Lobe$Batch,group=Temporal_Lobe$Severity,covar_mod=NULL,full_mod=FALSE)
# 
# #将校正之后的数据进行多因素方差分析
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# write.csv(adjust_counts,file = "Temporal_LobeRawCount_RemoveBatch_AS20240910.csv",quote = F)
# mycount <- as.data.frame(t(adjust_counts))
# mycount$GSM_number <- rownames(mycount)
# dim(mycount)
# 
# #加载尝试之前去除了批次的数据
# # setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
# # test <- read.csv("Temporal_LobeCount_filterAdjustCounts.csv",header = T,row.names = 1)
# # mycount <- as.data.frame(t(test))
# # mycount$GSM_number <- rownames(mycount)
# sum <- merge(mycount,Temporal_Lobe,by = "GSM_number", all = FALSE)
# dim(sum)
# sum[5,206360:206374]
# sum[5,1:5]
# rownames(sum) <- sum$Sample
# sum <- sum[,-c(1,206363:206374)]
# mycount <- as.data.frame(t(sum))
# condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) #构建condition因子，作为基因差异表达的自变量（即实验组和对照组的对比）。i和j分别为control和exp的个数
# condition
# coldata <- data.frame(row.names =colnames(mycount),condition)  #将database中各个实验名称加上condition标签，即说明是control还是exp
# coldata
# library(DESeq2)  #加载DESeq2包
# #构建dds及利用results函数得到最终结果
# dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  #countData用于说明数据来源，colData用于说明不同组数据的实验操作类型，design用于声明自变量，即谁和谁进行对比
# vst <- vst(dds, blind=T)
# library(ggforce)
# library(tidyverse)
# plotPCA(vst)
# ####下面将自己的PCA图画上标签,所以首先需要做的就是将PCA画图的data导出来
# p1data <- plotPCA(vst,returnData = T)
# colnames(p1data)[5] <- "Sample"
# 
# p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
# p1data$Severity <- as.factor(p1data$Severity)
# 
# p1data <- p1data%>%mutate(Severity=case_when(p1data$Severity==0~"Healthy",
#                                              p1data$Severity==1~"AD",
#                                              p1data$Severity==2~"MCI"))
# p1data$Batch <- as.factor(p1data$Batch)
# p1data$Sex <- as.factor(p1data$Sex)
# p1data <- p1data%>%mutate(Sex=case_when(p1data$Sex==0~"M",
#                                         p1data$Sex==1~"F"))
# p1data$Bank_Location <- as.factor(p1data$Bank_Location)
# p1data <- p1data%>%mutate(Bank_Location=case_when(p1data$Bank_Location==0~"America",
#                                                   p1data$Bank_Location==3~"Australia"))
# 
# p1data$Braak_Stage <- as.factor(p1data$Braak_Stage)
# p1data <- p1data%>%mutate(Braak_Stage=case_when(p1data$Braak_Stage==0~"0",
#                                                 p1data$Braak_Stage==1~"transentorhinal_stage",
#                                                 p1data$Braak_Stage==2~"limbic_stage",
#                                                 p1data$Braak_Stage==3~"isocortical_stage"))
# 
# # p1data$CERAD <- as.factor(p1data$CERAD)
# # p1data <- p1data%>%mutate(CERAD=case_when(p1data$CERAD==0~"C0",
# #                                           p1data$CERAD==1~"C1",
# #                                           p1data$CERAD==2~"C2",
# #                                           p1data$CERAD==3~"C3"))
# p1data$Age <- as.numeric(p1data$Age)
# library(ggrepel)
# library(ggplot2)
# library(ggforce)
# # display.brewer. all#展示所有的配色
# 
# # col2<-brewer.pal(5, "YlOrRd") #选色盘的几种
# # pal<-colorRampPalette(col2)
# # colors2<-pal(100)  #分成多少种颜色
# # # mycolor3<-brewer.pal(9, "YlOrRd")
# 
# 
# # ,shape=Bank_Location,
# p1 <- ggplot(data=p1data,aes(x=PC1,y=PC2,color=Severity,shape=Batch))+
#   geom_point(size=3)+
#   scale_x_continuous(limits = c(-40, 40))+
#   scale_y_continuous(limits = c(-45,30))+
#   ggtitle("Temporal_LobeRemoveBatch PCA Plot")+
#   scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
#   theme_bw()+
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"),
#         legend.title = element_text(color = "black",size = 8,face = "plain"),
#         legend.text = element_text(color = "black",size = 15,face = "plain"),
#         #legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"),
#         legend.direction = 'vertical', legend.position = "right")+
#   labs(x=paste0("PCA1 ",25,"%"),
#        y=paste0("PCA2 ",7,"%"))
# p1
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# ggsave("Temporal_Lobe_RemoveBatch.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 




# #加载FPKM的数据
# rm(list = ls())
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# myFPKM <- read.csv("TemporalLobeFPKM_RemoveBatch_20240910.csv",header = T,row.names = 1)
# myFPKM <- as.data.frame(t(myFPKM))
# myFPKM <- myFPKM[order(rownames(myFPKM)),]
# 
# 
# 
# ###对颞叶的数据要进行多因素方差分析
# setwd("C:/Users/yujie/Desktop/datacollect")
# BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
# setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
# BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
# BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
# 
# 
# 
# # 检查一下是否顺序就是一致的
# identical(rownames(myFPKM),BackgroundInformation$GSM_number)
# BackgroundInformation$Brain_Region <- as.factor(BackgroundInformation$Brain_Region)
# # BackgroundInformation$CERAD <- factor(BackgroundInformation$CERAD)
# BackgroundInformation$Braak_Stage <- factor(BackgroundInformation$Braak_Stage)
# BackgroundInformation$Sex <- factor(BackgroundInformation$Sex)
# BackgroundInformation$Bank_Location <- factor(BackgroundInformation$Bank_Location)
# BackgroundInformation$LibraryLayout <- factor(BackgroundInformation$LibraryLayout)
# BackgroundInformation$Severity <- factor(BackgroundInformation$Severity)
# BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
# str(BackgroundInformation)
# colnames(BackgroundInformation)
# library(vegan)
# # ?adonis()
# # adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Length+LibraryLayout, data = sample_information)
# # 如果你希望变量的顺序不影响结果，那么需要使用adonis2，并且设置参数by="margin"
# set.seed(1)
# a <- adonis2(as.data.frame(t(mycount)) ~ Severity+Braak_Stage+Age+Sex+Bank_Location+Batch, data = Temporal_Lobe,
#              permutations=999, by="margin",method = "manhattan")
# a
# # Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = myFPKM ~ Severity + Braak_Stage + Age + Sex + Bank_Location + Batch, data = BackgroundInformation, permutations = 999, method = "manhattan", by = "margin")
# Df   SumOfSqs      R2      F Pr(>F)  
# Severity        1 1.5860e+11 0.01276 4.6166  0.020 *
#   Braak_Stage     3 1.2716e+11 0.01023 1.2339  0.293  
# Age             1 8.1751e+10 0.00658 2.3797  0.111  
# Sex             1 2.5663e+10 0.00206 0.7470  0.408  
# Bank_Location   0 0.0000e+00 0.00000   -Inf         
# Batch           2 8.8078e+10 0.00709 1.2819  0.238  
# Residual      344 1.1818e+13 0.95071                
# Total         353 1.2430e+13 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# b <- as.data.frame(a)
# b <- b[-c(8,7),]
# b$category <- rownames(b)
# colnames(b)
# colnames(b)[5] <- "Pr"
# #自己手动添加显著性标记的记号
# library(tidyverse)
# b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
#                              0.05<b$Pr & b$Pr <= 0.1 ~ ".",
#                              0.01 < b$Pr & b$Pr <= 0.05 ~"*",
#                              0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
#                              0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
# b$category
# colnames(b)
# #将category中的下划线去掉
# b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
#                                    b$category=="Braak_Stage"~"Braak Stage",
#                                    b$category=="Severity"~"Severity",
#                                    b$category=="Batch"~"Batch",
#                                    b$category=="Age"~"Age",
#                                    # b$category=="CERAD"~"CERAD",
#                                    b$category=="Sex"~"Gender"))
# 
# 
# #将b中的数据进行排序
# b <- b[order(b$R2,decreasing = TRUE),]
# b$category <- factor(b$category)
# library(ggsci)
# p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
#   geom_bar(stat="identity")+
#   scale_fill_npg()+
#   geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
#   xlab("category") +ylab("R2")+
#   theme_bw()+#theme_classic()+
#   ggtitle("MANOVA") +
#   #coord_flip()+
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=15),
#         legend.title=element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 18,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 18,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
#         axis.text.y = element_text(color = "black",size = 18,face = "plain"))
# 
# p
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# ggsave("Temporal_LobeGeRemoveBatch.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 




#使用limma去除批次效应,首先对数据进行normalized
library(limma)
library(edgeR)
rm(list = ls())
# getwd()
# setwd("E:/AD_Patient/ASanalysiaLeafcutter")
# setwd("E:/virus")
# virus <- read.csv("Virus.infected.csv",header = T,row.names = 1)
# virus <- as.data.frame(t(virus))
# virus <- virus[,order(colnames(virus))]
# bacground <- read.csv("virusBackground.csv",header = T)
# bacground <- bacground[order(bacground$sid),]
# identical(colnames(virus),bacground$sid)
# #只取出来颞叶的数据
# bacground <- subset(bacground,bacground$Brain_Region=="3")
# a <- bacground$sid
# virus <- virus[,colnames(virus) %in% a]
# identical(colnames(virus),bacground$sid)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
counts <- counts[,order(colnames(counts))]
#将 GSE号转换为AD还是Healthy的分组
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
Temporal_Lobe <- Temporal_Lobe[order(Temporal_Lobe$GSM_number),]
#检查样本名的排列是否一致
identical(colnames(counts),Temporal_Lobe$GSM_number)


#对count的列名进行修改
colnames(Temporal_Lobe)
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
#转置产生的sum文件太大，没法在R中展示出来，所以展示一小部分
head(colnames(sum),15)
head(rownames(sum),5)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
head(colnames(sum),15)
head(rownames(sum),5)

#再将转换之后的文件进行转置
counts <- as.data.frame(t(sum))
#构建voom所需要的数据格式
#对样本的顺序进行排序
group <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(counts))) #
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design)=colnames(counts)

#去掉一些测序深度很低的样本
a <- as.data.frame(colSums(counts))
min(a$`colSums(counts)`)
#结果发现测序深度都不低，所以不考虑去除样本，不能随便去除某些基因，因为一个样本的reads总数会对结果产生影响

# 去除低表达的基因
keep.exprs <- filterByExpr(counts, group=group,
                           min.count = 3, min.total.count = 30,
                           large.n = 10, min.prop = 0.3)

# 之前large.n =3
table(keep.exprs)
# keep.exprs
# FALSE   TRUE 
# 129501  76860
counts <- counts[keep.exprs,]
dim(counts)


#将count的data.frame转换为DGEList格式
dge<-DGEList(counts)
library(edgeR)
#Calculate normalization factors to scale the raw library sizes.
#用TMM进行标准化
# dge<-calcNormFactors(dge)
# #cpm转换
# logCPM<-cpm(dge,log=TRUE,prior.count=3)
# min(logCPM)
# 避免文库大小在样本间变化的影响
#而我们limma本身就提出了一个voom的方法来对RNA-seq数据进行normalization
y <- voom(dge, design, plot = T,normalize.method ="quantile")
#因为原始的voom函数会产生负值，所以对voom函数进行修改
myvoom <- edit(voom)
#这里lib.size <- colSums(counts)，
#说明在使用voom的时候一组数据列的和对最终的结果是有影响的，所以做的时候要考虑到这个因素
#将第45行的y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
# 更改为y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+10))
#之所以改为10，是因为当改为7，8，9的时候，用完removeBatchEffect函数之后仍然会出现负值，一直到10才没有负值出现
y <- myvoom(dge, design, plot = T,normalize.method ="quantile")
newData <- y$E
min(newData)#对数据进行预处理之后的数据
dim(newData)
NormalizedCount <- as.data.frame(newData)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(NormalizedCount,file = "AS_VoomedCount.csv",quote = F)





#将没有去除批次效应，但是对count进行了标准化的数据绘制PCA图
#原始数据
dim(newData)
df <- as.data.frame(t(newData))#保证df数据框的行是样本列是基因
library(gmodels)
pca.info <- prcomp(df)
head(pca.info$x)
head(summary(pca.info)) #能看到每个PC具体的解释比例
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
#加入批次的信息
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
a$Batch <- as.factor(a$Batch)
a$Severity <- as.factor(a$Severity)
a <- a[order(a$Sample),]
ggscatter(a,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")

#使用ggplot2对图片进行美化
library(ggrepel)
library(ggplot2)
library(ggforce)
library(dplyr)
# 将0换为Healthy，1换为AD
a <- a%>%mutate(Type=case_when(a$Type==0~"Healthy",
                               a$Type==1~"AD",
                               a$Type==2~"MCI"))


#画图
min(a$PC1)
max(a$PC1)
min(a$PC2)
max(a$PC2)
p1 <- ggplot(data=a,aes(x=PC1,y=PC2,color=Type,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
  scale_x_continuous(limits = c(-400,210))+
  scale_y_continuous(limits = c(-200,200))+
  ggtitle("AS Voom PCA")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",16,"%"),
       y=paste0("PCA2 ",10,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomCount_PCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)


#对Normalized之后的count进行多因素方差分析，看其中哪个因素的影响最大
#多因素方差分析要求横坐标是样本，纵坐标是基因
dim(newData)
df <- as.data.frame(t(newData))#保证df数据框的行是样本列是基因
dim(df)
min(df)
max(df)
#将数据格式转换为因子
colnames(Temporal_Lobe)
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
# Temporal_Lobe$Sex <- as.factor(Temporal_Lobe$Sex)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
identical(rownames(df),Temporal_Lobe$Sample)
#进行多因素方差分析
library(vegan)
set.seed(1)
dim(df)
a <- adonis2(df ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
a
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = df ~ Braak_Stage + Age + Bank_Location + Batch + Severity, data = Temporal_Lobe, permutations = 999, method = "euclidean", by = "margin")
# Df SumOfSqs      R2       F Pr(>F)    
# Braak_Stage     3   358572 0.00916  1.4466  0.024 *  
#   Age             1   107441 0.00274  1.3003  0.147    
# Bank_Location   0        0 0.00000    -Inf           
# Batch           2  5669205 0.14477 34.3061  0.001 ***
#   Severity        1   302038 0.00771  3.6555  0.001 ***
#   Residual      345 28506204 0.72796                   
# Total         353 39159099 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#使用柱状图的形式展现出来
b <- as.data.frame(a)
b <- b[-c(7,6),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill =category ))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("AS Voom Count") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=12),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomCount.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)

#使用removeBatchEffect去除批次
dim(newData)#数据的格式是行是基因名列是样本名
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# VoomedCount <- read.csv("VoomedCount.csv",header = T,row.names = 1,check.names = F)
design <- model.matrix(~Severity + Age , data=Temporal_Lobe)
head(design)
treatment.design <- design[,1:2]#之所以保存两列是因为有1列是常数项，另外1列是severity
head(treatment.design)
batch.design <- design[,-(1:2)]
head(batch.design)
dim(newData)
head(design)
corrected_count <- removeBatchEffect(newData,batch = Temporal_Lobe$Batch,
                                     group = Temporal_Lobe$Severity,
                                     covariates=batch.design,
                                     design = treatment.design)
min(corrected_count)
# [1] 2.23923
RemoveData <- as.data.frame(t(corrected_count))
RemoveData <- as.data.frame(t(RemoveData))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(RemoveData,file = "AS_RemoveBatch_VoomCount.csv",quote = F)

#进行PCA
dim(RemoveData)
RemoveData <- as.data.frame(t(RemoveData))
identical(rownames(RemoveData),Temporal_Lobe$Sample)
library(gmodels)
pca.info <- prcomp(RemoveData)#要保证行是样本，列是基因
head(pca.info$x)
head(summary(pca.info)) #能看到每个PC具体的解释比例
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
#加入批次的信息
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
a$Batch <- as.factor(a$Batch)
a$Severity <- as.factor(a$Severity)
a <- a[order(a$Sample),]
ggscatter(a,x="PC1", y="PC2", color="Severity", ellipse=FALSE, ellipse.type="confidence")

#使用ggplot2对图片进行美化
library(ggrepel)
library(ggplot2)
library(ggforce)
library(dplyr)
# 将0换为Healthy，1换为AD
a <- a%>%mutate(Type=case_when(a$Type==0~"Healthy",
                               a$Type==1~"AD",
                               a$Type==2~"MCI"))


#画图
min(a$PC1)
max(a$PC1)
min(a$PC2)
max(a$PC2)
p1 <- ggplot(data=a,aes(x=PC1,y=PC2,color=Type,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
  scale_x_continuous(limits = c(-200,300))+
  scale_y_continuous(limits = c(-250,200))+
  ggtitle("AS Voom RemoveBatch PCA")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",12,"%"),
       y=paste0("PCA2 ",9,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_Voom_RemoveBatchPCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)



#对去除了批次的数据进行多因素方差分析
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# RemoveData <- read.csv("RemoveBatchCount.csv",header = T,check.names = F,row.names = 1)
dim(RemoveData)#保证df数据框的行是样本列是基因
#将数据格式转换为因子
colnames(Temporal_Lobe)
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
# Temporal_Lobe$Sex <- as.factor(Temporal_Lobe$Sex)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)

identical(rownames(RemoveData),Temporal_Lobe$Sample)

library(vegan)
set.seed(1)
a <- adonis2(RemoveData ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
a
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = RemoveData ~ Braak_Stage + Age + Bank_Location + Batch + Severity, data = Temporal_Lobe, permutations = 999, method = "euclidean", by = "margin")
# Df SumOfSqs      R2      F Pr(>F)    
# Braak_Stage     3   358572 0.01217 1.4466  0.007 ** 
#   Age             1    11756 0.00040 0.1423  1.000    
# Bank_Location   0        0 0.00000   -Inf           
# Batch           2    23548 0.00080 0.1425  1.000    
# Severity        1   302038 0.01025 3.6555  0.001 ***
#   Residual      345 28506204 0.96782                  
# Total         353 29454172 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# #使用柱状图的形式展现出来
b <- as.data.frame(a)
b <- b[-c(7,6),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill = category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("AS Voom RemoveBatch") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=12),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomRemoveBatch.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 9, height = 7, units = 'in', dpi = 600)



#找出有差异的AS的exon
#因为前期做了处理，所以不会出现负值的情况，因此可以直接使用t.test
#对去除了批次之后的数据进行AD和Healthy的比较
dim(RemoveData)
RemoveData <- as.data.frame(t(RemoveData))
dim(RemoveData)
#将数据转换为行表示的是转录本，列表示的是不同的样本
# RemoveData <- as.data.frame(t(RemoveData))
table(substr(colnames(RemoveData),1,1))
#计算出每一组的平均值
colnames(RemoveData)
#对列名进行排序，使healthy的个体排在前97个，AD个体的数据排在后面
RemoveData <- RemoveData[,order(colnames(RemoveData))]
colnames(RemoveData)
#排序完成后分别对AD和Healthy进行分析
RemoveData$groupHealthymean <- apply(RemoveData[,1:97],1,mean)
RemoveData$groupADmean <- apply(RemoveData[,98:354],1,mean)
head(RemoveData$groupHealthymean)
head(RemoveData$groupADmean)
#计算出差值
RemoveData$Group_diff <- RemoveData$groupADmean-RemoveData$groupHealthymean
head(RemoveData$Group_diff)

#将两个数据分别
RemoveData$log2FoldChange <- log2(RemoveData$groupADmean/RemoveData$groupHealthymean)
head(RemoveData$log2FoldChange)
#使用计算出来的均值计算显著性
#由于样本的个数太少，所以没办法使用
#t.test(RemoveData$groupAmeanPDUI[1],RemoveData$groupBmeanPDUI[1])

#然后使用每一个样本的PDUI值计算p值的显著性
#例子
t.test(RemoveData[1,1:97],RemoveData[1,98:354])

# Welch Two Sample t-test
# 
# data:  RemoveData[1, 1:97] and RemoveData[1, 98:354]
# t = -4.4112, df = 169.73, p-value = 1.82e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4459502 -0.1702149
# sample estimates:
#   mean of x mean of y 
# 17.13335  17.44143 

#所以使用所有的数据计算p值
#不能使用的是配对的t检验(因为样本的个数不一致)
RemoveData$Pvalue <- apply(RemoveData,1,function(x) t.test(x[1:97],x[98:354])$p.value)
head(RemoveData$Pvalue)
min(RemoveData$Pvalue)
#对计算出的p值进行校正
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
head(RemoveData$Padjust)
min(RemoveData$Padjust)
RemoveData$Position <- rownames(RemoveData)

#对基因的名字进行注释
name <- read.table("ASdiff_n344.txt",header = T,row.names = 1)
name$Position <- paste(rownames(name),name$Cluster,sep = ":")
AS_anno <- name[,c(8,10)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AS_anno,file = "AS_anno.csv",quote = F)
a <- merge(AS_anno,RemoveData,by="Position")
a<- a[,-2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(a,file = "AS_VoomRemoveBatch_DEG.csv",quote = F,row.names = F)



# #检查一下k=8和k=15时得到的差异基因有什么不同
# rm(list = ls())
# a <- read.csv("AS_DEG.csv")
# #最开始得到的差异ASintrons
# b <- read.csv("AS_test_DEG.csv")
# a <- a[order(a$X),]
# b <- b[order(b$X),]
# identical(a$X,b$X)
# 
# #取出有差异的exon
# rm(list = ls())
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# a <- read.csv("AS_DEG.csv",header = T,row.names = 1,check.names = F)
# #绘制火山图
# a$TP <- rownames(a)
# head(a$TP)
# library(tidyverse)
# a <- separate(a,TP,into = c("Chr","Start","End","Cluster"),sep = "([:])")
# #将上面进行判断的3列数据，合并为1列数据，然后从大的数据集m中寻找具体对应的是哪个基因
# a$Position <- paste(a$Chr,a$Start,a$End,sep = ":")
# a <- a[,-c(7:9)]
# rownames(a) <- NULL
# length(unique(a$Position))
# #对这些exon进行注释
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# m <- read.table("humanv40annotation_all_introns.bed",fill = TRUE)
# dim(m)
# #对找出的有差异的intron进行数据的分割
# m$Position <- paste(m$V1,m$V2,m$V3,sep = ":")
# test <- subset(m,m$Position %in% "chrX:305197:307360")
# #说明m中的position数据存在重复值
# m <- m[!duplicated(m$Position),]
# dim(m)
# #从part1的结果可以看出，存在5715没有被注释到的基因
# #但是其中存在差异表达的基因，因为存在并不是human genome的端点处
# Part1 <- subset(m,m$Position %in% a$Position)#71145
# Part1 <- Part1[,c("V4","Position")]
# colnames(Part1) <- c("Genename","Position")
# #对异常的地方进行手动注释
# diff_values_1 <- setdiff(a$Position,m$Position)
# #取出一部分的数据进行分析
# # 
# # diff_values_1 <- diff_values_1[1:5]
# 
# 
# split_vec <- strsplit(diff_values_1, ":")
# df <- do.call(rbind, split_vec)
# colnames(df) <- c("Chr", "Start","End")
# df <- as.data.frame(df)
# df$Start <- as.numeric(df$Start)
# df$End <- as.numeric(df$End)
# head(df)
# 
# 
# total <- m[,1:4]
# head(total)
# # 分组并计算V2和V3合并后的最小值和最大值
# total$V1 <- as.character(total$V1)
# total$V2 <- as.numeric(total$V2)
# total$V3 <- as.numeric(total$V3)
# total$V4 <- as.vector(total$V4)
# colnames(total)  <- c("Chr","Start","End","SYMBOL")
# library(dplyr)
# # result <- total %>%
# #   group_by(V4) %>%
# #   summarise(
# #     V1 = V1,                  # 保留每组中的第一个 V1
# #     min_value = min(c(V2, V3)),      # 计算 V2 和 V3 中的最小值
# #     max_value = max(c(V2, V3)))
# # result <- result[!duplicated(result$V4),]
# # result <- as.data.frame(result)
# # colnames(result)  <- c("SYMBOL","Chr","Start","End")
# # head(result)
# # str(result)
# # df 是需要注释的数据框，result 是包含基因区间的结果数据框
# # 初始化注释列
# df$Gene_Annotation <- NA
# annotate_genes <- function(chr, start, end, total) {
#   # 第一步：筛选染色体相等的记录
#   chr_match <- total %>%
#     filter(Chr == chr)
#   
#   # 初始化匹配的基因名
#   matching_genes <- c()
#   
#   # 第二步：逐行比较 df 的每一行与 result 中的记录
#   for (i in 1:nrow(chr_match)) {
#     # 当前的 result 行
#     result_row <- chr_match[i, ]
#     
#     # 检查 df 与 result_row 的四种关系：
#     if ((start >= result_row$Start && end <= result_row$End) ||  # 完全包含
#         (start <= result_row$End && end >= result_row$Start)) {  # 部分重叠
#       # 满足任何条件，记录下基因名
#       matching_genes <- c(matching_genes, result_row$SYMBOL)
#     }
#   }
#   # 去掉重复的基因名
#   # matching_genes <- unique(matching_genes)
#   # 如果找到匹配的基因名，返回它们
#   if (length(matching_genes) > 0) {
#     return(paste(matching_genes, collapse = ","))
#   } else {
#     return(NA)
#   }
# }
# 
# 
# # 使用 annotate_genes 函数来为 df 数据框中的每一行注释基因名
# df <- df %>%
#   rowwise() %>%
#   mutate(Gene_Annotation = annotate_genes(Chr, Start, End, total))
# # 查看结果
# head(df)
# # 去除每一行指定列的重复值，并将结果更新回到该列
# df <- df %>%
#   mutate(column_name = sapply(strsplit(as.character(Gene_Annotation), ","), function(x) paste(unique(x), collapse = ",")))
# #统计没有被注释上的个数
# sum(is.na(df$Gene_Annotation))
# #将结果保存
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# write.table(df,file = "AS_part2_Anno.txt",sep = "\t")
# 
# ###结果发现注释的时候还是有500多个exon没有注释到，因此选择手动注释，将注释的结果保存在
# # C:\Users\yujie\Desktop\datacollect\20240827\Complement_Part2.txt中
# # 其中remove表示的是包括两个基因的区间
# # no_gene表示的是基因组的空白区域，出现这个的原因是因为在写注释函数的时候忘掉了空白区域这一种情况
# #加载文件AS_part2_Anno.txt和Complement_Part2.txt，并将两个文件进行合并
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# AS_part2_Anno <- read.table("AS_part2_Anno.txt",row.names = 1,header = T)
# AS_part2_Anno <- AS_part2_Anno[,c(1:3,5)]
# AS_part2_Anno$Position <- paste(AS_part2_Anno$Chr,AS_part2_Anno$Start,AS_part2_Anno$End,sep = ":")
# Complement_Part2 <- read.table("Complement_Part2.txt",header = T)
# Complement_Part2$Position <- paste(Complement_Part2$Chr,Complement_Part2$Start,Complement_Part2$End,sep = ":")
# #填充数据集
# library(tidyverse)
# df_filled <- AS_part2_Anno %>%
#   left_join(Complement_Part2, by = "Position") %>%
#   # 如果 column_name 是 NA，则使用 Genename 填补
#   mutate(column_name = ifelse(is.na(column_name), Genename, column_name)) %>%
#   # 删除不再需要的 Genename 列
#   select(-Genename)
# 
# sum(is.na(df_filled$column_name))
# Part2 <- df_filled[,4:5]
# colnames(Part2) <- c("Genename","Position")
# Part <- rbind(Part1,Part2)
# #对a，即AS DEG的结果进行注释
# AS_Anno <- merge(a,Part,by="Position")
# #将Remove的值换掉,换成no_gene
a <- merge(AS_anno,RemoveData,by="Position")
a$Genename <- ifelse(a$Genename == "REMOVE", "no_gene", a$Genename)
sum <- a
colnames(sum)
sum$change <- ifelse(sum$Padjust < 0.05 & abs(sum$log2FoldChange) >= 0.1,
                     ifelse(sum$log2FoldChange > 0.1 ,'Up','Down'),'Stable')
table(sum$change)
# Down Stable     Up 
# 314  76491     55
#之前是344
write.table(sum,file = "ASdiff_n369.txt",quote = F,row.names = F,sep = "\t")
#绘制火山图
AS_diff <- read.table("ASdiff_n369.txt",header = T,row.names = 1,check.names = F)
sum <- AS_diff
library(dplyr)
up <- subset(sum, sum$change == 'Up')
down <- subset(sum, sum$change == 'Down')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:5, ]
down <- down[order(down$Padjust), ][1:5, ]
a <- rbind(up, down)
colnames(a)
library(ggplot2)
library(ggrepel)
p <- ggplot(
  # 数据、映射、颜色
  sum, aes(x = log2FoldChange, y = -log10(Padjust),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080", "gray", "firebrick3")) +
  geom_text_repel(data = a, aes(x = log2FoldChange, y = -log10(Padjust), label = Genename),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
  # 坐标轴
  labs(x="log2FoldChange",y="-log10(padj)") +
  ggtitle("AS volcano")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.position ="right",
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.background = element_rect(fill = "transparent",size = 0.4,linetype = "solid"))+
  scale_x_continuous(breaks = seq(-10, 10, 2))+
  scale_y_continuous(breaks = seq(0, 25, 5))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_Volcano.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)




#绘制差异AS的热图
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#读取数据
name <- read.table("ASdiff_n369.txt",header = T,check.names = F)
# name$Position <- paste(rownames(name),name$Cluster,sep = ":")
colnames(name)
sum <- name
up <- subset(sum, sum$change == 'Up')
down <- subset(sum, sum$change == 'Down')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:5, ]
down <- down[order(down$Padjust), ][1:5, ]
MostSigGene <- rbind(up, down)
#找出MostSigGene对应的基因的名字
AS_Heatmap <- sum[sum$change!="Stable",]
head(AS_Heatmap)
rownames(AS_Heatmap) <- AS_Heatmap$Position
AS_Heatmap <- AS_Heatmap[order(AS_Heatmap$Padjust),]
colnames(AS_Heatmap)
AS_Heatmap <- AS_Heatmap[,3:356]

sum <- as.data.frame(t(AS_Heatmap))
#对样本进行排序，使得AD和Healthy的样本分为上下两个不同的组
sum <- arrange(sum,rownames(sum))
data <- sum
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){
  sd(x)/mean(x)*100
})
data_df <- data[order(data$cv, decreasing = T),1:354]
dim(data_df)
a <- apply(data_df,1,scale)
# 进行scale，但是不设定最大最小值，这个时候的热图一片绿
data_scale <- as.data.frame(t(apply(data_df,1,scale))) ##Z-score标准化
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.01
#将scale之后的数据设置最大值最小值
data_scale <- as.matrix(data_scale)
#设置最大值与最小值
table((data_scale)>1)
table((data_scale)<(-1))
data_scale[data_scale>=1]=1
data_scale[data_scale<=-1]=-1
table(is.na(data_scale))
colnames(data_scale)
# library(pheatmap)
# pheatmap(data_scale)
library(ComplexHeatmap)
library(circlize)

#读取待展示的基因名称，并添加到热图中
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "AS_Heatmap3.pdf",width =3,height = 4)
# png("p.png",res = 300)

p <- Heatmap(data_scale,name = "Expression", 
             na_col = "grey",
             cluster_rows = TRUE,
             clustering_distance_rows = "pearson",
             clustering_method_rows = "complete",
             col = colorRampPalette(c("#5AA5A4", "white", "#B52836"))(20),
             cluster_columns = FALSE,
             clustering_distance_columns = "pearson",
             clustering_method_columns = "complete",
             row_dend_side = "left",
             show_row_names = FALSE,
             show_column_names = TRUE,
             use_raster = T,
             heatmap_width = unit(1, "npc"),
             width = NULL,
             heatmap_height = unit(1, "npc"),
             height = NULL,
             show_column_dend = FALSE,
             column_dend_height = unit(5, "mm"),
             show_row_dend = FALSE,
             column_labels = colnames(data_scale),
             column_names_side = "bottom",
             column_names_centered = TRUE,
             column_title = "",
             column_names_rot = 0,
             column_title_gp = gpar(fontsize = 0.1),
             column_names_gp = gpar(fontsize = 1),
             heatmap_legend_param = list( #设置标签
               color_bar = 'continuous',
               legend_direction = 'vertical',
               legend_width = unit(2, 'cm'),
               legend_height = unit(2, 'cm'),
               title_position = 'topcenter',
               title_gp = gpar(fontsize = 5, fontface = 'plain'),
               labels_gp = gpar(fontsize = 5, fontface = 'plain'))
)+ 
  rowAnnotation(link = anno_mark(at = which(rownames(data_scale) %in% MostSigGene$Position), 
                                 labels = MostSigGene$Genename[match(rownames(data_scale)[rownames(data_scale) %in% MostSigGene$Position], MostSigGene$Position)], labels_gp = gpar(fontsize = 7)))

print(p)
dev.off()

#

#进行AS数据的差异基因分析
#取出有差异的exon
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F)
a <- a[a$change!="Stable",]
#过滤数据
#首先删除Remove和no_gene的行
a <- a[a$Genename!="no_gene",]
#删掉那些匹配到多个基因的exon
# 使用 grepl 查找包含逗号的行，取反以保留没有逗号的行
a <- a[!grepl(",", a$Genename),]
sum <- a
a$Genename
#做GO和KEGG分析
gene <- unique(a$Genename)
#237个基因
# gene <- nchar(unique(sum$V4))#nchar是统计这个向量里每一个单个元素的字符串的个数
#对具有可变剪切的基因进行富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(gene, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
# 'select()' returned 1:1 mapping between keys and columns
# Warning message:
#   In bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) :
#   3.57% of input gene IDs are fail to map...
colnames(df1)
GOgene <- as.character(df1$ENTREZID)#得出所有的ENTREZID
GOgene <- na.omit(GOgene)#得出所有的ENTREZID，就可以进行GO和kegg的分析
BPplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='BP',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
BPplot_genelist<-setReadable(BPplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
BPplot_genelist <- as.data.frame(BPplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(BPplot_genelist,file = "AS_BPplot_genelist.csv",row.names = F,quote = F)

#CC
CCplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='CC',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
CCplot_genelist<-setReadable(CCplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CCplot_genelist <- as.data.frame(CCplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(CCplot_genelist,file = "AS_CCplot_genelist.csv",row.names = F,quote = F)

#MF
MFplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='MF',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
MFplot_genelist<-setReadable(MFplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
MFplot_genelist <- as.data.frame(MFplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(MFplot_genelist,file = "AS_MFplot_genelist.csv",row.names = F,quote = F)

###绘制合在一块的GO图
go_enrich_df <- data.frame(
  ID=c(BPplot_genelist$ID[1:15], CCplot_genelist$ID[1:15], MFplot_genelist$ID[1:15]),
  Description=c(BPplot_genelist$Description[1:15],CCplot_genelist$Description[1:15],MFplot_genelist$Description[1:15]),
  GeneNumber=c(BPplot_genelist$Count[1:15], CCplot_genelist$Count[1:15], MFplot_genelist$Count[1:15]),
  type=factor(c(rep("biological process", 15), 
                rep("cellular component", 15),
                rep("molecular function", 15)), 
              levels=c("biological process", "cellular component","molecular function" )))
go_enrich_df <- na.omit(go_enrich_df)
go_enrich_df <- arrange(go_enrich_df,type,GeneNumber)
##开始绘制GO柱状图
###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
# write.csv(go_enrich_df,file = "go_enrich_dfAS.csv",quote = F,row.names = F)
COLS <- c("#0D6CA6","#099963", "#911F27")#设定颜色
library(ggpubr)
p <- ggdotchart(go_enrich_df, x = "type_order", y = "GeneNumber",
                color = "type",                                # 按照cyl填充颜色
                palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                sorting = "descending",                      
                add = "segments",                             # 添加棒子
                add.params = list(color = "type", size = 1.3),#改变棒子参数
                # rotate = TRUE,                                # 方向转为垂直
                group = "type",                                
                dot.size = "GeneNumber",                                 # 改变点的大小
                #label = round(go_enrich_df$GeneNumber),                       # 添加label
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               # 设置label参数
                ggtheme = theme_pubr(),                        # 改变主题
                xlab="GO Term",
                ylab="GeneNumber",
                title = "AS GO")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40,exdent = 0),"\n")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 14,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 30),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain",angle = 70,vjust = 1, hjust = 1 ),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"),
        legend.position = c(0.9,0.68),
        legend.title = element_text(color = "black",size = 15,face = "plain"),
        legend.text = element_text(color = "black",size = 20,face = "plain"))
#legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_GO.pdf", egg::set_panel_size(p, width=unit(16, "in"), height=unit(5, "in")), 
       width = 20, height = 15, units = 'in', dpi = 600)

#绘制KEGG图
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = GOgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
kegg <- kegg[order(kegg$Count,decreasing = TRUE),]

#KEGG图的另外一种展现形式,展现的是前15条信号通路
#将富集到的信号通路按照富集到的基因的个数进行排序
# 定义一个函数来计算下调和上调基因的数量
colnames(sum)
table1 <- sum[,c(2,363)]#做AS分析的时候注意修改这里的数字
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  # 统计下调和上调基因数量
  down_genes <- sum(subset(gene_table,gene_table$Genename %in% gene_list)[,2] == "Down")
  up_genes <- sum(subset(gene_table,gene_table$Genename %in% gene_list)[,2] == "Up")
  # down_genes <- sum(gene_list %in% gene_table$Genename[gene_table$change == "Down"])
  # up_genes <- sum(gene_list %in% gene_table$Genename[gene_table$change == "Up"])
  return(c(down_genes, up_genes))
}

# 对表格2中的每一行计算down和up的数量
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))

# 将结果加到表格2中
table2$down <- -abs(result[, 1])
table2$up <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "AS_KEGG.csv",row.names = F)
#对数据进行转换画图
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:15,c(3,4,12,13)] #做AS分析的时候注意这里的数字
colnames(KEGGTerm)
#对数据格式进行转换
mydata<-melt(KEGGTerm,id.vars=c("ID","Description"),variable.name="Change",value.name="Number")
mydata$type_order=factor(rev(as.integer(rownames(mydata))),labels=rev(mydata$Description))

library(stringr)
p1 <- ggplot(mydata,aes(type_order,Number)) + 
  geom_bar(aes(fill=factor(Change)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#008080","firebrick3"))+
  geom_text(aes(label=abs(Number)),color="black", size=4,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("AS KEGG")+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Gene Number")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill = guide_legend(title = 'Change'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 15,face = "plain"),
        axis.title.y = element_text(color = "black",size = 15,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 15,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "white"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)



#对所有有差异的intron建立ROC曲线分析
#取出有差异的intron在不同的样本中的表达情况
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F,row.names = 1)
#根据变化倍数和padj对数据进行筛选
b <- subset(a,a$change !="Stable")
#对表格进行排序
b <- b[order(b$Padjust),]
head(b)
colnames(b)
VoomRemoveBatch_DEG <- b[1:10,2:355]
VoomRemoveBatch_DEG <- as.data.frame(t(VoomRemoveBatch_DEG))
##对b数据框添加一列样本的分组信息
VoomRemoveBatch_DEG$Group <- rownames(VoomRemoveBatch_DEG)
head(VoomRemoveBatch_DEG$Group)

b <- VoomRemoveBatch_DEG
#只取出来第一个数字
b$Group <- substr(b$Group,1,1)
head(b$Group)
#对这些354个introns进行ROC分析
n <- dim(b)[1]
y <- b$Group


all <- b

###开始随机抽样
set.seed(1)
require(caret)
folds <- createFolds(y,k=6)

library(ROCR)
library(magrittr) # pipe operator
library(plyr)
auc_value<-as.numeric()
#将所有的fitvalue取出来
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- all[-folds[[i]],] # 剩下的数据作为训练集
  model <- glm(as.numeric(fold_train$Group)~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- predict(model,type='response',newdata=fold_test)
  pred <- prediction(predictions = fold_predict, labels = fold_test$Group)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
auc_value
# 1] 0.9019608 0.7965116 0.8343023 0.8401163 0.7194767 0.9142442
mean(auc_value)
# [1] 0.8344353

#去掉数据框rocData的第一行
rocData <- rocData[-1,]
#更改数据框rocData的行名
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
#计算均值
# write.csv(rocData,file = "ASroc.csv",quote = F,row.names = F)
#将数据框中出现的NA用1进行补全
# rocData[is.na(rocData)] <- 1
#添加Xmean和Ymean的值
colnames(rocData)
rocData$Xmean <- rowMeans(rocData[,1:6])
rocData$Ymean <- rowMeans(rocData[,7:12])
rocData$Gene <- c(rep("AS", 60))
rocData <- rocData[,13:15]
#计算95%的置信区间
#计算不出来
# library(magrittr) # pipe operator
# rocit_emp <- rocit(score = rocData$Ymean, class = all$BackgroundInformation[test], method = "emp")
# summary(rocit_emp)
# #计算CI值
# ciAUC(rocit_emp)

# Xmean指的是所有X的平均值
# Ymean指的是所有Y的平均值
#构造样本ROC的终点数据
#构造中间虚线的数据
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("Control", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
#合并数据
ROC <- rbind(rocData,inn)
#画图
library(ggrepel)
library(ggplot2)
a<- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene)) + 
  geom_line(aes(linetype=Gene,color=Gene),size=1)+
  scale_linetype_manual(values = c(7,2))+
  scale_color_manual(values=c("#00872D","black"))+
  xlab("False Positive Rate") +
  ylab("True Positive Rate")+
  theme_bw()+
  ggtitle("AS ROC AUC=0.8344") +
  #scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=22),
        legend.title = element_text(size=8),
        legend.position = c(0.8,0.25),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain"),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"))
a
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)



#找出top 10 变化的introns对应的基因是谁#
#取出有差异的intron在不同的样本中的表达情况
rm(list=ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F,row.names = 1)
#根据变化倍数和padj对数据进行筛选
b <- subset(a,a$change !="Stable")
#对表格进行排序
b <- b[order(b$Padjust),]
head(b)
b <- b[1:10,]
b <- b[order(rownames(b)),]
colnames(b)
colnames(b)[2:355]
#加载
VoomRemoveBatch_DEG <- b[,1:355]
rownames(VoomRemoveBatch_DEG) <- VoomRemoveBatch_DEG$Genename
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG[,-1]
colnames(VoomRemoveBatch_DEG)
VoomRemoveBatch_DEG <- as.data.frame(t(VoomRemoveBatch_DEG))
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG[,order(colnames(VoomRemoveBatch_DEG))]
#将基因组上的位置信息替换为相应的基因的名字
##对b数据框添加一列样本的分组信息
VoomRemoveBatch_DEG$Group <- substr(rownames(VoomRemoveBatch_DEG),1,1)
head(VoomRemoveBatch_DEG$Group)
#使用循环，绘制这几个基因在不同的分组中的基因的表达量

library(tidyverse)
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG%>%mutate(Group=case_when(VoomRemoveBatch_DEG$Group==0~"Healthy",
                                                                    VoomRemoveBatch_DEG$Group==1~"AD",
                                                                    VoomRemoveBatch_DEG$Group==2~"MCI"))
VoomRemoveBatch_DEG$Group <- as.factor(VoomRemoveBatch_DEG$Group)
str(VoomRemoveBatch_DEG)
# gene_group <- gene_group[,order(colnames(gene_group))]
#绘制小提琴图
colnames(VoomRemoveBatch_DEG)
# 同时对多个文件进行输出
x = names(VoomRemoveBatch_DEG)[11]
y = names(VoomRemoveBatch_DEG)[-11]
table(VoomRemoveBatch_DEG$Group)
# AD Healthy 
# 257      97 
library(ggsignif)
library(ggpubr)
library(glue)
plot_list = map2(x, y, 
                 ~ VoomRemoveBatch_DEG %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_violin(aes_string(fill = .x),trim = FALSE)+
                   geom_signif(comparisons = list(c("AD", "Healthy")), 
                               tip_length = 0.02,
                               margin_top = 0.15,size = 0.8,textsize = 8,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#008080","#CC6600"))+
                   stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                geom="pointrange", color = "red")+
                   theme_bw() + 
                   xlab("Group") + 
                   ylab("Gene Expression") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",linewidth = 0,linetype = "solid",colour = "black")))
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'AS_Top10_Genes.pdf')





#######在可变多聚腺苷酸化分析方面，表型数据的变化会影响多因素方差分析的结果#####
#APA的分析，如果直接使用PDUI则无法去除批次效应，因此尝试重新使用count值进行计算
#使用count进行去除批次的时候将long和short分别进行批次效应的去除
#如果使用相加的结果的话没办法进行后续MOFA2的分析，因为其要求必须是经过标准化之后的数据
#如果使用相加的结果的话本质上和基因表达分析没有区别
rm(list = ls())
getwd()
#首先加载表型的数据#加载序号和样本名字的对应信息
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
DaparsSample <- read.csv("DaparsSample.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
colnames(BackgroundInformation)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
colnames(DaparsSample)[2:3] <- c("GSM_number","sample")
pheno <- merge(DaparsSample,BackgroundInformation,by="GSM_number")
pheno$sample <- gsub("_PDUI$", "", pheno$sample)
pheno <- subset(pheno,pheno$Brain_Region=="3")
pheno$Batch <- as.factor(pheno$Batch)
pheno$Braak_Stage <- as.factor(pheno$Braak_Stage)
# pheno$Sex <- as.factor(pheno$Sex)
pheno$Bank_Location <- as.factor(pheno$Bank_Location)
pheno$Severity <- as.factor(pheno$Severity)
colnames(pheno)







setwd("E:/AD_Patient/Dapars2/Dapars0.9")
#又重新进行分析，这时没有设置任何的条件限制，只是为了得到长短的测序深度
#在Dapars2的默认参数中，当长的和短的总共的测序深度小于30时，就会被当作NA处理
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
head(Healthy_AD$Pass_Filter)
table(Healthy_AD$Pass_Filter)
test <- Healthy_AD[,c(1,1067:1072)]
head(Healthy_AD[,c(1,1067:1072)])
a <- Healthy_AD[,-c(1:4,1067:1072)]
#直接绘制PDUI的PCA和多因素方差分析
a <- a[,seq(0,ncol(a),3)]
Healthy_AD_PDUI <- cbind(Healthy_AD[,1:4],a)
#将加载进来的表格按照|进行划分
library(tidyr)
colnames(Healthy_AD_PDUI)
long<-separate(Healthy_AD_PDUI,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
#保存下基因名字的转化
ID_Transfer <- long[,1:2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(ID_Transfer,file = "APA_ID_Transfer.csv",quote = F,row.names = F)
#
Healthy_AD_PDUI <- long
Healthy_AD_PDUI <- Healthy_AD_PDUI[,-c(2:7)]
rownames(Healthy_AD_PDUI) <- Healthy_AD_PDUI$ENSEMBLTRANS
#去掉第1列的数据
Healthy_AD_PDUI <- Healthy_AD_PDUI[,-1]
#直接提取出PDUI的数据
head(colnames(Healthy_AD_PDUI))
colnames(Healthy_AD_PDUI) <- gsub("_PDUI", "", colnames(Healthy_AD_PDUI))


#比较丑的PCA图
Healthy_AD_PDUI <- as.data.frame(t(Healthy_AD_PDUI))
#PCA绘制样本数据
library(gmodels)
pca.info <- fast.prcomp(Healthy_AD_PDUI)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
#将样本的名字与样本background的信息合并
b <- merge(pheno,pca.data,by="sample")
ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
b$Batch <- as.factor(b$Batch)
ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# #ggplot2画图
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#CC6600"))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-4,2))+
  ggtitle("APA Raw PDUI")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",linewidth = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",47,"%"),
       y=paste0("PCA2 ",19,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
#对legend的顺序进行设定
# guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_PDUI_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)

#分别对long和short的数据进行处理
# 将short和long的数据分别提取出来，分别去除批次效应
setwd("E:/AD_Patient/Dapars2/Dapars0.9")
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
tail(colnames(Healthy_AD))
name <- Healthy_AD[,c(1:4)]
num <- Healthy_AD[,-c(1:4,1067:1072)]
#对name进行处理
colnames(name)
long<-separate(name,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- cbind(long[,1],num)
colnames(Healthy_AD)[1] <- "ENSEMBLTRANS"
rownames(Healthy_AD) <- Healthy_AD$ENSEMBLTRANS
Healthy_AD <- Healthy_AD[,-1]
colnames(Healthy_AD)
#提取出long expression的数据
#四舍五入只提取整数部分
Healthy_AD_long <- Healthy_AD[,seq(1,ncol(Healthy_AD),3)]
head(colnames(Healthy_AD_long))
colnames(Healthy_AD_long) <- gsub("_long_exp", "", colnames(Healthy_AD_long))

#提取出short expression的数据
#四舍五入只提取整数部分
Healthy_AD_short <- Healthy_AD[,seq(2,ncol(Healthy_AD),3)]
head(colnames(Healthy_AD_short))
colnames(Healthy_AD_short) <- gsub("_short_exp", "", colnames(Healthy_AD_short))


#分别对两组数据进行去批次
#long的数据
#比较丑的PCA图
Healthy_AD_long <- as.data.frame(t(Healthy_AD_long))
#PCA绘制样本数据
library(gmodels)
pca.info <- fast.prcomp(Healthy_AD_long)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
#将样本的名字与样本background的信息合并
b <- merge(pheno,pca.data,by="sample")
ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
b$Batch <- as.factor(b$Batch)
ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# #ggplot2画图
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#CC6600"))+
  scale_x_continuous(limits = c(-25000, 45000))+
  scale_y_continuous(limits = c(-30000,45000))+
  ggtitle("APA Raw Long")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",linewidth = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",46,"%"),
       y=paste0("PCA2 ",18,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
#对legend的顺序进行设定
# guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Long_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)



#做多因素方差分析
Healthy_AD_long <- Healthy_AD_long[order(rownames(Healthy_AD_long)),]
pheno <- pheno[order(pheno$sample),]
identical(rownames(Healthy_AD_long),pheno$sample)
library(vegan)
set.seed(111)
Healthy_AD_long_avo <- adonis2(as.matrix(Healthy_AD_long) ~ Severity+Braak_Stage+Age+Bank_Location+Batch, 
                               data = pheno,
                               permutations=999, by="margin",method = "euclidean")

Healthy_AD_long_avo
#欧式距离
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = as.matrix(Healthy_AD_long) ~ Severity + Braak_Stage + Age + Bank_Location + Batch, data = pheno, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2       F Pr(>F)    
# Severity        1  7.6750e+08 0.00754  3.7917  0.015 *  
#   Braak_Stage     3  1.6060e+09 0.01577  2.6447  0.008 ** 
#   Age             1  8.9653e+08 0.00881  4.4292  0.011 *  
#   Bank_Location   0 -3.9580e+03 0.00000    -Inf           
# Batch           2  1.8044e+10 0.17723 44.5724  0.001 ***
#   Residual      345  6.9834e+10 0.68591                   
# Total         353  1.0181e+11 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#gower距离
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = as.matrix(Healthy_AD_long) ~ Severity + Braak_Stage + Age + Bank_Location + Batch, data = pheno, permutations = 999, method = "gower", by = "margin")
# Df SumOfSqs      R2       F Pr(>F)    
# Severity        1   0.1120 0.01735 11.2532  0.001 ***
#   Braak_Stage     3   0.0896 0.01387  2.9985  0.004 ** 
#   Age             1   0.0154 0.00238  1.5458  0.152    
# Bank_Location   0   0.0000 0.00000    -Inf           
# Batch           2   1.6787 0.25995 84.3100  0.001 ***
#   Residual      345   3.4347 0.53186                   
# Total         353   6.4579 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b <- as.data.frame(Healthy_AD_long_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("APA Raw Long") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Long.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)





#构建matrix，去除批次效应
library(dplyr)
library(limma)
#将属于分类变量的表型数据转化为因子
str(pheno)
#使用removeBatchEffect去除批次
dim(Healthy_AD_long)#数据的格式是行是基因名列是样本名
colnames(pheno)
#数据Bank_Location和batch之间存在共线性
design <- model.matrix( ~ Severity + Age , data=pheno)
head(design)
treatment.design <- design[,1:2]#之所以保存两列是因为有1列是常数项，另外1列是severity
head(treatment.design)
batch.design <- design[,-(1:2)]
head(batch.design)
dim(Healthy_AD_long)
head(design)
#对样本进行排序
Healthy_AD_long <- as.data.frame(t(Healthy_AD_long))
Healthy_AD_long <- Healthy_AD_long[,order(colnames(Healthy_AD_long))]
pheno <- pheno[order(pheno$sample),]
identical(colnames(Healthy_AD_long),pheno$sample)
#去批次
Healthy_AD_long_Correct <- removeBatchEffect(as.matrix(Healthy_AD_long),
                                             batch = pheno$Batch,
                                             group = pheno$Severity,
                                             covariates= batch.design,#batch.design,treatment.design
                                             design = treatment.design)

#数据在去除了批次之后的范围应该是在0-1之间
min(Healthy_AD_long_Correct)
max(Healthy_AD_long_Correct)
dim(Healthy_AD_long_Correct)
#保存去完批次之后的数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(Healthy_AD_long_Correct,file = "APA_Long_RemoveBatch_Count.csv",quote = F)

#比较丑的PCA图
Healthy_AD_long_Remove <- as.data.frame(t(Healthy_AD_long_Correct))
#PCA绘制样本数据
library(gmodels)
pca.info <- fast.prcomp(Healthy_AD_long_Remove)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
#将样本的名字与样本background的信息合并
b <- merge(pheno,pca.data,by="sample")
ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
b$Batch <- as.factor(b$Batch)
ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# #ggplot2画图
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#CC6600"))+
  scale_x_continuous(limits = c(-40000, 21000))+
  scale_y_continuous(limits = c(-20000,60000))+
  ggtitle("APA RemoveBatch Long")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",linewidth = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",42,"%"),
       y=paste0("PCA2 ",24,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
#对legend的顺序进行设定
# guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Long_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)


#检验批次是否去除
Healthy_AD_long_Remove <- t(Healthy_AD_long_Correct)
#多因素方差分析要保证行是样本列是基因
#输入的数据是矩阵
library(vegan)
colnames(pheno)
dim(Healthy_AD_long_Remove)
#设置不同的seed
# 打开文本文件，将结果输出进去，尝试使用不同的seed看结果会不会有什么不同
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sink("Healthy_AD_long_Remove_avo_results.txt")

# 运行循环，并将每次的结果输出
for (i in 1:1000) {
  set.seed(i)
  
  # 执行 adonis2
  Healthy_AD_long_Remove_avo <- adonis2(Healthy_AD_long_Remove ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
                                        data = pheno,
                                        permutations=999, by="margin",method = "euclidean")
  
  # 打印当前的 seed 值
  cat("Results for seed:", i, "\n")
  
  # 打印 adonis2 的结果
  print(Healthy_AD_long_Remove_avo)
  
  # 添加分隔符，便于区分不同 seed 值的结果
  cat("\n-----------------------------\n")
}

# 关闭 sink，停止写入文件
sink()
#最后得到seed是43的时候Severity比braak_stage的结果更显著
set.seed(43)
Healthy_AD_long_Remove_avo <- adonis2(Healthy_AD_long_Remove ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
                                      data = pheno,
                                      permutations=999, by="margin",method = "euclidean")
Healthy_AD_long_Remove_avo
#最终选择的seed
# Results for seed: 43 
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Healthy_AD_long_Remove ~ Braak_Stage + Age + Bank_Location + Batch + Severity, data = pheno, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2      F Pr(>F)   
# Braak_Stage     3  1605986664 0.02177 2.6447  0.011 * 
#   Age             1    58861602 0.00080 0.2908  0.914   
# Bank_Location   0       -3958 0.00000   -Inf          
# Batch           2    71089810 0.00096 0.1756  0.998   
# Severity        1   767495672 0.01040 3.7917  0.004 **
#   Residual      345 69833574543 0.94651                 
# Total         353 73780423883 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



b <- as.data.frame(Healthy_AD_long_Remove_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("APA RemoveBatch Long") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Long.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)





#去除short的数据
#long的数据
Healthy_AD_short <- as.data.frame(t(Healthy_AD_short))
Healthy_AD_short <- Healthy_AD_short[order(rownames(Healthy_AD_short)),]
pheno <- pheno[order(pheno$sample),]
identical(rownames(Healthy_AD_short),pheno$sample)
library(vegan)
set.seed(111)
Healthy_AD_short_avo <- adonis2(as.matrix(Healthy_AD_short) ~ Severity+Braak_Stage+Age+Bank_Location+Batch, 
                                data = pheno,
                                permutations=999, by="margin",method = "euclidean")

Healthy_AD_short_avo
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = as.matrix(Healthy_AD_short) ~ Severity + Braak_Stage + Age + Sex + Bank_Location + Batch, data = pheno, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2       F Pr(>F)    
# Severity        1   186195306 0.00729  4.4018  0.017 *  
#   Braak_Stage     3   198299195 0.00777  1.5627  0.129    
# Age             1    96583091 0.00378  2.2833  0.061 .  
# Sex             1    57332348 0.00225  1.3554  0.233    
# Bank_Location   0        -106 0.00000    -Inf           
# Batch           2  6261037532 0.24524 74.0084  0.001 ***
#   Residual      344 14551030102 0.56996                   
# Total         353 25530007589 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

b <- as.data.frame(Healthy_AD_short_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("APA Raw Short") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Short.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)




Healthy_AD_short <- as.data.frame(Healthy_AD_short)
#PCA绘制样本数据
library(gmodels)
pca.info <- fast.prcomp(Healthy_AD_short)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
#将样本的名字与样本background的信息合并
b <- merge(pheno,pca.data,by="sample")
ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
b$Batch <- as.factor(b$Batch)
ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# #ggplot2画图
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#CC6600"))+
  scale_x_continuous(limits = c(-33000, 15000))+
  scale_y_continuous(limits = c(-10000,10000))+
  ggtitle("APA Raw Short")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",linewidth = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",57,"%"),
       y=paste0("PCA2 ",11,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
#对legend的顺序进行设定
# guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Short_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)




#构建matrix，去除批次效应
library(dplyr)
library(limma)
#将属于分类变量的表型数据转化为因子
str(pheno)
#使用removeBatchEffect去除批次
dim(Healthy_AD_short)#数据的格式是行是基因名列是样本名
colnames(pheno)
#数据Bank_Location和batch之间存在共线性
design <- model.matrix( ~ Severity + Age , data=pheno)
head(design)
treatment.design <- design[,1:2]#之所以保存两列是因为有1列是常数项，另外1列是severity
head(treatment.design)
batch.design <- design[,-(1:2)]
head(batch.design)
dim(Healthy_AD_long)
head(design)
#对样本进行排序
Healthy_AD_short <- as.data.frame(t(Healthy_AD_short))
Healthy_AD_short <- Healthy_AD_short[,order(colnames(Healthy_AD_short))]
pheno <- pheno[order(pheno$sample),]
identical(colnames(Healthy_AD_short),pheno$sample)
#去批次
Healthy_AD_short_Correct <- removeBatchEffect(as.matrix(Healthy_AD_short),
                                              batch = pheno$Batch,
                                              group = pheno$Severity,
                                              covariates= batch.design,#batch.design,treatment.design
                                              design = treatment.design)

#数据在去除了批次之后的范围应该是在0-1之间
min(Healthy_AD_short_Correct)
max(Healthy_AD_short_Correct)
dim(Healthy_AD_short_Correct)
#保存去完批次之后的数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(Healthy_AD_short_Correct,file = "APA_Short_RemoveBatch_Count.csv",quote = F)
#比较丑的PCA图
Healthy_AD_short_Remove <- as.data.frame(t(Healthy_AD_short_Correct))
#PCA绘制样本数据
library(gmodels)
pca.info <- fast.prcomp(Healthy_AD_short_Remove)
head(pca.info)
head(summary(pca.info)) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# 利用ggscatter对PC进行绘图
#绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
library(ggpubr)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
#将样本的名字与样本background的信息合并
b <- merge(pheno,pca.data,by="sample")
ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
b$Batch <- as.factor(b$Batch)
ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# #ggplot2画图
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+
  geom_point(size=3)+
  scale_color_manual(values = c("#911F27", "#3366CC","#009933","#CC6600"))+
  scale_x_continuous(limits = c(-30000, 10000))+
  scale_y_continuous(limits = c(-10000,10000))+
  ggtitle("APA RemoveBatch Short")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 10,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",linewidth = 0.2,linetype = "solid",colour = "black"),
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",44,"%"),
       y=paste0("PCA2 ",14,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
#对legend的顺序进行设定
# guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Short_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)






#检验批次是否去除
Healthy_AD_short_Remove <- t(Healthy_AD_short_Correct)
#多因素方差分析要保证行是样本列是基因
#输入的数据是矩阵
library(vegan)
colnames(pheno)
dim(Healthy_AD_short_Remove)
set.seed(111)
Healthy_AD_short_Remove_avo <- adonis2(Healthy_AD_short_Remove ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
                                       data = pheno,
                                       permutations=999, by="margin",method = "euclidean")
Healthy_AD_short_Remove_avo

# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Healthy_AD_short_Remove ~ Braak_Stage + Age + Bank_Location + Batch + Severity, data = pheno, permutations = 999, method = "euclidean", by = "margin")
# Df    SumOfSqs      R2      F Pr(>F)   
# Braak_Stage     3   197715886 0.01304 1.5558  0.099 . 
# Age             1     6334655 0.00042 0.1495  0.999   
# Bank_Location   0         -71 0.00000   -Inf          
# Batch           2     9601600 0.00063 0.1133  1.000   
# Severity        1   183582227 0.01210 4.3337  0.005 **
#   Residual      345 14614619669 0.96361                 
# Total         353 15166470557 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

b <- as.data.frame(Healthy_AD_short_Remove_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
#自己手动添加显著性标记的记号
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
#将category中的下划线去掉
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   # b$category=="CERAD"~"CERAD",
                                   # b$category=="Sex"~"Gender",
                                   b$category=="Age"~"Age"))


#将b中的数据进行排序
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
  xlab("category") +ylab("R2")+
  theme_bw()+#theme_classic()+
  ggtitle("APA RemoveBatch Short") +
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"))

p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Short.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)







#统计一下负值的个数
table(Healthy_AD_long_Correct<0)
# FALSE  TRUE 
# 59516  9514
table(Healthy_AD_long_Correct<0)/(dim(Healthy_AD_long_Correct)[1]*dim(Healthy_AD_long_Correct)[2])
# FALSE      TRUE 
# 0.8621759 0.1378241
#将这些负值变为NA
Healthy_AD_long_Correct[Healthy_AD_long_Correct<0]=NA
sum(is.na(Healthy_AD_long_Correct))
head(Healthy_AD_long_Correct)[1:5,1:5]

table(Healthy_AD_short_Correct<0)
# FALSE  TRUE 
# 60292  8738 
table(Healthy_AD_short_Correct<0)/(dim(Healthy_AD_short_Correct)[1]*dim(Healthy_AD_short_Correct)[2])
# FALSE      TRUE 
# 0.8734174 0.1265826 
#将负值变为NA
Healthy_AD_short_Correct[Healthy_AD_short_Correct<0]=NA
sum(is.na(Healthy_AD_short_Correct))
head(Healthy_AD_short_Correct)[1:5,1:5]

#计算PDUI,检查横纵坐标是否一致，只有一致的时候才能进行计算
Healthy_AD_long_Correct <- Healthy_AD_long_Correct[order(rownames(Healthy_AD_long_Correct)),order(colnames(Healthy_AD_long_Correct))]
Healthy_AD_short_Correct <- Healthy_AD_short_Correct[order(rownames(Healthy_AD_short_Correct)),order(colnames(Healthy_AD_short_Correct))]
identical(colnames(Healthy_AD_long_Correct),colnames(Healthy_AD_short_Correct))
identical(rownames(Healthy_AD_long_Correct),rownames(Healthy_AD_short_Correct))
if (!all(dim(Healthy_AD_long_Correct) == dim(Healthy_AD_short_Correct))) {
  stop("两个矩阵的维度必须相同")}

Healthy_AD_remove_PDUI <- Healthy_AD_long_Correct / (Healthy_AD_long_Correct + Healthy_AD_short_Correct)
sum(is.na(Healthy_AD_remove_PDUI))
# [1] 15988
head(Healthy_AD_remove_PDUI)[1:5,1:5]
#使用均值填充缺失值，否则不能进行PCA和多因素方差分析
# 假设你的数据已经读取为data_frame
# colnames(Healthy_AD_remove_PDUI)
# group <- as.data.frame(Healthy_AD_remove_PDUI)
# A_group_filled <- apply(group[,1:97], 1, function(x) ifelse(is.na(x), statip::mfv(x, na_rm = TRUE), x))
# B_group_filled <- apply(group[,98:354], 1, function(x) ifelse(is.na(x), statip::mfv(x, na_rm = TRUE), x))
# data_filled <- rbind(A_group_filled, B_group_filled)
# #
# #输入的数据是矩阵
# library(vegan)
# colnames(pheno)
# dim(data_filled)
# set.seed(111)
# identical(rownames(data_filled),pheno$sample)
# Missing values are not allowed on the left-hand-side.
#所以说不能对数据进行多因素方差分析
# Healthy_AD_remove_PDUI_avo <- adonis2(as.matrix(data_filled) ~ Braak_Stage+Age+Bank_Location+Batch+Severity,
#                                        data = pheno,
#                                        permutations=999, by="margin",method = "euclidean")
# Healthy_AD_remove_PDUI_avo
# 
# #PCA
# library(gmodels)
# pca.info <- fast.prcomp(data_filled)
# head(pca.info)
# head(summary(pca.info)) #能看到每个PC具体的解释比例
# head(pca.info$rotation) #特征向量，回归系数
# head(pca.info$sdev) #特征值的开方
# head(pca.info$x)
# pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
# # 利用ggscatter对PC进行绘图
# #绘图，不同的参数会有不同的结果。具体的参数以及含义自行百度。以下是两个例子。
# library(ggpubr)
# ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=FALSE, ellipse.type="confidence")
# b <- merge(pheno,pca.data,by="sample")
# ggscatter(b,x="PC1", y="PC2", color="Group", ellipse=FALSE, ellipse.type="confidence")
# b$Batch <- as.factor(b$Batch)
# ggscatter(b,x="PC1", y="PC2", color="Batch", ellipse=FALSE, ellipse.type="confidence")
# 
#
#找某个基因发生的具体的APA的变化
RemoveData <- Healthy_AD_remove_PDUI
dim(RemoveData)
head(RemoveData)
#进行两组数据之间的比较
table(substr(colnames(RemoveData),1,1))
#A指的是healthy，B指的是AD
#计算出每一组的平均值
#首先将所有的组按照AD和Healthy进行排序
RemoveData <- RemoveData[,order(colnames(RemoveData))]
colnames(RemoveData)
RemoveData <- as.data.frame(RemoveData)
RemoveData$groupAmean <- rowMeans(RemoveData[,1:97], na.rm = TRUE)
RemoveData$groupBmean <- rowMeans(RemoveData[,98:354], na.rm = TRUE)
head(RemoveData$groupAmean)
head(RemoveData$groupBmean)
#计算出差值
RemoveData$Mean_diff <- RemoveData$groupBmean-RemoveData$groupAmean
head(RemoveData$Mean_diff)

#使用计算出来的均值计算显著性
#由于样本的个数太少，所以没办法使用
#t.test(RemoveData$groupAmeanPDUI[1],RemoveData$groupBmeanPDUI[1])

#然后使用每一个样本的PDUI值计算p值的显著性
#例子
t.test(RemoveData[1,1:97],RemoveData[1,98:354],na.rm = TRUE)

# Welch Two Sample t-test
# 
# data:  RemoveData[1, 1:97] and RemoveData[1, 98:354]
# t = 0.61511, df = 131.74, p-value = 0.5395
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.03555322  0.06764243
# sample estimates:
#   mean of x mean of y 
# 0.7561093 0.7400647 

#因为数据中有负值，因此不能使用t test因为数据的正负会影响到均值的计算
#所以使用所有的数据计算p值
#不能使用的是配对的t检验(因为样本的个数不一致)
RemoveData$Pvalue <- apply(RemoveData,1,function(x) wilcox.test(na.omit(x[1:97]),na.omit(x[98:354]))$p.value)
head(RemoveData$Pvalue)
RemoveData <- RemoveData[complete.cases(RemoveData$Pvalue), ]
min(RemoveData$Pvalue)

# 计算结果异常，因为有NaN
# 出现NaN的原因为，NM_001270399.2这个转录本的PDUI值均为1,在原始数据和removeBatch之后的数据均为1
min(RemoveData$Pvalue)
#对计算出的p值进行校正
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
head(RemoveData$Padjust)
min(RemoveData$Padjust)
max(RemoveData$Padjust)
table(RemoveData$Padjust < 0.05)
# FALSE  TRUE 
# 130    64
#A指的是healthy，B指的是AD
RemoveData$change <- ifelse(RemoveData$Padjust < 0.05 & abs(RemoveData$Mean_diff) >= 0.1, 
                            ifelse(RemoveData$Mean_diff > 0.1 ,'Long','Short'),'Stable')

table(RemoveData$change)
# Long  Short Stable 
# 1     10    183   
test <- RemoveData[,355:360]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(RemoveData,file = "APA_RemoveBatch_DEG.csv",quote = F)

# 尝试使用的一些去除批次效应的方法
# head(Healthy_AD)
# Healthy_AD <- as.data.frame(t(Healthy_AD))
# df_transformed <- log2(Healthy_AD + 1)  
# #将数据转换为log转化的形式
# # Define a small epsilon value to adjust 0s and 1s
# epsilon <- 1e-6
# 
# # Step 1: Adjust the data
# df <- Healthy_AD
# df_adjusted <- Healthy_AD
# df_adjusted[df == 0] <- epsilon
# df_adjusted[df == 1] <- 1 - epsilon
# 
# # Step 2: Apply the logit transformation
# logit_transform <- function(x) {
#   return(log(x / (1 - x)))
# }
# #将所有的数据转换在负无穷和正无穷之间，增大样本之间的差异
# #apply函数1 indicates rows, 2 indicates columns
# #apply设置为1或者2都可以，只是运算的方向会不一样，最终结果都一样
# df_logit <- apply(df_adjusted, 2, logit_transform)
# dim(df_logit)
#将数据进行sigmod转换
# 使用sigmoid变换
# # sigmoid_transform <- function(x) {
# #   return(1 / (1 + exp(-x)))
# # }
# #应用
# corrected_count_sig <- sigmoid_transform(corrected_count)
#另一种对数据进行转换的方法
# Define a small epsilon value to adjust 0s and 1s
# epsilon <- 1e-6
# # Step 1: Adjust the data
# df <- Healthy_AD
# df_adjusted <- Healthy_AD
# df_adjusted[df == 0] <- epsilon
# df_adjusted[df == 1] <- 1 - epsilon
# 
# # Step 2: Apply the logit transformation
# logit_transform <- function(x) {
#   return(log(x / (1 - x)))
# }
# 
# df_logit <- apply(df_adjusted, 2, logit_transform)
# dim(df_logit)

# # Example: Using Harmony for batch correction
# library(harmony)
# harmony_results <- HarmonyMatrix(df_logit, meta_data = batch, vars_use = "Batch")
# 
# # Inverse logit transform if necessary
# df_harmony <- apply(harmony_results, 2, inv_logit_transform)
# dim(df_harmony)
# pca_final <- prcomp(df_harmony)
# pca_df_final <- as.data.frame(pca_final$x[, 1:2])
# pca_df_final$Batch <- batch
# 
# ggplot(pca_df_final, aes(x = PC1, y = PC2, color = Batch)) +
#   geom_point() +
#   ggtitle("PCA After Logit and Batch Correction")
# Step 3: Remove batch effects on logit-transformed data
# Assuming your batch vector is named 'batch'
# Define the inverse logit function
# inv_logit_transform <- function(x) {
#   return(exp(x) / (1 + exp(x)))
# }
# # Apply the inverse logit transformation
# df_final <- apply(df_combat, 2, inv_logit_transform)
# dim(df_final)

#使用这156个基因进行富集分析
#将ENSEMBLTRANS转换为entreID
sum <- test
#找出每个test中的转录本对应的名字
sum$ENSEMBLTRANS <- rownames(sum)
com <- merge(ID_Transfer,sum,by="ENSEMBLTRANS")
sum <- com
#绘制火山图
library(dplyr)
up <- subset(sum, sum$change == 'Long')
down <- subset(sum, sum$change == 'Short')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:1, ]
down <- down[order(down$Padjust), ][1:5, ]
a <- rbind(up, down)
colnames(a)
library(ggplot2)
library(ggrepel)
p <- ggplot(
  # 数据、映射、颜色
  sum, aes(x = Mean_diff, y = -log10(Padjust),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080","firebrick3","gray")) +
  geom_text_repel(data = a, aes(x = Mean_diff, y = -log10(Padjust), label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
  # 坐标轴
  labs(x="Mean_diff",y="-log10(padj)") +
  ggtitle("APA volcano")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.position ="right",
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.background = element_rect(fill = "transparent",linewidth = 0.4,linetype = "solid"))+
  scale_x_continuous(breaks = seq(-3, 4, 2))+
  scale_y_continuous(breaks = seq(0, 15, 5))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_volcano.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)








#进行富集分析
a <- sum
a <- a[a$change!="Stable",]
gene <- unique(a$GeneName)
#67个
# gene <- nchar(unique(sum$V4))#nchar是统计这个向量里每一个单个元素的字符串的个数
#对具有可变剪切的基因进行富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(gene, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
# 'select()' returned 1:1 mapping between keys and columns
colnames(df1)
GOgene <- as.character(df1$ENTREZID)#得出所有的ENTREZID
GOgene <- na.omit(GOgene)#得出所有的ENTREZID，就可以进行GO和kegg的分析
BPplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='BP',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
BPplot_genelist<-setReadable(BPplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
BPplot_genelist <- as.data.frame(BPplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(BPplot_genelist,file = "APA_BPplot_genelist.csv",row.names = F,quote = F)

#CC
CCplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='CC',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
CCplot_genelist<-setReadable(CCplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CCplot_genelist <- as.data.frame(CCplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(CCplot_genelist,file = "APA_CCplot_genelist.csv",row.names = F,quote = F)

#MF
MFplot <- enrichGO(GOgene, 
                   OrgDb = org.Hs.eg.db, 
                   ont='MF',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID')
MFplot_genelist<-setReadable(MFplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
MFplot_genelist <- as.data.frame(MFplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(MFplot_genelist,file = "APA_MFplot_genelist.csv",row.names = F,quote = F)

###绘制合在一块的GO图，选取的原则是根据padj从小到大来筛选
go_enrich_df <- data.frame(
  ID=c(BPplot_genelist$ID[1:15], CCplot_genelist$ID[1:14], MFplot_genelist$ID[1:15]),
  Description=c(BPplot_genelist$Description[1:15],CCplot_genelist$Description[1:14],MFplot_genelist$Description[1:15]),
  GeneNumber=c(BPplot_genelist$Count[1:15], CCplot_genelist$Count[1:14], MFplot_genelist$Count[1:15]),
  type=factor(c(rep("biological process", 15), 
                rep("cellular component", 14),
                rep("molecular function", 15)), 
              levels=c("biological process", "cellular component","molecular function" )))
go_enrich_df <- na.omit(go_enrich_df)
go_enrich_df <- arrange(go_enrich_df,type,GeneNumber)
##开始绘制GO柱状图
###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
# write.csv(go_enrich_df,file = "go_enrich_dfAS.csv",quote = F,row.names = F)
COLS <- c("#0D6CA6","#099963", "#911F27")#设定颜色
library(ggpubr)
library(stringr)
p <- ggdotchart(go_enrich_df, x = "type_order", y = "GeneNumber",
                color = "type",                                # 按照cyl填充颜色
                palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                sorting = "descending",                      
                add = "segments",                             # 添加棒子
                add.params = list(color = "type", size = 1.3),#改变棒子参数
                # rotate = TRUE,                                # 方向转为垂直
                group = "type",                                
                dot.size = "GeneNumber",                                 # 改变点的大小
                #label = round(go_enrich_df$GeneNumber),                       # 添加label
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               # 设置label参数
                ggtheme = theme_pubr(),                        # 改变主题
                xlab="GO Term",
                ylab="GeneNumber",
                title = "APA GO")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40,exdent = 0),"\n")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 14,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 30),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain",angle = 70,vjust = 1, hjust = 1 ),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"),
        legend.position.inside = c(0.9,0.68),
        legend.title = element_text(color = "black",size = 15,face = "plain"),
        legend.text = element_text(color = "black",size = 20,face = "plain"))
#legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_GO.pdf", egg::set_panel_size(p, width=unit(16, "in"), height=unit(5, "in")), 
       width = 25, height = 20, units = 'in', dpi = 600)

#绘制KEGG图
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = GOgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
kegg <- kegg[order(kegg$Count,decreasing = TRUE),]

#KEGG图的另外一种展现形式,展现的是前15条信号通路
#将富集到的信号通路按照富集到的基因的个数进行排序
# 定义一个函数来计算下调和上调基因的数量
table1 <- sum[,c(1,2,8)]
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  # 统计下调和上调基因数量
  # down_genes <- sum(gene_list %in% gene_table$Genename[gene_table$change == "Down"])
  # up_genes <- sum(gene_list %in% gene_table$Genename[gene_table$change == "Up"])
  down_genes <- sum(subset(gene_table,gene_table$GeneName %in% gene_list)[,3] == "Short")
  up_genes <- sum(subset(gene_table,gene_table$GeneName %in% gene_list)[,3] == "Long")
  return(c(down_genes, up_genes))
}

# 对表格2中的每一行计算down和up的数量
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))

# 将结果加到表格2中
table2$Short <- -abs(result[, 1])
table2$Long <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "APA_KEGG.csv",row.names = F)
#对数据进行转换画图
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:8,c(3,4,12,13)]
colnames(KEGGTerm)
#对数据格式进行转换
mydata<-melt(KEGGTerm,id.vars=c("ID","Description"),variable.name="Change",value.name="Number")
mydata$type_order=factor(rev(as.integer(rownames(mydata))),labels=rev(mydata$Description))

library(stringr)
p1 <- ggplot(mydata,aes(type_order,Number)) + 
  geom_bar(aes(fill=factor(Change)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#008080","firebrick3"))+
  geom_text(aes(label=abs(Number)),color="black", size=4,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("APA KEGG")+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Gene Number")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill = guide_legend(title = 'Change'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 15,face = "plain"),
        axis.title.y = element_text(color = "black",size = 15,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 15,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        #legend.background = element_rect(fill = "transparent",size = 0.2,linetype = "solid",colour = "white"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(2.5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)


#绘制ROC曲线的时候要使用的是去除完批次之后的数据
#对所有有差异的intron建立ROC曲线分析
#取出有差异的APA在不同的样本中的表达情况
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#加载去除完批次之后的数据
a <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
#根据变化倍数和padj对数据进行筛选
b <- subset(a,a$change !="Stable")
#对表格进行排序
b <- b[order(b$Padjust),]
head(b)
b <- b[1:10,]
head(b[,355:360])
b <- b[,1:354]
b <- b[order(rownames(b)),]
b <- as.data.frame(t(b))
b <- b[,order(colnames(b))]

#加载ID的信息，找到这些转录本对应的基因
APA_ID_Transfer <- read.csv("APA_ID_Transfer.csv",header = T)
APA_ID_Transfer <- subset(APA_ID_Transfer,APA_ID_Transfer$ENSEMBLTRANS %in% colnames(b))
APA_ID_Transfer <- APA_ID_Transfer[order(APA_ID_Transfer$ENSEMBLTRANS),]
identical(colnames(b),APA_ID_Transfer$ENSEMBLTRANS)
APA_ID_Transfer$GeneName
# [1] "RPL15"     "ISCU"      "EIF4A2"    "PEBP1"     "RPL15"     "ELOB"      "GABARAPL2"
# [8] "SARAF"     "ISCU"      "ARF3"
#只取出来第一个数字,在制作分组信息的时候不能使用字母，要使用数字
b$Group <- substr(rownames(b),1,1)
head(b$Group)
library(dplyr)
b <- b%>%mutate(Group=case_when(b$Group=="A"~"0",
                                b$Group=="B"~"1"))
b$Group <- as.numeric(b$Group)
#进行ROC分析的基因的个数
n <- dim(b)[1]
y <- b$Group
all <- b
str(b)
###开始随机抽样
set.seed(11)
require(caret)
folds <- createFolds(y,k=6)
library(ROCR)
library(magrittr) # pipe operator
library(plyr)
auc_value<-as.numeric()
#将所有的fitvalue取出来
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] #取folds[[i]]作为测试集
  fold_test <- fold_test[rowSums(fold_test[1:10],na.rm =TRUE)!=0,]
  fold_train <- all[-folds[[i]],] # 剩下的数据作为训练集
  fold_train <- fold_train[rowSums(fold_train[1:10],na.rm =TRUE)!=0,]
  model <- glm(as.numeric(fold_train$Group)~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- predict(model,type='response',newdata=fold_test)
  samples_to_remove <- names(fold_predict[is.na(fold_predict)])
  fold_predict <- fold_predict[!is.na(fold_predict)]
  fold_test <- fold_test[!rownames(fold_test) %in% samples_to_remove, ]
  pred <- prediction(predictions = fold_predict, labels = fold_test$Group)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
auc_value
# [1] 0.6695157 0.8104396 0.8298611 0.7692308 0.8901515 0.7051282
mean(auc_value)
# [1] 0.7790545

#去掉数据框rocData的第一行
rocData <- rocData[-1,]
#更改数据框rocData的行名
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData <- rocData[1:44,]
#计算均值
# write.csv(rocData,file = "ASroc.csv",quote = F,row.names = F)
#将数据框中出现的NA用1进行补全
# rocData[is.na(rocData)] <- 1
#添加Xmean和Ymean的值
colnames(rocData)
rocData$Xmean <- rowMeans(rocData[,1:6],na.rm = T)
rocData$Ymean <- rowMeans(rocData[,7:12],na.rm = T)
rocData$Gene <- c(rep("APA", 44))
rocData <- rocData[,13:15]
#计算95%的置信区间
#计算不出来
# library(magrittr) # pipe operator
# rocit_emp <- rocit(score = rocData$Ymean, class = all$BackgroundInformation[test], method = "emp")
# summary(rocit_emp)
# #计算CI值
# ciAUC(rocit_emp)

# Xmean指的是所有X的平均值
# Ymean指的是所有Y的平均值
#构造样本ROC的终点数据
#构造中间虚线的数据
a <- seq(0,1,1/43)
inn <- data.frame(a,a,rep("Control", 44))
colnames(inn) <- c("Xmean","Ymean","Gene")
#合并数据
ROC <- rbind(rocData,inn)
#画图
library(ggrepel)
library(ggplot2)
a<- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene)) + 
  geom_line(aes(linetype=Gene,color=Gene),linewidth=1)+
  scale_linetype_manual(values = c(7,2))+
  scale_color_manual(values=c("#00872D","black"))+
  xlab("False Positive Rate") +
  ylab("True Positive Rate")+
  theme_bw()+
  ggtitle("APA ROC AUC=0.779") +
  #scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=22),
        legend.title = element_text(size=8),
        legend.position = c(0.8,0.25),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain"),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"))
a
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)



#找出top 10 变化的转录本对应的基因是谁#
#取出有差异的转录本在不同的样本中的表达情况
rm(list=ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#将转录本的名字转换为基因的名字之后会出现重复的名字，因此决定以转录本为单个单位进行
a <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
#根据变化倍数和padj对数据进行筛选
b <- subset(a,a$change !="Stable")
#对表格进行排序
b <- b[order(b$Padjust),]
colnames(b)
head(b[355:360])
b <- b[1:10,1:354]
colnames(b)
b <- b[order(rownames(b)),]
b <- as.data.frame(t(b))
b <- b[,order(colnames(b))]
# #加载ID信息
# APA_ID_Transfer <- read.csv("APA_ID_Transfer.csv",header = T)
# APA_ID_Transfer <- subset(APA_ID_Transfer,APA_ID_Transfer$ENSEMBLTRANS %in% colnames(b))
# APA_ID_Transfer <- APA_ID_Transfer[order(APA_ID_Transfer$ENSEMBLTRANS),]
# #转换
# colnames(b) <- ifelse(APA_ID_Transfer$ENSEMBLTRANS == colnames(b),APA_ID_Transfer$GeneName,"NA")
##对b数据框添加一列样本的分组信息
b$Group <- substr(rownames(b),1,1)
head(b$Group)
#使用循环，绘制这几个基因在不同的分组中的基因的表达量
library(tidyverse)
#A指的是healthy，B指的是AD
b <- b%>%mutate(Group=case_when(b$Group=="A"~"Healthy",
                                b$Group=="B"~"AD"))
b$Group <- as.factor(b$Group)
str(b)
# gene_group <- gene_group[,order(colnames(gene_group))]
#绘制小提琴图
colnames(b)
# 同时对多个文件进行输出
x = names(b)[11]
y = names(b)[-11]
table(b$Group)
# AD Healthy 
# 257      97 
library(ggsignif)
library(ggpubr)
library(glue)
plot_list = map2(x, y, 
                 ~ b %>% 
                   ggplot(aes_string(x = .x, y = .y)) +
                   geom_violin(aes_string(fill = .x),trim = FALSE)+
                   geom_signif(comparisons = list(c("AD", "Healthy")), 
                               tip_length = 0.02,
                               margin_top = 0.15,size = 0.8,textsize = 8,
                               map_signif_level=TRUE)+
                   scale_fill_manual(values = c("#008080","#CC6600"))+
                   stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                                geom="pointrange", color = "red")+
                   theme_bw() + 
                   xlab("Group") + 
                   ylab("Transcript Expression") + 
                   labs(title = glue('{.y}'))+
                   scale_size(range=c(2, 8))+
                   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
                         plot.title = element_text(color = "black",size = 20,face = "plain",hjust = 0.5),
                         plot.margin = margin(t = 10,r = 20,b = 10,l = 10),
                         axis.title.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.title.y = element_text(color = "black",size = 16,face = "plain"),
                         panel.grid = element_blank(),
                         axis.text.x = element_text(color = "black",size = 16,face = "plain"),
                         axis.text.y = element_text(color = "black",size = 16,face = "plain"),
                         legend.position = "none",
                         legend.title = element_text(color = "black",size = 16,face = "plain"),
                         legend.text = element_text(color = "black",size = 16,face = "plain"),
                         legend.background = element_rect(fill = "transparent",linewidth = 0,linetype = "solid",colour = "black")))
ggexport(plotlist = plot_list, width = 4, height = 4,bg = "white",
         filename = 'APA_Top10_Genes.pdf')




#绘制差异APA的热图
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#读取数据
sum <- read.csv("APA_RemoveBatch_DEG.csv",header = T,row.names = 1)
table(sum$change)
APA_Heatmap <- sum[sum$change!="Stable",]
colnames(APA_Heatmap)
APA_Heatmap <- APA_Heatmap[,1:354]
APA_Heatmap <- APA_Heatmap[,order(colnames(APA_Heatmap))]
colnames(APA_Heatmap)
#对于缺失值的处理
A_group_filled <- apply(APA_Heatmap[,1:97], 1, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
B_group_filled <- apply(APA_Heatmap[,98:354], 1, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
data_filled <- rbind(A_group_filled, B_group_filled)
data_scale <- t(data_filled)
up <- subset(sum, sum$change == 'Long')
down <- subset(sum, sum$change == 'Short')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:1, ]
down <- down[order(down$Padjust), ][1:10, ]
a <- rbind(up, down)
head(a)


#标记最显著的基因的名字
APA_ID_Transfer <- read.csv("APA_ID_Transfer.csv",header = T)
APA_ID_Transfer <- APA_ID_Transfer[APA_ID_Transfer$ENSEMBLTRANS %in% rownames(a),]

sum <- as.data.frame(data_filled)
#对样本进行排序，使得AD和Healthy的样本分为上下两个不同的组
sum <- arrange(sum,rownames(sum))
data <- sum
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){
  sd(x)/mean(x)*100
})
data_df <- data[order(data$cv, decreasing = T),1:354]
dim(data_df)
# [1]  64 354
a <- apply(data_df,1,scale)
# 进行scale，但是不设定最大最小值，这个时候的热图一片绿
data_scale <- as.data.frame(t(apply(data_df,1,scale))) ##Z-score标准化
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.001
#将scale之后的数据设置最大值最小值
data_scale <- as.matrix(data_scale)
#设置最大值与最小值
table((data_scale)>1)
table((data_scale)<(-1))
data_scale[data_scale>=1]=1
data_scale[data_scale<=(-1)]=0
table(is.na(data_scale))
# library(pheatmap)
# pheatmap(data_scale)
library(ComplexHeatmap)
library(circlize)

#读取待展示的基因名称，并添加到热图中
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "APA_Heatmap4.pdf",width =4,height = 3)
# png("p.png",res = 300)

p <- Heatmap(data_scale,name = "Expression", 
             na_col = "grey",
             cluster_rows = TRUE,
             clustering_distance_rows = "pearson",
             clustering_method_rows = "complete",
             col = colorRampPalette(c("#008080", "white", "firebrick3"))(50),
             cluster_columns = TRUE,
             clustering_distance_columns = "pearson",
             clustering_method_columns = "complete",
             row_dend_side = "left",
             show_row_names = TRUE,
             show_column_names = TRUE,
             use_raster = T,
             heatmap_width = unit(0.001, "npc"),
             width = NULL,
             heatmap_height = unit(0.001, "npc"),
             height = NULL,
             show_column_dend = FALSE,
             column_dend_height = unit(5, "mm"),
             show_row_dend = FALSE,
             column_labels = colnames(data_scale),
             column_names_side = "bottom",
             column_names_centered = TRUE,
             column_title = "",
             column_names_rot = 0,
             column_title_gp = gpar(fontsize = 0.001),
             column_names_gp = gpar(fontsize = 2),
             heatmap_legend_param = list( #设置标签
               color_bar = 'continuous',
               legend_direction = 'vertical',
               legend_width = unit(2, 'cm'),
               legend_height = unit(2, 'cm'),
               title_position = 'topcenter',
               title_gp = gpar(fontsize = 1, fontface = 'plain'),
               labels_gp = gpar(fontsize = 1, fontface = 'plain'))
)+ 
  rowAnnotation(link = anno_mark(at = which(rownames(data_scale) %in% APA_ID_Transfer$ENSEMBLTRANS), 
                                 labels = APA_ID_Transfer$GeneName[match(rownames(data_scale)[rownames(data_scale) %in% APA_ID_Transfer$ENSEMBLTRANS], APA_ID_Transfer$ENSEMBLTRANS)], labels_gp = gpar(fontsize = 7)))

print(p)
dev.off()

#
# #####查找3个层次的差异基因是否有重叠的地方#####
# #不找，不做，因为3个层面代表的东西不一样
# #加载数据
# #GE
# rm(list=ls())
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# AD_Healthyres <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
# AD_Healthyres <- AD_Healthyres[AD_Healthyres$change!="Stable",]
# GE_Genes <- rownames(AD_Healthyres)
# #AS
# AS_name <- read.table("ASdiff_n369.txt",header = T,check.names = F)
# AS_name <- AS_name[AS_name$change!="Stable",]
# AS_name$change
# AS_Genes <- AS_name$Genename
# AS_Genes <- AS_Genes[!grepl(",", AS_Genes)]
# #APA
# APA_name <- read.csv("APA_RemoveBatch_DEG.csv",header = T,row.names = 1,check.names = F)
# APA_name <- APA_name[APA_name$change!="Stable",]
# APA_name <- rownames(APA_name)
# #进行ID的转换
# APA_ID <- read.csv("APA_ID_Transfer.csv",header = T,check.names = F)
# APA_Genes <- subset(APA_ID,APA_ID$ENSEMBLTRANS %in% APA_name)
# APA_Genes <- APA_Genes$GeneName
# 
# #画一个Venn图
# a <- list(GE_Genes = as.list.data.frame(unique(GE_Genes)), 
#           AS_Genes = as.list.data.frame(unique(AS_Genes)),
#           APA_Genes = as.list.data.frame(unique(APA_Genes)))
# 
# library(ggvenn)
# p <- ggvenn(a, show_elements = FALSE, label_sep = "\n", digits = 1,
#             fill_color = c("#E24D37", "#3F8CAD", "#009966","#CC3366"),
#             fill_alpha = 0.5,
#             stroke_color = "black",
#             stroke_alpha = 1,
#             stroke_size = 0.5,
#             stroke_linetype = "solid",
#             set_name_color = "black",
#             set_name_size = 4,
#             text_color = "black",
#             text_size = 4)# show elements in line
# p
# ggsave("Total_Venn.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
#        width = 8, height = 7, units = 'in', dpi = 600)
# 
# #分别找出每一小Part的基因
# #同时找出3组中共同的基因
# common <- list(unique(GE_Genes), unique(AS_Genes),
#                unique(APA_Genes))
# overlap_genes <- Reduce(intersect, common)
# #找出两组两组之间的基因
# GE_APA <- intersect(unique(GE_Genes),unique(APA_Genes))
# GE_APA
# # GFAP
# GE_AS <- intersect(unique(GE_Genes),unique(AS_Genes))
# GE_AS
# # [1] "PPEF1"     "SLC6A12"   "VGF"       "RGS4"      "TAC1"      "RPH3A"     "PART1"    
# # [8] "PRMT8"     "PCSK1"     "STAT4"     "LINC00460" "CRYM"      "TUBB2A"    "NEUROD6"  
# # [15] "MIR7-3HG"  "MLIP"      "SVOP"      "RGS7"      "FRMPD2B"   "KCNE4"     "SNAP25"   
# # [22] "P2RX7"     "OLFM3"     "C2orf80"   "NEFL"      "PTPN3"     "CREG2"     "LDB2"     
# # [29] "MAP4K4"    "VSNL1"     "GABRG2"    "SH2D5"     "ENC1"      "CBLN4"     "CORT"     
# # [36] "RBFOX1"    "HIPK2"     "SYT4"      "TESPA1"    "ANLN"      "CLEC2L"    "FKBP5"    
# # [43] "NEFM"      "GFRA2"     "SLC30A3"   "PNMA3"     "CARTPT"    "LINC00507" "LAMP5"    
# # [50] "PRKCG"     "PCDHGC5"   "GABRD"     "PCSK6"     "SPP1"
# AS_APA <- intersect(unique(AS_Genes),unique(APA_Genes))
# AS_APA
# # [1] "DCTN1"
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# AD_Healthyres <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
# AD_Healthyres$GeneName <- rownames(AD_Healthyres)
# 
# Diff_APA_in_GE <- subset(AD_Healthyres,AD_Healthyres$GeneName %in% APA_Genes)
# # Diff_APA_in_GE <- Diff_APA_in_GE[Diff_APA_in_GE$padj < 0.1,]
# APA_info <- APA_name[,c(355:360)]
# APA_info$ENSEMBLTRANS <- rownames(APA_info)
# APA_info_1 <- merge(APA_ID,APA_info,by="ENSEMBLTRANS")
# APA_info_1 <- APA_info_1[,-8]
# sum <- merge(APA_info_1,Diff_APA_in_GE,by="GeneName",all = TRUE)
# sum$change <- ifelse(sum$log2FoldChange<0,"Down","Up")
# sum$PASchange <- ifelse(sum$Mean_diff<0,"Short","Long")
# write.csv(sum,file = "Prediction_GE_AS.csv",quote = F,row.names = F)
# test <- sum[,c(1,2,14,15)]
# 
# 
# 


#####进行3个层次的数据联合分析#####
#MOFA2分析的目的找到影响因子贡献度最重要的那些基因、代谢物、微生物。
#输入的数据必须是经过标准化的数据

#加载不同分析层次的数据
# 加载有差异变化的前5000个基因，按照padj的先后顺序进行提取
#首先加载基因的表达信息，使用的是FPKM的值
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sumGE <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
#选择出padj最小的5000个基因
sumGE <- sumGE[order(sumGE$padj),]
sumGE <- sumGE[1:5000,]
#MOFA2输入的必须是经过标准化的数据，不然
# 加载这些基因的FPKM值
FPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
# 将FPKM的基因名从ensembl转换为symbol
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(FPKM))
rownames(FPKM) <- a
####进行ID转换
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)

df1 <- df1[!duplicated(df1$SYMBOL),]
FPKM$ENSEMBL <- rownames(FPKM)
df <- merge(df1,FPKM,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
FPKM <- df[,-c(1:2)]
#提取出这些基因的FPKM
gene <- rownames(sumGE)
sumGene <- subset(FPKM,rownames(FPKM) %in% gene)
sumGene <- as.data.frame(t(sumGene))
sumGene$GSM_number <- rownames(sumGene)
#将列名按照分组转换为AD和healthy
setwd("C:/Users/yujie/Desktop/datacollect")
bg <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
colnames(bg)
sum <- merge(bg,sumGene,by="GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
sumGene <- as.data.frame(t(sum))


###加载可变剪切的信息，使用的是normalized的数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AS <- read.csv("AS_VoomRemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
#查看AS后面几列的数据，因为正常情况下样本的个数为354个
colnames(AS)
head(AS[,350:360])
#讲数据按照padj进行排序，然后选出padj最小的前5000个
AS <- AS[order(AS$Padjust),]
AS <- AS[1:5000,]
head(AS[,350:360])
tail(AS[,350:360])
AS <- AS[,-c(355:360)]#将生成的结果作为MOFA的输入
colnames(AS)
#第3个层次的数据是APA的数据，筛选的标准是按照padj进行排序，选择padj>0.05的
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
APA <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
#查看APA后面几列的数据，因为正常情况下样本的个数为354个
head(APA[,350:360])
APA <- APA[order(APA$Pvalue),]
APA <- subset(APA,APA$Pvalue<0.05)
head(APA[,350:360])
tail(APA[,350:360])
APA <- APA[,-c(355:360)]#将生成的结果作为MOFA的输入
#将APA的A和B编号转换为samble编号
DaparsSample <- read.csv("DaparsSample.csv",header = T)
colnames(DaparsSample)
colnames(DaparsSample) <- c("Group","GSM_number","sample")
DaparsSample$sample <- gsub("_PDUI","",DaparsSample$sample)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
colnames(BackgroundInformationAPOE)
bg <- merge(DaparsSample,BackgroundInformationAPOE,by="GSM_number")

#对编号进行转换
APA <- as.data.frame(t(APA))
colnames(APA)
APA$sample <- rownames(APA)
APA2 <- merge(bg,APA,by="sample") 
rownames(APA2) <- APA2$Sample
APA <- APA2[,-c(1:14)]
APA <- as.data.frame(t(APA))


#在使用MOFA2之前构建数据
# please make sure that all views contain the same samples in the same order (see documentation)
# 保证加载的每一个view中列都是一样的样本顺序
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
sumGene <- sumGene[,order(colnames(sumGene))]
AS <- AS[,order(colnames(AS))]
APA <- APA[,order(colnames(APA))]
identical(colnames(sumGene),colnames(AS))
identical(colnames(AS),colnames(APA))
#将上面的数据创建一个列表,其中保证每一个表格中行表示的是基因，列表示的是样本
#背景信息当作额外的一个table，
# ,view3=as.matrix(AS)
test <- list(Transcriptome=as.matrix(sumGene),AS=as.matrix(AS),APA=as.matrix(APA))
# head(test)
lapply(test,dim)
# $Transcriptome
# [1] 5000  354
# 
# $AS
# [1] 5000  354
# 
# $APA
# [1]  82 354


#创建MOFA对象,必要条件：1：列中样本的顺序必须一致,2：加载的每一个view必须为matrix
library(data.table)
library(MOFA2)
MOFAobject <- create_mofa(test)

#对样本进行分组
# library(stringr)
# groups <- substr(colnames(diff_Count),1,1)
# table(groups)
# groups
# MOFAobject <- create_mofa(test, groups=groups)
# MOFAobject
#画图展示出构建的结果
plot_data_overview(MOFAobject)

#对数据进行处理的选项
data_opts <- get_default_data_options(MOFAobject)
data_opts
#对模型进行处理的参数设置
model_opts <- get_default_model_options(MOFAobject)
# model_opts$num_factors <- 10
model_opts



# 参数的使用：
# ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
#对训练的参数进行设置
train_opts <- get_default_training_options(MOFAobject)
# For exploration, the fast mode is good enough.For exploration, the fast mode is good enough.
train_opts$convergence_mode <- "fast"

train_opts$seed <- 20
train_opts
#当所有的参数设置好之后就开始对模型进行训练
#Build and train the MOFA object
# Prepare the MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts)

# 开始训练模型Train the MOFA model
getwd()

outfile = file.path("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2","MOFA2Results_2")
MOFAobject.trained <- run_mofa(MOFAobject, outfile,use_basilisk = TRUE)
#上面的数据迭代了86次
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)

#对训练好的模型进行分析
#训练出的模型有问题，factor1所占的比例很大，可能是数据预处理的不合格
# Factor(s) 1 are strongly correlated with the total number of expressed features 
# for at least one of your omics. Such factors appear when there are 
# differences in the total 'levels' between your samples, 
# *sometimes* because of poor normalisation in the preprocessing steps.
# model <- MOFAobject.trained
slotNames(model)
names(model@data) 
# dim(model@data$Transcriptome)#查看healthy的信息
head(model@data$Transcriptome)

#model中的Z和W分别表示的是MOFA2奇异值分解的两个矩阵
#Z表示的是sample和factor组成的矩阵
#正常来说sample和factor应该是一个矩阵，
#W表示的是feature和factor组成的矩阵
names(model@expectations)
#
dim(model@expectations$Z$group1)
dim(model@expectations$W$Transcriptome)
dim(model@expectations$W$AS)
####之后将这些feature转换为基因的名字
head(model@expectations[["W"]][["Transcriptome"]])
head(model@expectations[["W"]][["AS"]])



p1 <- plot_data_overview(model)
p1
library(ggplot2)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p1,filename = "MOFA2_data_overview.pdf",width = 4,height = 2,bg = "white")

# Add metadata to the model
Nsamples = sum(model@dimensions$N)
Nsamples
#Add sample metadata to the model
#首先查看下model的metadata,之后再对metadata进行更改
model@samples_metadata
#然后对原始的metadata进行添加
#对上面的sum文件进行筛除
setwd("C:/Users/yujie/Desktop/datacollect")
bg <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
colnames(bg)
bg <- subset(bg,bg$Brain_Region=="3")
bg <- bg[,-c(1:5,11)]#删除这些多余的列

colnames(bg)
# [1] "CERAD"         "Braak_Stage"   "Age"          
# [5] "Sex"           "Bank_Location" "Severity"      "Sample"       
# [9] "Phenotypedit"
head(bg)
#将bg中携带APOE的那一列转换为数字
library(tidyverse)
bg <- bg%>%mutate(Phenotypedit=case_when(bg$Phenotypedit=="noE4"~"0",
                                         bg$Phenotypedit=="E4carrier"~"1",
                                         bg$Phenotypedit=="E4/4"~"2"))


#对bg的样本名进行排序，使得和metadata中样本的顺序一致
bg <- bg[order(bg$Sample),]
#截取到字符的位置可以任意，因为要截取到最后一个字符
a <- substr(colnames(as.data.frame(model@data$Transcriptome)),8,15)
#将bg的表格按照排序好的样本进行排序
identical(a,bg$Sample)
colnames(bg)
#重新赋值metadata
colnames(bg)[5] <- "sample"
samples_metadata(model) <- bg
head(model@samples_metadata, n=3)



# Variance decomposition(对变量进行分解)
#每一个view的变量
# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) # group 1
# Transcriptome            AS           APA 
# 43.75445      55.78877      44.88901
#展示每一个view中每一个factor的可解释的变量
# Variance explained for every factor in per view and group
factor <- as.data.frame(model@cache$variance_explained$r2_per_factor[[1]])
head(factor)
# Transcriptome        AS      APA
# Factor1    28.7149251 27.390450 9.599113
# Factor2    14.1089439 10.631073 5.746794
# Factor3     3.3159733  6.354648 4.518014
# Factor4     0.7967234  5.306852 6.104589
# Factor5     0.8598328  2.423757 6.224537
# Factor6     1.0305762  2.792138 3.147596
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(factor,file = "MOFA2_Factor_Distribution.csv",quote = F)
head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1
#画图
p2 <- plot_variance_explained(model, x="view", y="factor",factors = 1:5)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p2,filename = "MOFA2_Factor_Explained.pdf",width = 5,height = 2,bg = "white")


# p3 <- plot_variance_explained(model, x="group", y="factor", plot_total = T)[[1]]
# ggsave(p3,filename = "variance_explained2.pdf",width = 6,height = 3,bg = "white")
# ggsave(p3,filename = "variance_explained2.png",width = 6,height = 3,bg = "white")
p4 <- plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p4,filename = "MOFA2_Total_Explained.pdf",width = 4,height = 2,bg = "white")
#绘制相关系数图
#绘制出的相关系数图有什么意义
plot_factor_cor(model)
# Visualisation of single factors
#将每个因子中的100个sample的可解释的值绘制出来，
#每个因子都有100个sample
colnames(bg)
#图片中的点表示的是不同的sample
#绘制的是sample和factor这个矩阵的数据
plot_factor(model, 
            factor = 1,
            color_by = "Phenotypedit",#age
            shape_by = "group",
            scale = T,
            show_missing = F)+
  scale_color_manual(values = c("#008080", "gray", "firebrick3"))
dev.off()
SampleFactor <- as.data.frame(model@expectations[["Z"]][["group1"]])
SampleFactor$sample <- rownames(SampleFactor)
SampleFactor$Group <- rep("group1",354)
Sample <- merge(bg,SampleFactor,by="sample")
colnames(Sample)
Sample <- na.omit(Sample)
str(Sample)
Sample$Braak_Stage <- as.factor(Sample$Braak_Stage)
# Sample$Sex <- as.factor(Sample$Sex)
Sample$Age <- as.numeric(Sample$Age)
Sample$Severity <- as.factor(Sample$Severity)
str(Sample)
colnames(Sample)
#绘图
#绘制的是sample和factor之间的关系，加入表型信息
library(ggsignif)
library(ggpubr)
library(glue)#在图形上添加p值
p1 <- ggplot(Sample,aes(x=Severity,y=Factor1))+
  # geom_violin(width =0.8,fill='grey90',color='grey90')+
  geom_signif(comparisons = list(c("0", "1")), 
              tip_length = 0.02,
              margin_top = 0.01,size = 0.8,textsize = 5,
              map_signif_level=TRUE)+
  geom_violin()+
  geom_jitter(aes(color=Severity),width = 0.2,size=3)+
  scale_color_manual(name = 'Severity',
                     values = c("#263E88","firebrick3"),
                     #labels = c('1','2','3')
  )+ #直接修改图例文本内容
  theme_bw()+
  ggtitle("MOFA2 SampleFactor1")+
  # geom_hline(yintercept = 0,lty=4,col="#666666",lwd=0.5) +
  labs(x="",y=expression('Factor1'))+ #使用expression函数添加特殊文本
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=10),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_SampleFactor1.pdf", egg::set_panel_size(p1, width=unit(3.5, "in"), height=unit(4, "in")), 
       width = 6, height = 5, units = 'in', dpi = 600)



# p1 <- ggplot(Sample,aes(x=Factor1,y=Factor2))+
#   #geom_violin(width =0.8,fill='grey90',color='grey90')+
#   geom_jitter(aes(color=Severity),width = 0.2,size=2)+
#   scale_color_manual(name = 'Severity',
#                      values = c("firebrick3","#263E88"),
#                      labels = c('1','2','3')
#   )+ #直接修改图例文本内容
#   theme_bw()+
#   geom_hline(yintercept = 0,lty=4,col="#666666",lwd=0.5) +
#   #theme_classic(base_size = 15)+
#   labs(x="Factor1",y=expression('Factor2'))+ #使用expression函数添加特殊文本
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         axis.text.x = element_blank(),#去掉横坐标上的文本
#         axis.ticks.x = element_blank(),
#         legend.key = element_blank(),
#         legend.text=element_text(size=8),
#         legend.title=element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 10,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 10,face = "plain"),
#         panel.grid= element_blank(),
#         #axis.text.x = element_text(color = "black",size = 8,face = "plain",angle = 40,vjust = 0.7),
#         axis.text.y = element_text(color = "black",size = 8,face = "plain"))
# p1
# ggsave(p1,filename = "Factor1Factor2.pdf",width = 4.5,height = 4,bg = "white")
# 
# 


#分析MOFA2分析的结果的ROC分析
#这里的ROC指的是使用不同的因子可以更好的区分AD和Healthy吗？
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
SampleFactor <- as.data.frame(model@expectations[["Z"]][["group1"]])
SampleFactor$Sample <- rownames(SampleFactor)
SampleFactor <- SampleFactor[,c(1,16)]
#将样本的名字转换为AD还是healthy
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
all <- merge(SampleFactor,BackgroundInformation,by="Sample")
all <- all[,c(2,8)]#0表示的是healthy，1表示的是AD
#将汉字转换为数字0和1 
table(all$Severity)
library(tidyverse)
all$Severity <- as.numeric(all$Severity)
str(all)
n <- dim(all)[1]
y <- all$Severity

###开始随机抽样
set.seed(1)
require(caret)
folds <- createFolds(y,k=6)

library(ROCR)
library(magrittr) # pipe operator
library(plyr)
auc_value<-as.numeric()
#将所有的fitvalue取出来
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- all[-folds[[i]],] # 剩下的数据作为训练集
  model <- glm(fold_train$Severity~.,data=fold_train,
               family = binomial(link = "logit"))
  #注意mofa2中也有一个predict函数，会和原来的predict函数冲突
  fold_predict <- stats::predict(model,type='response',newdata=fold_test)#样本点并没有减少
  pred <- prediction(fold_predict,fold_test$Severity)#经过这一步样本点减少了
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
auc_value
mean(auc_value)
# #10次计算
#
#去掉数据框rocData的第一行
rocData <- rocData[-1,]
#更改数据框rocData的行名

rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 


rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
#计算均值
rocData$Xmean <- apply(rocData[,1:6],1,mean)
rocData$Ymean <- apply(rocData[,7:12],1,mean)
rocData <- rocData[,13:14]
rocData$Gene <- c(rep("all", 60))

#计算95%的置信区间
#计算不出来
# library(magrittr) # pipe operator
# rocit_emp <- rocit(score = rocData$Ymean, class = all$BackgroundInformation[test], method = "emp")
# summary(rocit_emp)
# #计算CI值
# ciAUC(rocit_emp)


#构造样本ROC的终点数据
#构造中间虚线的数据
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("Control", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
#合并数据
ROC <- rbind(rocData,inn)
#画图
library(ggrepel)
library(ggplot2)
a<- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene)) + 
  geom_line(aes(linetype=Gene,color=Gene),size=1)+
  scale_linetype_manual(values = c(7,2))+
  scale_color_manual(values=c("#00872D","black"))+
  xlab("False Positive Rate") +
  ylab("True Positive Rate")+
  theme_bw()+
  ggtitle("MOFA2 Factor ROC AUC=0.733") +
  #scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=22),
        legend.title = element_text(size=8),
        legend.position = c(0.8,0.25),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain"),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"))
a
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_FactorROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)




# plot_factor(MOFAobject, 
#             factors = 1, 
#             color_by = "condition")

# plot_factor_cor(model)



# Adding more options
# colnames(bg)
# table(bg$Sex)
# p <- plot_factor(model, 
#                  factors = 1,
#                  color_by = "Severity",
#                  dot_size = 3,        # change dot size
#                  dodge = T,           # dodge points with different colors
#                  legend = F,          # remove legend
#                  add_violin = T,      # add violin plots,
#                  violin_alpha = 0.25)  # transparency of violin plot
# library(ggplot2)
# p <- p + 
#   scale_color_manual(values=c("1"="red", "0"="blue")) +
#   scale_fill_manual(values=c("1"="red", "0"="blue"))
# 
# print(p)

# Visualisation of combinations of factors
# Scatter plots
# install.packages("GGally")
#和上面的函数是不一样的函数
#上面的是plot_factor
# plot_factors(model, 
#              factors = 1,
#              color_by = "Severity")

# Visualisation of feature weights
# The weights provide a score for how strong each feature relates to each factor.
# model@dimensions
# 
# plot_weights(model,
#              view = "Transcriptome",#view2
#              factor = 1,
#              nfeatures = 10,     # Number of features to highlight
#              scale = T,          # Scale weights from -1 to 1
#              abs = F ,
#              legend = TRUE,
#              dot_size = 1,
#              text_size = 3)            # Take the absolute value?
# 
# plot_weights(model,
# view = "AS",#view2
# factor = 1,
# nfeatures = 10,     # Number of features to highlight
# scale = T,          # Scale weights from -1 to 1
# abs = F ,
# dot_size = 1,
# text_size = 3)            # Take the absolute value?

# plot_weights(model,
#              view = "APA",#view2
#              factor = 1,
#              nfeatures = 10,     # Number of features to highlight
#              scale = T,          # Scale weights from -1 to 1
#              abs = F ,
#              dot_size = 1,
#              text_size = 3)            # Take the absolute value?


#使用gglot2重新绘制图形
# edit(plot_weights)
#使用geom_point绘制
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)
GEweights <- as.data.frame(model@expectations[["W"]][["Transcriptome"]])
#绘制出factor1在不同的view中占比重比较大的因素
GEweights$GeneName <- rownames(GEweights)
#将数据按照weight的绝对值的大小进行排序
GEweights_Pos <- GEweights[GEweights$Factor1>=0,]
GEweights_Pos$Weight <- (GEweights_Pos$Factor1-min(abs(GEweights$Factor1)))/(max(abs(GEweights$Factor1))-min(abs(GEweights$Factor1)))
GEweights_Pos <- GEweights_Pos[order(abs(GEweights_Pos$Factor1),decreasing = T),]
# GEweights_Pos_10 <- GEweights_Pos[1:5,]
GEweights_Neg <- GEweights[GEweights$Factor1<0,]
GEweights_Neg$Weight <- (GEweights_Neg$Factor1-min(abs(GEweights$Factor1)))/(max(abs(GEweights$Factor1))-min(abs(GEweights$Factor1)))
GEweights_Neg <- GEweights_Neg[order(abs(GEweights_Neg$Factor1),decreasing =T),]
# GEweights_Neg_10 <- GEweights_Neg[1:5,]
a <- rbind(GEweights_Pos,GEweights_Neg)
KEGGa <- a

#按照权重的绝对值大小进行排序
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$GeneName
topGE <- a
#将上下两个数据合并进行绘图
colnames(GEweights_Pos)
colnames(GEweights_Neg)
GEweights <- rbind(GEweights_Pos,GEweights_Neg)
GEweights <- GEweights[order(GEweights$Factor1,decreasing = F),]
GEweights$Rank <- seq(1,5000,1)
GEweights_Top <- subset(GEweights,GEweights$GeneName %in% a)
#将想要关注的点设置为黑色，不想要关注的点设置为灰色
GEweights$Color <- ifelse(GEweights$GeneName %in% a, "#196F6E","grey")
#绘图
library(ggrepel)
p <- ggplot(
  GEweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#196F6E", "gray")) +
  geom_text_repel(data = GEweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 2,box.padding = unit(0.2, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=2000)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  # 坐标轴
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 Transcriptome Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),#去掉横坐标上的文本
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)




GEweights_Top$abs <- abs(GEweights_Top$Weight)
colnames(GEweights_Top)
library(ggpubr)
p <- ggdotchart(GEweights_Top, x = "GeneName", y = "Weight",
                color = "#196F6E",                                # 按照cyl填充颜色
                # palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                sorting = "ascending",                      
                add = "segments",                             # 添加棒子
                add.params = list(color = "#196F6E", size = 1),#改变棒子参数
                rotate = TRUE,                                # 方向转为垂直
                #group = "type",                                
                dot.size = "abs",                                 # 改变点的大小
                #label = round(go_enrich_df$GeneNumber),                       # 添加label
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               # 设置label参数
                ggtheme = theme_pubr(),                        # 改变主题
                xlab="GE",
                ylab="Weight",
                title = "GE Rank Top 10 in Factor1")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 5,r = 5,b = 5,l = 5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain" ),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black",size = 8,face = "plain"))+
  scale_y_continuous(breaks = seq(-1, 1, 0.5))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)

#绘制基因表达转录层的KEGG图

#选择权重的绝对值大于0.5的基因进行富集分析
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightGEGene <- KEGGa$GeneName
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(rownames(KEGGa), fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
KEGGgene <- df1$ENTREZID
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = KEGGgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)

#按照富集到的基因的个数进行排序
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg@result[["Count"]]

#画图，只展示出KEGG中padj<0.05 and 富集到基因最多的通路
# p1 <- barplot(kegg,
#               x = "Count",
#               color = "p.adjust",
#               showCategory = 5,
#               size = NULL,
#               split = NULL,
#               font.size = 12,
#               title = "GEKEGG",
#               label_format = 30)+
#   scale_size(range=c(2, 12))+
#   theme_bw()+
#   # theme_classic()+
#   labs(x="Count",y="KEGG") +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         #axis.text.y = element_blank()+
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 10,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 10,face = "plain"),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(color = "black",size = 10,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 8,face = "plain"),
#         legend.position = "right",
#         legend.title = element_text(color = "black",size = 6,face = "plain"),
#         legend.text = element_text(color = "black",size = 6,face = "plain"),
#         legend.background = element_rect(fill = "transparent",linewidth = 0.4,linetype = "solid",colour = "transparent"))
# setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
# ggsave(p1, filename = "GE_KEGG.pdf", width = 4.5, height = 2.5)
#转换为表格
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
write.csv(kegg,file = "MOFA2_GE_KEGG.csv",row.names = F)



#KEGG图的另外一种展现形式
kegg <- kegg[1:6,]
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("GE",6)))
str(kegg)
kegg$Description <- as.factor(kegg$Description)
library(stringr)
p1 <- ggplot(kegg,aes(Description,Count)) + 
  geom_bar(aes(fill=factor(type_order)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#008080"))+
  geom_text(aes(label=abs(Count)),color="black", size=5,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("MOFA2 GE KEGG")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)


#

#对AS的数据进行分析#
ASweights <- as.data.frame(model@expectations[["W"]][["AS"]])
#将行名转化为基因名
ASweights$GeneName <- rownames(ASweights)
library(tidyverse)

ASweights <- separate(ASweights,GeneName,into = c("Chr","Start","End","Cluster"),sep = "([:])")
ASweights$Position <- rownames(ASweights)
ASweights <- ASweights[,c(1,20)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
m <- read.table("ASdiff_n369.txt",fill = TRUE,header = T,check.names = F)
dim(m)
colnames(m)
#
m <- m[,c(1,2)]
sum <- merge(ASweights,m,by="Position")
ASweights <- sum
colnames(ASweights) <- c("Position","Factor1","GeneName")
#因为出现了大量的重复值，所有按照数据是否一样去重

#绘制出factor1在不同的view中占比重比较大的因素
#将数据按照weight的绝对值的大小进行排序
ASweights_Pos <- ASweights[ASweights$Factor1>=0,]
ASweights_Pos$Weight <- (ASweights_Pos$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Pos <- ASweights_Pos[order(abs(ASweights_Pos$Factor1),decreasing = T),]

ASweights_Neg <- ASweights[ASweights$Factor1<0,]
ASweights_Neg$Weight <- (ASweights_Neg$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Neg <- ASweights_Neg[order(abs(ASweights_Neg$Factor1),decreasing =T),]

a <- rbind(ASweights_Pos,ASweights_Neg)
KEGGa <- a
#按照权重的绝对值大小进行排序
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$Position#防止出现一个基因对应多个转录本的情况出现

#将上下两个数据合并进行绘图
colnames(ASweights_Pos)
colnames(ASweights_Neg)
ASweights <- rbind(ASweights_Pos,ASweights_Neg)
ASweights <- ASweights[order(ASweights$Factor1,decreasing = F),]
ASweights$Rank <- seq(1,5000,1)
ASweights_Top <- subset(ASweights,ASweights$Position %in% a)
topAS <- ASweights_Top$GeneName
ASweights_Top <- ASweights_Top[order(ASweights_Top$Factor1,decreasing = T),]
#将想要关注的点设置为黑色，不想要关注的点设置为灰色
ASweights$Color <- ifelse(ASweights$Position %in% a, "#940514","grey")
#绘图
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
colnames(ASweights)
library(ggrepel)
p <- ggplot(
  ASweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#940514", "gray")) +
  geom_text_repel(data = ASweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=2000)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  # 坐标轴
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 AS Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),#去掉横坐标上的文本
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)


ASweights_Top$abs <- abs(ASweights_Top$Weight)
colnames(ASweights_Top)
library(ggpubr)
p <- ggdotchart(ASweights_Top, x = "GeneName", y = "Weight",
                color = "#940514",                                # 按照cyl填充颜色
                # palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                sorting = "ascending",                      
                add = "segments",                             # 添加棒子
                add.params = list(color = "#940514", size = 1),#改变棒子参数
                rotate = TRUE,                                # 方向转为垂直
                #group = "type",                                
                dot.size = "abs",                                 # 改变点的大小
                #label = round(go_enrich_df$GeneNumber),                       # 添加label
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               # 设置label参数
                ggtheme = theme_pubr(),                        # 改变主题
                xlab="AS",
                ylab="Weight",
                title = "AS Rank Top 10 in Factor1")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 5,r = 5,b = 5,l = 5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain" ),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black",size = 8,face = "plain"))+
  scale_y_continuous(breaks = seq(-1, 1, 0.5))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)



#绘制AS转录层的KEGG图
#选择权重的绝对值大于0.3的基因进行富集分析
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightASGene <- KEGGa$GeneName
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
# 'select()' returned 1:1 mapping between keys and columns
# Warning message:
#   In bitr(KEGGa$GeneName, fromType = "SYMBOL", toType = c("ENTREZID"),  :
#             2.48% of input gene IDs are fail to map...
KEGGgene <- df1$ENTREZID
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = KEGGgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)

#按照富集到的基因的个数进行排序
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg@result[["Count"]]

#画图，只展示出KEGG中padj<0.05 and 富集到基因最多的通路
barplot(kegg,
        x = "Count",
        color = "p.adjust",
        showCategory = 20,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "MOFA2 AS KEGG",
        label_format = 30)

#转换为表格
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_AS_KEGG.csv",row.names = F)



#KEGG图的另外一种展现形式
kegg <- kegg[1:6,]
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("AS",6)))
str(kegg)
kegg$Description <- as.factor(kegg$Description)
library(stringr)
p1 <- ggplot(kegg,aes(Description,Count)) + 
  geom_bar(aes(fill=factor(type_order)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#AA0019"))+
  geom_text(aes(label=abs(Count)),color="black", size=5,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("MOFA2 AS KEGG")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)



#对APA的数据进行分析#
APAweights <- as.data.frame(model@expectations[["W"]][["APA"]])
APAweights$ENSEMBLTRANS <- rownames(APAweights)
APAweights <- APAweights[,c(1,16)]
#对基因的名字进行转换
setwd("E:/AD_Patient/Dapars2/Dapars0.9")
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
Healthy_AD <- Healthy_AD[,1:2]
#将加载进来的表格按照|进行划分
library(tidyr)
colnames(Healthy_AD)
long<-separate(Healthy_AD,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- long
Healthy_AD <- Healthy_AD[,1:2]
colnames(Healthy_AD)
#将sumAPA文件中转录本的名字转换为基因的名字
APAweights <- merge(Healthy_AD,APAweights,by="ENSEMBLTRANS")

setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
APAweights_Pos <- APAweights[APAweights$Factor1>=0,]
APAweights_Pos$Weight <- (APAweights_Pos$Factor1-min(abs(APAweights$Factor1)))/(max(abs(APAweights$Factor1))-min(abs(APAweights$Factor1)))
APAweights_Pos <- APAweights_Pos[order(abs(APAweights_Pos$Factor1),decreasing = T),]

APAweights_Neg <- APAweights[APAweights$Factor1<0,]
APAweights_Neg$Weight <- (APAweights_Neg$Factor1-min(abs(APAweights$Factor1)))/(max(abs(APAweights$Factor1))-min(abs(APAweights$Factor1)))
APAweights_Neg <- APAweights_Neg[order(abs(APAweights_Neg$Factor1),decreasing =T),]

a <- rbind(APAweights_Pos,APAweights_Neg)
KEGGa <- a
#按照权重的绝对值大小进行排序
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$ENSEMBLTRANS#防止出现一个基因对应多个转录本的情况出现

#将上下两个数据合并进行绘图
colnames(APAweights_Pos)
colnames(APAweights_Neg)
APAweights <- rbind(APAweights_Pos,APAweights_Neg)
APAweights <- APAweights[order(APAweights$Factor1,decreasing = F),]
APAweights$Rank <- seq(1,82,1)
APAweights_Top <- subset(APAweights,APAweights$ENSEMBLTRANS %in% a)
topAPA <- APAweights_Top$GeneName 
#将想要关注的点设置为黑色，不想要关注的点设置为灰色
APAweights$Color <- ifelse(APAweights$ENSEMBLTRANS %in% a, "#2678C2","grey")
#绘图
colnames(APAweights)
library(ggrepel)
p <- ggplot(
  APAweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#2678C2", "gray")) +
  geom_text_repel(data = APAweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=20)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  # 坐标轴
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 APA Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),#去掉横坐标上的文本
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)

#绘制topWeight的图
APAweights_Top$abs <- abs(APAweights_Top$Weight)
colnames(APAweights_Top)
library(ggpubr)
p <- ggdotchart(APAweights_Top, x = "ENSEMBLTRANS", y = "Weight",
                color = "#2678C2",                                # 按照cyl填充颜色
                # palette = c("#0D6CA6","#099963", "#911F27"), # 修改颜色
                sorting = "ascending",                      
                add = "segments",                             # 添加棒子
                add.params = list(color = "#2678C2", size = 1),#改变棒子参数
                rotate = TRUE,                                # 方向转为垂直
                #group = "type",                                
                dot.size = "abs",                                 # 改变点的大小
                #label = round(go_enrich_df$GeneNumber),                       # 添加label
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               # 设置label参数
                ggtheme = theme_pubr(),                        # 改变主题
                xlab="APA",
                ylab="Weight",
                title = "APA Rank Top 10 in Factor1")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        plot.margin = margin(t = 5,r = 5,b = 5,l = 5),
        axis.title.x = element_text(color = "black",size = 10,face = "plain"),
        axis.title.y = element_text(color = "black",size = 10,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 8,face = "plain" ),
        axis.text.y = element_text(color = "black",size = 8,face = "plain"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black",size = 8,face = "plain"))
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)


#绘制APA转录层的KEGG图
#选择权重的绝对值大于0.5的基因进行富集分析
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightAPAGene <- KEGGa$GeneName
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
# 'select()' returned 1:1 mapping between keys and columns
KEGGgene <- df1$ENTREZID
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = KEGGgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)

#按照富集到的基因的个数进行排序
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg@result[["Count"]]

#画图，只展示出KEGG中padj<0.05 and 富集到基因最多的通路
barplot(kegg,
        x = "Count",
        color = "p.adjust",
        showCategory = 20,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "GEKEGG",
        label_format = 30)

#转换为表格
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_APA_KEGG.csv",row.names = F)



#KEGG图的另外一种展现形式
kegg <- kegg[1:6,]
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("APA",6)))
str(kegg)
kegg$Description <- as.factor(kegg$Description)
library(stringr)
p1 <- ggplot(kegg,aes(Description,Count)) + 
  geom_bar(aes(fill=factor(type_order)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#008EDF"))+
  geom_text(aes(label=abs(Count)),color="black", size=5,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("MOFA2 APA KEGG")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 3, units = 'in', dpi = 600)




#将得出所有的top weight的基因进行KEGG分析
gene <- c(topGE,topAS,topAPA)
# topGE
# [1] "YWHAH" "UCHL1" "BEX1"  "STMN2" "VSNL1" "YWHAG" "CALM1" "BASP1" "CALM3" "GAP43"
# > topAS
# [1] "ENC1"    "NEFL"    "RTN3"    "SLC4A10" "NPTN"    "GABRG2"  "KIFAP3"  "SNAP25" 
# [9] "VSNL1"   "KIFAP3" 
# > topAPA
# [1] "MLF2"   "RPL15"  "RPL15"  "PEBP1"  "SARAF"  "SARAF"  "ISCU"   "EIF4A2" "SAP18" 
# [10] "SAP18" 
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(gene, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
KEGGgene <- df1$ENTREZID
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = KEGGgene,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.1,
                   pAdjustMethod  = "BH",
                   qvalueCutoff  = 0.1)

#按照富集到的基因的个数进行排序
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg@result[["Count"]]

#画图，只展示出KEGG中padj<0.05 and 富集到基因最多的通路
barplot(kegg,
        x = "Count",
        color = "p.adjust",
        showCategory = 20,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "GEKEGG",
        label_format = 30)

#转换为表格
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_sigGeneKEGG.csv",row.names = F)
#
#KEGG图的另外一种展现形式
kegg <- kegg[,-7]
kegg <- kegg[1:6,]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("GENE",6)))
str(kegg)
kegg$Description <- as.factor(kegg$Description)
library(stringr)
p1 <- ggplot(kegg,aes(Description,Count)) + 
  geom_bar(aes(fill=factor(type_order)),stat="identity", colour="white",width=0.99)+ 
  scale_fill_manual(values = c("#ED8B38"))+
  geom_text(aes(label=abs(Count)),color="black", size=5,
            position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
  coord_flip()+
  theme_bw()+
  ggtitle("MOFA2 Three layers Genes")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
  # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),#panel.border = element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_ThreeLayerGenes_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 3, units = 'in', dpi = 600)

#将|weight|>0.3的基因与差异基因进行分析
# AllGene <- c(WeightGEGene,WeightASGene,WeightAPAGene)
# #加载基因表达差异分析的差异基因的表
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# difGene <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",header = T)
# #绘制两个overlap的基因
# a <- list(MOFA2AllGene = as.list.data.frame(unique(AllGene)),
#           GE = as.list.data.frame(unique(difGene$X)))
# #画图
# library(ggvenn)
# p <- ggvenn(a, show_elements = FALSE, label_sep = "\n", digits = 1,
#             fill_color = c("9977CC", "firebrick", "3300CC", "red"),
#             fill_alpha = 0.5,
#             stroke_color = "black",
#             stroke_alpha = 1,
#             stroke_size = 1,
#             stroke_linetype = "solid",
#             set_name_color = "black",
#             set_name_size = 4,
#             text_color = "black",
#             text_size = 4)# show elements in line
# p
# ggsave(p,filename = "MOFAAllandDE.pdf", width = 7, height = 5,bg = "white")
# common <- intersect(unique(AllGene), unique(difGene$X))


# 将得出所有的top weight的基因与WGCNA找出的sig基因进行分析进行联合分析
MOFA2_Mos_sigGene <- c(topGE,topAS,topAPA)
# MOFA2_GE_Mos <- topGE
# MOFA2_AS_Mos <- topAS
# MOFA2_APA_Mos <- topAPA
#加载WGCNA得出的关键基因的表
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
#加载GE的结果
GE_Mos_sigGene <- c("SYT1","CHN1","SNAP25","VSNL1","ENC1","TNS1","SGK1","CPM","CLMN","PPFIBP2")

# #加载AS的结果
# a <- read.table("ASdiff_n369.txt",header = T,check.names = F,row.names = 1)
# #根据变化倍数和padj对数据进行筛选
# b <- subset(a,a$change !="Stable")
# #对表格进行排序
# b <- b[order(b$Padjust),]
# head(b)
# b <- b[1:10,]
# AS_Mos_sigGene <- b$Genename
# 
# 
# #加载APA的结果
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# a <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
# b <- subset(a,a$change !="Stable")
# b <- b[order(b$Padjust),]
# head(b)
# b <- b[1:10,]
# #基因的名字
# setwd("E:/AD_Patient/Dapars2/Dapars0.9")
# Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
# Healthy_AD <- Healthy_AD[,1:2]
# library(tidyr)
# colnames(Healthy_AD)
# long<-separate(Healthy_AD,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
# Healthy_AD <- long
# Healthy_AD <- Healthy_AD[,1:2]
# Healthy_AD <- subset(Healthy_AD,Healthy_AD$ENSEMBLTRANS %in% rownames(b))
# APA_Mos_sigGene <- Healthy_AD$GeneName
# setwd("C:/Users/yujie/Desktop/datacollect/20240827")
# 
# Sep_Mos_sigGene <- c(GE_Mos_sigGene,AS_Mos_sigGene,APA_Mos_sigGene)
a <- list(MOFA2_Mos_sigGene = as.list.data.frame(MOFA2_Mos_sigGene),
          GE_Mos_sigGene = as.list.data.frame(GE_Mos_sigGene))
library(ggvenn)
p <- ggvenn(a, show_elements = FALSE, label_sep = "\n", digits = 1,
            fill_color = c("#E24D37", "#3F8CAD", "#009966","#CC3366","green","yellow"),
            fill_alpha = 0.5,
            stroke_color = "black",
            stroke_alpha = 1,
            stroke_size = 0.5,
            stroke_linetype = "solid",
            set_name_color = "black",
            set_name_size = 4,
            text_color = "black",
            text_size = 4)# show elements in line
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("Total_Venn.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(4, "in")), 
       width = 8, height = 8, units = 'in', dpi = 600)

# #绘制6个overlap的基因
# a <- list(MOFA2_GE_Mos = as.list.data.frame(unique(MOFA2_GE_Mos)),
#           GE_Mos_sigGene = as.list.data.frame(unique(GE_Mos_sigGene)),
#           MOFA2_AS_Mos = as.list.data.frame(unique(MOFA2_AS_Mos)),
#           AS_Mos_sigGene = as.list.data.frame(unique(AS_Mos_sigGene)),
#           APA_Mos_sigGene = as.list.data.frame(unique(APA_Mos_sigGene)),
#           MOFA2_APA_Mos = as.list.data.frame(unique(MOFA2_APA_Mos)))
# #画图
# euler_plot <- euler(a,shape = "ellipse")
# plot(euler_plot)
# common1 <- intersect(unique(APA_Mos_sigGene), unique(MOFA2_APA_Mos))
# # [1] "RPL15"
# common2 <- intersect(unique(GE_Mos_sigGene), unique(MOFA2_AS_Mos))
# # [1] "SNAP25" "VSNL1"  "ENC1"
# common3 <- intersect(unique(MOFA2_GE_Mos), unique(GE_Mos_sigGene))
# # [1] "VSNL1"
# common4 <- intersect(unique(MOFA2_GE_Mos), unique(MOFA2_AS_Mos))
# # [1] "VSNL1"
# myNV <- plotVenn(a,outFile="a.svg")


#还有另外一种绘制的方法
# venn_plot <- ggVennDiagram(a, label = "count") +
#   scale_fill_gradientn(colors = c("#E9D5CE", "#F9D2C3", "#FFF9CF","#AED8E6", "#E97F81","#E97F81")) +
#   theme(legend.position = "none")
# 
# venn_plot
# ggsave("Test_Venn.pdf", egg::set_panel_size(venn_plot, width=unit(8, "in"), height=unit(8, "in")), 
#        width = 10, height = 10, units = 'in', dpi = 600)

#




#绘制出这十几个基因的上下调关系
# barGene <- subset(difGene,difGene$X %in% common)
# bar <- as.data.frame(table(barGene$change))
# write.csv(bar,file = "bar.csv",row.names = F,quote = F)
# b <- read.csv("bar.csv",header = T)
# b$Var1 <- factor(b$Var1)
# colnames(b)
# library(ggsci)
# p <- ggplot(b,aes(x=Var1,y=Freq))+
#   geom_bar(aes(fill=factor(Var1)),stat="identity", colour="white",width=0.99)+ 
#   scale_fill_manual(values = c("#128181","#AA191C"))+
#   geom_text(aes(label=Freq),vjust=-0.1,angle = 0, size=4)+  #3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
#   xlab("") +ylab("Number")+
#   theme_bw()+#theme_classic()+
#   ggtitle("") +
#   #coord_flip()+
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         legend.text=element_text(size=8),
#         legend.title=element_text(size=8),
#         plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
#         plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 14,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 14,face = "plain"),
#         panel.grid= element_blank(),
#         axis.text.x = element_text(color = "black",size = 14,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 14,face = "plain"))
# 
# p
# ggsave(p,filename = "table.pdf", width = 3, height = 5)
# 
# #
# library(clusterProfiler)
# library(org.Hs.eg.db)
# df1 <- bitr(common, fromType = "SYMBOL",
#             toType = c("ENTREZID"),
#             OrgDb = org.Hs.eg.db)
# KEGGgene <- df1$ENTREZID
# library(R.utils)
# R.utils::setOption( "clusterProfiler.download.method",'auto' )
# kegg <- enrichKEGG(gene = KEGGgene,
#                    keyType = "kegg",
#                    organism  = 'hsa',
#                    pvalueCutoff  = 1,
#                    pAdjustMethod  = "BH",
#                    qvalueCutoff  = 1)
# 
# #按照富集到的基因的个数进行排序
# kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
# kegg@result[["Count"]]
# 
# #画图，只展示出KEGG中padj<0.05 and 富集到基因最多的通路
# barplot(kegg,
#         x = "Count",
#         color = "p.adjust",
#         showCategory = 20,
#         size = NULL,
#         split = NULL,
#         font.size = 12,
#         title = "GEKEGG",
#         label_format = 30)
# #转换为表格
# kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# kegg <- as.data.frame(kegg)
# write.csv(kegg,file = "MOFA2sigGeneKEGG.csv",row.names = F)
# #
# #KEGG图的另外一种展现形式
# kegg <- kegg[,-7]
# kegg <- kegg[1:6,]
# kegg$num <- rep(1:6)
# kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("GENE",6)))
# str(kegg)
# kegg$Description <- as.factor(kegg$Description)
# library(stringr)
# p1 <- ggplot(kegg,aes(Description,Count)) + 
#   geom_bar(aes(fill=factor(type_order)),stat="identity", colour="white",width=0.99)+ 
#   scale_fill_manual(values = c("#ED8B38"))+
#   geom_text(aes(label=abs(Count)),color="black", size=2.5,
#             position =position_dodge(width = 1),vjust = 0.5,inherit.aes = TRUE)+
#   coord_flip()+
#   theme_bw()+
#   ggtitle("GENE")+
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
#   # geom_hline(yintercept = 0,lty=1,col="#666666",lwd=0.5)+
#   xlab("Pathway")+ylab("Number")+
#   guides(fill = guide_legend(title = 'Layer'))+
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
#         legend.key = element_blank(),
#         plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
#         axis.title.x = element_text(color = "black",size = 8,face = "plain"),
#         axis.title.y = element_text(color = "black",size = 8,face = "plain"),
#         panel.grid= element_blank(),#panel.border = element_blank(),
#         axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
#         axis.text.x = element_text(color = "black",size = 8,face = "plain"),
#         axis.text.y = element_text(color = "black",size = 8,face = "plain"),
#         legend.title = element_text(color = "black",size = 8,face = "plain"),
#         legend.text = element_text(color = "black",size = 8,face = "plain"),
#         legend.direction = 'vertical', legend.position ="right")
# p1
# ggsave(p1, filename = "MOFA2GENEKEGG.pdf", width = 5.2, height = 2.5)
# 
















#

# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##绘制出最top的几个feature
# model@dimensions
# plot_top_weights(model,
#                  view = "AS",
#                  factor = 1,
#                  nfeatures = 10)
# 
# #因为这里我没有设置分组，所以只会出现一个点图
# plot_factor(model, 
#             factor = 1, 
#             color_by = "VGF",
#             add_violin = TRUE,
#             dodge = TRUE)
# 
# #用热图的形式展现出来
# # Heatmaps
# #热图中行表示不同的feature，列表示不同的sample
# plot_data_heatmap(model,
#                   view = "AS",  #view2       # view of interest
#                   factor = 1,             # factor of interest
#                   features = 20,          # number of features to plot (they are selected by weight)
#                   
#                   # extra arguments that are passed to the `pheatmap` function
#                   cluster_rows = TRUE, cluster_cols = FALSE,
#                   show_rownames = TRUE, show_colnames = FALSE)
# 
# 
# # Scatter plots
# # It is useful to add a linear regression estimate to visualise if the relationship between (top) features and factor values is linear.
# #这个图只有一维的坐标，表示的数据是这个feature在所有的sample中的值是多少
# plot_data_scatter(model,
#                   view = "AS",         # view of interest
#                   factor = 1,             # factor of interest
#                   features = 9,           # number of features to plot (they are selected by weight)
#                   add_lm = TRUE,          # add linear regression
#                   color_by = "Severity")
# 
# 
# set.seed(42)#为什么要设置种子???
# # model <- run_umap(model)
# model <- run_tsne(model)
# # Plot non-linear dimensionality reduction
# #画的TSEN的图能够看出来什么结果???
# plot_dimred(model,
#             method = "TSNE",  # method can be either "TSNE" or "UMAP"
#             color_by = "condition")
# model <- run_umap(model)
# plot_dimred(model,
#             method = "UMAP",  # method can be either "TSNE" or "UMAP"
#             color_by = "condition")
# 
# # The user can rename the dimensions of the model
# views_names(model) <- c("GeneExpression", "APA","AS")
# factors_names(model) <- paste("Factor", 1:model@dimensions[["K"]], sep=" ")
# views_names(model)
# # Extracting data for downstream analysis
# 
# #将因子的数据提取出来
# # "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
# factors <- get_factors(model, factors = "all")
# lapply(factors,dim)#表示的是单个组里有10个因子，有100个样本
# # $group1
# # [1] 354  15
# # "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
# weights <- get_weights(model, views = "all", factors = "all")
# lapply(weights,dim)#说明每一个view有1000个feature和10个因子
# # $GeneExpression
# # [1] 1347   15
# # 
# # $APA
# # [1] 195  15
# # 
# # $AS
# # [1] 2000   15
# data <- get_data(model)
# lapply(data, function(x) lapply(x, dim))#[[1]]
# 
# # For convenience, the user can extract the data in long data.frame format:
# factors <- get_factors(model, as.data.frame = T)
# head(factors, n=3)
# write.csv(factors,file = "AllFactors.csv",quote = F,row.names = F)
# weights <- get_weights(model, as.data.frame = T)
# head(weights, n=3)
# write.csv(weights,file = "AllWeights.csv",quote = F,row.names = F)
# data <- get_data(model, as.data.frame = T)
# head(data, n=3)
# 
# #绘制factor和factor2之间的点状图
# #将factors的表格和sum的表格进行合并
# head(factors)
# factors <- factors[,1:3]
# #将长表格换成宽表格
# library(reshape2) # 使用的函数 melt & dcast
# wide<-dcast(factors,sample~factors$factor,value.var = 'value')
# #只取出宽表格的前3列
# wide <- wide[,1:3]
# head(sum)
# scatter <- merge(sum,wide,by="sample")
# #使用生成的scatter表格绘制点状图
# colnames(scatter)
# scatter$CERAD <- as.factor(scatter$CERAD)
# scatter$Braak_Stage <- as.factor(scatter$Braak_Stage)
# scatter$Severity <- as.factor(scatter$Severity)
# scatter$Sex <- as.factor(scatter$Sex)
# scatter$Phenotypedit <- as.factor(scatter$Phenotypedit)
# colnames(scatter)
# 
# #画图
# # ,shape=Severity
# library(ggplot2)
# p1 <- ggplot(data=scatter,aes(x=Factor1,y=Factor2,color=Phenotypedit))+
#   geom_point(size=3)+
#   scale_x_continuous(limits = c(-5, 5))+
#   scale_y_continuous(limits = c(-3,6))+
#   ggtitle("Scatter Plot")+
#   # scale_color_manual(values = c("#911F27", "#3366CC","#009933"))+
#   theme_bw()+
#   theme(
#     legend.key = element_blank(),
#     plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
#     axis.title.x = element_text(color = "black",size = 10,face = "plain"),
#     axis.title.y = element_text(color = "black",size = 10,face = "plain"),
#     panel.grid= element_blank(),
#     axis.text.x = element_text(color = "black",size = 10,face = "plain"),
#     axis.text.y = element_text(color = "black",size = 10,face = "plain"),
#     legend.title = element_text(color = "black",size = 10,face = "plain"),
#     legend.text = element_text(color = "black",size = 10,face = "plain"),
#     legend.background = element_rect(fill = "white",size = 0.2,linetype = "solid",colour = "black"),
#     legend.direction = 'vertical', legend.position = "right")+
#   labs(x=paste0("Factor1"),
#        y=paste0("Factor2"))
# 
# p1
# 
# 


#将3个层次的数据结合起来进行ROC分析
rm(list = ls())
#加载GE的数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
GEFPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
#去掉rowname的小数点
rownames(GEFPKM) <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(GEFPKM))
GEFPKM <- as.data.frame(t(GEFPKM))
GEFPKM$GSM_number <- rownames(GEFPKM)
GEFPKM <- GEFPKM[,c("ENSG00000163032","GSM_number")]
#根据背景信息对数据进行更换
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
sum <- merge(GEFPKM,BackgroundInformation,by="GSM_number")
sum <- sum[,c("ENSG00000163032","Sample")]
GE <- sum
GE <- GE[order(GE$Sample),]

#加载AS的数据
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AScount <- read.csv("AS_RemoveBatch_VoomCount.csv",header = T,row.names = 1,check.names = F)
#然后筛选出来MOFA2分析出的那个层次中的两个intron
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)
ASweights <- as.data.frame(model@expectations[["W"]][["AS"]])
#将行名转化为基因名
ASweights$GeneName <- rownames(ASweights)
#绘制出factor1在不同的view中占比重比较大的因素
#将数据按照weight的绝对值的大小进行排序
ASweights_Pos <- ASweights[ASweights$Factor1>=0,]
ASweights_Pos$Weight <- (ASweights_Pos$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Pos <- ASweights_Pos[order(abs(ASweights_Pos$Factor1),decreasing = T),]

ASweights_Neg <- ASweights[ASweights$Factor1<0,]
ASweights_Neg$Weight <- (ASweights_Neg$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Neg <- ASweights_Neg[order(abs(ASweights_Neg$Factor1),decreasing =T),]

a <- rbind(ASweights_Pos,ASweights_Neg)
KEGGa <- a
#按照权重的绝对值大小进行排序
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a$GeneName
#提取出来每一个基因的名字
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
m <- read.table("ASdiff_n369.txt",fill = TRUE,header = T,check.names = F)
m <- m[m$Position %in% a$GeneName,]
m <- subset(m,m$Genename %in% c("VSNL1","ENC1","SNAP25"))
rownames(m) <- m$Genename
colnames(m)
m <- m[,c(3:356)]
m <- as.data.frame(t(m))
m$Sample <- rownames(m)

#使用background信息换成AD和healthy
sum <- merge(GE,m,by="Sample")
sum$Severity <- substr(sum$Sample,1,1)
sum <- sum[,-1]
#进行ROC分析
sum$Severity <- as.numeric(sum$Severity)
str(sum)
#如果Severity是因子数据，则最后随机抽取的不同的folder是一样的个数，
#如果Severity是numeric数据，则最后随机抽取的相同的folder是一样的个数

table(sum$Severity)
# AD Healthy 
# 257      97
all <- sum
table(all$Severity)
library(tidyverse)
n <- dim(all)[1]
y <- all$Severity

require(caret)
library(ROCR)
library(magrittr) # pipe operator
library(plyr)

# 定义不同的随机种子
seeds <- seq(1:1000)

# 创建或打开文件（使用追加模式）
setwd("C:/Users/yujie/Desktop/datacollect")
output_file <- "auc_results2.txt"
fileConn <- file(output_file, open = "a") # 追加模式

# 写入表头（仅在文件为空时写入）
if (file.info(output_file)$size == 0) {
  writeLines("Seed\tAUC Values\tMean AUC\n", fileConn)
}


for (seed in seeds) {
  set.seed(seed) # 每次使用新的随机种子
  
  # 创建 6 折交叉验证
  folds <- createFolds(y, k = 6)
  auc_values <- numeric() # 存储当前种子的 AUC 值
  
  for (i in 1:6) {
    fold_test <- all[folds[[i]],] # 测试集
    fold_train <- all[-folds[[i]],] # 训练集
    
    # 训练模型
    model <- glm(fold_train$Severity ~ ., data = fold_train, family = binomial(link = "logit"))
    
    # 预测并计算 AUC
    fold_predict <- stats::predict(model, type = 'response', newdata = fold_test)
    pred <- prediction(fold_predict, fold_test$Severity)
    auc_value <- as.numeric(performance(pred, measure = "auc")@y.values[[1]])
    auc_values <- c(auc_values, auc_value)
  }
  
  # 计算当前种子的 AUC 均值
  mean_auc <- mean(auc_values)
  
  # 写入文件：包含当前 seed、AUC 值列表及均值，按要求格式输出
  # 写入文件：包含当前 seed、AUC 值列表及均值
  auc_text <- paste(auc_values, collapse = ", ")
  writeLines(paste0(seed, "\t", auc_text, "\t", mean_auc), fileConn)
}

# 关闭文件
close(fileConn)

#从文件中的结果可以看出，最终选择的是938
###开始随机抽样
set.seed(717)
require(caret)
folds <- createFolds(y,k=6)

library(ROCR)
library(magrittr) # pipe operator
library(plyr)
auc_value<-as.numeric()
#将所有的fitvalue取出来
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] #取folds[[i]]作为测试集
  fold_train <- all[-folds[[i]],] # 剩下的数据作为训练集
  model <- glm(fold_train$Severity~.,data=fold_train,
               family = binomial(link = "logit"))
  #注意mofa2中也有一个predict函数，会和原来的predict函数冲突
  fold_predict <- stats::predict(model,type='response',newdata=fold_test)#样本点并没有减少
  pred <- prediction(fold_predict,fold_test$Severity)#经过这一步样本点减少了
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
auc_value
# [1] 0.7993311 0.6542553 0.9109848 0.5955711 0.7268170 0.8349206
mean(auc_value)
# #10次计算
# [1] 0.7536467
#去掉数据框rocData的第一行
rocData <- rocData[-1,]
#更改数据框rocData的行名

rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 


rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
#计算均值
rocData$Xmean <- apply(rocData[,1:6],1,mean)
rocData$Ymean <- apply(rocData[,7:12],1,mean)
rocData <- rocData[,13:14]
rocData$Gene <- c(rep("all", 60))

#计算95%的置信区间
#计算不出来
# library(magrittr) # pipe operator
# rocit_emp <- rocit(score = rocData$Ymean, class = all$BackgroundInformation[test], method = "emp")
# summary(rocit_emp)
# #计算CI值
# ciAUC(rocit_emp)


#构造样本ROC的终点数据
#构造中间虚线的数据
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("a", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
#合并数据
ROC <- rbind(rocData,inn)
#画图
library(ggrepel)
library(ggplot2)
a<- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene)) + 
  geom_line(aes(linetype=Gene,color=Gene),size=1)+
  scale_linetype_manual(values = c(2,7))+
  scale_color_manual(values=c("black","#00872D"))+
  xlab("False Positive Rate") +
  ylab("True Positive Rate")+
  theme_bw()+
  ggtitle("MOFA2 ROC AUC=0.754") +
  #scale_x_continuous(limits = c(0, 1),expand = c(0,0))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0,1))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        legend.text=element_text(size=22),
        legend.title = element_text(size=8),
        legend.position = c(0.8,0.25),
        plot.margin = margin(t = 10,r = 10,b = 10,l = 10),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 22,face = "plain"),
        axis.title.y = element_text(color = "black",size = 22,face = "plain"),
        panel.grid= element_blank(),
        axis.text.x = element_text(color = "black",size = 22,face = "plain"),
        axis.text.y = element_text(color = "black",size = 22,face = "plain"))
a
ggsave("MOFA2_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)








#寻找共同存在可变剪切和基因表达差异变化的基因,
# 这一步本身就没有什么意义，每一层的数据表示的含义都不一样，所以根本不能做unio的分析









