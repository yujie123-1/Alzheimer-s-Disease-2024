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
sapply(sum, function(x) sum(is.na(x)))
str(sum)
sum$Braak_Stage <- as.factor(sum$Braak_Stage)
sum$Brain_Region <- as.factor(sum$Brain_Region)
sum$Severity <- as.factor(sum$Severity)
sum$Bank_Location <- as.factor(sum$Bank_Location)
table(sum$Brain_Region)
set.seed(1)
sumKnn <- VIM::kNN(sum,trace = T,k = 15)
sumKnn$Age <- round(sumKnn$Age)
sumKnn$Age <- as.integer(sumKnn$Age)
which(sumKnn$Brain_Region != sum$Brain_Region)
sum <- read.csv("sum.csv",header = T,fill = T)
data <- cbind(sum[,c(2:5)],sumKnn[,1:5])
summary(data)
anyNA(data)
setwd("C:/Users/yujie/Desktop/datacollect")
write.csv(data,file = "sumkNN20240827_removeSex.csv",row.names = F,quote = F)
sum_compelete <- sum_compelete[order(sum_compelete$GSM_number),]
sum_compelete <- sum_compelete[!duplicated(sum_compelete$GSM_number),]
sum_compelete <- sum_compelete[,-c(2:5,7,10)]
sum_compelete$Braak_Stage <- as.factor(sum_compelete$Braak_Stage)
sum_compelete$Brain_Region <- as.factor(sum_compelete$Brain_Region)
sum_compelete$Severity <- as.factor(sum_compelete$Severity)
sum_compelete$Bank_Location <- as.factor(sum_compelete$Bank_Location)
sum_compelete$Age <- as.numeric(sum_compelete$Age)
colnames(sum_compelete)
set.seed(1)  
missing_ratios <- list(Braak_Stage = 0.2, Age = 0.12 ) 
missing_ratios <- list(Braak_Stage = 0.05, Age = 0.05) 
missing_ratios <- list(Braak_Stage = 0.09, Age = 0.09 ) 
missing_ratios <- list(Braak_Stage = 0.1, Age = 0.1) 
missing_ratios <- list(Braak_Stage = 0.12, Age = 0.12) 
missing_ratios <- list(Braak_Stage = 0.15, Age = 0.15) 
missing_ratios <- list(Braak_Stage = 0.2, Age = 0.2 ) 
missing_ratios <- list(Braak_Stage = 0.25, Age = 0.25)
missing_ratios <- list(Braak_Stage = 0.30, Age = 0.30) 
df <- sum_compelete
for (col in names(missing_ratios)) {
  num_values <- nrow(df)  
  num_missing <- round(missing_ratios[[col]] * num_values)  
  set.seed(1)
  missing_indices <- sample(1:num_values, num_missing)
  df[missing_indices, col] <- NA
}
k_values <- c(5,8,10,15, 20,25, 30,35,40,45,50)
k_values <- 15
results <- data.frame(k = integer(),
                      BA_accuracy = numeric(),
                      mse = numeric(),
                      mae = numeric(),
                      stringsAsFactors = FALSE)
na_indices <- which(is.na(df), arr.ind = TRUE)
for (k in k_values) {
  set.seed(1)
  df_Knn <- VIM::kNN(df, trace = F, k = k)
  df_Knn$Age <- round(df_Knn$Age, 0)
  na_values_from_Knn <- df_Knn[na_indices]
  Braak_Stage_Filled <- na_values_from_Knn[na_indices[, "col"] == 3]
  Age_Filled <- na_values_from_Knn[na_indices[, "col"] == 4]
  na_values_from_sum <- sum_compelete[na_indices]
  Braak_Stage_original <- na_values_from_sum[na_indices[, "col"] == 3]
  Age_original <- na_values_from_sum[na_indices[, "col"] == 4]
  BA_accuracy <- table(Braak_Stage_original == Braak_Stage_Filled)[2] / length(Braak_Stage_original)
  Age_Filled <- as.numeric(Age_Filled)
  Age_original <- as.numeric(Age_original)
  mse_age_filled <- mean((Age_Filled - Age_original)^2, na.rm = TRUE)
  mae_age_filled <- mean(abs(Age_Filled - Age_original), na.rm = TRUE)
  
  results <- rbind(results, data.frame(k = k,
                                       BA_accuracy = BA_accuracy,
                                       mse = mse_age_filled,
                                       mae = mae_age_filled))
  
}
print(results)

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
gene_Length <- read.table(file = "HumanGeneLength.txt",header = T)
count <- read.csv(file="ADPatientCount.csv",header = T)
merge <- merge(count,gene_Length,by = 'Geneid') 
dim(merge)
count <- merge[,1:(dim(merge)[2]-1)]
gene_num <- dim(merge)[1]
sample_num <- dim(merge)[2]-2 
i <- 2
repeat{
  mapped_reads <- sum(merge[1:gene_num,i])
  FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])
  FPKM <- pmax(FPKM,0)
  count = data.frame(count[1:gene_num,],FPKM)
  i <- i + 1
  if(i > sample_num+1){
    break
  }
}
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(FPKM,"ADPatientFPKM.csv",row.names = FALSE, quote = FALSE)
setwd("C:/Users/yujie/Desktop/datacollect")
SampleInformtion <- read.csv("sumkNN20240827_removeSex.csv",header = T)
version <- read.csv("BackgroundInformationAPOE.csv",header = T)
colnames(version)
merged_file <- merge(SampleInformtion, version[, c("GSM_number", "Sample","Batch","Phenotypedit")], by = "GSM_number")
write.csv(merged_file,file = "BackgroundInformationkNN20240827_k_15.csv",row.names = F,quote = F)
setwd("C:/Users/yujie/Desktop/datacollect")
old_version <- read.csv("BackgroundInformationAPOE.csv",header = T)
old_version <- old_version[order(old_version$GSM_number),]
colnames(old_version)
new_version <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
new_version <- new_version[order(new_version$GSM_number),]
colnames(new_version)
identical(old_version$GSM_number,new_version$GSM_number)
length(which(old_version$Braak_Stage != new_version$Braak_Stage))
length(which(old_version$Age != new_version$Age))
length(which(old_version$Bank_Location != new_version$Bank_Location))
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
FPKM[1:5,1:5]
dim(FPKM)
FPKM <- as.matrix(FPKM)
FPKM[FPKM<=1] <- 1
FPKM <- as.data.frame(FPKM)
FPKM <- na.omit(FPKM)
FPKM <- log10(FPKM)
head(colnames(FPKM))
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
a <- summary(pca.info) 
b <- as.data.frame(a[["importance"]])
head(pca.info$rotation) 
head(pca.info$sdev) 
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
colnames(pca.data)[1] <- "Sample"
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="confidence")
b <- merge(SampleInformtion,pca.data,by="Sample")
b$Brain_Region <- as.factor(b$Brain_Region)
library(tidyverse)
pca.data <- pca.data%>%mutate(Type=case_when(pca.data$Type==0~"Healthy",
                                             pca.data$Type==1~"AD",
                                             pca.data$Type==2~"MCI"))
str(pca.data$Type)
pca.data$Type <- as.factor(pca.data$Type)
library(ggforce)
p1 <- ggplot(data=b,aes(x=PC1,y=PC2,color=Type))+
  geom_point(size=3)+
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",25,"%"),
       y=paste0("PCA2 ",23,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AllSamplePCAFPKM.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
FPKM <- as.data.frame(t(FPKM))
FPKM <- FPKM[order(rownames(FPKM)),]
setwd("C:/Users/yujie/Desktop/datacollect")
sumknnImputation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sumknnImputation <- sumknnImputation[order(sumknnImputation$GSM_number),]
identical(rownames(FPKM),sumknnImputation$GSM_number)
sumknnImputation <- sumknnImputation[,-c(1,2)]
colnames(sumknnImputation)
sumknnImputation$Brain_Region <- as.factor(sumknnImputation$Brain_Region)
sumknnImputation$Batch <- factor(sumknnImputation$Batch)
sumknnImputation$Braak_Stage <- factor(sumknnImputation$Braak_Stage)
sumknnImputation$Bank_Location <- factor(sumknnImputation$Bank_Location)
sumknnImputation$LibraryLayout <- factor(sumknnImputation$LibraryLayout)
sumknnImputation$Severity <- factor(sumknnImputation$Severity)
library(vegan)
set.seed(1)
a <- adonis2(FPKM ~ Severity+Braak_Stage+Brain_Region+Age+Batch+Bank_Location, 
             data = sumknnImputation,
             permutations=999, by="margin",method = "gower")
b <- as.data.frame(a)
b <- na.omit(b)
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- b%>%mutate(category=case_when(b$category=="Brain_Region"~"Region",
                                   b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=0.2,angle = 0, size=4)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("All sample MANOVA") +
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

rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sample <- read.csv("SampleDistribution.csv",header = T)
colnames(sample)
p <- ggplot(sample,aes(x=reorder(Region,-Number),y=Number,fill=Region))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Number),vjust=-0.2,angle = 0, size=3)+  
  theme_bw()+
  ggtitle("All Sample") +
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
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
FPKM <- FPKM[,colnames(FPKM) %in% BackgroundInformation$GSM_number]
FPKM <- as.data.frame(t(FPKM))
FPKM$GSM_number <- rownames(FPKM)
sum <- merge(FPKM,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
library(gmodels)
pca.info <- fast.prcomp(sum)
head(pca.info)
head(summary(pca.info)) 
head(pca.info$rotation) 
head(pca.info$sdev) 
head(pca.info$x)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
library(ggpubr)
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
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
myFPKM <- read.csv("ADPatientFPKM.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- subset(myFPKM,rownames(myFPKM) %in% BackgroundInformation$GSM_number)
myFPKM <- myFPKM[order(rownames(myFPKM)),]
setwd("C:/Users/yujie/Desktop/datacollect")
sumKnn <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
sumKnn <- subset(sumKnn,sumKnn$GSM_number %in% BackgroundInformation$GSM_number)
sumKnn <- sumKnn[order(sumKnn$GSM_number),]
identical(rownames(myFPKM),sumKnn$GSM_number)
sample_information <- sumKnn
colnames(sample_information)
sample_information$Brain_Region <- as.factor(sample_information$Brain_Region)
sample_information$Braak_Stage <- factor(sample_information$Braak_Stage)
sample_information$Bank_Location <- factor(sample_information$Bank_Location)
sample_information$LibraryLayout <- factor(sample_information$LibraryLayout)
sample_information$Severity <- factor(sample_information$Severity)
sample_information$Batch <- factor(sample_information$Batch)
library(vegan)
set.seed(1)
a <- adonis2(myFPKM ~ Severity+Braak_Stage+Age+Batch+Bank_Location, data = sample_information,
             permutations=999, by="margin",method = "euclidean")
a
b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age",
                                   b$category=="Sex"~"Gender"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=0.3,angle = 0, size=4)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("GE Raw MANOVA") +
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
rm(list = ls())
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount$GSM_number <- rownames(mycount)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) 
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  
coldata
library(DESeq2)  
max(mycount)
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  
vst <- vst(dds, blind=T)
head(assay(vst), 3)
library(ggforce)
plotPCA(vst)
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
colnames(p1data)
p1data <- p1data%>%mutate(Severity=case_when(p1data$Severity==0~"Healthy",
                                             p1data$Severity==1~"AD",
                                             p1data$Severity==2~"MCI"))
p1data$Severity <- as.factor(p1data$Severity)
p1data$Batch <- as.factor(p1data$Batch)
p1data$Bank_Location <- as.factor(p1data$Bank_Location)
library(ggrepel)
library(ggplot2)
library(ggforce)
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",32,"%"),
       y=paste0("PCA2 ",17,"%"))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_Raw_PCACount.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
BackgroundInformation$Severity <- as.factor(BackgroundInformation$Severity)
BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount <- subset(mycount,rownames(mycount) %in% BackgroundInformation$GSM_number)
mycount <- mycount[order(rownames(mycount)),]
identical(rownames(mycount),BackgroundInformation$GSM_number)
mycount <- as.matrix(t(mycount))
library(sva)
adjust_counts <- ComBat_seq(mycount, batch=BackgroundInformation$Batch,group=BackgroundInformation$Severity,covar_mod=NULL,full_mod=FALSE)
head(adjust_counts)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(adjust_counts,file = "GE_Count_RemoveBatch.csv",quote = F)
mycount <- as.data.frame(t(adjust_counts))
mycount$GSM_number <- rownames(mycount)
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
sum[5,67940:67951]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67951)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) 
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  
library(DESeq2)  
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition) 
vst <- vst(dds, blind=T)
library(ggforce)
plotPCA(vst)
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"
p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
p1data$Severity <- as.factor(p1data$Severity)
p1data <- p1data%>%mutate(Severity=case_when(p1data$Severity==0~"Healthy",
                                             p1data$Severity==1~"AD",
                                             p1data$Severity==2~"MCI"))
p1data$Batch <- as.factor(p1data$Batch)
p1data$Sex <- as.factor(p1data$Sex)
p1data$Bank_Location <- as.factor(p1data$Bank_Location)
p1data <- p1data%>%mutate(Bank_Location=case_when(p1data$Bank_Location==0~"America",
                                                  p1data$Bank_Location==3~"Australia"))
p1data$Braak_Stage <- as.factor(p1data$Braak_Stage)
p1data <- p1data%>%mutate(Braak_Stage=case_when(p1data$Braak_Stage==0~"0",
                                                p1data$Braak_Stage==1~"transentorhinal_stage",
                                                p1data$Braak_Stage==2~"limbic_stage",
                                                p1data$Braak_Stage==3~"isocortical_stage"))
p1data$Age <- as.numeric(p1data$Age)
library(ggrepel)
library(ggplot2)
library(ggforce)
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
        legend.direction = 'vertical', legend.position = "right")+
  labs(x=paste0("PCA1 ",16,"%"),
       y=paste0("PCA2 ",10,"%"))
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_RemoveBatch_PCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
dim(adjust_counts)
count <- as.data.frame(adjust_counts)
count$Geneid <- rownames(count)
gene_Length <- read.table(file = "HumanGeneLength.txt",header = T)
merge <- merge(count,gene_Length,by = 'Geneid') 
dim(merge)
count <- merge[,1:(dim(merge)[2]-1)]
gene_num <- dim(merge)[1]
sample_num <- dim(merge)[2]-2 
i <- 2
repeat{
  mapped_reads <- sum(merge[1:gene_num,i])
  FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])
  FPKM <- pmax(FPKM,0)
  count = data.frame(count[1:gene_num,],FPKM)
  i <- i + 1
  if(i > sample_num+1){
    break
  }
}
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(FPKM,"GE_RemoveBatch_FPKM.csv",row.names = FALSE, quote = FALSE)
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myFPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- myFPKM[order(rownames(myFPKM)),]
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/STARFeatureCount/result/ReAnalysis")
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
identical(rownames(myFPKM),BackgroundInformation$GSM_number)
BackgroundInformation$Brain_Region <- as.factor(BackgroundInformation$Brain_Region)
BackgroundInformation$Braak_Stage <- factor(BackgroundInformation$Braak_Stage)
BackgroundInformation$Bank_Location <- factor(BackgroundInformation$Bank_Location)
BackgroundInformation$LibraryLayout <- factor(BackgroundInformation$LibraryLayout)
BackgroundInformation$Severity <- factor(BackgroundInformation$Severity)
BackgroundInformation$Batch <- as.factor(BackgroundInformation$Batch)
str(BackgroundInformation)
colnames(BackgroundInformation)
library(vegan)
set.seed(1)
a <- adonis2(myFPKM ~ Severity+Braak_Stage+Age+Bank_Location+Batch, data = BackgroundInformation,
             permutations=999, by="margin",method = "euclidean")
a
b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b$category
colnames(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("GE Remove Batch MANOVA") +
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
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(mycount))
rownames(mycount) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
mycount$ENSEMBL <- rownames(mycount)
df <- merge(df1,mycount,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
mycount <- df[,-c(1:2)]
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) 
condition
coldata <- data.frame(row.names =colnames(mycount),condition)  
coldata
library(DESeq2)  
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  
library(edgeR)
keep <- rowSums(cpm(mycount)>0.5) >= 20
dds <- dds[keep,]
dds_norm <- DESeq(dds)
dds_norm
resultsNames(dds_norm)
res1 <- results(dds_norm,contrast = c("condition","1","0"))
head(res1)
res1 <- res1[order(res1$padj,decreasing = F),]
head(res1)
sum(res1$padj<0.01, na.rm = TRUE)
summary(res1,alpha = 0.01)
res1 <-na.omit(res1)
res1 <- res1[!duplicated(rownames(res1)),]
dim(res1)
res1 <- as.data.frame(res1)
AD_Healthyres <- res1
AD_Healthyres$change <- ifelse(AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >= 0.585, 
                               ifelse(AD_Healthyres$log2FoldChange > 0.585 ,'Up','Down'),'Stable')
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres,file = "GE_RemoveBatch_DEseq2.csv")
AD_Healthyres_def <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585 )
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_def,file = "GE_RemoveBatch_DEseq2_Significant.csv")
df1 <- bitr(rownames(AD_Healthyres_def), fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
GOgene=as.character(df1$ENTREZID)
GOgene <- na.omit(GOgene)
AD_Healthyres_GOplot <- enrichGO(GOgene, 
                                 OrgDb = org.Hs.eg.db, 
                                 ont='All',
                                 pAdjustMethod = 'BH',
                                 pvalueCutoff = 0.1, 
                                 qvalueCutoff = 0.1,
                                 keyType = 'ENTREZID')
AD_Healthyres_GOplot_genelist<-setReadable(AD_Healthyres_GOplot, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_GOplot_genelist <- as.data.frame(AD_Healthyres_GOplot_genelist)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_GOplot_genelist, file='GE_GO_genelist.csv')
AD_Healthyres_GOBP <- enrichGO(GOgene, 
                               OrgDb = org.Hs.eg.db, 
                               ont='BP',
                               pAdjustMethod = 'BH',
                               pvalueCutoff = 0.1, 
                               qvalueCutoff = 0.1,
                               keyType = 'ENTREZID')
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
  labs(x="Count",y="GO BP") +
  theme(
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
                               qvalueCutoff = 0.1,
                               keyType = 'ENTREZID')
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
  labs(x="Count",y="GO CC") +
  theme(
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
  labs(x="Count",y="GO MF") +
  theme(
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
library(ggpubr)
p <- ggarrange(p1, p2,p3, ncol = 3, nrow = 1)+
  theme(plot.margin = margin(t = 20,r = 20,b = 20,l = 20))
go_enrich_df <- data.frame(
  ID=c(AD_Healthyres_GOBP$ID[1:15], AD_Healthyres_GOCC$ID[1:15], AD_Healthyres_GOMF$ID[1:15]),
  Description=c(AD_Healthyres_GOBP$Description[1:15],AD_Healthyres_GOCC$Description[1:15],AD_Healthyres_GOMF$Description[1:15]),
  GeneNumber=c(AD_Healthyres_GOBP$Count[1:15], AD_Healthyres_GOCC$Count[1:15], AD_Healthyres_GOMF$Count[1:15]),
  type=factor(c(rep("biological process", 15), 
                rep("cellular component", 15),
                rep("molecular function", 15)), 
              levels=c("biological process", "cellular component","molecular function" )))
go_enrich_df <- arrange(go_enrich_df,type,GeneNumber)
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#0D6CA6","#099963", "#911F27")
library(ggpubr)
p4 <- ggdotchart(go_enrich_df, x = "type_order", y = "GeneNumber",
                 color = "type",                                
                 palette = c("#0D6CA6","#099963", "#911F27"), 
                 sorting = "descending",                      
                 add = "segments",                            
                 add.params = list(color = "type", size = 1.3),
                 group = "type",                                
                 dot.size = "GeneNumber",                                
                 font.label = list(color = "black", size = 7,
                                   vjust = 0.5),              
                 ggtheme = theme_pubr(),                      
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
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
AD_Healthyres_kegg <- enrichKEGG(gene = GOgene,
                                 keyType = "kegg",
                                 organism  = 'hsa',
                                 pvalueCutoff  = 0.1,
                                 pAdjustMethod  = "BH",
                                 qvalueCutoff  = 0.1)
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
  labs(x="Count",y="KEGG") +
  theme(
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
AD_Healthyres_kegg<-setReadable(AD_Healthyres_kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
AD_Healthyres_kegg <- as.data.frame(AD_Healthyres_kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AD_Healthyres_kegg,file = "AD_Healthyres_kegg_genelist0.05and0.585.csv",row.names = F)
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
table1 <- AD_Healthyres_def[,c(1,2,7)]
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  down_genes <- sum(gene_list %in% rownames(AD_Healthyres_def)[AD_Healthyres_def$change == "Down"])
  up_genes <- sum(gene_list %in% rownames(AD_Healthyres_def)[AD_Healthyres_def$change == "Up"])
  return(c(down_genes, up_genes))
}
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))
table2$down <- -abs(result[, 1])
table2$up <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "GE_KEGG.csv",row.names = F)
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:15,c(3,4,12,13)] 
colnames(KEGGTerm)
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
  xlab("Pathway")+ylab("Gene Number")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  guides(fill = guide_legend(title = 'Change'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 15,face = "plain"),
        axis.title.y = element_text(color = "black",size = 15,face = "plain"),
        panel.grid= element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 15,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
p1
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("GE_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)
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
  AD_Healthyres, aes(x = log2FoldChange, y = -log10(padj),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080", "gray", "firebrick3")) +
  geom_text_repel(data = rbind(up, down), aes(x = log2FoldChange, y = -log10(padj), label = rownames(a)),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
  labs(x="log2FoldChange",y="-log10(padj)") +
  ggtitle("AD_Healthy volcano")+
  theme_bw()+
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
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AD_Healthyres <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
table(AD_Healthyres$padj<0.05 & abs(AD_Healthyres$log2FoldChange)>0.585 )    
AD_Healthyres <- AD_Healthyres[order(AD_Healthyres$padj,decreasing = F),]
head(AD_Healthyres)
choose_gene <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585)
dim(choose_gene)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
MostSigGene <- read.csv("GE_MostSignificantGene.csv",header = T,row.names = 1)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myfpkm <- read.csv(file = "GE_RemoveBatch_FPKM.csv",header = T ,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
AD_Healthyres_heatmap <- subset(myfpkm, rownames(myfpkm) %in% rownames(choose_gene))
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AD_Healthyres_heatmap <- as.data.frame(t(AD_Healthyres_heatmap))
AD_Healthyres_heatmap$GSM_number <- rownames(AD_Healthyres_heatmap)
sum <- merge(AD_Healthyres_heatmap,BackgroundInformation,by = "GSM_number",all= FALSE)
sum[5,999:1010]
sum[5,1:5]
rownames(sum) <- sum$Sample
sum <- arrange(sum,rownames(sum))
data <- sum[-c(1,1000:1010)]
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){
  sd(x)/mean(x)*100
})
data_df <- data[order(data$cv, decreasing = T),1:354]
dim(data_df)
a <- apply(data_df,1,scale)
data_scale <- as.data.frame(t(apply(data_df,1,scale))) ##Z-score标准化
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.01
data_scale <- as.matrix(data_scale)
table((data_scale)>2)
table((data_scale)<(-2))
data_scale[data_scale>=2]=2
data_scale[data_scale<=-2]=-2
table(is.na(data_scale))
library(ComplexHeatmap)
library(circlize)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "GE_Heatmap.pdf",width =3,height = 4)
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
             heatmap_legend_param = list( 
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
rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sigGene <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",row.names = 1)
library(stringr)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
BackgroundInformationAPOE <- BackgroundInformationAPOE[BackgroundInformationAPOE$Brain_Region==3,]
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% rownames(sigGene))
sig_FPKM <- as.data.frame(t(sig_FPKM))
library(WGCNA)
gsg <- goodSamplesGenes(sig_FPKM, verbose=3)
gsg
BackgroundInformation <- subset(BackgroundInformationAPOE,BackgroundInformationAPOE$GSM_number %in% rownames(sig_FPKM))
rownames(BackgroundInformation) <- BackgroundInformation$GSM_number
colnames(BackgroundInformation)
BackgroundInformation <- BackgroundInformation[,-c(1:6,10:11)]
library(dplyr)
colnames(BackgroundInformation)
BackgroundInformation <- BackgroundInformation%>%mutate(Phenotypedit=case_when(BackgroundInformation$Phenotypedit=="noE4"~"0",
                                                                               BackgroundInformation$Phenotypedit=="E4carrier"~"1",
                                                                               BackgroundInformation$Phenotypedit=="E4/4"~"2"))
BackgroundInformation$Phenotypedit <- as.numeric(BackgroundInformation$Phenotypedit)
BackgroundInformation$Braak_Stage <- as.numeric(BackgroundInformation$Braak_Stage)
BackgroundInformation$Bank_Location <- as.numeric(BackgroundInformation$Bank_Location)
str(BackgroundInformation)
sampleTree <- hclust(dist(sig_FPKM), method="average")
sizeGrWindow(2, 2)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="",
     xlab="", cex.lab=1.5, 
     cex.axis=1.5, cex.main=2)
sampleTree2 <- hclust(dist(sig_FPKM), method="average")
traitColors <- numbers2colors(BackgroundInformation, signed=F)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=colnames(BackgroundInformation), 
                    main="Sample dendrogram and trait heatmap")
powers <- c(c(1:10), seq(from=12, to=50, by=2))
library("doParallel")
sft <- pickSoftThreshold(sig_FPKM, powerVector=powers, verbose=5, networkType="unsigned")
sft$powerEstimate
sizeGrWindow(2, 2)
par(mfrow=c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels=powers, cex=cex1, col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels=powers, cex=cex1, col="red")
sft$powerEstimate <- 6
cor <- WGCNA::cor
net <- blockwiseModules(sig_FPKM, power = sft$powerEstimate,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,randomSeed=47,
                        saveTOMFileBase = "ADPatientTOM20240910", 
                        verbose = 3)
table(net$colors)
cor<-stats::cor
sizeGrWindow(2, 2)
mergedColors <- labels2colors(net$colors)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "GE_WGCNA_Module_colors.pdf",width =5,height = 3,bg = "white")
p1 <- plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          cex.rowText = 0.9,cex.colorLabels = 0.9, cex.dendroLabels = 0.9)
print(p1)
dev.off()
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
save(MEs, moduleLabels, moduleColors, geneTree,
     file="GE_networkConstruction-auto.RData")
library(WGCNA)
options(stringsAsFactors=F)

allowWGCNAThreads()
nGenes = ncol(sig_FPKM)
nSamples = nrow(sig_FPKM)
a <- moduleEigengenes(sig_FPKM, moduleColors)
table(moduleColors)
MEs0 = moduleEigengenes(sig_FPKM, moduleColors)$eigengenes
identical(rownames(MEs0),rownames(BackgroundInformation))
moduleTraitCor = cor(MEs0, BackgroundInformation,use = 'p',method = "spearman")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(2,2)
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
APOE <- as.data.frame(BackgroundInformation$APOE)
colnames(APOE) = "APOE"
modNames <- substring(names(MEs0), 3)
identical(rownames(MEs0),rownames(sig_FPKM))
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) = paste("p.MM", modNames, sep="")
identical(rownames(sig_FPKM),rownames(BackgroundInformation))
geneTraitSignificance = as.data.frame(cor(sig_FPKM, APOE, use = "p",method = "spearman"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) = paste("GS.", names(APOE), sep="");
colnames(GSPvalue) = paste("p.GS.", names(APOE), sep="")
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
allLLIDs <- colnames(sig_FPKM)
module = "blue"
modGenes = (moduleColors==module)
table(modGenes)
modLLIDs = allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]
m$MMblue <- abs(m$MMblue)
m$GS.APOE <- abs(m$GS.APOE)
m1 <- m
m <- subset(m,m$p.MMblue<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$GS.APOE) > 0.3)
blue_APOE <- n$Row.names
write.csv(n,file = "GE_nblueAPOE.csv")
colnames(n)
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMblue) , y = abs(GS.APOE))) +
  geom_point(size=4,colour="blue") +
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
module = "turquoise"
modGenes = (moduleColors==module)
modLLIDs = allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]
m$MMturquoise <- abs(m$MMturquoise)
m$GS.APOE <- abs(m$GS.APOE)
m1 <- m
m <- subset(m,m$p.MMturquoise<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$GS.APOE) > 0.3)
colnames(n)
turquoise_APOE <- n$Row.names
write.csv(n,file = "GE_nturquoiseAPOE.csv")
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMturquoise) , y = abs(GS.APOE))) +
  geom_point(size=4,colour="turquoise") +
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
Braak_Stage <- as.data.frame(BackgroundInformation$`Braak Stage`)
colnames(Braak_Stage) = "Braak Stage"
modNames <- substring(names(MEs0), 3)
identical(rownames(MEs0),rownames(sig_FPKM))
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) = paste("p.MM", modNames, sep="")
identical(rownames(sig_FPKM),rownames(BackgroundInformation))
geneTraitSignificance = as.data.frame(cor(sig_FPKM, Braak_Stage, use = "p",method = "spearman"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) = paste("GS.", names(Braak_Stage), sep="");
colnames(GSPvalue) = paste("p.GS.", names(Braak_Stage), sep="")
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
module = "blue"
modGenes = (moduleColors==module)
modLLIDs = allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]
m$MMblue <- abs(m$MMblue)
m$`GS.Braak Stage` <- abs(m$`GS.Braak Stage`)
m1 <- m
m <- subset(m,m$p.MMblue<0.05 & m$`p.GS.Braak Stage` <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$`GS.Braak Stage`) > 0.3)
write.csv(n,file = "GE_nblueBraak_Stage.csv")
colnames(n)
blue_BraakStage <- n$Row.names
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMblue) , y = abs(`GS.Braak Stage`))) +
  geom_point(size=4,colour="blue") +
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
module = "turquoise"
modGenes = (moduleColors==module)
modLLIDs = allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]
m$MMturquoise <- abs(m$MMturquoise)
m$`GS.Braak Stage` <- abs(m$`GS.Braak Stage`)
m1 <- m
m <- subset(m,m$p.MMturquoise<0.05 & m$`p.GS.Braak Stage` <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$`GS.Braak Stage`) > 0.3)
colnames(n)
turquoise_BraakStage <- n$Row.names
write.csv(n,file = "GE_nturquoiseBraak_Stage.csv")
library(ggpubr)
library(ggrepel)
p1 <- ggplot(
  m1, aes(x =abs(MMturquoise) , y = abs(`GS.Braak Stage`))) +
  geom_point(size=4,colour="turquoise") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
turquoise_mark_gene <- append(turquoise_APOE,turquoise_BraakStage)
gene <- turquoise_mark_gene
gene <- gene[!duplicated(gene)]
ADvsHealthy_DEseq2 <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]
gene <- ADvsHealthy_DEseq2$X
turquoise_gene <- ADvsHealthy_DEseq2$X
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)
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
x = names(gene_group)[6]
y = names(gene_group)[-6]
table(gene_group$Severity)
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
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,12)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
x = names(gene_group)[6]
y = names(gene_group)[-6]
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
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,17)]
gene_group <- na.omit(gene_group)
colnames(gene_group)
gene_group <- gene_group%>%mutate(Phenotypedit=case_when(gene_group$Phenotypedit=="noE4"~"0",
                                                         gene_group$Phenotypedit=="E4carrier"~"1",
                                                         gene_group$Phenotypedit=="E4/4"~"2"))
gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
str(gene_group)
colnames(gene_group)
table(gene_group$Phenotypedit)
x = names(gene_group)[6]
y = names(gene_group)[-6]
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
blue_mark_gene <- append(blue_APOE,blue_BraakStage)
gene <- blue_mark_gene 
gene <- gene[!duplicated(gene)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ADvsHealthy_DEseq2 <- read.csv("GE_RemoveBatch_DEseq2_Significant.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]
gene <- ADvsHealthy_DEseq2$X
blue_gene <- ADvsHealthy_DEseq2$X
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)
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
colnames(gene_group)
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
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,12)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
str(gene_group)
table(gene_group$Braak_Stage)
colnames(gene_group)
x = names(gene_group)[6]
y = names(gene_group)[-6]
table(gene_group$Braak_Stage)
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
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,17)]
gene_group <- na.omit(gene_group)
colnames(gene_group)
gene_group <- gene_group%>%mutate(Phenotypedit=case_when(gene_group$Phenotypedit=="noE4"~"0",
                                                         gene_group$Phenotypedit=="E4carrier"~"1",
                                                         gene_group$Phenotypedit=="E4/4"~"2"))
gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
str(gene_group)
colnames(gene_group)
table(gene_group$Phenotypedit)
x = names(gene_group)[6]
y = names(gene_group)[-6]
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
gene <- c("SYT1","CHN1","SNAP25","VSNL1","ENC1","TNS1","SGK1","CPM","CLMN","PPFIBP2")
myfpkm <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(a, fromType = "ENSEMBL",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:11,16)]
gene_group <- gene_group%>%mutate(Severity=case_when(gene_group$Severity==0~"Healthy",
                                                     gene_group$Severity==1~"AD",
                                                     gene_group$Severity==2~"MCI"))
gene_group$Severity <- as.factor(gene_group$Severity)
b <- gene_group
colnames(b)[11] <- "Group"
library(dplyr)
b <- b%>%mutate(Group=case_when(b$Group=="Healthy"~"0",
                                b$Group=="AD"~"1"))
n <- dim(b)[1]
y <- b$Group
all <- b
set.seed(11)
require(caret)
folds <- createFolds(y,k=6)
library(ROCR)
library(magrittr) 
library(plyr)
auc_value<-as.numeric()
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] 
  fold_train <- all[-folds[[i]],] 
  model <- glm(as.numeric(fold_train$Group)~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- predict(model,type='response',newdata=fold_test)
  pred <- prediction(predictions = fold_predict, labels = fold_test$Group)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
mean(auc_value)
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
colnames(rocData)
rocData$Xmean <- rowMeans(rocData[,1:6])
rocData$Ymean <- rowMeans(rocData[,7:12])
rocData$Gene <- c(rep("GE", 60))
rocData <- rocData[,13:15]
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("Control", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
ROC <- rbind(rocData,inn)
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
