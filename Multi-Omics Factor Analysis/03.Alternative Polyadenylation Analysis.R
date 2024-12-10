rm(list = ls())
library(tidyr)
library(gmodels)
library(ggpubr)
library(vegan)
library(tidyverse)
library(ggsci)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
DaparsSample <- read.csv("DaparsSample.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
colnames(DaparsSample)[2:3] <- c("GSM_number","sample")
pheno <- merge(DaparsSample,BackgroundInformation,by="GSM_number")
pheno$sample <- gsub("_PDUI$", "", pheno$sample)
pheno <- subset(pheno,pheno$Brain_Region=="3")
pheno$Batch <- as.factor(pheno$Batch)
pheno$Braak_Stage <- as.factor(pheno$Braak_Stage)
pheno$Bank_Location <- as.factor(pheno$Bank_Location)
pheno$Severity <- as.factor(pheno$Severity)
setwd("E:/AD_Patient/Dapars2/Dapars0.9")
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
test <- Healthy_AD[,c(1,1067:1072)]
a <- Healthy_AD[,-c(1:4,1067:1072)]
a <- a[,seq(0,ncol(a),3)]
Healthy_AD_PDUI <- cbind(Healthy_AD[,1:4],a)
long<-separate(Healthy_AD_PDUI,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
ID_Transfer <- long[,1:2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(ID_Transfer,file = "APA_ID_Transfer.csv",quote = F,row.names = F)
Healthy_AD_PDUI <- long
Healthy_AD_PDUI <- Healthy_AD_PDUI[,-c(2:7)]
rownames(Healthy_AD_PDUI) <- Healthy_AD_PDUI$ENSEMBLTRANS
Healthy_AD_PDUI <- Healthy_AD_PDUI[,-1]
colnames(Healthy_AD_PDUI) <- gsub("_PDUI", "", colnames(Healthy_AD_PDUI))
Healthy_AD_PDUI <- as.data.frame(t(Healthy_AD_PDUI))
pca.info <- fast.prcomp(Healthy_AD_PDUI)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
b <- merge(pheno,pca.data,by="sample")
b$Batch <- as.factor(b$Batch)
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",47,"%"),
       y=paste0("PCA2 ",19,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_PDUI_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)
setwd("E:/AD_Patient/Dapars2/Dapars0.9")
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
name <- Healthy_AD[,c(1:4)]
num <- Healthy_AD[,-c(1:4,1067:1072)]
long<-separate(name,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- cbind(long[,1],num)
colnames(Healthy_AD)[1] <- "ENSEMBLTRANS"
rownames(Healthy_AD) <- Healthy_AD$ENSEMBLTRANS
Healthy_AD <- Healthy_AD[,-1]
Healthy_AD_long <- Healthy_AD[,seq(1,ncol(Healthy_AD),3)]
colnames(Healthy_AD_long) <- gsub("_long_exp", "", colnames(Healthy_AD_long))
Healthy_AD_short <- Healthy_AD[,seq(2,ncol(Healthy_AD),3)]
colnames(Healthy_AD_short) <- gsub("_short_exp", "", colnames(Healthy_AD_short))
Healthy_AD_long <- as.data.frame(t(Healthy_AD_long))
pca.info <- fast.prcomp(Healthy_AD_long)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
b <- merge(pheno,pca.data,by="sample")
b$Batch <- as.factor(b$Batch)
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",46,"%"),
       y=paste0("PCA2 ",18,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Long_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_long <- Healthy_AD_long[order(rownames(Healthy_AD_long)),]
pheno <- pheno[order(pheno$sample),]
identical(rownames(Healthy_AD_long),pheno$sample)
Healthy_AD_long_avo <- adonis2(as.matrix(Healthy_AD_long) ~ Severity+Braak_Stage+Age+Bank_Location+Batch, 
                               data = pheno,
                               permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(Healthy_AD_long_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("APA Raw Long") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Long.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
design <- model.matrix( ~ Severity + Age , data=pheno)
treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]
Healthy_AD_long <- as.data.frame(t(Healthy_AD_long))
Healthy_AD_long <- Healthy_AD_long[,order(colnames(Healthy_AD_long))]
pheno <- pheno[order(pheno$sample),]
identical(colnames(Healthy_AD_long),pheno$sample)
Healthy_AD_long_Correct <- removeBatchEffect(as.matrix(Healthy_AD_long),
                                             batch = pheno$Batch,
                                             group = pheno$Severity,
                                             covariates= batch.design,
                                             design = treatment.design)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(Healthy_AD_long_Correct,file = "APA_Long_RemoveBatch_Count.csv",quote = F)
Healthy_AD_long_Remove <- as.data.frame(t(Healthy_AD_long_Correct))
pca.info <- fast.prcomp(Healthy_AD_long_Remove)
b <- merge(pheno,pca.data,by="sample")
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",42,"%"),
       y=paste0("PCA2 ",24,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Long_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_long_Remove <- t(Healthy_AD_long_Correct)
Healthy_AD_long_Remove_avo <- adonis2(Healthy_AD_long_Remove ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
                                      data = pheno,
                                      permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(Healthy_AD_long_Remove_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+ 
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("APA RemoveBatch Long") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Long.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_short <- as.data.frame(t(Healthy_AD_short))
Healthy_AD_short <- Healthy_AD_short[order(rownames(Healthy_AD_short)),]
pheno <- pheno[order(pheno$sample),]
identical(rownames(Healthy_AD_short),pheno$sample)
Healthy_AD_short_avo <- adonis2(as.matrix(Healthy_AD_short) ~ Severity+Braak_Stage+Age+Bank_Location+Batch, 
                                data = pheno,
                                permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(Healthy_AD_short_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("APA Raw Short") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Short.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_short <- as.data.frame(Healthy_AD_short)
pca.info <- fast.prcomp(Healthy_AD_short)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
b <- merge(pheno,pca.data,by="sample")
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",57,"%"),
       y=paste0("PCA2 ",11,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_Raw_Short_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)
design <- model.matrix( ~ Severity + Age , data=pheno)
treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]
Healthy_AD_short <- as.data.frame(t(Healthy_AD_short))
Healthy_AD_short <- Healthy_AD_short[,order(colnames(Healthy_AD_short))]
pheno <- pheno[order(pheno$sample),]
identical(colnames(Healthy_AD_short),pheno$sample)
Healthy_AD_short_Correct <- removeBatchEffect(as.matrix(Healthy_AD_short),
                                              batch = pheno$Batch,
                                              group = pheno$Severity,
                                              covariates= batch.design,
                                              design = treatment.design)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(Healthy_AD_short_Correct,file = "APA_Short_RemoveBatch_Count.csv",quote = F)
Healthy_AD_short_Remove <- as.data.frame(t(Healthy_AD_short_Correct))
pca.info <- fast.prcomp(Healthy_AD_short_Remove)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
b <- merge(pheno,pca.data,by="sample")
b$Batch <- as.factor(b$Batch)
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",44,"%"),
       y=paste0("PCA2 ",14,"%"))+
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Short_PCA.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")),
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_short_Remove <- t(Healthy_AD_short_Correct)
Healthy_AD_short_Remove_avo <- adonis2(Healthy_AD_short_Remove ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
                                       data = pheno,
                                       permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(Healthy_AD_short_Remove_avo)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- na.omit(b)
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("APA RemoveBatch Short") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_RemoveBatch_Short.pdf", egg::set_panel_size(p, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
Healthy_AD_long_Correct[Healthy_AD_long_Correct<0]=NA
Healthy_AD_short_Correct[Healthy_AD_short_Correct<0]=NA
Healthy_AD_long_Correct <- Healthy_AD_long_Correct[order(rownames(Healthy_AD_long_Correct)),order(colnames(Healthy_AD_long_Correct))]
Healthy_AD_short_Correct <- Healthy_AD_short_Correct[order(rownames(Healthy_AD_short_Correct)),order(colnames(Healthy_AD_short_Correct))]
identical(colnames(Healthy_AD_long_Correct),colnames(Healthy_AD_short_Correct))
identical(rownames(Healthy_AD_long_Correct),rownames(Healthy_AD_short_Correct))
if (!all(dim(Healthy_AD_long_Correct) == dim(Healthy_AD_short_Correct))) {
  stop("两个矩阵的维度必须相同")}
Healthy_AD_remove_PDUI <- Healthy_AD_long_Correct / (Healthy_AD_long_Correct + Healthy_AD_short_Correct)
RemoveData <- Healthy_AD_remove_PDUI
RemoveData <- RemoveData[,order(colnames(RemoveData))]
RemoveData <- as.data.frame(RemoveData)
RemoveData$groupAmean <- rowMeans(RemoveData[,1:97], na.rm = TRUE)
RemoveData$groupBmean <- rowMeans(RemoveData[,98:354], na.rm = TRUE)
RemoveData$Mean_diff <- RemoveData$groupBmean-RemoveData$groupAmean
RemoveData$Pvalue <- apply(RemoveData,1,function(x) wilcox.test(na.omit(x[1:97]),na.omit(x[98:354]))$p.value)
RemoveData <- RemoveData[complete.cases(RemoveData$Pvalue), ]
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
RemoveData$change <- ifelse(RemoveData$Padjust < 0.05 & abs(RemoveData$Mean_diff) >= 0.1, 
                            ifelse(RemoveData$Mean_diff > 0.1 ,'Long','Short'),'Stable')
test <- RemoveData[,355:360]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(RemoveData,file = "APA_RemoveBatch_DEG.csv",quote = F)
sum <- test
sum$ENSEMBLTRANS <- rownames(sum)
com <- merge(ID_Transfer,sum,by="ENSEMBLTRANS")
sum <- com
up <- subset(sum, sum$change == 'Long')
down <- subset(sum, sum$change == 'Short')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:1, ]
down <- down[order(down$Padjust), ][1:5, ]
a <- rbind(up, down)
p <- ggplot(
  sum, aes(x = Mean_diff, y = -log10(Padjust),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080","firebrick3","gray")) +
  geom_text_repel(data = a, aes(x = Mean_diff, y = -log10(Padjust), label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_volcano.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)
a <- sum
a <- a[a$change!="Stable",]
gene <- unique(a$GeneName)
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(gene, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
GOgene <- as.character(df1$ENTREZID)
GOgene <- na.omit(GOgene)
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
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#0D6CA6","#099963", "#911F27")#设定颜色
library(ggpubr)
library(stringr)
p <- ggdotchart(go_enrich_df, x = "type_order", y = "GeneNumber",
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_GO.pdf", egg::set_panel_size(p, width=unit(16, "in"), height=unit(5, "in")), 
       width = 25, height = 20, units = 'in', dpi = 600)
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
table1 <- sum[,c(1,2,8)]
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  down_genes <- sum(subset(gene_table,gene_table$GeneName %in% gene_list)[,3] == "Short")
  up_genes <- sum(subset(gene_table,gene_table$GeneName %in% gene_list)[,3] == "Long")
  return(c(down_genes, up_genes))
}
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))
table2$Short <- -abs(result[, 1])
table2$Long <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "APA_KEGG.csv",row.names = F)
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:8,c(3,4,12,13)]
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(2.5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
b <- subset(a,a$change !="Stable")
b <- b[order(b$Padjust),]
b <- b[1:10,]
head(b[,355:360])
b <- b[,1:354]
b <- b[order(rownames(b)),]
b <- as.data.frame(t(b))
b <- b[,order(colnames(b))]
APA_ID_Transfer <- read.csv("APA_ID_Transfer.csv",header = T)
APA_ID_Transfer <- subset(APA_ID_Transfer,APA_ID_Transfer$ENSEMBLTRANS %in% colnames(b))
APA_ID_Transfer <- APA_ID_Transfer[order(APA_ID_Transfer$ENSEMBLTRANS),]
identical(colnames(b),APA_ID_Transfer$ENSEMBLTRANS)
b$Group <- substr(rownames(b),1,1)
library(dplyr)
b <- b%>%mutate(Group=case_when(b$Group=="A"~"0",
                                b$Group=="B"~"1"))
b$Group <- as.numeric(b$Group)
n <- dim(b)[1]
y <- b$Group
all <- b
folds <- createFolds(y,k=6)
library(ROCR)
library(magrittr)
library(plyr)
auc_value<-as.numeric()
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] 
  fold_test <- fold_test[rowSums(fold_test[1:10],na.rm =TRUE)!=0,]
  fold_train <- all[-folds[[i]],] 
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
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData <- rocData[1:44,]
rocData$Xmean <- rowMeans(rocData[,1:6],na.rm = T)
rocData$Ymean <- rowMeans(rocData[,7:12],na.rm = T)
rocData$Gene <- c(rep("APA", 44))
rocData <- rocData[,13:15]
a <- seq(0,1,1/43)
inn <- data.frame(a,a,rep("Control", 44))
colnames(inn) <- c("Xmean","Ymean","Gene")
ROC <- rbind(rocData,inn)
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("APA_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
b <- subset(a,a$change !="Stable")
b <- b[order(b$Padjust),]
b <- b[1:10,1:354]
b <- b[order(rownames(b)),]
b <- as.data.frame(t(b))
b <- b[,order(colnames(b))]
b$Group <- substr(rownames(b),1,1)
library(tidyverse)
b <- b%>%mutate(Group=case_when(b$Group=="A"~"Healthy",
                                b$Group=="B"~"AD"))
b$Group <- as.factor(b$Group)
x = names(b)[11]
y = names(b)[-11]
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sum <- read.csv("APA_RemoveBatch_DEG.csv",header = T,row.names = 1)
APA_Heatmap <- sum[sum$change!="Stable",]
APA_Heatmap <- APA_Heatmap[,1:354]
APA_Heatmap <- APA_Heatmap[,order(colnames(APA_Heatmap))]
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
APA_ID_Transfer <- read.csv("APA_ID_Transfer.csv",header = T)
APA_ID_Transfer <- APA_ID_Transfer[APA_ID_Transfer$ENSEMBLTRANS %in% rownames(a),]
sum <- as.data.frame(data_filled)
sum <- arrange(sum,rownames(sum))
data <- sum
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){
  sd(x)/mean(x)*100
})
data_df <- data[order(data$cv, decreasing = T),1:354]
a <- apply(data_df,1,scale)
data_scale <- as.data.frame(t(apply(data_df,1,scale))) 
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.001
data_scale <- as.matrix(data_scale)
table((data_scale)>1)
table((data_scale)<(-1))
data_scale[data_scale>=1]=1
data_scale[data_scale<=(-1)]=0
library(ComplexHeatmap)
library(circlize)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "APA_Heatmap.pdf",width =4,height = 3)
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
             heatmap_legend_param = list(
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

