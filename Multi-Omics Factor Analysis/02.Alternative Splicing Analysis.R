rm(list = ls())
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
a <- load("E:/AD_Patient/ASanalysiaLeafcutter/testtemporal_lobeADpatient.Rdata")
a <- colnames(counts)
b <- substr(a,1,10)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(counts,file = "RawIntroRead.csv",quote = F)
exons_table <- as.data.frame(exons_table)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(exons_table,file = "exons_table.csv",quote = F)
intron <- introns[,1:7]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(intron,file = "Intron_To_ClusterNames.csv",quote = F)
rm(list = ls())
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
a <- as.data.frame(rowSums(counts))
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
counts <- sum
counts <- as.data.frame(t(counts))
condition <- factor(substring(colnames(counts), 1,1))
coldata <- data.frame(row.names =colnames(counts),condition) 
library(DESeq2)  
dds <- DESeqDataSetFromMatrix(countData = counts,colData = coldata,design = ~ condition) 
vst <- vst(dds, blind=T)
library(ggforce)
library(ggrepel)
plotPCA(vst,intgroup = "condition")
p1data <- plotPCA(vst,returnData = T,intgroup = "condition")
colnames(p1data)[5] <- "Sample"
pca_data <- merge(p1data,Temporal_Lobe,by="Sample")
pca_data$condition <- as.factor(pca_data$condition)
pca_data$Batch <- as.factor(pca_data$Batch)
p1 <- ggplot(data=pca_data,aes(x=PC1,y=PC2,color=condition,shape =Batch))+
  geom_point(size=3)+
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_RawCount_PCA.pdf", egg::set_panel_size(p1, width=unit(4.5, "in"), height=unit(4, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
df <- as.data.frame(t(counts))
Temporal_Lobe$Brain_Region <- as.factor(Temporal_Lobe$Brain_Region)
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
library(vegan)
a <- adonis2(df ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(6,7),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
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
                                   b$category=="Sex"~"Gender"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2, fill= category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+  
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("AS Raw Count") +
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
library(limma)
library(edgeR)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
counts <- counts[,order(colnames(counts))]
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("E:/AD_Patient/ASanalysiaLeafcutter")
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
Temporal_Lobe <- Temporal_Lobe[order(Temporal_Lobe$GSM_number),]
identical(colnames(counts),Temporal_Lobe$GSM_number)
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
counts <- as.data.frame(t(sum))
group <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(counts))) #
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design)=colnames(counts)
a <- as.data.frame(colSums(counts))
keep.exprs <- filterByExpr(counts, group=group,
                           min.count = 3, min.total.count = 30,
                           large.n = 10, min.prop = 0.3)
counts <- counts[keep.exprs,]
dge<-DGEList(counts)
library(edgeR)
y <- voom(dge, design, plot = T,normalize.method ="quantile")
myvoom <- edit(voom)
y <- myvoom(dge, design, plot = T,normalize.method ="quantile")
newData <- y$E
NormalizedCount <- as.data.frame(newData)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(NormalizedCount,file = "AS_VoomedCount.csv",quote = F)
df <- as.data.frame(t(newData))
library(gmodels)
pca.info <- prcomp(df)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
library(ggpubr)
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
a$Batch <- as.factor(a$Batch)
a$Severity <- as.factor(a$Severity)
a <- a[order(a$Sample),]
library(ggrepel)
library(ggplot2)
library(ggforce)
library(dplyr)
a <- a%>%mutate(Type=case_when(a$Type==0~"Healthy",
                               a$Type==1~"AD",
                               a$Type==2~"MCI"))
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",16,"%"),
       y=paste0("PCA2 ",10,"%"))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomCount_PCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
df <- as.data.frame(t(newData))
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
identical(rownames(df),Temporal_Lobe$Sample)
library(vegan)
a <- adonis2(df ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(7,6),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill =category ))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+ 
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("AS Voom Count") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomCount.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
design <- model.matrix(~Severity + Age , data=Temporal_Lobe)
treatment.design <- design[,1:2]
batch.design <- design[,-(1:2)]
corrected_count <- removeBatchEffect(newData,batch = Temporal_Lobe$Batch,
                                     group = Temporal_Lobe$Severity,
                                     covariates=batch.design,
                                     design = treatment.design)
RemoveData <- as.data.frame(t(corrected_count))
RemoveData <- as.data.frame(t(RemoveData))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(RemoveData,file = "AS_RemoveBatch_VoomCount.csv",quote = F)
RemoveData <- as.data.frame(t(RemoveData))
identical(rownames(RemoveData),Temporal_Lobe$Sample)
library(gmodels)
pca.info <- prcomp(RemoveData)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
library(ggpubr)
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
a$Batch <- as.factor(a$Batch)
a$Severity <- as.factor(a$Severity)
a <- a[order(a$Sample),]
library(ggrepel)
library(ggplot2)
library(ggforce)
library(dplyr)
a <- a%>%mutate(Type=case_when(a$Type==0~"Healthy",
                               a$Type==1~"AD",
                               a$Type==2~"MCI"))
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
        legend.direction = 'vertical', legend.position ="right")+
  labs(x=paste0("PCA1 ",12,"%"),
       y=paste0("PCA2 ",9,"%"))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_Voom_RemoveBatchPCA.pdf", egg::set_panel_size(p1, width=unit(5, "in"), height=unit(5, "in")), 
       width = 8, height = 7, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
Temporal_Lobe$Braak_Stage <- as.factor(Temporal_Lobe$Braak_Stage)
Temporal_Lobe$Bank_Location <- as.factor(Temporal_Lobe$Bank_Location)
Temporal_Lobe$Batch <- as.factor(Temporal_Lobe$Batch)
Temporal_Lobe$Severity <- as.factor(Temporal_Lobe$Severity)
identical(rownames(RemoveData),Temporal_Lobe$Sample)
library(vegan)
a <- adonis2(RemoveData ~ Braak_Stage+Age+Bank_Location+Batch+Severity, 
             data = Temporal_Lobe,
             permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(7,6),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
library(tidyverse)
b <- b%>%mutate(Pr=case_when(b$Pr>0.1~"",
                             0.05<b$Pr & b$Pr <= 0.1 ~ ".",
                             0.01 < b$Pr & b$Pr <= 0.05 ~"*",
                             0.001 < b$Pr & b$Pr <= 0.01 ~ "**",
                             0.0001 <= b$Pr & b$Pr <= 0.001 ~"***"))
b <- b%>%mutate(category=case_when(b$category=="Bank_Location"~"Bank Location",
                                   b$category=="Braak_Stage"~"Braak Stage",
                                   b$category=="Severity"~"Severity",
                                   b$category=="Batch"~"Batch",
                                   b$category=="Age"~"Age"))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
library(ggsci)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill = category))+
  geom_bar(stat="identity")+
  scale_fill_npg()+
  geom_text(aes(label=Pr),vjust=-0.1,angle = 0, size=5)+ 
  xlab("category") +ylab("R2")+
  theme_bw()+
  ggtitle("AS Voom RemoveBatch") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_VoomRemoveBatch.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 9, height = 7, units = 'in', dpi = 600)
RemoveData <- as.data.frame(t(RemoveData))
RemoveData <- RemoveData[,order(colnames(RemoveData))]
RemoveData$groupHealthymean <- apply(RemoveData[,1:97],1,mean)
RemoveData$groupADmean <- apply(RemoveData[,98:354],1,mean)
RemoveData$Group_diff <- RemoveData$groupADmean-RemoveData$groupHealthymean
RemoveData$log2FoldChange <- log2(RemoveData$groupADmean/RemoveData$groupHealthymean)
RemoveData$Pvalue <- apply(RemoveData,1,function(x) t.test(x[1:97],x[98:354])$p.value)
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
RemoveData$Position <- rownames(RemoveData)
name <- read.table("ASdiff_n344.txt",header = T,row.names = 1)
name$Position <- paste(rownames(name),name$Cluster,sep = ":")
AS_anno <- name[,c(8,10)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(AS_anno,file = "AS_anno.csv",quote = F)
a <- merge(AS_anno,RemoveData,by="Position")
a<- a[,-2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(a,file = "AS_VoomRemoveBatch_DEG.csv",quote = F,row.names = F)
a <- merge(AS_anno,RemoveData,by="Position")
a$Genename <- ifelse(a$Genename == "REMOVE", "no_gene", a$Genename)
sum <- a
sum$change <- ifelse(sum$Padjust < 0.05 & abs(sum$log2FoldChange) >= 0.1,
                     ifelse(sum$log2FoldChange > 0.1 ,'Up','Down'),'Stable')
write.table(sum,file = "ASdiff_n369.txt",quote = F,row.names = F,sep = "\t")
AS_diff <- read.table("ASdiff_n369.txt",header = T,row.names = 1,check.names = F)
sum <- AS_diff
library(dplyr)
up <- subset(sum, sum$change == 'Up')
down <- subset(sum, sum$change == 'Down')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:5, ]
down <- down[order(down$Padjust), ][1:5, ]
a <- rbind(up, down)
library(ggplot2)
library(ggrepel)
p <- ggplot(
  sum, aes(x = log2FoldChange, y = -log10(Padjust),colour=change)) +
  geom_point(aes(color = change), size=3) +
  scale_color_manual(values = c("#008080", "gray", "firebrick3")) +
  geom_text_repel(data = a, aes(x = log2FoldChange, y = -log10(Padjust), label = Genename),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=100)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="#666666",lwd=0.5) +
  geom_hline(yintercept = -(log10(0.05)),lty=4,col="#666666",lwd=0.5) +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_Volcano.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(5, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
name <- read.table("ASdiff_n369.txt",header = T,check.names = F)
sum <- name
up <- subset(sum, sum$change == 'Up')
down <- subset(sum, sum$change == 'Down')
a <- rbind(up, down)
up <- up[order(up$Padjust), ][1:5, ]
down <- down[order(down$Padjust), ][1:5, ]
MostSigGene <- rbind(up, down)
AS_Heatmap <- sum[sum$change!="Stable",]
rownames(AS_Heatmap) <- AS_Heatmap$Position
AS_Heatmap <- AS_Heatmap[order(AS_Heatmap$Padjust),]
colnames(AS_Heatmap)
AS_Heatmap <- AS_Heatmap[,3:356]
sum <- as.data.frame(t(AS_Heatmap))
sum <- arrange(sum,rownames(sum))
data <- sum
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
table((data_scale)>1)
table((data_scale)<(-1))
data_scale[data_scale>=1]=1
data_scale[data_scale<=-1]=-1
library(ComplexHeatmap)
library(circlize)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
pdf(file = "AS_Heatmap3.pdf",width =3,height = 4)
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
             heatmap_legend_param = list(
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F)
a <- a[a$change!="Stable",]
a <- a[a$Genename!="no_gene",]
a <- a[!grepl(",", a$Genename),]
sum <- a
gene <- unique(a$Genename)
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
write.csv(BPplot_genelist,file = "AS_BPplot_genelist.csv",row.names = F,quote = F)
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
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#0D6CA6","#099963", "#911F27")
library(ggpubr)
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_GO.pdf", egg::set_panel_size(p, width=unit(16, "in"), height=unit(5, "in")), 
       width = 20, height = 15, units = 'in', dpi = 600)
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
table1 <- sum[,c(2,363)]
table2 <- kegg
calculate_up_down <- function(genes, gene_table) {
  gene_list <- unlist(strsplit(genes, "/"))
  down_genes <- sum(subset(gene_table,gene_table$Genename %in% gene_list)[,2] == "Down")
  up_genes <- sum(subset(gene_table,gene_table$Genename %in% gene_list)[,2] == "Up")
  return(c(down_genes, up_genes))
}
result <- t(apply(table2, 1, function(row) {
  calculate_up_down(row['geneID'], table1)
}))
table2$down <- -abs(result[, 1])
table2$up <- result[, 2]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
write.csv(table2,file = "AS_KEGG.csv",row.names = F)
library(reshape2)
library(knitr)
KEGGTerm <- table2[1:15,c(3,4,12,13)]
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
        legend.direction = 'vertical', legend.position ="right")
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
ggsave("AS_KEGG.pdf", egg::set_panel_size(p1, width=unit(2.5, "in"), height=unit(5, "in")), 
       width = 9, height = 6, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F,row.names = 1)
b <- subset(a,a$change !="Stable")
b <- b[order(b$Padjust),]
VoomRemoveBatch_DEG <- b[1:10,2:355]
VoomRemoveBatch_DEG <- as.data.frame(t(VoomRemoveBatch_DEG))
VoomRemoveBatch_DEG$Group <- rownames(VoomRemoveBatch_DEG)
b <- VoomRemoveBatch_DEG
b$Group <- substr(b$Group,1,1)
n <- dim(b)[1]
y <- b$Group
all <- b
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
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData$Xmean <- rowMeans(rocData[,1:6])
rocData$Ymean <- rowMeans(rocData[,7:12])
rocData$Gene <- c(rep("AS", 60))
rocData <- rocData[,13:15]
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("Control", 60))
colnames(inn) <- c("Xmean","Ymean","Gene")
ROC <- rbind(rocData,inn)
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
ggsave("AS_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
a <- read.table("ASdiff_n369.txt",header = T,check.names = F,row.names = 1)
b <- subset(a,a$change !="Stable")
b <- b[order(b$Padjust),]
b <- b[1:10,]
b <- b[order(rownames(b)),]
VoomRemoveBatch_DEG <- b[,1:355]
rownames(VoomRemoveBatch_DEG) <- VoomRemoveBatch_DEG$Genename
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG[,-1]
VoomRemoveBatch_DEG <- as.data.frame(t(VoomRemoveBatch_DEG))
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG[,order(colnames(VoomRemoveBatch_DEG))]
VoomRemoveBatch_DEG$Group <- substr(rownames(VoomRemoveBatch_DEG),1,1)
library(tidyverse)
VoomRemoveBatch_DEG <- VoomRemoveBatch_DEG%>%mutate(Group=case_when(VoomRemoveBatch_DEG$Group==0~"Healthy",
                                                                    VoomRemoveBatch_DEG$Group==1~"AD",
                                                                    VoomRemoveBatch_DEG$Group==2~"MCI"))
VoomRemoveBatch_DEG$Group <- as.factor(VoomRemoveBatch_DEG$Group)
x = names(VoomRemoveBatch_DEG)[11]
y = names(VoomRemoveBatch_DEG)[-11]
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

