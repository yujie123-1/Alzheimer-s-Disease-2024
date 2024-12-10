rm(list = ls())
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sumGE <- read.csv("GE_RemoveBatch_DEseq2.csv",header = T,row.names = 1)
sumGE <- sumGE[order(sumGE$padj),]
sumGE <- sumGE[1:5000,]
FPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(FPKM))
rownames(FPKM) <- a
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
gene <- rownames(sumGE)
sumGene <- subset(FPKM,rownames(FPKM) %in% gene)
sumGene <- as.data.frame(t(sumGene))
sumGene$GSM_number <- rownames(sumGene)
setwd("C:/Users/yujie/Desktop/datacollect")
bg <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
sum <- merge(bg,sumGene,by="GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:12)]
sumGene <- as.data.frame(t(sum))
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AS <- read.csv("AS_VoomRemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
AS <- AS[order(AS$Padjust),]
AS <- AS[1:5000,]
AS <- AS[,-c(355:360)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
APA <- read.csv("APA_RemoveBatch_DEG.csv",header = T,check.names = F,row.names = 1)
APA <- APA[order(APA$Pvalue),]
APA <- subset(APA,APA$Pvalue<0.05)
APA <- APA[,-c(355:360)]
DaparsSample <- read.csv("DaparsSample.csv",header = T)
colnames(DaparsSample) <- c("Group","GSM_number","sample")
DaparsSample$sample <- gsub("_PDUI","",DaparsSample$sample)
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformationAPOE <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
bg <- merge(DaparsSample,BackgroundInformationAPOE,by="GSM_number")
APA <- as.data.frame(t(APA))
APA$sample <- rownames(APA)
APA2 <- merge(bg,APA,by="sample") 
rownames(APA2) <- APA2$Sample
APA <- APA2[,-c(1:14)]
APA <- as.data.frame(t(APA))
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
sumGene <- sumGene[,order(colnames(sumGene))]
AS <- AS[,order(colnames(AS))]
APA <- APA[,order(colnames(APA))]
identical(colnames(sumGene),colnames(AS))
identical(colnames(AS),colnames(APA))
test <- list(Transcriptome=as.matrix(sumGene),AS=as.matrix(AS),APA=as.matrix(APA))
library(data.table)
library(MOFA2)
MOFAobject <- create_mofa(test)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 20
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts)
outfile = file.path("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2","MOFA2Results_2")
MOFAobject.trained <- run_mofa(MOFAobject, outfile,use_basilisk = TRUE)
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)
p1 <- plot_data_overview(model)
library(ggplot2)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p1,filename = "MOFA2_data_overview.pdf",width = 4,height = 2,bg = "white")
Nsamples = sum(model@dimensions$N)
setwd("C:/Users/yujie/Desktop/datacollect")
bg <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
bg <- subset(bg,bg$Brain_Region=="3")
bg <- bg[,-c(1:5,11)]
library(tidyverse)
bg <- bg%>%mutate(Phenotypedit=case_when(bg$Phenotypedit=="noE4"~"0",
                                         bg$Phenotypedit=="E4carrier"~"1",
                                         bg$Phenotypedit=="E4/4"~"2"))
bg <- bg[order(bg$Sample),]
a <- substr(colnames(as.data.frame(model@data$Transcriptome)),8,15)
colnames(bg)[5] <- "sample"
samples_metadata(model) <- bg
factor <- as.data.frame(model@cache$variance_explained$r2_per_factor[[1]])
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(factor,file = "MOFA2_Factor_Distribution.csv",quote = F)
p2 <- plot_variance_explained(model, x="view", y="factor",factors = 1:5)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p2,filename = "MOFA2_Factor_Explained.pdf",width = 5,height = 2,bg = "white")
p4 <- plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave(p4,filename = "MOFA2_Total_Explained.pdf",width = 4,height = 2,bg = "white")
SampleFactor <- as.data.frame(model@expectations[["Z"]][["group1"]])
SampleFactor$sample <- rownames(SampleFactor)
SampleFactor$Group <- rep("group1",354)
Sample <- merge(bg,SampleFactor,by="sample")
Sample <- na.omit(Sample)
Sample$Braak_Stage <- as.factor(Sample$Braak_Stage)
Sample$Age <- as.numeric(Sample$Age)
Sample$Severity <- as.factor(Sample$Severity)
library(ggsignif)
library(ggpubr)
library(glue)
p1 <- ggplot(Sample,aes(x=Severity,y=Factor1))+
  geom_signif(comparisons = list(c("0", "1")), 
              tip_length = 0.02,
              margin_top = 0.01,size = 0.8,textsize = 5,
              map_signif_level=TRUE)+
  geom_violin()+
  geom_jitter(aes(color=Severity),width = 0.2,size=3)+
  scale_color_manual(name = 'Severity',
                     values = c("#263E88","firebrick3"))+
  theme_bw()+
  ggtitle("MOFA2 SampleFactor1")+
  labs(x="",y=expression('Factor1'))+ 
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_SampleFactor1.pdf", egg::set_panel_size(p1, width=unit(3.5, "in"), height=unit(4, "in")), 
       width = 6, height = 5, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
SampleFactor <- as.data.frame(model@expectations[["Z"]][["group1"]])
SampleFactor$Sample <- rownames(SampleFactor)
SampleFactor <- SampleFactor[,c(1,16)]
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
all <- merge(SampleFactor,BackgroundInformation,by="Sample")
all <- all[,c(2,8)]
library(tidyverse)
all$Severity <- as.numeric(all$Severity)
n <- dim(all)[1]
y <- all$Severity
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
  model <- glm(fold_train$Severity~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- stats::predict(model,type='response',newdata=fold_test)
  pred <- prediction(fold_predict,fold_test$Severity)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData$Xmean <- apply(rocData[,1:6],1,mean)
rocData$Ymean <- apply(rocData[,7:12],1,mean)
rocData <- rocData[,13:14]
rocData$Gene <- c(rep("all", 60))
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
  ggtitle("MOFA2 Factor ROC AUC=0.733") +
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_FactorROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)
GEweights <- as.data.frame(model@expectations[["W"]][["Transcriptome"]])
GEweights$GeneName <- rownames(GEweights)
GEweights_Pos <- GEweights[GEweights$Factor1>=0,]
GEweights_Pos$Weight <- (GEweights_Pos$Factor1-min(abs(GEweights$Factor1)))/(max(abs(GEweights$Factor1))-min(abs(GEweights$Factor1)))
GEweights_Pos <- GEweights_Pos[order(abs(GEweights_Pos$Factor1),decreasing = T),]
GEweights_Neg <- GEweights[GEweights$Factor1<0,]
GEweights_Neg$Weight <- (GEweights_Neg$Factor1-min(abs(GEweights$Factor1)))/(max(abs(GEweights$Factor1))-min(abs(GEweights$Factor1)))
GEweights_Neg <- GEweights_Neg[order(abs(GEweights_Neg$Factor1),decreasing =T),]
a <- rbind(GEweights_Pos,GEweights_Neg)
KEGGa <- a
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$GeneName
topGE <- a
GEweights <- rbind(GEweights_Pos,GEweights_Neg)
GEweights <- GEweights[order(GEweights$Factor1,decreasing = F),]
GEweights$Rank <- seq(1,5000,1)
GEweights_Top <- subset(GEweights,GEweights$GeneName %in% a)
GEweights$Color <- ifelse(GEweights$GeneName %in% a, "#196F6E","grey")
library(ggrepel)
p <- ggplot(
  GEweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#196F6E", "gray")) +
  geom_text_repel(data = GEweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 2,box.padding = unit(0.2, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=2000)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 Transcriptome Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
GEweights_Top$abs <- abs(GEweights_Top$Weight)
library(ggpubr)
p <- ggdotchart(GEweights_Top, x = "GeneName", y = "Weight",
                color = "#196F6E",                                
                sorting = "ascending",                      
                add = "segments",                             
                add.params = list(color = "#196F6E", size = 1),
                rotate = TRUE,                               
                dot.size = "abs",                                
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               
                ggtheme = theme_pubr(),                     
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
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
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
write.csv(kegg,file = "MOFA2_GE_KEGG.csv",row.names = F)
kegg <- kegg[1:6,]
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("GE",6)))
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
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_GE_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)
ASweights <- as.data.frame(model@expectations[["W"]][["AS"]])
ASweights$GeneName <- rownames(ASweights)
library(tidyverse)
ASweights <- separate(ASweights,GeneName,into = c("Chr","Start","End","Cluster"),sep = "([:])")
ASweights$Position <- rownames(ASweights)
ASweights <- ASweights[,c(1,20)]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
m <- read.table("ASdiff_n369.txt",fill = TRUE,header = T,check.names = F)
m <- m[,c(1,2)]
sum <- merge(ASweights,m,by="Position")
ASweights <- sum
colnames(ASweights) <- c("Position","Factor1","GeneName")
ASweights_Pos <- ASweights[ASweights$Factor1>=0,]
ASweights_Pos$Weight <- (ASweights_Pos$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Pos <- ASweights_Pos[order(abs(ASweights_Pos$Factor1),decreasing = T),]
ASweights_Neg <- ASweights[ASweights$Factor1<0,]
ASweights_Neg$Weight <- (ASweights_Neg$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Neg <- ASweights_Neg[order(abs(ASweights_Neg$Factor1),decreasing =T),]
a <- rbind(ASweights_Pos,ASweights_Neg)
KEGGa <- a
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$Position
ASweights <- rbind(ASweights_Pos,ASweights_Neg)
ASweights <- ASweights[order(ASweights$Factor1,decreasing = F),]
ASweights$Rank <- seq(1,5000,1)
ASweights_Top <- subset(ASweights,ASweights$Position %in% a)
topAS <- ASweights_Top$GeneName
ASweights_Top <- ASweights_Top[order(ASweights_Top$Factor1,decreasing = T),]
ASweights$Color <- ifelse(ASweights$Position %in% a, "#940514","grey")
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
library(ggrepel)
p <- ggplot(
  ASweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#940514", "gray")) +
  geom_text_repel(data = ASweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=2000)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 AS Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
ASweights_Top$abs <- abs(ASweights_Top$Weight)
library(ggpubr)
p <- ggdotchart(ASweights_Top, x = "GeneName", y = "Weight",
                color = "#940514",                              
                sorting = "ascending",                      
                add = "segments",                           
                add.params = list(color = "#940514", size = 1),
                rotate = TRUE,                              
                dot.size = "abs",                               
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),               
                ggtheme = theme_pubr(),                       
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightASGene <- KEGGa$GeneName
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",
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
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_AS_KEGG.csv",row.names = F)
kegg <- kegg[1:6,]
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("AS",6)))
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
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_AS_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 4, units = 'in', dpi = 600)
APAweights <- as.data.frame(model@expectations[["W"]][["APA"]])
APAweights$ENSEMBLTRANS <- rownames(APAweights)
APAweights <- APAweights[,c(1,16)]
setwd("E:/AD_Patient/Dapars2/Dapars0.9")
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
Healthy_AD <- Healthy_AD[,1:2]
library(tidyr)
long<-separate(Healthy_AD,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- long
Healthy_AD <- Healthy_AD[,1:2]
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
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
a <- a$ENSEMBLTRANS
APAweights <- rbind(APAweights_Pos,APAweights_Neg)
APAweights <- APAweights[order(APAweights$Factor1,decreasing = F),]
APAweights$Rank <- seq(1,82,1)
APAweights_Top <- subset(APAweights,APAweights$ENSEMBLTRANS %in% a)
topAPA <- APAweights_Top$GeneName 
APAweights$Color <- ifelse(APAweights$ENSEMBLTRANS %in% a, "#2678C2","grey")
library(ggrepel)
p <- ggplot(
  APAweights, aes(x = Weight, y = Rank,colour=Color)) +
  geom_point(aes(color = Color), size=3) +
  scale_color_manual(values = c("#2678C2", "gray")) +
  geom_text_repel(data = APAweights_Top, aes(x = Weight, y = Rank, label = GeneName),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE,
                  color = "black",max.overlaps=20)+
  geom_vline(xintercept=c(-1,0,1),lty=7,col="#666666",lwd=0.2) +
  labs(x="Weight",y="Rank") +
  ggtitle("MOFA2 APA Rank")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(color = "black",size = 12,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"))+
  scale_x_continuous(breaks = seq(-1, 1, 1))
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_Rank.pdf", egg::set_panel_size(p, width=unit(3.5, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
APAweights_Top$abs <- abs(APAweights_Top$Weight)
library(ggpubr)
p <- ggdotchart(APAweights_Top, x = "ENSEMBLTRANS", y = "Weight",
                color = "#2678C2",                               
                sorting = "ascending",                      
                add = "segments",                            
                add.params = list(color = "#2678C2", size = 1),
                rotate = TRUE,                             
                dot.size = "abs",                               
                font.label = list(color = "black", size = 7,
                                  vjust = 0.5),              
                ggtheme = theme_pubr(),                      
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
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_Top.pdf", egg::set_panel_size(p, width=unit(2, "in"), height=unit(2.5, "in")), 
       width = 6, height = 6, units = 'in', dpi = 600)
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightAPAGene <- KEGGa$GeneName
library(clusterProfiler)
library(org.Hs.eg.db)
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",
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
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_APA_KEGG.csv",row.names = F)
kegg <- kegg[,-7]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("APA",6)))
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
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_APA_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 3, units = 'in', dpi = 600)
gene <- c(topGE,topAS,topAPA)
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
kegg@result <- kegg@result[order(kegg@result[["Count"]],decreasing = T),]
kegg<-setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg <- as.data.frame(kegg)
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
write.csv(kegg,file = "MOFA2_sigGeneKEGG.csv",row.names = F)
kegg <- kegg[,-7]
kegg <- kegg[1:6,]
kegg$num <- rep(1:6)
kegg$type_order=factor(rev(as.character(kegg$Description)),labels=rev(rep("GENE",6)))
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
  xlab("Pathway")+ylab("Number")+
  guides(fill = guide_legend(title = 'Layer'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.8, linetype="solid"),
        legend.key = element_blank(),
        plot.title = element_text(color = "black",size = 10,face = "plain",hjust = 0.5),
        axis.title.x = element_text(color = "black",size = 18,face = "plain"),
        axis.title.y = element_text(color = "black",size = 18,face = "plain"),
        panel.grid= element_blank(),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 18,face = "plain"),
        axis.text.y = element_text(color = "black",size = 18,face = "plain"),
        legend.title = element_text(color = "black",size = 8,face = "plain"),
        legend.text = element_text(color = "black",size = 15,face = "plain"),
        legend.direction = 'vertical', legend.position ="right")
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("MOFA2_ThreeLayerGenes_KEGG.pdf", egg::set_panel_size(p1, width=unit(3, "in"), height=unit(2, "in")), 
       width = 8, height = 3, units = 'in', dpi = 600)
MOFA2_Mos_sigGene <- c(topGE,topAS,topAPA)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
GE_Mos_sigGene <- c("SYT1","CHN1","SNAP25","VSNL1","ENC1","TNS1","SGK1","CPM","CLMN","PPFIBP2")
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
            text_size = 4)
p
setwd("C:/Users/yujie/Desktop/datacollect/20240827/MOFA2")
ggsave("Total_Venn.pdf", egg::set_panel_size(p, width=unit(5, "in"), height=unit(4, "in")), 
       width = 8, height = 8, units = 'in', dpi = 600)
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
GEFPKM <- read.csv("GE_RemoveBatch_FPKM.csv",header = T,row.names = 1)
rownames(GEFPKM) <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(GEFPKM))
GEFPKM <- as.data.frame(t(GEFPKM))
GEFPKM$GSM_number <- rownames(GEFPKM)
GEFPKM <- GEFPKM[,c("ENSG00000163032","GSM_number")]
setwd("C:/Users/yujie/Desktop/datacollect")
BackgroundInformation <- read.csv("BackgroundInformationkNN20240827_k_15.csv",header = T)
sum <- merge(GEFPKM,BackgroundInformation,by="GSM_number")
sum <- sum[,c("ENSG00000163032","Sample")]
GE <- sum
GE <- GE[order(GE$Sample),]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
AScount <- read.csv("AS_RemoveBatch_VoomCount.csv",header = T,row.names = 1,check.names = F)
library(MOFA2)
filepath <- "C:/Users/yujie/Desktop/datacollect/20240827/MOFA2/MOFA2Results_2"
model <- load_model(filepath)
ASweights <- as.data.frame(model@expectations[["W"]][["AS"]])
ASweights$GeneName <- rownames(ASweights)
ASweights_Pos <- ASweights[ASweights$Factor1>=0,]
ASweights_Pos$Weight <- (ASweights_Pos$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Pos <- ASweights_Pos[order(abs(ASweights_Pos$Factor1),decreasing = T),]
ASweights_Neg <- ASweights[ASweights$Factor1<0,]
ASweights_Neg$Weight <- (ASweights_Neg$Factor1-min(abs(ASweights$Factor1)))/(max(abs(ASweights$Factor1))-min(abs(ASweights$Factor1)))
ASweights_Neg <- ASweights_Neg[order(abs(ASweights_Neg$Factor1),decreasing =T),]
a <- rbind(ASweights_Pos,ASweights_Neg)
KEGGa <- a
a <- a[order(abs(a$Weight),decreasing = T),]
a <- a[1:10,]
setwd("C:/Users/yujie/Desktop/datacollect/20240827")
m <- read.table("ASdiff_n369.txt",fill = TRUE,header = T,check.names = F)
m <- m[m$Position %in% a$GeneName,]
m <- subset(m,m$Genename %in% c("VSNL1","ENC1","SNAP25"))
rownames(m) <- m$Genename
m <- m[,c(3:356)]
m <- as.data.frame(t(m))
m$Sample <- rownames(m)
sum <- merge(GE,m,by="Sample")
sum$Severity <- substr(sum$Sample,1,1)
sum <- sum[,-1]
sum$Severity <- as.numeric(sum$Severity)
all <- sum
library(tidyverse)
n <- dim(all)[1]
y <- all$Severity
folds <- createFolds(y,k=6)
library(ROCR)
library(magrittr) 
library(plyr)
auc_value<-as.numeric()
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){
  fold_test <- all[folds[[i]],] 
  fold_train <- all[-folds[[i]],] 
  model <- glm(fold_train$Severity~.,data=fold_train,
               family = binomial(link = "logit"))
  fold_predict <- stats::predict(model,type='response',newdata=fold_test)
  pred <- prediction(fold_predict,fold_test$Severity)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))
}
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3",
                       "X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData$Xmean <- apply(rocData[,1:6],1,mean)
rocData$Ymean <- apply(rocData[,7:12],1,mean)
rocData <- rocData[,13:14]
rocData$Gene <- c(rep("all", 60))
a <- seq(0,1,1/59)
inn <- data.frame(a,a,rep("a", 60))
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
  ggtitle("MOFA2 ROC AUC=0.754") +
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
ggsave("MOFA2_ROC.pdf", egg::set_panel_size(a, width=unit(5, "in"), height=unit(4, "in")), 
       width = 7, height = 6, units = 'in', dpi = 600)

