a <- read.csv("sum.csv",header = T,fill = TRUE)
library(ggsci)
b <- a[,c("Brain_Region","Severity")]
library(dplyr)
d <- b |> group_by(Severity,Brain_Region) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d%>%mutate(Severity=case_when(d$Severity==0~"Healthy",
                                   d$Severity==1~"AD",
                                   d$Severity==2~"MCI"))
d <- d%>%mutate(Brain_Region=case_when(d$Brain_Region==0~"Hippo",
                                       d$Brain_Region==1~"PC",
                                       d$Brain_Region==2~"Blood",
                                       d$Brain_Region==3~"TL",
                                       d$Brain_Region==4~"TIPC",
                                       d$Brain_Region==5~"PL",
                                       d$Brain_Region==6~"OB",
                                       d$Brain_Region==7~"CB"))
d$Severity <- as.factor(d$Severity)
d$Brain_Region <- as.factor(d$Brain_Region)
library(ggplot2)
p <- ggplot(d,aes(x=Severity,y=freq,fill=Brain_Region))+geom_bar(stat="identity")
b <- a[,c("Bank_Location","Severity")]
d <- b |> group_by(Severity,Bank_Location) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d%>%mutate(Bank_Location=case_when(d$Bank_Location==0~"America",
                                        d$Bank_Location==1~"Britain",
                                        d$Bank_Location==2~"China",
                                        d$Bank_Location==3~"Australia",
                                        d$Bank_Location==4~"Italy",
                                        d$Bank_Location==5~"Japan",
                                        d$Bank_Location==6~"Colombia"))
d$Severity <- as.factor(d$Severity)
d$Bank_Location <- as.factor(d$Bank_Location)
p <- ggplot(d,aes(x=Severity,y=freq,fill=Bank_Location))+geom_bar(stat="identity")
b <- a[,c("Braak_Stage","Severity")]
d <- b |> group_by(Severity,Braak_Stage) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d%>%mutate(Braak_Stage=case_when(d$Braak_Stage==0~"0",
                                      d$Braak_Stage==1~"1~2",
                                      d$Braak_Stage==2~"3~4",
                                      d$Braak_Stage==3~"5~6"))
d$Severity <- as.factor(d$Severity)
d$Braak_Stage <- as.factor(d$Braak_Stage)
p <- ggplot(d,aes(x=Severity,y=freq,fill=Braak_Stage))+geom_bar(stat="identity")
ggsave("p3.pdf", width = 8, height = 7)
b <- a[,c("CERAD","Severity")]
d <- b |> group_by(Severity,CERAD) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d%>%mutate(CERAD=case_when(d$CERAD==0~"0",
                                d$CERAD==1~"A",
                                d$CERAD==2~"B",
                                d$CERAD==3~"C"))
d$Severity <- as.factor(d$Severity)
d$CERAD <- as.factor(d$CERAD)
p <- ggplot(d,aes(x=Severity,y=freq,fill=CERAD))+geom_bar(stat="identity")
a <- read.csv("BackgroundInformationAPOE.csv",header = T,fill = TRUE)
b <- a[,c("Phenotypedit","Severity")]
d <- b |> group_by(Severity,Phenotypedit) |> summarise(freq=n()) 
d <- as.data.frame(d)
d$Severity <- as.factor(d$Severity)
d$Phenotypedit <- as.factor(d$Phenotypedit)
colnames(d)[2] <- "APOE"
p <- ggplot(d,aes(x=Severity,y=freq,fill=APOE))+geom_bar(stat="identity")
a <- read.csv("sum.csv",header = T,fill = TRUE)
b <- a[,c("Age","Severity")]
d <- b |> group_by(Severity,Age) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d %>% mutate(Age=case_when(40 < d$Age & d$Age <= 50 ~"40~50",
                                50 < d$Age & d$Age <= 60 ~"50~60",
                                60 < d$Age & d$Age <= 70 ~"60~70",
                                70 < d$Age & d$Age <= 80 ~"70~80",
                                80 < d$Age & d$Age <= 90 ~"80~90",
                                90 < d$Age & d$Age <= 100 ~"90~100",
                               d$Age >= 100 ~">100"))
m <- d |> group_by(Severity,Age) |> summarise(sum = sum(freq))
m <- as.data.frame(m)
m$Severity <- as.factor(m$Severity)
m$Age <- as.factor(m$Age)
m <- arrange(m,Severity,Age)
p <- ggplot(m,aes(x=Severity,y=sum,fill=Age))+geom_bar(stat="identity")
b <- a[,c("Sex","Severity")]
d <- b |> group_by(Severity,Sex) |> summarise(freq=n()) 
d <- as.data.frame(d)
d <- d%>%mutate(Sex=case_when(d$Sex==1~"F",
                              d$Sex==0~"M"))
d$Severity <- as.factor(d$Severity)
d$Sex <- as.factor(d$Sex)
p <- ggplot(d,aes(x=Severity,y=freq,fill=Sex))+geom_bar(stat="identity")






FPKM <- read.csv("ADPatientFPKM.csv",header = T ,row.names = 1)
FPKM <- as.matrix(FPKM)
FPKM[FPKM<=1] <- 1
FPKM <- as.data.frame(FPKM)
FPKM <- na.omit(FPKM)
SampleInformtion <- read.csv("BackgroundInformation.csv",header = T)
FPKM <- as.data.frame(t(FPKM))
FPKM$GSM_number <- rownames(FPKM)
sum <- merge(FPKM,SampleInformtion,by = "GSM_number", all = FALSE)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67949)]
pca.info <- prcomp(sum)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
pca.data <- pca.data%>%mutate(Type=case_when(pca.data$Type==0~"Healthy",
                                             pca.data$Type==1~"AD",
                                             pca.data$Type==2~"MCI"))
pca.data$Type <- as.factor(pca.data$Type)
p <- ggplot(data=pca.data,aes(x=PC1,y=PC2,color=Type))
FPKM <- as.data.frame(t(FPKM))
FPKM <- FPKM[order(rownames(FPKM)),]
sumknnImputation <- read.csv("sumknnImputation.csv",header = T)
sumknnImputation <- sumknnImputation[order(sumknnImputation$GSM_number),]
sumknnImputation <- sumknnImputation[,-7]
a <- adonis2(FPKM ~ Severity+CERAD+Braak_Stage+Brain_Region+Age+Sex+Bank_Location+Batch, data = sumknnImputation,permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- na.omit(b)
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2,fill=category))+geom_bar(stat="identity")








myFPKM <- read.csv("ADPatientFPKM.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- subset(myFPKM,rownames(myFPKM) %in% BackgroundInformation$GSM_number)
myFPKM <- myFPKM[order(rownames(myFPKM)),]
sumKnn <- read.csv("sumknnImputation.csv",header = T)
sumKnn <- subset(sumKnn,sumKnn$GSM_number %in% BackgroundInformation$GSM_number)
sumKnn <- sumKnn[order(sumKnn$GSM_number),]
a <- adonis2(myFPKM ~ Severity+CERAD+Braak_Stage+Age+Sex+Batch+Bank_Location, data = sumKnn, permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(8,9),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr<=0.001~"***",b$Pr>0.001~""))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount$GSM_number <- rownames(mycount)
BackgroundInformation <- read.csv("BackgroundInformation.csv",header = T)
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67949)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount)))
coldata <- data.frame(row.names =colnames(mycount),condition) 
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)
vst <- vst(dds, blind=T)
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"
p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
p <- ggplot(data=p1data,aes(x=PC1,y=PC2,color=Severity,shape =Batch))









BackgroundInformation <- read.csv("BackgroundInformation.csv",header = T )
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
mycount <- read.csv("ADPatientCount.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount <- subset(mycount,rownames(mycount) %in% BackgroundInformation$GSM_number)
mycount <- mycount[order(rownames(mycount)),]
mycount <- as.matrix(t(mycount))
adjust_counts <- ComBat_seq(mycount, batch=BackgroundInformation$Batch,group=BackgroundInformation$Severity,covar_mod=NULL)
write.csv(adjust_counts,file = "Temporal_LobeRawCount_RemoveBatch20220831.csv",quote = F)
mycount <- as.data.frame(t(adjust_counts))
mycount$GSM_number <- rownames(mycount)
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67949)]
mycount <- as.data.frame(t(sum))
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) 
coldata <- data.frame(row.names =colnames(mycount),condition)  
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition) 
vst <- vst(dds, blind=T)
p1data <- plotPCA(vst,returnData = T)
colnames(p1data)[5] <- "Sample"
p1data <- merge(p1data,BackgroundInformation,by= "Sample",ALL=FALSE)
p <- ggplot(data=p1data,aes(x=PC1,y=PC2,color=Severity,shape=Batch))
myFPKM <- read.csv("TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T,row.names = 1)
myFPKM <- as.data.frame(t(myFPKM))
myFPKM <- myFPKM[order(rownames(myFPKM)),]
BackgroundInformation <- read.csv("sumknnImputation.csv",header = T)
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
a <- adonis2(myFPKM ~ Severity+CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch, data = BackgroundInformation,permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(8,9),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b%>%mutate(Pr=case_when(b$Pr<=0.1~".",b$Pr>0.1~""))
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")








mycount <- read.csv("Temporal_LobeRawCount_RemoveBatch20220831.csv",header = T,row.names = 1)
mycount <- as.data.frame(t(mycount))
mycount$GSM_number <- rownames(mycount)
BackgroundInformation <- read.csv("BackgroundInformation.csv",header = T)
BackgroundInformation <- BackgroundInformation[BackgroundInformation$Brain_Region=="3",]
BackgroundInformation <- BackgroundInformation[order(BackgroundInformation$GSM_number),]
sum <- merge(mycount,BackgroundInformation,by = "GSM_number", all = FALSE)
rownames(sum) <- sum$Sample
sum <- sum[,-c(1,67941:67949)]
mycount <- as.data.frame(t(sum))
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(mycount))
rownames(mycount) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
mycount$ENSEMBL <- rownames(mycount)
df <- merge(df1,mycount,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
mycount <- df[,-c(1:2)]
condition <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(mycount))) 
coldata <- data.frame(row.names =colnames(mycount),condition)  
dds <- DESeqDataSetFromMatrix(countData = mycount,colData = coldata,design = ~condition)  
keep <- rowSums(cpm(mycount)>0.5) >= 20
dds <- dds[keep,]
dds_norm <- DESeq(dds)
res1 <- results(dds_norm,contrast = c("condition","1","0"))
res1 <- res1[order(res1$padj,decreasing = F),]
res1 <- na.omit(res1)
res1 <- res1[!duplicated(rownames(res1)),]
res1 <- as.data.frame(res1)
AD_Healthyres <- res1
AD_Healthyres$change <- ifelse(AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >= 0.585, 
                               ifelse(AD_Healthyres$log2FoldChange > 0.585 ,'Up','Down'),'Stable')
write.csv(AD_Healthyres,file = "Temporal_Lobe_ADvsHealthy_DEseq2.csv")
AD_Healthyres_def <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585)
write.csv(AD_Healthyres_def,file = "Temporal_Lobe_ADvsHealthy_DEseq2_Significant.csv")
df1 <- bitr(rownames(AD_Healthyres_def),fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
GOgene <- as.character(df1$ENTREZID)
GOgene <- na.omit(GOgene)
AD_Healthyres_GOBP <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(AD_Healthyres_GOBP,x = "Count",color = "p.adjust",showCategory = 15)
AD_Healthyres_GOCC <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(AD_Healthyres_GOCC, x = "Count", color = "p.adjust",showCategory = 15)
AD_Healthyres_GOMF <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(AD_Healthyres_GOMF,x = "Count",color = "p.adjust",showCategory = 15)
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
AD_Healthyres_kegg <- enrichKEGG(gene = GOgene,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.1,pAdjustMethod  = "BH",qvalueCutoff  = 0.1)
p <- barplot(AD_Healthyres_kegg,x = "Count",color = "p.adjust",showCategory = 20)
p <- ggplot(AD_Healthyres, aes(x = log2FoldChange, y = -log10(padj),colour=change)) 
AD_Healthyres <- AD_Healthyres[order(AD_Healthyres$padj,decreasing = F),]
choose_gene <-subset(AD_Healthyres, AD_Healthyres$padj < 0.05 & abs(AD_Healthyres$log2FoldChange) >0.585)
myfpkm <- read.csv(file = "TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T ,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
AD_Healthyres_heatmap <- subset(myfpkm, rownames(myfpkm) %in% rownames(choose_gene))
BackgroundInformation <- read.csv("BackgroundInformation.csv",header = T )
AD_Healthyres_heatmap <- as.data.frame(t(AD_Healthyres_heatmap))
AD_Healthyres_heatmap$GSM_number <- rownames(AD_Healthyres_heatmap)
sum <- merge(AD_Healthyres_heatmap,BackgroundInformation,by = "GSM_number",all= FALSE)
rownames(sum) <- sum$Sample
sum <- arrange(sum,rownames(sum))
data <- sum[-c(1,1343:1351)]
data <- as.data.frame(t(data))
data$cv <- apply(data, 1, function(x){sd(x)/mean(x)*100})
data_df <- data[order(data$cv, decreasing = T),1:354]
data_scale <- as.data.frame(t(apply(data_df,1,scale)))
names(data_scale) <- names(data_df)
data_scale[is.na(data_scale)] <- min(data_scale,na.rm = T)*0.01
data_scale <- as.matrix(data_scale)
data_scale[data_scale>=2] <- 2
data_scale[data_scale<=-2] <- -2
pdf(file = "Heatmap.pdf",width =3,height = 4)
p <- Heatmap(data_scale)
print(p)
dev.off()








BackgroundInformationAPOE <- read.csv("BackgroundInformationAPOE.csv",header = T)
gene <- read.csv("Temporal_Lobe_ADvsHealthy_DEseq2_Significant.csv",header = T,row.names = 1)
myfpkm <- read.csv("TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% rownames(gene))
sig_FPKM <- as.data.frame(t(sig_FPKM))
gsg <- goodSamplesGenes(sig_FPKM, verbose=3)
BackgroundInformation <- subset(BackgroundInformationAPOE,BackgroundInformationAPOE$GSM_number %in% rownames(sig_FPKM))
rownames(BackgroundInformation) <- BackgroundInformation$GSM_number
BackgroundInformation <- BackgroundInformation[,-c(1,2,8,9,10)]
powers <- c(c(1:10), seq(from=12, to=50, by=2))
sft <- pickSoftThreshold(sig_FPKM, powerVector=powers, verbose=5, networkType="unsigned")
sft$powerEstimate <- 9
cor <- WGCNA::cor
net <- blockwiseModules(sig_FPKM, power = sft$powerEstimate,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "ADPatientTOM", 
                        verbose = 3)
cor<-stats::cor
mergedColors <- labels2colors(net$colors)
pdf(file = "Module_colors.pdf",width =5,height = 3)
p <- plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          cex.rowText = 0.9,cex.colorLabels = 0.9, cex.dendroLabels = 0.9)
print(p)
dev.off()
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
nGenes <- ncol(sig_FPKM)
nSamples <- nrow(sig_FPKM)
a <- moduleEigengenes(sig_FPKM, moduleColors)
MEs0 <- moduleEigengenes(sig_FPKM, moduleColors)$eigengenes
moduleTraitCor <- cor(MEs0, BackgroundInformation,use = 'p',method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(2,2)
textMatrix <-  paste(signif(moduleTraitCor, 2),"\n(",signif(moduleTraitPvalue, 1),")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file = "Module_trait_relationships.pdf",width =5,height = 2.5)
par(mar = c(5, 7, 3, 3))
p <- labeledHeatmap(Matrix = moduleTraitCor,
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
                     main = paste("Module-trait relationships"))
print(p)
dev.off()








APOE <- as.data.frame(BackgroundInformation$APOE)
colnames(APOE) = "APOE"
modNames <- substring(names(MEs0), 3)
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) <- paste("MM", modNames, sep="")
colnames(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <- as.data.frame(cor(sig_FPKM, APOE, use = "p",method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- paste("GS.", names(APOE), sep="");
colnames(GSPvalue) <- paste("p.GS.", names(APOE), sep="")
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
allLLIDs <- colnames(sig_FPKM)
module <- "blue"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]
m$MMblue <- abs(m$MMblue)
m$GS.APOE <- abs(m$GS.APOE)
m <- subset(m,m$p.MMblue<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$GS.APOE) > 0.3)
p <- ggplot(m, aes(x =abs(MMblue) , y = abs(GS.APOE))) +geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")
module <- "turquoise"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]
m$MMturquoise <- abs(m$MMturquoise)
m$GS.APOE <- abs(m$GS.APOE)
m <- subset(m,m$p.MMturquoise<0.05 & m$p.GS.APOE <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$GS.APOE) > 0.3)
p <- ggplot(m,aes(x =abs(MMturquoise),y = abs(GS.APOE)))+geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")
Braak_Stage <- as.data.frame(BackgroundInformation$Braak Stage)
modNames <- substring(names(MEs0), 3)
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) <- paste("MM", modNames, sep="")
colnames(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <- as.data.frame(cor(sig_FPKM, Braak_Stage, use = "p",method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- paste("GS.", names(Braak_Stage), sep="");
colnames(GSPvalue) <- paste("p.GS.", names(Braak_Stage), sep="")
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
module <- "blue"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]
m$MMblue <- abs(m$MMblue)
m$GS.Braak Stage <- abs(m$GS.Braak Stage)
m <- subset(m,m$p.MMblue<0.05 & m$p.GS.Braak Stage <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$`GS.Braak Stage`) > 0.3)
p <- ggplot(m, aes(x =abs(MMblue) , y = abs(GS.Braak Stage))) +geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")
module <- "turquoise"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]
m$MMturquoise <- abs(m$MMturquoise)
m$GS.Braak Stage <- abs(m$GS.Braak Stage)
m <- subset(m,m$p.MMturquoise<0.05 & m$p.GS.Braak Stage <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$GS.Braak Stage) > 0.3)
p <- ggplot(m, aes(x =abs(MMturquoise) , y = abs(GS.Braak Stage))) +geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")
CERAD <- as.data.frame(BackgroundInformation$CERAD)
colnames(CERAD) = "CERAD"
modNames <- substring(names(MEs0), 3)
geneModuleMembership <- as.data.frame(cor(sig_FPKM, MEs0, use = "p",method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <- as.data.frame(cor(sig_FPKM, CERAD, use = "p",method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- paste("GS.", names(CERAD), sep="");
colnames(GSPvalue) <- paste("p.GS.", names(CERAD), sep="")
geneModuleMembership <- merge(geneModuleMembership,MMPvalue,by="row.names")
geneTraitSignificance <- merge(geneTraitSignificance,GSPvalue,by="row.names")
module <- "blue"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,2,5,8,9)]
m$MMblue <- abs(m$MMblue)
m$GS.CERAD <- abs(m$GS.CERAD)
m <- subset(m,m$p.MMblue<0.05 & m$p.GS.CERAD <0.05)
n <- subset(m,abs(m$MMblue > 0.8) & abs(m$GS.CERAD) > 0.3)
p <- ggplot(m, aes(x =abs(MMblue) , y = abs(GS.CERAD))) +geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")
module <- "turquoise"
modGenes <- (moduleColors==module)
modLLIDs <- allLLIDs[modGenes]
b <- subset(geneModuleMembership,geneModuleMembership$Row.names %in% modLLIDs)
d <- subset(geneTraitSignificance,geneTraitSignificance$Row.names %in% modLLIDs)
m <- merge(b,d,by="Row.names")
m <- m[,c(1,4,7,8,9)]
m$MMturquoise <- abs(m$MMturquoise)
m$GS.CERAD <- abs(m$GS.CERAD)
m <- subset(m,m$p.MMturquoise<0.05 & m$p.GS.CERAD <0.05)
n <- subset(m,abs(m$MMturquoise > 0.8) & abs(m$GS.CERAD) > 0.3)
p <- ggplot(m, aes(x =abs(MMturquoise) , y = abs(GS.CERAD))) +geom_smooth(method = 'lm',se = F)+stat_cor(data=m, method = "spearman")








gene <- read.table("turquoise_mark_gene.txt",header = T)
gene <- gene[!duplicated(gene),]
ADvsHealthy_DEseq2 <- read.csv("Temporal_Lobe_ADvsHealthy_DEseq2.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]
gene <- ADvsHealthy_DEseq2$X
myfpkm <- read.csv("TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)
BackgroundInformationAPOE <- read.csv("BackgroundInformationAPOE.csv",header = T)
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,13)]
gene_group$Severity <- as.factor(gene_group$Severity)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_violin(aes_string(fill = .x),trim = FALSE)+geom_signif(comparisons = list(c("AD", "Healthy")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,9)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),c("1","2"),c("1","3"),c("2","3")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,8)]
gene_group$CERAD <- as.factor(gene_group$CERAD)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),c("1","2"),c("1","3"),c("2","3")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,16)]
gene_group <- na.omit(gene_group)
gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("1","2")))
gene <- read.table("blue_mark_gene.txt",header = T)
gene <- gene[!duplicated(gene),]
ADvsHealthy_DEseq2 <- read.csv("Temporal_Lobe_ADvsHealthy_DEseq2.csv",header = T)
ADvsHealthy_DEseq2 <- subset(ADvsHealthy_DEseq2,ADvsHealthy_DEseq2$X %in% gene)
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[order(ADvsHealthy_DEseq2$baseMean,decreasing = T),]
ADvsHealthy_DEseq2 <- ADvsHealthy_DEseq2[1:5,]
gene <- ADvsHealthy_DEseq2$X
myfpkm <- read.csv("TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(myfpkm))
rownames(myfpkm) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
myfpkm$ENSEMBL <- rownames(myfpkm)
df <- merge(df1,myfpkm,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
myfpkm <- df[,-c(1:2)]
sig_FPKM <- subset(myfpkm,rownames(myfpkm) %in% gene)
sig_FPKM <- as.data.frame(t(sig_FPKM))
sig_FPKM$GSM_number <- rownames(sig_FPKM)
BackgroundInformationAPOE <- read.csv("BackgroundInformationAPOE.csv",header = T)
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,13)]
gene_group$Severity <- as.factor(gene_group$Severity)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_violin(aes_string(fill = .x),trim = FALSE)+geom_signif(comparisons = list(c("AD", "Healthy")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,9)]
gene_group$Braak_Stage <- as.factor(gene_group$Braak_Stage)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y,  ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),c("1","2"),c("1","3"),c("2","3")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,8)]
gene_group$CERAD <- as.factor(gene_group$CERAD)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y, ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("0","3"),c("1","2"),c("1","3"),c("2","3")))
gene_group <- merge(sig_FPKM,BackgroundInformationAPOE,by = "GSM_number")
gene_group <- gene_group[,c(2:6,16)]
gene_group <- na.omit(gene_group)
gene_group$Phenotypedit <- as.factor(gene_group$Phenotypedit)
x <- names(gene_group)[6]
y <- names(gene_group)[-6]
plot_list <- map2(x, y,  ~ gene_group %>% ggplot(aes_string(x = .x, y = .y)) +geom_boxplot(aes_string(fill = .x),outlier.shape = NA)+geom_signif(comparisons = list(c("0", "1"),c("0","2"),c("1","2")))
sessionInfo()
