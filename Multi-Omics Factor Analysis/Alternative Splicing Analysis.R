counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
BackgroundInformation <- read.csv("BackgroundInformationAPOE.csv",header = T)
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:11)]
counts <- sum
pca.info <- prcomp(counts)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
a <- a[order(a$Sample),]
p1 <- ggplot(data=a,aes(x=PC1,y=PC2,color=Type,shape =Batch))+geom_point()
df <- counts
a <- adonis2(df ~ CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch+Severity,data = Temporal_Lobe,permutations=9999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(9,10),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- na.omit(b)
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")+geom_text(aes(label=Pr))








counts <- read.csv("RawIntroRead.csv",header = T,row.names = 1)
counts <- counts[,order(colnames(counts))]
BackgroundInformation <- read.csv("BackgroundInformationAPOE.csv",header = T)
Temporal_Lobe <- subset(BackgroundInformation,BackgroundInformation$Brain_Region==3)
Temporal_Lobe <- Temporal_Lobe[order(Temporal_Lobe$GSM_number),]
counts <- as.data.frame(t(counts))
counts$GSM_number <- rownames(counts)
sum <- merge(Temporal_Lobe,counts,by = "GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:11)]
counts <- as.data.frame(t(sum))
group <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(counts))) #
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design) <- colnames(counts)
keep.exprs <- filterByExpr(counts, group=group,min.count = 3,min.total.count = 3,large.n = 3, min.prop = 0.3)
counts <- counts[keep.exprs,]
dge<-DGEList(counts)
y <- voom(dge, design, plot = T,normalize.method ="quantile")
newData <- y$E
NormalizedCount <- as.data.frame(newData)
df <- as.data.frame(t(newData))
pca.info <- prcomp(df)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
p1 <- ggplot(data=a,aes(x=PC1,y=PC2,color=Type,shape =Batch))+geom_point()
df <- as.data.frame(t(newData))
a <- adonis2(df ~ CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch+Severity,data = Temporal_Lobe,permutations=9999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(9,8),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")+geom_text(aes(label=Pr)) 
batch <- as.factor(Temporal_Lobe$Batch)
group <- factor(gsub(pattern = '\\_\\d*',replacement = '',x=colnames(counts))) #
design <- model.matrix(~group)
colnames(design) <- levels(group)
rownames(design) <- colnames(counts)
RemoveData <- removeBatchEffect(newData,batch = batch,design = design)
RemoveData <- as.data.frame(t(RemoveData))
pca.info <- prcomp(RemoveData)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
colnames(pca.data)[1] <- 'Sample'
a <- merge(Temporal_Lobe,pca.data,by="Sample")
p <- ggplot(data=a,aes(x=PC1,y=PC2,color=Type,shape =Batch))+geom_point()
a <- adonis2(RemoveData ~ CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch+Severity,data = Temporal_Lobe,permutations=9999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(9,8),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b[order(b$R2,decreasing = TRUE),]
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")+geom_text(aes(label=Pr))
RemoveData <- RemoveData[,order(colnames(RemoveData))]
RemoveData$groupHealthymean <- apply(RemoveData[,1:97],1,mean)
RemoveData$groupADmean <- apply(RemoveData[,98:354],1,mean)
RemoveData$Group_diff <- RemoveData$groupADmean-RemoveData$groupHealthymean
RemoveData$log2FoldChange <- log2(RemoveData$groupADmean/RemoveData$groupHealthymean)
RemoveData$Pvalue <- apply(RemoveData,1,function(x) t.test(x[1:97],x[98:354])$p.value)
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
write.csv(RemoveData,file = "VoomRemoveBatch.csv",quote = F)





                           

a <- read.csv("VoomRemoveBatchDiff.csv",header = T,row.names = 1)
a <- a[,355:360]
a$TP <- rownames(a)
a <- separate(a,TP,into = c("Chr","Start","End","Cluster"),sep = "([:])")
a$Position <- paste(a$Chr,a$Start,a$End,sep = ":")
a <- a[,-c(7:10)]
rownames(a) <- NULL
m <- read.table("humanv40annotation_all_introns.bed",fill = TRUE)
m$Position <- paste(m$V1,m$V2,m$V3,sep = ":")
m <- m[,c("V4","Position")]
m <- subset(m,m$Position %in% a$Position)
sum <- merge(a,m,by="Position")
sum <- sum[!duplicated(sum),]
sum$change <- ifelse(sum$Padjust < 0.05 & abs(sum$log2FoldChange) >= 0.1,ifelse(sum$log2FoldChange > 0.1 ,'Up','Down'),'Stable')
p <- ggplot(sum, aes(x = log2FoldChange, y = -log10(Padjust),colour=change)) +geom_point(aes(color = change)) 
a <- read.csv("VoomRemoveBatchDiff.csv",header = T,row.names = 1,check.names = F)
b <- subset(a,a$Padjust < 0.05 & abs(a$log2FoldChange) > 0.1)
b <- b[order(b$Padjust),]
b <- b[1:10,]
b <- b[,-c(355:360)]
b <- as.data.frame(t(b))
b$Group <- rownames(b)
b$Group <- substr(b$Group,1,1)
n <- dim(b)[1]
y <- b$Group
all <- b
folds <- createFolds(y,k=6)
auc_value<-as.numeric()
rocData <- data.frame(matrix(ncol = length(folds[[1]]), nrow = 1))
for(i in 1:6){ 
  fold_test <- all[folds[[i]],] 
  fold_train <- all[-folds[[i]],] 
  model <- glm(as.numeric(fold_train$Group)~.,data=fold_train,family = binomial(link = "logit"))
  fold_predict <- predict(model,type='response',newdata=fold_test)
  pred <- prediction(predictions = fold_predict, labels = fold_test$Group)
  roc <- performance(pred,"tpr","fpr")
  auc_value<- append(auc_value,as.numeric(performance(pred, measure = "auc")@y.values[[1]]))
  rocData <- rbind.fill(rocData,data.frame(t(roc@x.values[[1]])),data.frame(t(roc@y.values[[1]])))}
rocData <- rocData[-1,]
rownames(rocData) <- c("X1","Y1","X2","Y2","X3","Y3","X4","Y4","X5","Y5","X6","Y6") 
rocData <- rocData[order(rownames(rocData)),]
rocData <- as.data.frame(t(rocData))
rocData <- read.csv("ASroc.csv",header = T)
rocData <- rocData[,13:14]
rocData$Gene <- c(rep("all", 61))
a <- seq(0,1,1/60)
inn <- data.frame(a,a,rep("a", 61))
colnames(inn) <- c("Xmean","Ymean","Gene")
ROC <- rbind(rocData,inn)
a <- ggplot(data = ROC, mapping = aes(x = Xmean, y = Ymean,group=Gene))+geom_line(aes(linetype=Gene,color=Gene))
a <- read.csv("VoomRemoveBatchDiff.csv",header = T,row.names = 1)
a <- a[,355:360]
b <- subset(a,a$Padjust < 0.05 & abs(a$log2FoldChange) > 0.1)
b$TP <- rownames(b)
b <- separate(b,TP,into = c("Chr","Start","End","Cluster"),sep = "([:])")
b$Position <- paste(b$Chr,b$Start,b$End,sep = ":")
b <- b[,-c(7:10)]
m <- read.table("humanv40annotation_all_introns.bed",fill = TRUE)
m <- m[,-10]
m$Position <- paste(m$V1,m$V2,m$V3,sep = ":")
sum <- merge(b,m,by="Position")
sum$change <- ifelse(sum$Padjust < 0.05 & abs(sum$log2FoldChange) >= 0.1,ifelse(sum$log2FoldChange > 0.1 ,'Up','Down'),'Stable')
gene <- unique(sum$V4)
df1 <- bitr(gene, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
GOgene <- as.character(df1$ENTREZID)
GOgene <- na.omit(GOgene)
BPplot <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(BPplot,x = "Count",color = "p.adjust",showCategory = 15)
CCplot <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(CCplot,x = "Count",color = "p.adjust",showCategory = 15)
MFplot <- enrichGO(GOgene,OrgDb = org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
p <- barplot(MFplot,x = "Count",color = "p.adjust",showCategory = 15)
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = GOgene,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
p <- barplot(AD_Healthyres_kegg,x = "Count",color = "p.adjust",showCategory = 20)
sessionInfo()
