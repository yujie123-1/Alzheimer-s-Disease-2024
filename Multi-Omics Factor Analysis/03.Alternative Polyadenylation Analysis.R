Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
a <- Healthy_AD[,-c(1:4,1067:1072)]
a <- a[,seq(0,ncol(a),3)]
Healthy_AD <- cbind(Healthy_AD[,1:4],a)
long<-separate(Healthy_AD,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- long
Healthy_AD <- Healthy_AD[,-c(2:7)]
rownames(Healthy_AD) <- Healthy_AD$ENSEMBLTRANS
Healthy_AD <- Healthy_AD[,-1]
Healthy_AD <- as.data.frame(t(Healthy_AD))
pca.info <- fast.prcomp(Healthy_AD)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=substr(rownames(pca.info$x),1,1), pca.info$x)
DaparsSample <- read.csv("DaparsSample.csv",header = T)
BackgroundInformation <- read.csv("BackgroundInformationAPOE.csv",header = T)
colnames(DaparsSample)[2:3] <- c("GSM_number","sample")
a <- merge(DaparsSample,BackgroundInformation,by="GSM_number")
b <- merge(a,pca.data,by="sample")
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+geom_point()
BackgroundInformation <- subset(BackgroundInformation,BackgroundInformation$Brain_Region=="3")
a <- adonis2(as.matrix(Healthy_AD) ~ Severity+CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch,data = BackgroundInformation,permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(8,9),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")+geom_text(aes(label=Pr))
back <- merge(DaparsSample,BackgroundInformation,by="GSM_number")
back <- arrange(back,back$sample)
Healthy_AD <- arrange(Healthy_AD,rownames(Healthy_AD))
group <- factor(substr(rownames(Healthy_AD),1,1))
design <- model.matrix(~group)
colnames(design) <- levels(group)
Healthy_AD <- as.data.frame(t(Healthy_AD))
RemoveData <- removeBatchEffect(Healthy_AD,batch = batch,design = design)
RemoveData <- as.data.frame(t(RemoveData))
RemoveData <- as.matrix(RemoveData)
a <- adonis2(as.matrix(RemoveData) ~ Severity+CERAD+Braak_Stage+Age+Sex+Bank_Location+Batch,data = back,permutations=999, by="margin",method = "euclidean")
b <- as.data.frame(a)
b <- b[-c(8,9),]
b$category <- rownames(b)
colnames(b)[5] <- "Pr"
b <- b[order(b$R2,decreasing = TRUE),]
b$category <- factor(b$category)
p <- ggplot(b,aes(x=reorder(category,-R2),y=R2))+geom_bar(stat="identity")+geom_text(aes(label=Pr))
RemoveData <- as.data.frame(RemoveData)
pca.info <- fast.prcomp(RemoveData)
pca.data <- data.frame(sample = rownames(pca.info$x), Type=gsub(pattern = '\\_\\d*',replacement = '',x=rownames(pca.info$x)), pca.info$x)
b <- merge(back,pca.data,by="sample")
p <- ggplot(data=b,aes(x=PC1,y=PC2,color=Group,shape =Batch))+geom_point()
RemoveData <- as.data.frame(t(RemoveData))
RemoveData <- RemoveData[,order(colnames(RemoveData))]
RemoveData$groupAmeanPDUI <- apply(RemoveData[,1:97],1,mean)
RemoveData$groupBmeanPDUI <- apply(RemoveData[,98:354],1,mean)
RemoveData$PDUI_diff <- RemoveData$groupBmeanPDUI-RemoveData$groupAmeanPDUI
RemoveData$log2FoldChangePDUI <- log2(RemoveData$groupBmeanPDUI/RemoveData$groupAmeanPDUI)
RemoveData$Pvalue <- apply(RemoveData,1,function(x) wilcox.test(x[1:97],x[98:354])$p.value)
RemoveData <- na.omit(RemoveData)
RemoveData$Padjust <- p.adjust(RemoveData$Pvalue,method = "BH")
sessionInfo()
