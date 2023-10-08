sumGE <- read.csv("Temporal_Lobe_ADvsHealthy_DEseq2.csv",header = T,row.names = 1)
sumGE <- sumGE[order(sumGE$padj),]
sumGE <- sumGE[1:5000,]
FPKM <- read.csv("TemporalLobeFPKM_RemoveBatch_20220831.csv",header = T,row.names = 1)
a <- gsub(pattern = '\\.\\d*',replacement = '',x=rownames(FPKM))
rownames(FPKM) <- a
df1 <- bitr(a, fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
df1 <- df1[!duplicated(df1$SYMBOL),]
FPKM$ENSEMBL <- rownames(FPKM)
df <- merge(df1,FPKM,by="ENSEMBL",all=FALSE)
rownames(df) <- df$SYMBOL
FPKM <- df[,-c(1:2)]
gene <- rownames(sumGE)
sumGene <- subset(FPKM,rownames(FPKM) %in% gene)
sumGene <- as.data.frame(t(sumGene))
sumGene$GSM_number <- rownames(sumGene)
bg <- read.csv("BackgroundInformationAPOE.csv",header = T)
sum <- merge(bg,sumGene,by="GSM_number")
rownames(sum) <- sum$Sample
sum <- sum[,-c(1:11)]
sumGene <- as.data.frame(t(sum))
AS <- read.csv("VoomRemoveBatchDiff.csv",header = T,check.names = F,row.names = 1)
AS <- AS[order(AS$Padjust),]
AS <- AS[1:5000,]
AS <- AS[,-c(355:360)]
APA <- read.csv("RemoveDataDifference.csv",header = T,check.names = F,row.names = 1)
APA <- APA[order(APA$Padjust),]
APA <- subset(APA,APA$Padjust<0.05)
APA <- APA[,-c(355:360)]
DaparsSample <- read.csv("DaparsSample.csv",header = T)
colnames(DaparsSample) <- c("Group","GSM_number","sample")
BackgroundInformationAPOE <- read.csv("BackgroundInformationAPOE.csv",header = T)
bg <- merge(DaparsSample,BackgroundInformationAPOE,by="GSM_number")
APA <- as.data.frame(t(APA))
APA$sample <- rownames(APA)
APA2 <- merge(bg,APA,by="sample") 
APA <- APA2[,-c(1:10,12:13)]
rownames(APA) <- APA$Sample
APA <- APA[,-1]
APA <- as.data.frame(t(APA))








sumGene <- sumGene[,order(colnames(sumGene))]
AS <- AS[,order(colnames(AS))]
APA <- APA[,order(colnames(APA))]
test <- list(Transcriptome=as.matrix(sumGene),AS=as.matrix(AS),APA=as.matrix(APA))
library(MOFA2)
MOFAobject <- create_mofa(test)
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)
MOFAobject <- prepare_mofa(object = MOFAobject,data_options = data_opts,model_options = model_opts,training_options = train_opts)
outfile <- file.path("AAA","BBB")
MOFAobject.trained <- run_mofa(MOFAobject, outfile,use_basilisk = TRUE)
filepath <- "AAA/BBB"
model <- load_model(filepath)
p <- plot_data_overview(model)
Nsamples <- sum(model@dimensions$N)
bg <- read.csv("BackgroundInformationAPOE.csv",header = T)
bg <- subset(bg,bg$Brain_Region=="3")
bg <- bg[,-c(1,2,10)]
bg <- bg[order(bg$Sample),]
colnames(bg)[7] <- "sample"
samples_metadata(model) <- bg
SampleFactor <- as.data.frame(model@expectations[["Z"]][["group1"]])
SampleFactor$sample <- rownames(SampleFactor)
SampleFactor$Group <- rep("group1",354)
Sample <- merge(bg,SampleFactor,by="sample")
Sample <- na.omit(Sample)
p <- ggplot(Sample,aes(x=Severity,y=Factor1))+geom_signif(comparisons = list(c("0", "1")))+geom_violin()+geom_jitter(aes(color=Severity))
factorSample <- as.data.frame(model@expectations[["Z"]][["group1"]])
factorSample <- factorSample[,1:5]
bg <- read.csv("BackgroundInformationAPOE.csv",header = T)
bg <- subset(bg,bg$Brain_Region=="3")
rownames(bg) <- bg$Sample
bg <- bg[,-c(1,2,8,9,10)]
factorSample <- factorSample[order(rownames(factorSample),decreasing = T),]
bg <- bg[order(rownames(bg),decreasing = T),]
nSamples <- dim(bg)[1]
bg$Phenotypedit <- as.integer(bg$Phenotypedit)
library(WGCNA)
moduleTraitCor <- cor(factorSample, bg,use = 'p',method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(2,2)
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue,1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file = "Module_trait_relationships.pdf",width =4,height = 3)
p <- labeledHeatmap(Matrix = moduleTraitCor,xLabels = colnames(bg),yLabels = colnames(factorSample),ySymbols = colnames(factorSample),textMatrix = textMatrix)
print(p)
dev.off()







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
p <- ggplot(GEweights, aes(x = Weight, y = Rank,colour=Color))+geom_point(aes(color = Color))+geom_text_repel(data = GEweights_Top, aes(x = Weight, y = Rank, label = GeneName))
GEweights_Top$abs <- abs(GEweights_Top$Weight)
p <- ggdotchart(GEweights_Top, x = "GeneName", y = "Weight",sorting = "ascending")
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightGEGene <- KEGGa$GeneName
df1 <- bitr(rownames(KEGGa), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
KEGGgene <- df1$ENTREZID
library(R.utils)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene = KEGGgene,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
p <- barplot(kegg,x = "Count",color = "p.adjust",showCategory = 5)
ASweights <- as.data.frame(model@expectations[["W"]][["AS"]])
ASweights$GeneName <- rownames(ASweights)
ASweights <- separate(ASweights,GeneName,into = c("Chr","Start","End","Cluster"),sep = "([:])")
ASweights$Position <- paste(ASweights$Chr,ASweights$Start,ASweights$End,sep = ":")
ASweights <- ASweights[,c(1,20)]
m <- read.table("humanv40annotation_all_introns.bed",fill = TRUE)
m <- m[,-10]
m$Position <- paste(m$V1,m$V2,m$V3,sep = ":")
sum <- merge(ASweights,m,by="Position")
sum <- sum[,c(2,6,9)]
ASweights <- sum
colnames(ASweights) <- c("Factor1","GeneName","Transcript")
ASweights <- ASweights[!duplicated(ASweights$Factor1),]
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
a <- a$Transcript
ASweights <- rbind(ASweights_Pos,ASweights_Neg)
ASweights <- ASweights[order(ASweights$Factor1,decreasing = F),]
ASweights$Rank <- seq(1,4779,1)
ASweights_Top <- subset(ASweights,ASweights$Transcript %in% a)
topAS <- ASweights_Top$GeneName
ASweights_Top <- ASweights_Top[order(ASweights_Top$Factor1,decreasing = T),]
ASweights_Top <- ASweights_Top[!duplicated(ASweights_Top$GeneName),]
ASweights$Color <- ifelse(ASweights$Transcript %in% a, "#940514","grey")
p <- ggplot(ASweights, aes(x = Weight, y = Rank,colour=Color))+geom_point(aes(color = Color))+geom_text_repel(data = ASweights_Top, aes(x = Weight, y = Rank, label = GeneName))
ASweights_Top$abs <- abs(ASweights_Top$Weight)
p <- ggdotchart(ASweights_Top, x = "GeneName", y = "Weight",sorting = "ascending")
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightASGene <- KEGGa$GeneName
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
KEGGgene <- df1$ENTREZID
kegg <- enrichKEGG(gene = KEGGgene,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
p <- barplot(kegg,x = "Count",color = "p.adjust",showCategory = 20)
APAweights <- as.data.frame(model@expectations[["W"]][["APA"]])
APAweights$ENSEMBLTRANS <- rownames(APAweights)
APAweights <- APAweights[,c(1,16)]
Healthy_AD <- read.table("Dapars2_All_Prediction_Results.txt",header = T)
Healthy_AD <- Healthy_AD[,1:2]
long<-separate(Healthy_AD,Gene,into = c("ENSEMBLTRANS","GeneName","Chr","Strand"),sep = "([|])")
Healthy_AD <- long
Healthy_AD <- Healthy_AD[,1:2]
APAweights <- merge(Healthy_AD,APAweights,by="ENSEMBLTRANS")
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
APAweights$Rank <- seq(1,36,1)
APAweights_Top <- subset(APAweights,APAweights$ENSEMBLTRANS %in% a)
topAPA <- APAweights_Top$GeneName 
APAweights$Color <- ifelse(APAweights$ENSEMBLTRANS %in% a, "#2678C2","grey")
p <- ggplot(APAweights, aes(x = Weight, y = Rank,colour=Color))+geom_point(aes(color = Color))+geom_text_repel(data = APAweights_Top, aes(x = Weight, y = Rank, label = GeneName))
APAweights_Top$abs <- abs(APAweights_Top$Weight)
p <- ggdotchart(APAweights_Top, x = "ENSEMBLTRANS", y = "Weight",sorting = "ascending")
KEGGa <- subset(KEGGa,abs(KEGGa$Weight)>= 0.3)
WeightAPAGene <- KEGGa$GeneName
df1 <- bitr(KEGGa$GeneName, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
KEGGgene <- df1$ENTREZID
kegg <- enrichKEGG(gene = KEGGgene,keyType = "kegg",organism  = 'hsa',pvalueCutoff  = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
p <- barplot(kegg,x = "Count",color = "p.adjust",showCategory = 20)
sigGene <- read.csv("sigGene.csv",header = T)
weightGene <- read.csv("weightGene.csv",header = T)
a <- list(MOFA2 = as.list.data.frame(unique(weightGene$weightGene)),WGCNA = as.list.data.frame(unique(sigGene$sigGene)))
p <- ggvenn(a, show_elements = FALSE, label_sep = "\n")
sessionInfo()
