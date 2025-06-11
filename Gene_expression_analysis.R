### All microarray datasets were analyzed separately using the following codes with variation in the datasets GEO accession ID ###

install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
library(forcats)
library(stringr)
library(ggplot2)
library(ggrepel)
library(readr)
library(tidyr)
library(survminer)
library(GEOquery)
library(limma)
library(pheatmap)

### Differential expression analysis of microarray datasets using Limma ###

my_id <- "GSE113513"
gse <- getGEO(my_id)
length(gse)
gse <- gse[[1]]
gse
pData(gse)
fData(gse)
exprs(gse)
summary(exprs(gse))
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)
library(dplyr)
sampleInfo <- pData(gse)
sampleInfo
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1)
sampleInfo <- rename(sampleInfo,patient = source_name_ch1, group=characteristics_ch1)
sampleInfo
library(pheatmap)
## argument use="c" stops an error if there are any missing data points
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 
rownames(sampleInfo)
colnames(corMatrix)
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)   
library(ggplot2)
library(ggrepel)
library(PCAtools)
pca <- prcomp(t(exprs(gse)))
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", characteristics_ch1))) + geom_point() + geom_text_repel()
features <- read.csv("GPL15207-17536.csv")
View(features)
library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write_csv(full_output, path="gse_full_output1.csv")

features <- select(features,Gene_Symbol,Gene_Title,Entrez_Gene)
library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
colnames(design) <- c("Tumor","Normal")
design
summary(exprs(gse))
#cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
#is_expressed <- exprs(gse) > cutoff
#View(is_expressed)
### Identify genes expressed in more than 2 samples ###
#keep <- rowSums(is_expressed) > 2
## check how many genes are removed / retained.
#table(keep)
## subset to just those expressed genes
#gse <- gse[keep,]
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)
contrasts <- makeContrasts(Normal - Tumor, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable(fit2, coef=1)
decideTests(fit2)
table(decideTests(fit2))
aw <- arrayWeights(exprs(gse),design)
aw
fit <- lmFit(exprs(gse), design,
             weights = aw)
contrasts <- makeContrasts(Normal - Tumor, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
anno <- read.csv("GPL15207-17536.csv")
anno
anno <- select(anno,Gene_Symbol,Gene_Title,Entrez_Gene)
fit2$genes <- anno
topTable(fit2)
full_results <- topTable(fit2, number=Inf)
full_results
full_results <- tibble::rownames_to_column(full_results,"ID")
full_results
library(ggplot2)
ggplot(full_results,aes(x = logFC, y= B)) + geom_point()
write_csv(full_results, path="gse_full_resultsanno1.csv")
