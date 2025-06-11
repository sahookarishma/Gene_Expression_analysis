
### Estimation of Upregulated and Downregulated genes###
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

p_cutoff <- 0.01
fc_cutoff <- 1


library(readr)
filter(full_results, P.Value < 0.01, abs(logFC) > 1) %>%
  write_csv(path="filtered_de_results.csv")
full_results

### Visualization ###
library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 50

full_results %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene_Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col='black')
library(pheatmap)
topN <- 20
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)
gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Gene.Symbol) 
gene_matrix <- exprs(gse)[ids_of_interest,]
View(gene_matrix)
pheatmap(gene_matrix,
         labels_row = gene_names)
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row")
### Visualization using volcano plots for differential Upregulated and Downregulated genes###
library(ggplot2)
library(ggrepel)
library(dplyr)

de_genes <-read.csv('gse_full_resultsanno1.csv')
head(de_genes)
colnames(de_genes)

de_genes$diffexpressed <- "NO"
head(de_genes)

dim(de_genes)
de_genes$diffexpressed[de_genes$logFC>0.1 & de_genes$P.Value <0.05] <-"UP"

de_genes$diffexpressed[de_genes$logFC<0.1 & de_genes$P.Value <0.05] <-"DOWN"

head(de_genes)
de_genes$delabel<- NA

ggplot(data=de_genes, aes(x=logFC,y=-log10(P.Value),col=diffexpressed,label=delabel))+ 
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c('blue','black','red'))+
  geom_vline(xintercept=c(-0.8, 0.8), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(text=element_text(size=20))

Trt_up <- de_genes %>% filter(de_genes$diffexpressed == "UP")
Trt_down <- de_genes %>% filter(de_genes$diffexpressed == "DOWN")
lt1 <- Trt_up
lt2 <- Trt_down
write.csv(file = "UPREGULATED_GSE113513.csv", lt1)
write.csv(file = "DOWNREGULATED_GSE113513.csv", lt2)

library(EnhancedVolcano)
EnhancedVolcano(de_genes, lab = de_genes$Symbol, x = "logFC", y="P.Value")
EnhancedVolcano(de_genes, lab = de_genes$Gene_Symbol, x = "logFC", y="P.Value", border = "full", borderWidth = 1.5, borderColour = "black", gridlines.major = FALSE, gridlines.minor = FALSE, title = "CRC versus normal")
View(de_genes)
write.csv(file = "de_genes.csv", de_genes)
anno_genes <-read.csv('de_genes.csv')
