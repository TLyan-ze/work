#更换路径
getwd()
setwd("c:/Users/47423/Desktop/")
#合并表达矩阵
options(stringsAsFactors = FALSE)
c1<-read.table("2.count",sep = "\t",col.names = c("gene_id","c1"))
k1<-read.table("3.count",sep = "\t",col.names = c("gene_id","k1"))
k2<-read.table("4.count",sep = "\t",col.names = c("gene_id","k2"))
c2<-read.table("5.count",sep = "\t",col.names = c("gene_id","c2"))
raw_count <- merge(merge(c1, c2, by="gene_id"), merge(k1, k2, by="gene_id"))
head(raw_count)
raw_count_filt <- raw_count[-1:-5,]
head(raw_count_filt)
ENSEMBL <- gsub("\\.\\d*", "", raw_count_filt$gene_id)
head(ENSEMBL)
row.names(raw_count_filt) <- ENSEMBL
head(raw_count_filt)
raw_count_filt <- cbind(ENSEMBL,raw_count_filt)
head(raw_count_filt)
colnames(raw_count_filt)[1] <- c("ensembl_gene_id")
#基因进行注释，获取gene_symbol
library('biomaRt')
library('curl')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
my_ensembl_gene_id <- row.names(raw_count_filt)
options(timeout = 4000000)
hg_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"chromosome_name", "start_position","end_position", "band"), filters= 'ensembl_gene_id', values = my_ensembl_gene_id, mart = ensembl)
head(hg_symbols)
#整合合并并后的表达矩阵和注释文件
readcount <-merge(raw_count_filt,hg_symbols, by="ensembl_gene_id")
head(readcount)
#输出count文件
write.csv(readcount,file = 'readcount_all,csv')
readcount<-raw_count_filt[ ,-1:-2]
write.csv(readcount,file='readcount.csv')
head(readcount)
# DEseq2
mycounts <- read.csv("readcount.csv")
head(mycounts)
#把多的x删除
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
#构建condition
condition <- factor(c(rep("control",2),rep("treat",2)),levels =c("control","treat"))
condition
colData <-data.frame(row.names = colnames(mycounts),condition)
colData
#构建dds对象，开始DESeq
library(DESeq2)
dds <- DESeqDataSetFromMatrix(mycounts,colData,design= ~ condition)
dds <- DESeq(dds)
dds
#查看结果
res = results(dds,contrast = c("condition","control","treat"))
res = res[order(res$pvalue),]
head(res)
summary(res)
q#输出结果
write.csv(res,file="all_results.csv")
table(res$padj<0.05)
#提取差异表达genes
diff_geng_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_geng_deseq2)
head(diff_geng_deseq2)
write.csv(diff_geng_deseq2,file= "DEG_C_V_T.csv")
#注释差异表达基因
library('biomaRt')
library('curl')
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
my_ensembl_gene_id <- row.names(diff_geng_deseq2)
hg_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), filters= 'ensembl_gene_id', values = my_ensembl_gene_id, mart = ensembl)
head(hg_symbols)
#合并数据
ensembl_gene_id<-rownames(diff_geng_deseq2)
diff_geng_deseq2<-cbind(ensembl_gene_id,diff_geng_deseq2)
colnames(diff_geng_deseq2)[1]<-c("ensembl_gene_id")
diff_name<-merge(diff_geng_deseq2,hg_symbols,by="ensembl_gene_id")
head(diff_name)
write.csv(diff_name,file= "DEG_zhu.csv")
#画图
library('ggplot2')
library('pheatmap')
library('vioplot')
library('dplyr')
#主成分分析图
rld <- vst(dds,blind =FALSE)
head(rld)
plotPCA(rld,intgroup = 'condition')
#火山图
for_volcano <- data.frame('log2FoldChange'=res$log2FoldChange,'padj'=res$padj,'trend'=rep('no',length(res$log2FoldChange)))
up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 2),which(for_volcano$padj < 0.05))
down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -2),which(for_volcano$padj < 0.05))
for_volcano[up_sig_indices,'trend'] <- 'up'
for_volcano[down_sig_indices,'trend'] <- 'down'
for_volcano$trend <- as.factor(for_volcano$trend)
for_volcano$padj <- -log10(for_volcano$padj)
p <- ggplot(for_volcano,aes(x=log2FoldChange,y=padj,colour=trend))+
  geom_point(size=I(0.7))+
  scale_color_manual(values = c('no'='blue','up'='red','down'='black'))+
  geom_vline(xintercept = c(2,-2),lty=2,size=I(0.4),colour = 'grey11')+
  geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2,size=I(0.4),colour = 'grey11')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'))+
  labs(x='log2FoldChange',y='-log10Pvalue')
p
#绘制差异表达基因热图
library('pheatmap')
pheatmap(res)
#富集分析
library(clusterProfiler)
library('AnnotationHub')
library('topGO')
library('enrichplot')
ah <- AnnotationHub()
grs <- query(ah,"human")
grs
ah[ah$species == 'Homo sapiens' & ah$rdataclass == 'OrgDb']
org.hs <- ah[['AH79577']]
deseq2.sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
ego <- enrichGO(gene = row.names(deseq2.sig),OrgDb = org.hs,keyType = "ENSEMBL",ont = "ALL",pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2)
head(ego)
dim(ego)
dim(ego[ego$ONTOLOGY=='BP',])
dim(ego[ego$ONTOLOGY=='CC',])
dim(ego[ego$ONTOLOGY=='MF',])
dotplot(ego,font.size=5)
dotplot(ego,showCategory=30)
barplot(ego,showCategory=20, width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,)
ego.BP <- enrichGO(gene = row.names(deseq2.sig), OrgDb = org.hs, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENSEMBL',readable=TRUE)
plotGOgraph(ego.BP)
ego.Bp2 <-simplify(ego.BP,cutoff=0.7,by="p.adjust",select_fun=min)#去掉，调整cutoff
dim(ego.Bp2)
#热图
heatplot(ego.Bp2)
write.csv(ego.BP,file='ego_BP.csv')
library(UpSetR)
library(ggupset)
upsetplot(ego.BP)
goplot(ego.BP)
heatplot(ego.BP, showCategory = 10, foldChange = NULL)
emapplot(ego.BP, showCategory = 30)
cnetplot(ego.BP,showCategory = 5,categorySize="pvalue",circular =TRUE,colorEdge =TRUE)
#GSEA分析
genelist <- deseq2.sig$log2FoldChange
head(genelist)
names(genelist) <- rownames(deseq2.sig)
genelist <- sort(genelist, decreasing = TRUE)
gsemf <- gseGO(genelist,
               OrgDb = org.hs,
               keyType = "ENSEMBL",
               ont="MF")
head(gsemf)
gseaplot(gsemf, geneSetID="GO:0005201")
#gseGO(geneList, ont = "BP", OrgDb, keyType = "ENSEMBL", exponent = 1,nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,pAdjustMethod = "BH", verbose = TRUE, seed = FALSE, by = "fgsea")
#KEGG分析
library(pathview)
gene_list <- mapIds(org.hs, keys = row.names(deseq2.sig),
                    column = "ENTREZID", keytype = "ENSEMBL" )
head(gene_list)
kk <- enrichKEGG(gene_list, organism="hsa",
                 keyType = "ncbi-geneid",
                 pvalueCutoff=0.05, pAdjustMethod="BH",
                 qvalueCutoff=0.1)
head(summary(kk))
barplot(kk, showCategory = 10)
dotplot(kk)
emapplot(kk,showCategory = 30)
cnetplot(kk, showCategory = 5,categorySize="pvalue",circular =TRUE,colorEdge =TRUE)
browseKEGG(kk, "hsa04510")
#DO分析
library('DOSE')
do <- enrichDO(gene = gene_list,ont = "DO",pvalueCutoff = 0.05,qvalueCutoff = 0.5)
barplot(do, showCategory = 10)
### Reactome pathway analysis
library(ReactomePA)
pp <- enrichPathway(gene         = gene,
                    organism     = 'human',
                    pvalueCutoff = 0.05)
head(pp)[,1:6]
pp2 <- gsePathway(geneList     = geneList,
                  organism     = 'human',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
head(pp2)[,1:6]
