######################################################################################################################
#####' STAT 6358 REAL DATA ANALYSIS - COMPARE CAMP WITH DESEQ2, ALDEx2, and edgeR ####################################
######################################################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(camp)
library(edgeR)
library(DESeq2)
library(ALDEx2)

## set up the airway data,
se <- airway
XX <- assay(se)
coldata <- colData(se)
grp <- coldata$dex # as factor
grp.vec <- as.vector(grp) # as vector

#### CAMP ####
camp.pval <- camp(t(XX), grp.vec)
table(p.adjust(camp.pval, method = 'BH') < 0.05) 
camp.DEGs <- rownames(XX[which(p.adjust(camp.pval, method = 'BH') < 0.05),])

#### edgeR ####
yy <- DGEList(counts = XX, group = grp.vec)
keep <- filterByExpr(yy)
yy <- yy[keep,]
yy <- calcNormFactors(yy)
yy <- estimateDisp(yy)
fit <- glmQLFit(yy, design = model.matrix(~ grp.vec))
edgeR.res <- glmQLFTest(fit)
edgeR.pval <- p.adjust(edgeR.res$table$PValue, method = 'BH')
table(edgeR.pval < 0.05)
edgeR.DEGs <- rownames(subset(edgeR.res$table, p.adjust(PValue, method = 'BH') < 0.05))

#### DESeq2 ####
dds <- DESeqDataSetFromMatrix(countData = XX, colData = coldata, 
                              design = ~ cell + dex)
dds <- DESeq(dds)
deseq.res <- results(dds)
deseq.pvals <- deseq.res$padj
table(deseq.pvals < 0.05)
deseq.DEGs <- rownames(subset(deseq.res, padj < 0.05))

#### ALDEx2 ####
aldex.res <- aldex(XX, grp.vec)
aldex.pval <- aldex.res$we.eBH 
table(aldex.pval < 0.05) 
aldex.DEGs <- rownames(subset(aldex.res, we.eBH < 0.05))



