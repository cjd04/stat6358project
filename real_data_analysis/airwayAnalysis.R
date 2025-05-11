######################################################################################################################
#####' STAT 6358 REAL DATA ANALYSIS - COMPARE CAMP WITH DESEQ2, ALDEx2, and edgeR ####################################
######################################################################################################################
# setwd("~/Desktop/SMU/Class stuff/2025 Spring/6358 High Throughput Data Analysis/stat6358project/real_data_analysis")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(camp)
library(edgeR)
library(DESeq2)
library(ALDEx2)
library(airway)
data(airway)


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


#### Real Data analysis 2 ####

BiocManager::install("recount")
library('recount') # browseVignettes('recount')
# browse_study('SRP009615')
# https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=SRP009615
url <- download_study('SRP009615') # https://jhubiostatistics.shinyapps.io/recount/
load(file.path('SRP009615', 'rse_gene.Rdata'))
rse <- scale_counts(rse_gene)
colData(rse)


## below code taken from https://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html
## Find a project of interest
project_info <- abstract_search("GSE32465")

## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

## Load the data
load(file.path(project_info$project, "rse_gene.Rdata"))

## Browse the project at SRA
browse_study(project_info$project)

geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Note that the information for this study is a little inconsistent, so we
## have to fix it.
geochar <- do.call(rbind, lapply(geochar, function(x) {
  if ("cells" %in% colnames(x)) {
    colnames(x)[colnames(x) == "cells"] <- "cell.line"
    return(x)
  } else {
    return(x)
  }
}))

## We can now define some sample information to use
sample_info <- data.frame(
  run = colData(rse_gene)$run,
  group = ifelse(grepl("uninduced", colData(rse_gene)$title), "uninduced", "induced"),
  gene_target = sapply(colData(rse_gene)$title, function(x) {
    strsplit(strsplit(
      x,
      "targeting "
    )[[1]][2], " gene")[[1]][1]
  }),
  cell.line = geochar$cell.line
)

## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

## Add sample information for DE analysis
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target

#### DESEQ2 ####
## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse, ~ gene_target + group)
## Perform DE analysis
dds <- DESeq(dds, test = "LRT", reduced = ~gene_target, fitType = "local")
res <- results(dds)
deseq.pvals <- res$padj
table(deseq.pvals < 0.05)
deseq.DEGs <- rownames(subset(res, padj < 0.05))


#### CAMP ####
camp.pval <- camp(t(assay(rse)), rse$group)
table(p.adjust(camp.pval, method = 'BH') < 0.05) 
camp.DEGs <- rownames(XX[which(p.adjust(camp.pval, method = 'BH') < 0.05),])


#### ALDEx2 ####
aldex.res <- aldex(assay(rse), rse$group)
aldex.pval <- aldex.res$we.eBH 
table(aldex.pval < 0.05) 
aldex.DEGs <- rownames(subset(aldex.res, we.eBH < 0.05))
