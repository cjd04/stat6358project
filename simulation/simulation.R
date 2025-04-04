#########################################################################################################################
#####' STAT 6358 SIMULATION FOR CAMP METHOD COMPARISON FOR RNA-SEQ DATASETS. COMPARE TO DESEQ2, ALDEx2, and edgeR #######
#####' only lines 148 onwards are needed for simulation results. everything prior to l148 is testing individually #######
#########################################################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("compcodeR")

library(compcodeR)
library(camp)

# setwd("~/Desktop/SMU/Class stuff/2025 Spring/STAT6358 High Throughput Data Analysis/stat6358project/simulation")

#### Generate RNA-seq data using compcodeR package ####
dat <- generateSyntheticData(dataset = "simulated.dat",
                             n.vars = 12500, # initial number of genes in dataset
                             samples.per.cond = 5, # number of samples in each of the two conditions (c1 and c2)
                             n.diffexp = 1250, # number of genes to be DE between the two conditions
                             fraction.upregulated = 0.5, # the fraction of DEGs that is upregulated in c2 compared to c1
                             repl.id = 1, # a replicate ID for the specific simulation instance
                             output.file = "mydata_5spc.rds" # the path to the file where the data object should be saved.
)
## extract count matrix
XX <- dat@count.matrix #; head(XX)

## extract true DEGs, named by row. needed for simulation to examine method performance
true.DEGs <- rownames(subset(dat@variable.annotations, differential.expression==1))

## experimental conditions / groups
grp <- rep(c("Group1", "Group2"), each = 5) # condition

#### useful function for power calculations -------------------------------------------------------- start ####
get_pwr_tbl <- function(pp, true.DEGs, ct.mat, methodname){ 
  #' function that takes in p vector, true DEGs, and count matrix.
  #' returns a matrix of power and type 1 error.
  p.adj <- p.adjust(pp, method = "BH")
  sim.DEGs <- rownames(ct.mat)[p.adj < 0.05]
  tp <- sum(sim.DEGs %in% true.DEGs) # which DEGs are truly detected. true positive
  fp <- length(sim.DEGs) - tp # false positive. found a gene to be DE but the truth is it wasn't
  fn <- sum(!true.DEGs %in% sim.DEGs) # false negative. categorised truly DEGs as non-DE
  
  ## create table
  t1e <- fp / (length(rownames(ct.mat)) - length(true.DEGs))
  pwr <- tp / length(true.DEGs)
  pwr.mat <- matrix(c(pwr, t1e), ncol = 2, 
                    dimnames = list(methodname, c('power','Type1 error')))
  pwr.mat
}
### useful function -------------------------------------------------------------------------------- end

#### CAMP ####
camp.sim1 <- camp(t(XX), grp) # vector of p values

## get power
camp.mat <- get_pwr_tbl(camp.sim1, true.DEGs, XX, 'CAMP')
#' low power good t1e, suggests CAMP is a conservative approach.


#### edgeR ####
library(edgeR)
yy <- DGEList(counts = XX, group = grp) # raw count matrix and groups/condtions

## filter out low-expressed genes. 
keep <- filterByExpr(yy) 
yy <- yy[keep,]
#' this step removes genes that are not expressed enough (across samples) to 
#' extract any meaningful information. the default setting is based on minimum CPM
#' where the cutoff/minimum level is computed by edgeR based on the library size.
#' For simulation, it would be good to know where it cuts off. On the other hand,
#' we should assume the default setting is probably the best, since in real data 
#' analysis we wouldn't know what to tune it to (maybe an field-expert would, but 
#' that's entirely situational and a bit too granular for the purpose of this project).

## normalise (TMM)
yy <- calcNormFactors(yy)
#' default method is fine. it essentially scales the counts to account for composition
#' bias to prevent highly expressed genes from dominating the analysis

## Estimate dispersion
yy <- estimateDisp(yy)
#' edgeR assumes negative binomial because of the dispersion parameter.more ideal than 
#' Poisson because biological data has extra variability, leading to overdispersion.
#' NB is better equipped to handle this extra variance.
#' gene w high dispersion => larger variance between samples in same group => harder to claim DEG
#' gene w low dispersion => smaller variance, small changes in expression could be meaningful!

## plot the biological coefficient of variation (BCV)
plotBCV(yy)
#' this is the sqrt of the dispersion estimate (y) against the average log CPM (x)
#' each black dot is a single gene. 
#' low x (CPM) usually means more variability
#' large y => high variability => hard to detect DE

## check for batch effects / clustering. 
plotMDS(yy)
#' ideally, we see samples from the same group clustering together, indicating a 
#' clear difference between conditions/groups => improve DEG detection.

## perform DE testing
fit <- glmQLFit(yy, design = model.matrix(~grp))
edgeR.res <- glmQLFTest(fit)
#' for each gene separately, we test the null hypothesis H0: the gene is not DE, and that 
#' any observed differences we see in expression are due to biol. variability and tech. error vs.
#' H1: the gene is DE, and the observed differences in expression are from a true change in 
#' expression due to the group / experimental condition.

## get p vals and power 
edgeR.pvals <- edgeR.res$table$PValue
edgeR.mat <- get_pwr_tbl(edgeR.pvals, true.DEGs, XX, 'edgeR')



#### DESeq2 ####
library(DESeq2)
## create object for DESeq2
coldata <- data.frame(
  condition = factor(grp),
  row.names = colnames(XX)  # rows of coldata must match cols of count matrix (i.e. sample 1, 2, 3, ...)
)
dds <- DESeqDataSetFromMatrix(countData = XX, colData = coldata, design = ~ condition)

## visualisation steps
plotDispEsts(dds)
plotMA(res, main="MA Plot (DESeq2)")

## DE testing via DESeq2
dds <- DESeq(dds) # does all the dispersion, normalisation stuff in this function.
res <- results(dds) # extract results.
de.seq.pval <- res$pvalue
## get power 
deseq.mat <- get_pwr_tbl(de.seq.pval, true.DEGs, XX, "DESeq2")



#### ALDEx2 ####
library(ALDEx2)
## from our set up prevoiusly, we can basically run aldex2 off the bat
aldex.res <- aldex(XX, grp)
aldex.pval <- aldex.res$we.ep # extract pvalues

## get power 
aldex.mat <- get_pwr_tbl(aldex.pval, true.DEGs, XX, "ALDEx2")

rbind(camp.mat, edgeR.mat, deseq.mat, aldex.mat)



####### Simulate under different scenarios ########
#' assess power and type 1 error rate for detecting differentially expressed genes in RNA-seq
#' data for different methods. Two methods are traditional reference based, the other two
#' methods assume the data are compositional.


simTest_comparison <- function(n.genes = 12500, n.per.grp = 5, n.trueDEGs = 0.1){
  dat <- generateSyntheticData(dataset = "simulated.dat",
                               n.vars = n.genes, # initial number of genes in dataset
                               samples.per.cond = n.per.grp, # number of samples in each of the two conditions (c1 and c2)
                               n.diffexp = n.trueDEGs * n.genes, # number of genes to be DE between the two conditions
                               fraction.upregulated = 0.5, # the fraction of DEGs that is upregulated in c2 compared to c1
                               repl.id = 1, # a replicate ID for the specific simulation instance
                               output.file = "mydata_5spc.rds" # the path to the file where the data object should be saved.
  )
  ## extract count matrix
  XX <- dat@count.matrix 
  ## extract true DEGs, named by row.
  true.DEGs <- rownames(subset(dat@variable.annotations, differential.expression==1))
  ## experimental conditions / groups
  grp <- rep(c("Group1", "Group2"), each = n.per.grp)
  
  
  #### CAMP ####
  #print("working on CAMP now")
  camp.sim1 <- camp(t(XX), grp) # vector of p values
  camp.mat <- get_pwr_tbl(camp.sim1, true.DEGs, XX, 'CAMP') # pwr matrix
  
  #### edgeR ####
  #print("working on edgeR now")
  yy <- DGEList(counts = XX, group = grp) # raw count matrix and groups/condtions
  ## filter out low-expressed genes. 
  keep <- filterByExpr(yy) 
  yy <- yy[keep,]
  yy <- calcNormFactors(yy)
  ## Estimate dispersion
  yy <- estimateDisp(yy)
  ## perform DE testing
  fit <- glmQLFit(yy, design = model.matrix(~grp))
  edgeR.res <- glmQLFTest(fit)
  ## get p vals and power 
  edgeR.pvals <- edgeR.res$table$PValue
  edgeR.mat <- get_pwr_tbl(edgeR.pvals, true.DEGs, XX, 'edgeR')
  
  #### DESeq2 ####
  #print("working on DESeq2 now")
  ## create object for DESeq2
  coldata <- data.frame(
    condition = factor(grp),
    row.names = colnames(XX)  # rows of coldata must match cols of count matrix (i.e. sample 1, 2, 3, ...)
  )
  dds <- DESeqDataSetFromMatrix(countData = XX, colData = coldata, design = ~ condition)
  ## DE testing via DESeq2
  dds <- DESeq(dds) # does all the dispersion, normalisation stuff in this function.
  res <- results(dds) # extract results.
  de.seq.pval <- res$pvalue
  ## get power 
  deseq.mat <- get_pwr_tbl(de.seq.pval, true.DEGs, XX, "DESeq2")
  
  #### ALDEx2 ####
  #print("working on ALDEx2 now")
  aldex.res <- aldex(XX, grp, mc.samples = 64) # default MC samples is 128.
  aldex.pval <- aldex.res$we.ep # extract pvalues
  ## get power 
  aldex.mat <- get_pwr_tbl(aldex.pval, true.DEGs, XX, "ALDEx2")
  
  return(rbind(camp.mat, edgeR.mat, deseq.mat, aldex.mat))
}

#### simulation settings ####
opt <- list(n.genes = c(10000, 15000, 20000), 
            n.per.grp = c(10, 25, 50), # previously tried c(10,50,100), too large for bigger sim on local. 
                                       #' if necessary we can change this to be larger and run on the cluster
            n.trueDEGs = c(0.05, 0.1))#, # express as a proportion of total genes 
#inflate0 = 1) # we ignore inflate0 for now. The goal is to manually adjust sparsity 
sim.scenarios <- expand.grid(opt)

## Loop, run all code (between sim.results to close(..) at the same time
seed <- 4
set.seed(seed)
Sys.time()

run_sim <- function(i){
  cat('starting simulation ', i, ' of ', nrow(sim.scenarios), '...')
  sim.res <- do.call(simTest_comparison, sim.scenarios[i,])
  return(sim.res)
}
sim.results <- sapply(1:nrow(sim.scenarios),
                      function(i) run_sim(i), simplify = FALSE)
saveRDS(sim.results, "sim_results.rds")
Sys.time() # about 45 mins to run under 10k,15k,20k, 10,25,50


## extract results into a nice table
sim.results.df <- do.call(rbind, sim.results) %>% 
  as.data.frame() %>% 
  mutate(method = rep(c("CAMP", "edgeR", "DESeq2", "ALDEx2"), times = length(sim.results)),
         n.genes = rep(sim.scenarios[1:length(sim.results), 'n.genes'], each = 4),
         n.per.grp = rep(sim.scenarios[1:length(sim.results), 'n.per.grp'], each = 4),
         n.trueDEGs = rep(sim.scenarios[1:length(sim.results), 'n.trueDEGs'], each = 4))
rownames(sim.results.df) <- NULL  # Remove row names
sim.results.df <- sim.results.df %>%
  select(method, everything()) # rearrange columns
print(sim.results.df)
saveRDS(sim.results.df, "sim_results.rds")

library(knitr)
kable(sim.results.df, digits = 3)

library(ggplot2)

ggplot(sim.results.df, aes(x = as.factor(n.per.grp), y = power, 
                           color = method, group = method, linetype = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ n.genes + n.trueDEGs) +
  labs(title = "Power Trends Across Sample Sizes",
       x = "Sample Size per Group",
       y = "Power") +
  theme_minimal() + 
  theme(legend.position = "top")

ggplot(sim.results.df, aes(x = as.factor(n.per.grp), y = "Type1 error", 
                           color = method, group = method, linetype = method)) +
  geom_line() +
  geom_point() +
  facet_grid(~ n.genes + n.trueDEGs) +
  labs(title = "Type 1 error Across Sample Sizes",
       x = "Sample Size per Group",
       y = "Type 1 error") +
  theme_minimal()










