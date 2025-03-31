## Can we treat zeros in RNA-seq data as censored values?

RNA-seq data is very sparse and an ongoing challenge in this field is the treatment of zeros.
CAMP (Chan and Li 2024) is a recently developed method that treats zeros as censored, enabling the use of survival analysis methods for differential abundance analysis (DAA). Although CAMP was developed for microbiome data, the authors note that the framework can be applied to RNA-seq data due to its sparse nature.

Essentially, this project aims to examine how results differ when zeros are treated as censored values compared to more traditional RNA-seq methods.
Some questions we would like to answer

  1) Does treating zeros as censored change the differentially expressed genes (DEGs) results meaningfully?
  2) Do the unique DEGs detected by the CAMP framework make sense biologically?
  3) Under what scenarios is there a lot of/no overlap between CAMP and standard methods (i.e., when do they agree/disagree)?
 


## Abstract 
RNA-seq data is sparse and an ongoing challenge in this field is the choice of zero-handling method when the goal is to test for differential expression. 
CAMP (@chan2024zero) is a recently developed method that treats zero counts as censored values, facilitating the use of survival analysis methods for differential abundance analysis (DAA). 
The goal of this project is to assess the feasibility of treating zero counts as censored values in RNA-seq data, in contrast to traditional methods which remove lowly-expressed genes or add a pseudocount to zeros. 
In particular, we aim to investigate the impact on differential analysis of RNA-seq data when zeros are treated as censored by exploring the following interesting questions. 
1. Does treating zeros as censored change the differentially expressed genes (DEGs) results meaningfully? 
2. Do the unique DEGs detected by the CAMP framework make sense biologically? 
3. Under what conditions (for example, sample size or sparsity) is there a large overlap between CAMP and more-traditional RNA-seq methods, and equivalently, when do they disagree? 
By answering these questions, we aim to determine if treating zeros as censored is a valid approach that could, in-turn, open up new ideas to existing methodology by including the idea that a zero-count is partial information. 
Although CAMP was developed for microbiome data, the authors suggest that it could extend to RNA-seq data due to the similarity in sparsity.

This project will be uploaded to <https://github.com/cjd04/stat6358project>. Ideally, it will include the simulation and real data analysis code, the simulated datasets, references, and whatever else that might be relevant (if space/memory limits allow, particularly for the datasets).
 
 