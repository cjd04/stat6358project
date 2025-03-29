## Can we treat zeros in RNA-seq data as censored values?

RNA-seq data is very sparse and an ongoing challenge in this field is the treatment of zeros.
CAMP (Chan and Li 2024) is a recently developed method that treats zeros as censored, enabling the use of survival analysis methods for differential abundance analysis (DAA). Although CAMP was developed for microbiome data, the authors note that the framework can be applied to RNA-seq data due to its sparse nature.

Essentially, this project aims to examine how results differ when zeros are treated as censored values compared to more traditional RNA-seq methods.
Some questions we would like to answer

  1) Does treating zeros as censored change the differentially expressed genes (DEGs) results meaningfully?
  2) Do the unique DEGs detected by the CAMP framework make sense biologically?
  3) Under what scenarios is there a lot of/no overlap between CAMP and standard methods (i.e., when do they agree/disagree)?
 
## References                                           

  - CAMP paper
  - Silverman Paper
  - other references