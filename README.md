# LEMUR

**LEMUR** is a tool that performs generalized linear mixed model based differential expression analysis on single-cell data with multiple samples. Unlike existing DE methods, we suggest to use raw UMI counts without pre-normalization on data and any normalization technique in the model. We provided two models, Poisson-glmm and Binomial-glmm, modeling UMI counts and zero proportion respectively.

## Installation

**LEMUR** can be installed from github directly as follows:
```{r}
devtools::install_github("C-HW/LEMUR")
```

## Read Data
Read in an example data set from Kang et al. 2018. This dataset comprised 10X droplet-based scRNA-seq PBCM data from eight Lupus patients obtained before and after 6h-treatment with IFN-β. 
```{r}
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "Kang")
sce = eh[["EH2259"]]

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
```
After removing undetected and lowly expressed genes, the dataset consisted of 29065 cells and 7661 genes. 

For this vignette, we focus on a cell type, B cells.
```{r}
celltype = "B cells"
sub_sce = sce[, sce$cell %in% celltype]
```

## DE analysis

```{r}
pois_glmm_df = poisson_glmm_DE(sub_sce, cellgroups = sub_sce$stim, repgroups = sub_sce$ind)
  Kang_pois_glmm_df[[key]]$hits = identifyhits(Kang_pois_glmm_df[[key]]$BH, Kang_pois_glmm_df[[key]]$log2FC, subgroupsce@metadata$log2mean, subgroupsce@metadata$log2meandiff, newcriteria = T, log2FCcutoff = 1)
```

## Author
Chih-Hsuan Wu (U Chicago)    
Bug report, comments or questions please send to chihhsuan@uchicago.edu.
