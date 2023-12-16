# LEMUR

**LEMUR** is a tool that performs generalized linear mixed model based differential expression analysis on single-cell data with multiple samples. Unlike existing DE methods, we suggest to use raw UMI counts without pre-normalization on data and any normalization technique in the model. We provided two models, Poisson-glmm and Binomial-glmm, modeling UMI counts and zero proportion respectively.

In this vignette, we will demonstrate how to install and use this package.

## Installation

**LEMUR** can be installed from github directly as follows:
```{r}
devtools::install_github("C-HW/LEMUR")
```

## Read Data

We provided a subset of a `SingleCellExperiment` object from Kang et al. 2018 in the package. The original dataset comprised 10X droplet-based scRNA-seq PBCM data from eight Lupus patients obtained before and after 6h-treatment with IFN-β. The subset `Bcells_sce` contains only B cells with 1488 cells in control group and 1392 cells in stimulated group. 

The original dataset is available via `ExperimentHub`.

```{r}
# Download data from ExperimentHub
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

```{r}
table(sce$cell, sce$stim)
```

For this vignette, we focus on a cell type, B cells. These cells were contributed from eight different individuals.
```{r}
celltype = "B cells"
Bcells_sce = sce[, sce$cell %in% celltype]
table(Bcells_sce$ind)
```

`Bcells_sce` is also accessible with `data(Bcells_sce)`.


## DE analysis

For this dataset, the goal of DE analysis is to determine the DEGs between two different conditions (groups). The function `poisson_glmm_DE` and `binomial_glmm_DE` are designed to perform DE analysis.

```{r}
pois_glmm_df = poisson_glmm_DE(Bcells_sce, cellgroups = Bcells_sce$stim, repgroups = Bcells_sce$ind)
binomial_glmm_df = binomial_glmm_DE(Bcells_sce, cellgroups = Bcells_sce$stim, repgroups = Bcells_sce$ind)
```

## Determine DEGs

To determine DEGs, we implemented `identifyDEGs` to perform either conventional criteria or our new criteria. The conventional one selects genes satisfying its adjusted p-value passes 0.05 (default) and its absolute value of log2 fold change passes log2(1.5) (default).

```{r}
# conventional criteria
pois_glmm_df$conv_DEGs = identifyDEGs(pois_glmm_df$BH, pois_glmm_df$log2FC, log2FCcutoff = 1)
```

We proposed new criteria that are based on the convention plus the gene mean and the difference in mean. If the log2 gene mean in two groups is lower than a certain value (-2.25 by default) and the log2 mean difference is smaller than a threshold (-1 by default), the gene would not be considered as a DEG. These can also be used as a filter before any DE analysis to speed up the computation. Both criteria are adjustable, depending on the dataset's performance and characteristics.

```{r}
# new criteria
pois_glmm_df$new_DEGs = identifyDEGs(pois_glmm_df$BH, pois_glmm_df$log2FC, pois_glmm_df$log2mean, pois_glmm_df$log2meandiff, newcriteria = T, log2FCcutoff = 1)
```

The volcano plots below demonstrate the DEGs from different criteria.
```{r}
library(ggplot2)
pvalcutoff = 0.05
log2FCcutoff = 1
ggplot(pois_glmm_df, aes(x = log2FC, y = -log10(BH), colour = conv_DEGs)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_hline(yintercept=-log10(pvalcutoff),linetype=2) +
    geom_vline(xintercept=log2FCcutoff,linetype=2) + 
    geom_vline(xintercept=-log2FCcutoff,linetype=2) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    scale_color_manual(values = c("gray", "blue")) +
    xlim(c(-3,3)) + 
    ylim(c(0,50)) + 
    xlab("Log2 Fold Change") + 
    ylab("-Log10 (adj.pval)")
    
ggplot(pois_glmm_df, aes(x = log2FC, y = -log10(BH), colour = new_DEGs)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_hline(yintercept=-log10(pvalcutoff),linetype=2) +
    geom_vline(xintercept=log2FCcutoff,linetype=2) + 
    geom_vline(xintercept=-log2FCcutoff,linetype=2) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    scale_color_manual(values = c("gray", "blue")) +
    xlim(c(-3,3)) + 
    ylim(c(0,50)) + 
    xlab("Log2 Fold Change") + 
    ylab("-Log10 (adj.pval)")
```

We can compare the DEGs from different criteria via heatmaps. It shows that our new criteria largely exclude out the noise.
```{r}
# conventional criteria
library(pheatmap)
sort_log2FC = sort(abs(pois_glmm_df$log2FC[pois_glmm_df$conv_DEGs]), index.return = T, decreasing = T)
sorted_DEGs = pois_glmm_df$genes[which(pois_glmm_df$conv_DEGs)][sort_log2FC$ix]
mat = counts(Bcells_sce)[sorted_DEGs,]

#row annotation
annotation_log2FC = data.frame(abs_log2FC = sort_log2FC$x)
rownames(annotation_log2FC) = rownames(mat)

#col annotation
annotation_df = data.frame(Donor = paste0("D",Bcells_sce$ind), Group = Bcells_sce$stim)
rownames(annotation_df) = colnames(mat)
annotation_df = annotation_df[with(annotation_df, order(Group, Donor)), ]

pheatmap(mat[,rownames(annotation_df)], 
         main = paste0("Poisson-glmm DEGs\nin UMI counts (", nrow(mat)," DEGs)"),
         annotation_col = annotation_df, 
         annotation_row = annotation_log2FC,
         cluster_rows=F, cluster_cols=F, 
         show_colnames = F, show_rownames = F,
         annotation_names_row = F,
         color=colorRampPalette(c("navy", "white", "red"))(10), breaks = seq(0,3, length.out = 10)
         )
```

```{r}
sort_log2FC = sort(abs(pois_glmm_df$log2FC[pois_glmm_df$new_DEGs]), index.return = T, decreasing = T)
sorted_DEGs = pois_glmm_df$genes[which(pois_glmm_df$new_DEGs)][sort_log2FC$ix]
mat = counts(Bcells_sce)[sorted_DEGs,]

#row annotation
annotation_log2FC = data.frame(abs_log2FC = sort_log2FC$x)
rownames(annotation_log2FC) = rownames(mat)

#col annotation
annotation_df = data.frame(Donor = paste0("D",Bcells_sce$ind), Group = Bcells_sce$stim)
rownames(annotation_df) = colnames(mat)
annotation_df = annotation_df[with(annotation_df, order(Group, Donor)), ]

pheatmap(mat[,rownames(annotation_df)], 
         main = paste0("Poisson-glmm DEGs\nin UMI counts (", nrow(mat)," DEGs)"),
         annotation_col = annotation_df, 
         annotation_row = annotation_log2FC,
         cluster_rows=F, cluster_cols=F, 
         show_colnames = F, show_rownames = F,
         annotation_names_row = F,
         color=colorRampPalette(c("navy", "white", "red"))(10), breaks = seq(0,3, length.out = 10)
         )
```
## Author
Chih-Hsuan Wu (U Chicago)    
Bug report, comments or questions please send to chihhsuan@uchicago.edu.
