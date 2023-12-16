
#' Two-sample t-test for gene expressions between two groups
#'
#' @param counts gene expression matrix with gene rows and cell columns
#' @param cellgroup1 the index of cells from group1
#' @param cellgroup2 the index of cells from group2
#' @return A dataframe with t-score, p-value, and BH p-value for each gene in a row
#' @examples
#' data(Bcells_sce)
#' simple_mean_DE(counts(Bcells_sce), Bcells_sce$stim == "ctrl", Bcells_sce$stim == "stim")
#' @export
simple_mean_DE = function(counts, cellgroup1, cellgroup2){
  tmpcount1 = counts[,cellgroup1]
  tmpcount2 = counts[,cellgroup2]
  if(length(cellgroup1) == 1){
    tmpmean1 = mean(tmpcount1)
  }else{
    tmpmean1 = Matrix::rowMeans(tmpcount1)
  }
  if(length(cellgroup2) == 1){
    tmpmean2 = mean(tmpcount2)
  }else{
    tmpmean2 = Matrix::rowMeans(tmpcount2)
  }
  df = data.frame(genes = rownames(counts))
  df$t = sapply(1:nrow(counts), FUN = function(i){stats::t.test(tmpcount1[i,], tmpcount2[i,])$statistic})
  df$pval = 2*stats::pt(-abs(df$t), df = length(cellgroup1) + length(cellgroup2) -1)
  df$BHpval = stats::p.adjust(df$pval, method = "BH")
  return(df)
}

#' An implementation Poisson glmm DE method.
#'
#' @param sce a SingleCellExperiment object with raw count matrix.
#' @param cellgroups a vector of group labels for cells. The cells consist of two different groups.
#' @param repgroups a vector of donor labels for cells.
#' @param freq_expressed a threshold for gene detection rate.
#' @return A dataframe of Poisson glmm DE results with each row for a gene.
#' @examples
#' data(Bcells_sce)
#' poisson_glmm_DE(Bcells_sce, Bcells_sce$stim, Bcells_sce$ind)
#' @export
poisson_glmm_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  countdf = data.frame(cellgroups = as.factor(cellgroups),
                       repgroups = as.factor(repgroups))
  pval = rep(NA,nrow(sce@assays@data$counts))
  mu = rep(NA,nrow(sce@assays@data$counts))
  beta_cellgroups = rep(NA,nrow(sce@assays@data$counts))
  status = rep("done", nrow(sce@assays@data$counts))
  res_square = rep(NA,nrow(sce@assays@data$counts))
  # REvariation = rep(NA,nrow(sce@assays@data$counts))
  # FEvariation = rep(NA,nrow(sce@assays@data$counts))
  # RESvariation = rep(NA,nrow(sce@assays@data$counts))

  for(i in 1:nrow(sce@assays@data$counts)){
    countdf$count = round(pmax(sce@assays@data$counts[i,],0))
    if (mean(countdf$count!=0, na.rm = TRUE) <= freq_expressed) {
      if(mean(countdf$count!=0, na.rm = TRUE) == 0){
        status[i] = "zero mean"
        next
      }else{
        status[i] = "lowly expressed"
        next
      }
    }
    gm = tryCatch(summary(MASS::glmmPQL(count~cellgroups, random = ~1|repgroups, family = stats::poisson, data = countdf, verbose = FALSE)),
                  error = function(e){NULL})
    if (is.null(gm)){
      status[i] = "not converge"
      next
    }
    gm_null = tryCatch(summary(MASS::glmmPQL(count~1, random = ~1|repgroups, family = stats::poisson, data = countdf, verbose = FALSE)),
                  error = function(e){NULL})
    pval[i] = gm$tTable[2 ,"p-value"]
    res_square[i] = gm$sigma^2
    mu[i] = gm$coefficients$fixed[1]
    beta_cellgroups[i] = gm$coefficients$fixed[2]
    # rsquared = tryCatch(r.squaredGLMM(gm, gm_null, pj2014 = T), error = function(e){NULL})
    # if (!is.null(rsquared)){
    #   REvariation[i] = rsquared[1,2] - rsquared[1,1]
    #   FEvariation[i] = rsquared[1,1]
    #   RESvariation[i] = 1-rsquared[1,2]
    # }
  }
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$mu = mu
  df$beta_cellgroups = beta_cellgroups
  df$log2FC = log2(exp(beta_cellgroups))
  df$sigma_square = res_square
  df$pval = pval
  df$status = status
  df$BH = stats::p.adjust(df$pval, method = "BH")
  # df$REvariation = REvariation
  # df$FEvariation = FEvariation
  # df$RESvariation = RESvariation
  return(df)
}

#' An implementation of Binomial glmm DE method.
#'
#' @param sce a SingleCellExperiment object with raw count matrix.
#' @param cellgroups a vector of group labels for cells. The cells consist of two different groups.
#' @param repgroups a vector of donor labels for cells.
#' @param freq_expressed a threshold for gene detection rate.
#' @return A dataframe of Binomial glmm DE results with each row for a gene.
#' @examples
#' data(Bcells_sce)
#' binomial_glmm_DE(Bcells_sce, Bcells_sce$stim, Bcells_sce$ind)
#' @export
binomial_glmm_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  countdf = data.frame(cellgroups = as.factor(cellgroups), repgroups = as.factor(repgroups))
  pval = rep(NA,nrow(sce@assays@data$counts))
  beta_cellgroups = rep(NA,nrow(sce@assays@data$counts))
  mu = rep(NA,nrow(sce@assays@data$counts))
  status = rep("done", nrow(sce@assays@data$counts))
  res_square = rep(NA,nrow(sce@assays@data$counts))
  for(i in 1:nrow(sce@assays@data$counts)){
    countdf$count = 1*(sce@assays@data$counts[i,]>0)
    if (mean(countdf$count) <= freq_expressed) {
      if(mean(countdf$count) == 0){
        status[i] = "zero mean"
        next
      }else{
        status[i] = "lowly expressed"
        next
      }
    }
    genemean_bygroup = stats::aggregate(count ~ cellgroups, data = countdf, FUN = mean)
    if (genemean_bygroup$count[1] == genemean_bygroup$count[2]){
      status[i] = "no difference between groups"
      next
    }
    gm = tryCatch(summary(MASS::glmmPQL(count~cellgroups, random = ~1|repgroups, family = stats::binomial, data = countdf, verbose = FALSE, niter = 50)),
                  error = function(e){NULL})
    if (is.null(gm)){
      status[i] = "not converge"
      next
    }
    pval[i] = gm$tTable[2 ,"p-value"]
    beta_cellgroups[i] = gm$coefficients$fixed[2]
    mu[i] = gm$coefficients$fixed[1]
    res_square[i] = gm$sigma^2
  }
  df = data.frame(genes = rownames(sce@assays@data$counts))
  df$beta_cellgroups = beta_cellgroups
  df$log2FC = log2(exp(beta_cellgroups))
  df$mu = mu
  df$sigma_square = res_square
  df$pval = pval
  df$status = status
  df$BH = stats::p.adjust(df$pval, method = "BH")
  return(df)
}

# pseudobulk_deseq2 = function(sce, cellgroups, repgroups){
#   raw_count = melt(t(round(pmax(sce@assays@data$counts,0))), value.name = "count")
#   colnames(raw_count) = c("cell", "gene", "count")
#   raw_count$cellgroups = rep(as.factor(cellgroups), nrow(raw_count)/length(cellgroups))
#   raw_count$repgroups = rep(as.factor(repgroups), nrow(raw_count)/length(repgroups))
#   pseudobulk_count = aggregate(count ~ repgroups + cellgroups + gene , data = raw_count, FUN = sum)
#   coldata = unique(pseudobulk_count[c("repgroups","cellgroups")])
#   pseudobulk_count$cell_rep = paste0(pseudobulk_count$cellgroups, "_", pseudobulk_count$repgroups)
#   pseudobulk_count = reshape(pseudobulk_count[,c("cell_rep","gene", "count")], idvar = "gene", timevar = c("cell_rep"), direction = "wide")
#   rownames(pseudobulk_count) = pseudobulk_count[,1]
#   pseudobulk_count[,1] = NULL
#   dds = DESeqDataSetFromMatrix(countData = pseudobulk_count,
#                                colData = coldata,
#                                design= ~ repgroups + cellgroups)
#   dds = DESeq(dds)
#   df = results(dds)[c("pvalue", "padj", "baseMean", "log2FoldChange")]
#   colnames(df) = c("pval", "BH", "basemean", "log2FC")
#   return(df)
# }

#' An implementation of MAST DE method.
#'
#' @param sce a SingleCellExperiment object with raw count matrix. The CPM counts are computed in the funciton.
#' @param cellgroups a vector of group labels for cells. The cells consist of two different groups.
#' @param repgroups a vector of donor labels for cells.
#' @param freq_expressed a threshold for gene detection rate.
#' @return A dataframe of MAST DE results with each row for a gene.
#' @examples
#' data(Bcells_sce)
#' MAST_DE(Bcells_sce, Bcells_sce$stim, Bcells_sce$ind)
#' @export
MAST_DE = function(sce, cellgroups, repgroups, freq_expressed = 0.05){
  sca = MAST::FromMatrix(log2(edgeR::cpm(round(pmax(sce@assays@data$counts,0)))+1),
                      data.frame(cellgroups = as.factor(cellgroups), repgroups = as.factor(repgroups)),
                      data.frame(gene = rownames(sce@assays@data$counts)))
  expressed_genes = MAST::freq(sca) > freq_expressed
  sca = sca[expressed_genes,]
  SummarizedExperiment::colData(sca)$cdr = colSums(SummarizedExperiment::assay(sca)>0)/dim(sca)[1]
  zlmCellgroups = MAST::zlm( ~ cellgroups + cdr + repgroups, sca)
  lrt = MAST::lrTest(zlmCellgroups, "cellgroups")
  df = data.frame(genes = rownames(sca))
  df$coef_cellgroups = zlmCellgroups@coefC[,2]
  df$log2FC = log2(exp(df$coef_cellgroups))
  df = merge(df, data.frame(genes = names(lrt[,3,3]), pval = lrt[,3,3]), by = "genes", all.x = TRUE)
  df$BH = stats::p.adjust(df$pval, method = "BH")
  return(df)
}
#' Identify DEGs for a list of genes after performing DE analysis.
#'
#' An old framework select genes with adjusted p-values smaller than
#' a threshold and absolute log2 fold change greater than a threshold.
#' A new framework filters out genes with small average log2 gene means,
#' but genes showing large difference in mean would be considered as a
#' candidate for DEGs.
#'
#' @param adj_pval a vector of adjusted p-values obtained from a DE analysis
#' @param log2FC a vector of log2 fold change obtained from a DE analysis
#' @param log2mean a vector of log2(genemean1*genemean2)/2 with genemean1 and genemean2 representing the gene mean from raw counts
#' @param log2meandiff a vector of log2(abs(genemean1-genemean2)) with genemean1 and genemean2 representing the gene mean from raw counts
#' @param pvalcutoff the p-value threshold to determine DEGs
#' @param log2FCcutoff the log2 fold change threshold to determine DEGs
#' @param log2meancutoff the log2 mean threshold to determine DEGs
#' @param log2meandiffcutoff the log2 difference of gene mean threshold to determine DEGs
#' @param newcriteria logical. Whether the gene mean and difference of mean should be included in the criteria
#' @return A logical vector indicating DEGs
#' @examples
#' identifyDEGs(runif(1000), rnorm(1000))
#' @export
identifyDEGs = function(adj_pval, log2FC, log2mean = NA, log2meandiff = -Inf,
                        pvalcutoff = 0.05, log2FCcutoff = log2(1.5),
                        log2meancutoff = -2.25, log2meandiffcutoff = -1, newcriteria = F){
  if(newcriteria){
    DEGs = adj_pval<pvalcutoff & abs(log2FC)>log2FCcutoff & (log2mean > log2meancutoff | log2meandiff > log2meandiffcutoff)
  }else{
    DEGs = adj_pval<pvalcutoff & abs(log2FC)>log2FCcutoff
  }
  DEGs = ifelse(is.na(adj_pval), NA, DEGs)
  return(DEGs)
}

