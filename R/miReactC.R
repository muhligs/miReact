#' @import magrittr
#' @import sparseMatrixStats
NULL

#' @title Prepare data
#' @description Prepare data for miReactC calculations
#' @param sco List containing raw counts (exp) and run parameters (runparameters)
#' @param seqs.path Path to motif counts (default="./hs.seqXmot.counts.utr3_mrs_7mer.rds")
#' @param seqlist.path Path to motif counts (default="./hs.seqXmot.counts.utr3_mrs_7mer.rds")
#' @param patterns.path Path to motif counts (default="./hs.seqXmot.counts.utr3_mrs_7mer.rds")
#' @param pval.path Path to motif counts (default="./hs.seqXmot.counts.utr3_mrs_7mer.rds")
#' @param counts.path Path to motif counts (default="./hs.seqXmot.counts.utr3_mrs_7mer.rds")
#' @param verbose Show progress (default=T)
#' @return A list with RNA count data, sequences, motif p-values, motif counts, and indexes
#' @export
prepareData <- function(sco,
                        seqs.path = "./seqs/hs.utr3.seqs.rds", 
                        seqlist.path = "./seqs/hs.utr3.seqlist.rds", 
                        patterns.path = "./motif.models/patterns.7mer.Rdata", 
                        pval.path = "./motif.models/hs.seqXmot.utr3_mrs_7mer.rds", 
                        counts.path = "./motif.models/hs.seqXmot.counts.utr3_mrs_7mer.rds", 
                        verbose = T) {
  cat(format(Sys.time()),"Loading data .")
  # addSeqs.R
  seqs <- readRDS(seqs.path)
  
  sco$geneID <- "gid"
  rownames(sco$exp) <- sub("[.].*","",rownames(sco$exp))  
  idx <- rownames(sco$exp) %in% seqs[[sco$geneID]]
  sco$exp <- sco$exp[idx,]
  sco$seqs <- seqs[match(rownames(sco$exp),seqs[[sco$geneID]]),c("sequence","gid","tid","gsym","nchar")]
  sco$seqlist <- readRDS(seqlist.path) %>% .[match(rownames(sco$exp),seqs[[sco$geneID]])]
  names(sco$seqlist) <- rownames(sco$exp)
  if(verbose) cat(".")
  
  # addMotifs.R
  load(patterns.path)
  sco$motifModels <- patterns
  names(sco$motifModels) <- Regmex::all.mers(sco$runparameters$motifs)
  sco$pval.mat <- readRDS(pval.path) %>% .[,sco$seqs$tid]
  if(verbose) cat(". ")
  sco$motifCounts <- readRDS(counts.path) %>% .[,sco$seqs$tid]
  sco$motifCounts[sco$motifCounts>2] <- 2
  
  # addMotifExpression.R
  if(verbose) cat("done!\n",format(Sys.time()),"Calculating medians ... ")
  sco$medianExp <- sparseMatrixStats::rowMedians(sco$exp) # Much faster
  
  if(verbose) {
    cat("done!\n",format(Sys.time()),"Ordering ... \n")
    sco$var <- pbapply::pbapply(sco$exp, 2, function(x) order(x-sco$medianExp, decreasing = TRUE), cl=1)
  } else {
    sco$var <- apply(sco$exp, 2, function(x) order(x-sco$medianExp, decreasing = TRUE))
  }
  
  if(verbose) cat(format(Sys.time()),"All done!")
  
  return(sco)
}

setChunks <- function(chunks, n.cores) {
  lapply(1:chunks, function(x) {
    c(((x-1)*n.cores+1):(x*n.cores))
  })
}

#' @title Calculate miRNA activity using C++
#' @description Calculates miRNA activity using C++
#' @param sco A list containing indexes (var), motif p-values (pval.mat), motif counts (motifCounts), and alpha (runparameters$alpha)
#' @param verbose Show progress (default=T)
#' @return A matrix with the estimated miRNA activities
#' @export
miReactC <- function(sco, n.cores = 1, verbose = T) {
  samples <- ncol(sco$var)
  if(n.cores > samples) n.cores <- samples
  chunks = samples/n.cores
  
  if(chunks %% 1 != 0) {
    tmp <- round(chunks, 0)
    chunks <- ifelse(tmp <= chunks, tmp, tmp - 1)
    chunks %<>% setChunks(n.cores)
    chunks %<>% append(list((chunks %>% unlist() %>% length() + 1):samples))
  } else {
    chunks %<>% setChunks(n.cores)
  }
  
  if(verbose) cat(format(Sys.time()),"Analyzing with",n.cores,"cores in",length(chunks),"chunks ... \n")
  
  res <- sccore:::plapply(chunks, function(c) {
    wcmodCPP(sco$var[,c], sco$pval.mat, sco$motifCounts, sco$runparameters$alpha)
  }, n.cores=1, progress=verbose) %>% do.call(cbind, .)
  
  dimnames(res) <- list(names(sco$motifModels), colnames(sco$var))
  
  res[res == Inf] <- max(res[is.finite(res)])+0.01
  res[res == (-Inf)] <- min(res[is.finite(res)])-0.01
  
  res <- res * -1
  
  if(verbose) cat("\n",format(Sys.time()),"All done!")
  
  return(res)
}