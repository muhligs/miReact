####################################################################################################################
# first define the sequences. addSeqs is fine as is.
addSeqs <- function(sco,species="hs",seq.type="utr3"){
  if(!"runparameters" %in% names(sco)) sco$runparameters <- list()
  sco$runparameters <- append(sco$runparameters,list(species=species, seq.type=seq.type))
  seqs <- readRDS(paste0("./seqs/",sco$runparameters$species,".",sco$runparameters$seq.type,".seqs.rds"))
  cat("genes in exp matrix: ",	dim(sco$exp)[1],"\n")
  cat("Samples/cells in exp matrix: ",	dim(sco$exp)[2],"\n")
  # we will streamline so that only gene IDs can be used, and a function can convert between gid, tid and gsym.
  sco$geneID <- "gid"
  cat("Chosen identifier: ",	switch(EXPR = sco$geneID, "gid" = "Ensembl gene ID","tid" = "Ensembl transcript ID", "gsym" = "gene symbol"),".\n") #
  # remove trailing [.]xx from gid/tid
  if(sco$geneID %in% c("gid","tid")) rownames(sco$exp) <- sub("[.].*","",rownames(sco$exp))  
  # report stats
  cat("genes not found: ",	sum(!rownames(sco$exp) %in% seqs[[sco$geneID]]),"\n") #
  cat("genes found: ",	sum(rownames(sco$exp) %in% seqs[[sco$geneID]]),"\n") #
  if(sum(rownames(sco$exp) %in% seqs[[sco$geneID]])<100)stop("Gene IDs should be Ensemble Gene IDs")
  # remove genes without seq
  cat("reducing expression matrix to contain genes with sequence data.\n")
  idx <- rownames(sco$exp) %in% seqs[[sco$geneID]]
  sco$exp <- sco$exp[idx,]
  cat("expression matrix now contains ",dim(sco$exp)[1]," genes with ",length(unique(rownames(sco$exp)))," names.\n")
  cat("adding $seqs to output object.\n")
  sco$seqs <- seqs[match(rownames(sco$exp),seqs[[sco$geneID]]),c("sequence","gid","tid","gsym","nchar")]
  # make a regmex style sequence object
  cat("adding $seqslist (sequence stats) to output object.\n")
  seqlist <- readRDS(paste0("./seqs/",sco$runparameters$species,".",sco$runparameters$seq.type,".seqlist.rds"))
  sco$seqlist <- seqlist[match(rownames(sco$exp),seqs[[sco$geneID]])]
  names(sco$seqlist) <- rownames(sco$exp)
  cat("Finished.\n")
  return(sco)
}
####################################################################################################################
#