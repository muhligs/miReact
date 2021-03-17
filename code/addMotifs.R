####################################################################################################################
# add the motifs and calculate motif probabilities for the sequences.
addMotifs <- function(sco,motifs=7){
  #### checks and balances ####################
  if(!"runparameters" %in% names(sco)) stop("Run addSeqs first")
  sco$runparameters <- append(sco$runparameters,list(motifs=motifs))
  if(!runparameters$tarbaserun){
  require(Regmex)
  #################################
  # adding the motifs:
  cat("adding motif models.\n")
  if(!is.numeric(motifs)|!length(motifs)==1) stop("Motifs needs to be a number representing a k-mer from 5 to 8")
    file <- paste0("./motif.models/patterns.",motifs,"mer.Rdata")
    if(file.exists(file)) load(file) else stop(paste0("Please prepare motif models for ",moitfs,"-mers."))
    sco$motifModels <- patterns
    names(sco$motifModels) <- Regmex::all.mers(motifs)
  cat(length(sco$motifModels)," motif models added.\n\n")
  ################################# 
  # check if seqXmot calculations were already done (True for 5,6,7 mers)
    file <- paste0("./motif.models/",sco$runparameters$species,".seqXmot.",sco$runparameters$seq.type,"_mrs_",sco$runparameters$motifs,"mer.rds")
    if(!file.exists(file)) stop(paste0("Please prepare sequence specific probabilities for ",moitfs,"-mers in species ",sco$runparameters$species,"."))
      cat("Sequence specific probabilities were pre-calculated. Loading...\n")
      pval.mat <- readRDS(file)
      # load("./motif.models/utr3.seqs.Rdata") # load match file (not needed because colnames are tid per construction)
      pval.mat <- pval.mat[,sco$seqs$tid]
      sco$pval.mat <- pval.mat
  #################################
  # count number of motifs in sequences?
    cat("Making motif counts...\n")
    ################################# # in case motifs are standard k-mers, re-use if possible (5,6,7mers)
      if(is.numeric(motifs))ksize <- motifs
      file <- paste0("./motif.models/",sco$runparameters$species,".seqXmot.counts.",sco$runparameters$seq.type,"_mrs_",sco$runparameters$motifs,"mer.rds")
      if(!file.exists(file)) stop(paste0("Please prepare motif counts for ",moitfs,"-mers in species ",sco$runparameters$species,"."))
      cat("Premade. Loading...\n")
      counts <- readRDS(file)
      sco$motifCounts <- counts[,sco$seqs$tid]
      sco$motifCounts[sco$motifCounts>2] <- 2 # added on 2020-03-03 to avoid single gene bias for high motif containing genes. 
      cat("Finished.\n")
  return(sco)
  } else {
    print("tarbaserun")
        print(getwd())
    tar <- readRDS(tar, file="./data/tarbase.rds")
    if(runparameters$species=="mm") tar <- tar[tar$species=="Mus musculus"&tar$up_down %in% c(NA,"DOWN"),]
    if(runparameters$species=="hs") tar <- tar[tar$species=="Homo sapiens"&tar$up_down %in% c("DOWN"),]
        
    sco$motifModels <- unique(tar$mirna)
        
    names(sco$motifModels) <- unique(tar$mirna)
        
    pval.mat <- matrix(0.01,nrow=length(unique(tar$mirna)),ncol=nrow(sco$exp))
        
    colnames(pval.mat) <- sco$seqs$tid
        
    rownames(pval.mat) <- unique(tar$mirna)
        
    sco$pval.mat <- pval.mat
        
    gnames <- sco$seqs$gsym
        
    motifCounts <- t(as.matrix(as.data.frame(lapply(unique(tar$mirna),function(x) as.numeric(gnames %in% tar$geneName[tar$mirna==x])))))
    #print(motifCounts[1:5,1:5])
    colnames(motifCounts) <- sco$seqs$tid
    rownames(motifCounts) <- unique(tar$mirna)
    sco$motifCounts <- motifCounts
    return(sco)
    }
  }
  ####################################################################################################################
