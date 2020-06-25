#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#library(Regmex)
#dir <- as.character(args[1]) # dir
expfile <- as.character(args[1])# expfile
alpha <- as.numeric(args[2])# expfile
cat("efile",expfile,"alpha",alpha,"\n")
#cat(args[1])
#cat("\n")
#cat(args[2])
cat("load expression",format(Sys.time()),"\n")
#setwd(paste0("./",dir))
if(file.exists(sub("exp_","res_",expfile))){
  system(paste0("rm ",expfile))
  stop("Output already produced. Input removed.")
}
load(expfile) 
cat("load motif p1om p values",format(Sys.time()),"\n")
load("pval.mat.Rdata")
cat("load motif counts",format(Sys.time()),"\n")
load("counts.Rdata")
cat(format(Sys.time()),"\n")
var <- get(sub("[.]Rdata","",expfile))
nmots <- dim(pval.mat)[1]
samples <- 1:dim(var)[2]
ngenes <- dim(var)[1]
source("../code/wcmod.p.R")
#res.mat <- matrix(NA,ncol=dim(mlp)[1],nrow=10)
cat("ws size ",sum(sort( sapply(ls(),function(x){object.size(get(x))})))/(1024*1024*1024)," gb ",date())
cat("\n")
cat("Running ",dim(var)[2],"samples",format(Sys.time()),"\n")
#mexp <- lapply(samples,function(s)unlist(lapply(1:nmots,function(x)bbfun(mlp[x,var[,s]],mlpVariance[x],ngenes))))
mexp <- lapply(samples,function(s)unlist(lapply(1:nmots,function(x)wcmod.p(pval.mat[x,var[,s]],counts[x,var[,s]],alpha))))
cat("Saving results",sub("exp_","res_",expfile),format(Sys.time()),"\n")
save(mexp,file=sub("exp_","res_",expfile))
system(paste0("rm ",expfile))
cat(format(Sys.time()),"\n")
cat("Finished.\n")
# end

