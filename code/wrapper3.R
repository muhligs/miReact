#!/usr/bin/env Rscript
cat(getwd(),"\n")
load("runparameters.Rdata")
# make  mail message
sink("mailmessage.txt")
cat(paste(names(runparameters)[names(runparameters)!="motifs"],runparameters[names(runparameters)!="motifs"]),sep="\n")
cat("runparameters$motifs: 7-mers\n")
sink()

# mail if address is defined
if(!is.null(runparameters$mail)) system(paste0("mail -s 'Started ",runparameters$wd,"' ",runparameters$mail," < mailmessage.txt"))

cat("\nLoading expression data...\n")
cat(format(Sys.time()),"\n")
suppressWarnings(if(class(try(load("e.Rdata"),silent=T))=="try-error") e <- readRDS("e.Rdata"))
setwd("..") # we were in runparameters$wd (oDate folder) now move up to opd.
source("./code/addSeqs.R")
 source("./code/addMotifs.R")
 source("./code/addMotifExpression.R")
sco <- list(exp=e)
# add sequence objects
cat("\nrunning addSeqs...\n")
cat(format(Sys.time()),"\n")
sco <- addSeqs(sco,species=runparameters$species, seq.type=runparameters$seq.type)
cat("Done\n")
cat(format(Sys.time()),"\n")
# add motif objects
cat("\nrunning addMotifs...\n")
cat(format(Sys.time()),"\n")
sco <- addMotifs(sco,motifs=runparameters$motifs)
cat("Done\n")
cat(format(Sys.time()),"\n")

# add motif activity objects
cat("\nrunning addMotifExpression...\n")
cat(format(Sys.time()),"\n")
sco <- addMotifExpression(sco)
cat("Done\n")
cat(format(Sys.time()),"\n")

params <- unlist(lapply(sco$runparameters,unlist))
sco$me <- sco$me*-1 # miRNA setting, activity is inverse to motif activity
if(runparameters$out.meonly) sco <- sco$me
setwd(runparameters$wd) # move into wd again
cat("saving results...\n",runparameters$out.file,"\n")
saveRDS(sco,file=runparameters$out.file)
system("rm e.Rdata")
cat("Finished!\n")
cat(format(Sys.time()),"\n")
if(!is.null(runparameters$mail)) system(paste0("mail -s 'Finished ",runparameters$wd,"' ",runparameters$mail," < mailmessage.txt"))
