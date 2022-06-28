####################################################################################################################
addMotifExpression <- function(sco){
  cat("Starting addMotifExpression...\n")
  # here we need to feed the variance and the number of genes.
  sco$runparameters$alpha <- 1e-10
  sco$ngenes <- dim(sco$exp)[1]
  sco$nmotifs <- length(sco$motifModels)
  sco$nsamples <- dim(sco$exp)[2]
  cat("calculating expected (median) expression values for genes\n") # Faster since data should be sparse: sco$medianExp <- sco$exp %>% Matrix::drop0() %>% sparseMatrixStats::rowMedians()
  sco$medianExp <- apply(sco$exp,1,median) 
  cat("calculating fold change rank values for each sample\n")
  var <- apply(sco$exp,2,function(x)order(x-sco$medianExp,decreasing = TRUE))
   ################################# 
    # we need a loop that will control the distribution of these calculations for fewer samples at a time.
    # we should do it so that the samples are the distributed entity, since they are potentially large.
      byparameter <- min(20,sco$nsamples) # number of samples in a job
      v.start <- seq(1,length(colnames(sco$exp)),by = byparameter) # make expression chunk intervals
      v.end <- unique(c(seq(byparameter,length(colnames(sco$exp)),by = byparameter),length(colnames(sco$exp))))
      fileapp <- paste0(v.start,"-",v.end) # file append taxt
      gwd <- getwd() # store wd for later return
      dir <- sub(" ","_",as.character(Sys.time())) # create working directory with time stamp name
      cat("Creating temporary working directory ./",dir,"\n",sep="")
        system(paste0("mkdir ",dir))
        setwd(paste0("./",dir)) # move to temp wd
        cat("saving temporary data set...\n")
        pval.mat <- sco$pval.mat
        save(pval.mat, file="pval.mat.Rdata") # we need to save these because the source files are much larger, containing all transcripts.
        counts <- sco$motifCounts
        save(counts, file="counts.Rdata")
        cat("Submitting ",length(v.start)," jobs...")
        modulus <- 100 # step interval for reporting.
        for(i in seq_along(v.start)){ # loop through the intervals for the expression chunks
        expfile <- paste0("exp_",fileapp[i]) # generate exp chunk file name
        assign(expfile,as.matrix(var[,v.start[i]:v.end[i]])) # assign data. as matrix needed if single samples apply (else vector produced).
        save(list=expfile, file=paste0(expfile,".Rdata")) # save in wd for worker to load
        # jobstarts

#################################
# prepare SBATCH script
sink("sbatch.script")
cat("#!/bin/sh","#SBATCH --nodes 1","#SBATCH -c 1","#SBATCH --time 4:00:00","#SBATCH --mem 16g","#SBATCH --job-name Rscript",sep="\n")
cat("\n[ $SLURM_PROCID -ne 0 ] && exit 0\n")
cat("cd ",getwd(),"\n")
cat("Rscript --vanilla ../code/me_mrs_script.R ",expfile,".Rdata ",sco$runparameters$alpha,"\n",sep="")
sink()
#################################

# run script
waste <- system("sbatch sbatch.script", ignore.stdout = T, ignore.stderr = T)

        #waste <- system(paste("qx --no-scratch  --time=4:00:00 --mem=16g ","\"Rscript --vanilla ../code/me_mrs_script.R ",expfile,".Rdata ",sco$runparameters$alpha,"\"\n",sep=""),intern=T)
        if(i %% modulus == 0)cat("\rjob ",i," Samples ",v.start[i]," to ",v.end[i]," ") # stepwise report submission process.
      } # Interval loop end 
        
      #################################
      # monitoring of job processes
      nfold <- -1 
      cat("\rAll jobs submitted.\n")
      cat("Jobs started...s0.\n",sep="")
      cat("Jobs completed...c0.\n",sep="")
      # monitor results until all are complete
      lfc <- length(list.files(pattern = "res*"))  # jobs completed
      lfs <- length(list.files(pattern = "[.]out$")) # jobs started
      oji <- c(0,0) # job inspector
      linel <- 0 # a variable that holds the position of a write out cursor.
      nproc <- length(v.start) # total number of jobs
      while(lfc<nproc){ 
        Sys.sleep(60) # take a nap.
        if(!identical(c(lfs,lfc),oji)){ # if something new happened
          oji <- c(lfs,lfc) # update job inspector
          cat(".s",lfs,"c",lfc,sep="") # write progress out
          linel <- linel+nchar(as.character(lfs))+nchar(as.character(lfc))+3 # update cursor position
          if(linel>80){linel <- 0;cat("\n")} # add a new line if needed.
        }
        lfc <- length(list.files(pattern = "res*"))  # update jobs completed
        lfs <- length(list.files(pattern = "[.]out$")) # update jobs started
      }
      #################################      

      #################################      
      # collect results and cleanup.
      lf <- list.files(pattern = "res_")
      cat("\nLoading ",length(lf)," result files... \n")
      lf <- lf[order(as.numeric(sub("-.*","",sub("res_","",lf))),decreasing = FALSE)] # define the order for loading results
      for(i in seq_along(lf)){ # load results into mexp.res
        load(lf[i])
        if(i==1) mexp.res <- mexp else mexp.res <- append(mexp.res,mexp)
      }
      mexp <- mexp.res
      
      #################################      
      # rearrange data
      cat("re-arranging data...\n")
      sco$me <- do.call(rbind,mexp)
      rownames(sco$me) <- colnames(sco$exp) # add gene IDs
      sco$me <- t(sco$me)
      if(is.null(rownames(sco$me))) rownames(sco$me) <- names(sco$motifModels) # add motifs if needed 
      cat("Truncating overflow values.\n") # truncate INF values
      cat(sum(!is.finite(sco$me))," values truncated. ",sum(sco$me==Inf)," Inf values and ",sum(sco$me==(-Inf))," -Inf values.\n")
      sco$me[sco$me==Inf] <- max(sco$me[is.finite(sco$me)])+0.01
      sco$me[sco$me==(-Inf)] <- min(sco$me[is.finite(sco$me)])-0.01
      #################################      
      
      #################################      
      # cleaning up
      cat("Cleaning up\n")
      system("rm res_*") 
      system("rm Rscript*") 
      system("rm pval.mat.Rdata");system("rm counts.Rdata")
      system("rm sbatch.script")
      setwd(gwd) # jump to old wd      
      system(paste0("rmdir ",dir)) # remove temp wd
      #################################      

      return(sco) # return results
}
####################################################################################################################
