mireact <- function(exp, motifs=7, 
	species="hs", seq.type="utr3", out.file = "res.data", out.meonly = T,	verbose = F , install.dir = "~/temp/miRNAsc", mail = NULL, memory="128g", wt="48:00:00", tarbaserun=FALSE){
	#################################

  ##### runparameters ############################
  runparameters <- as.list(environment())
  runparameters$exp <- ifelse(class(exp)=="character",exp,class(exp))
  if(verbose){
  p.runparameters <- runparameters
  p.runparameters$motifs <- "7-mers"
  print(p.runparameters)
  }
  #################################

		
  #### input checks #############################
  # input checks
  # 
  owd <- getwd() 
  setwd(install.dir) # jump to install directory
  if(!file.exists("./code/lambdaProbdist.R"))stop("Please set the install directory: install.dir='path/to/install/directory'") # check we are in the directory where data and code lives.
  # check for sbatch installation
  var <- try(system("command -v sbatch",intern=T),silent = T)
  if(grepl("Error",var))stop("sbatch is needed. Please install")
  # check for Regmex and expm
  require(Regmex)
  require(expm)
  #################################
	

	#### run directory #############################
	# create run directory and save run parameter object in it.
	dir <- paste0("o",sub(" ","",as.character(Sys.time())))
  wd <- paste0(getwd(),"/",dir)
  cat("creating directory: ",wd,"\n")
  system(paste0("mkdir ",wd))
  setwd(wd)
  runparameters$wd <- wd
	#################################

	#### start wrapper #############################
    # place expression matrix in directory
    save(runparameters, file="runparameters.Rdata")
  
		if(is.matrix(exp))save(exp,file="e.Rdata") else system(paste0("ln -s ",exp," ",wd,"/e.Rdata"))

#################################
# prepare SBATCH script
sink("sbatch.script")
cat("#!/bin/sh","#SBATCH --nodes 1","#SBATCH -c 1",paste0("#SBATCH --time ",wt),paste0("#SBATCH --mem ",memory),"#SBATCH --job-name Rscript",sep="\n")
cat("\n[ $SLURM_PROCID -ne 0 ] && exit 0\n")
cat("cd ",wd,"\n")
cat("Rscript --vanilla ../code/wrapper3.R\n")
sink()
#################################

# run script
waste <- system("sbatch sbatch.script", ignore.stdout = T, ignore.stderr = T)
setwd(owd)
return(paste0("Follow progress in ",wd,"/Rscript-[jobid].out"))
}

	

