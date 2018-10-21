#' Distributed parallel version of bigcor at https://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
#' 
#' @param x                 matrix. Matrix to compute the column correlation matrix from
#' @param nblocks           integer. Number of partitions to separate data
#' @param cluster           list. Cluster definition.
#'                          e.g.: list(config=list(type="raw", def=list(user = "[userid]", nodes = c("[address]"), n.cores = c([number-of-processors]), verbose = T)))
#' @param method            string. Correlation metric to compute
#' @param verbose           logical. If TRUE, information is printed in the console when running. 
#' @param monitor.progress  logical. If TRUE, show progress bar of evolution of computation.
#'
big.cor <- function(x, nblocks = 10, cluster = NULL, method = "pearson", verbose = FALSE, monitor.progress = TRUE, ...) {
  
  library(ff, quietly = TRUE)
  NCOL <- ncol(x)
  
  block.size<-NCOL/nblocks
  if(block.size < 2000) {
    warning(paste0("Blocksize is too low (",block.size,"), approximation could be bad. Try to decrease the number of blocks."))
  } else {
    print(paste0("Blocksize: ", block.size))
  }
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  print("Computing the correlation values...")
  if(is.null(cluster)) {
    print("... using sequential version.")
    for (i in 1:nrow(COMBS)) {
      COMB <- COMBS[i, ]
      G1 <- SPLIT[[COMB[1]]]
      G2 <- SPLIT[[COMB[2]]]
      if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
      flush.console()
      COR <- cor(x[, G1], x[, G2], method = method, ...)
      corMAT[G1, G2] <- COR
      corMAT[G2, G1] <- t(COR)
      COR <- NULL
    }
  } else {
    print("... using parallel version.")
    source("https://raw.githubusercontent.com/mase5/r-utils/master/mnp.R")
    cor.l<-mnp(l = 1:nrow(COMBS), f = function(i) {
      COMB <- COMBS[i, ]
      G1 <- SPLIT[[COMB[1]]]
      G2 <- SPLIT[[COMB[2]]]
      if (verbose) {
        cat("\n (Block", COMB[1], "with Block", COMB[2], ")\n")
      }
      flush.console()
      COR <- cor(x[, G1], x[, G2], method = method, ...)
      l<-list()
      l[[1]]<-list("G1"=G1, "G2"=G2, "COR"=COR)
      return (l)
    }, combine = c, cluster = cluster, cluster.type = "PSOCK", cluster.keep.open = F, monitor.progress = monitor.progress, verbose = verbose, env = environment())
    
    print("Updating the correlation matrix...")
    for(i in 1:length(cor.l)) {
      corMAT[cor.l[[i]][["G1"]], cor.l[[i]][["G2"]]]<-cor.l[[i]][["COR"]]
      corMAT[cor.l[[i]][["G2"]], cor.l[[i]][["G1"]]]<-t(cor.l[[i]][["COR"]])
    }
  }
  gc()
  return(corMAT)
}