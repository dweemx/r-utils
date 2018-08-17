#'
#' MNP: mulit-node parallelism in R at its best
#' Author: Maxime De Waegeneer
#' 
#'@example # Define a cluster with 1 node [address] and 2 cores.
#'         cluster<-list(config=list(type="raw"
#'                                 , def=list(user = "[username]"
#'                                 , nodes = c("[address]")
#'                                 , n.cores = c(2)
#'                                 , verbose = T)))
#'         out<-mnp(l = c(1:100), f = function(x) {
#'            return (sqrt(t))
#'         }, combine = sum, cluster = cluster, monitor.progress = monitorProgress)
#'

#'@name open_PSC
#'@description    Create PSOCK cluster with the given nodes and the respective number of cores.
#'@param user     User name to connect to the cluster.
#'@param nodes    List of server names or IP addresses of the nodes.
#'@param n.cores  List of number cores for each node. 
#'@param verbose  Display some debug information.
#'@param out.file.path File path where output will be written to.
open_PSC<-function(user, nodes, n.cores, verbose = F, out.file.path = "") {
  library(doSNOW)
  library(doRNG)
  library(foreach)
  if(length(nodes) != length(n.cores)) {
    stop("Number of nodes should be the same length as the list of number of cores/node.")
  }
  machine.addresses<-lapply(X = seq_along(nodes), FUN = function(i) {
    list(host=nodes[i], user=user, ncore=n.cores[i])
  })
  
  # Set the first node as master node
  primary<-machine.addresses[[1]]$host
  
  spec<-lapply(machine.addresses, 
               function(machine) {
                 rep(list(list(host=machine$host,
                               user=machine$user)),
                     machine$ncore)
               }
  )
  spec<-unlist(spec,recursive=FALSE)
  
  if(verbose) {
    message("Building the cluster...")
  }
  cluster<-parallel::makeCluster(type='PSOCK',
                                 master=primary,
                                 spec=spec,
                                 outfile=out.file.path)
  # Register the cluster
  if(verbose) {
    message("\nRegistering the workers in doPar backend...")
  }
  registerDoSNOW(cluster)
  if(verbose) {
    message(paste("Specifications of the cluster:", foreach::getDoParWorkers(), "cores", "spread over", length(machine.addresses), "nodes.\n"))
  }
  
  psc.config<-list("master"=primary
                   , "machine.addresses"=machine.addresses
                   , "n.cores"=n.cores
                   , "cluster"=cluster)
  invisible(psc.config)
}

#'@name close_PSC
#'@description    Close the given PSOCK cluster.
#'@param cluster  Cluster object returned by open_PSC.
#'@param verbose  Display some debug information.
close_PSC<-function(cluster, verbose = F) {
  # stop cluster and remove clients
  if(verbose) {
    message("\nStop cluster and remove clients.")
  }
  stopCluster(cluster$bin$cluster)
  
  # insert serial backend, otherwise error in repetetive tasks
  if(verbose) {
    message("\nInsert serial backend.")
  }
  registerDoSEQ()
  
  # clean up a bit. (https://github.com/tobigithub/R-parallel/blob/gh-pages/R/code-setups/Install-doSNOW-parallel-DeLuxe.R)
  invisible(gc)
  remove(cluster)
}

#'@name start_mnp
#'@description    Start the multi-node parallelism cluster.
#'@param cluster  Cluster object returned by open_PSC.
start_mnp<-function(cluster) {
  if(cluster$config$type == "raw") {
    # Build the cluster
    message("Building the PS cluster...")
    cl <- open_PSC(user = cluster$config$def$user, nodes = cluster$config$def$nodes, n.cores = cluster$config$def$n.cores, verbose = cluster$config$def$verbose)
    cluster$bin<-cl
    return (cluster)
  } else if(cluster$config$type == "psc") {
    message("Please make sure that the PSC is closed after all computations are done!")
  } else {
    stop("The given cluster object is invalid!")
  }
}

#'@name stop_mnp
#'@description              Stop the multi-node parallelism cluster.
#'@param cluster            Cluster object returned by open_PSC.
#'@param monitor.progress   Whether to display the progress bar.
#'@param pb                 Progress bar object.
stop_mnp<-function(cluster, monitor.progress, pb) {
  # Close the monitor progress
  if(monitor.progress) {
    close(pb)
  }
  
  if(cluster$config$type == "raw") {
    # Close the PS cluster
    message("Closing the PS cluster...")
    close_PSC(cluster = cluster, verbose = T)
  }
}

#'@name mnp
#'@description            Run the given function f for each element in the given list l in the given multi-node cluster.
#'@param l                List of elements which the given function f will be applied on.
#'@param f                Function to apply on each of the elements in the given list l.
#'@param combine          Function to combine all the results from calling the function on each of the elements.
#'@param cluster          List defining the configuration of the cluster.
#'@param monitor.progress Whether to show a progress bar.
#'@param ...              
mnp<-function(l, f, combine, cluster, monitor.progress, ...) {
  if(monitor.progress) {
    # Monitoring the progress
    pb <- txtProgressBar(min=1, max=length(l), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  } else {
    opts <- NULL
  }
  
  library(doSNOW)
  library(doRNG)
  library(foreach)
  library(doParallel)
  
  cl<-start_mnp(cluster = cluster)

  out <- tryCatch({
    suppressPackageStartupMessages(out <- doRNG::"%dorng%"(foreach::foreach(x=l, .combine=combine, .options.snow = opts), f(x)))
    return (out)
  }, error=function(cond) {
    message(cond)
    stop_mnp(cluster = cl, monitor.progress = monitor.progress, pb = pb)
    return(NA)
  },
  warning=function(cond) {
    message(cond)
  },
  finally={
    stop_mnp(cluster = cl, monitor.progress = monitor.progress, pb = pb)
  })
  return (out)
}