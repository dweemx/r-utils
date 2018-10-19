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
#'            return (sqrt(l))
#'         }, combine = sum, cluster = cluster, monitor.progress = monitorProgress)
#'

library(doSNOW)
library(doRNG)
library(foreach)
library(doParallel)

#'@name open
#'@description  Create PSOCK cluster with the given nodes and the respective number of cores.
#'@param cluster.type   Type of the cluster to 
#'@param user           User name to connect to the cluster.
#'@param nodes          List of server names or IP addresses of the nodes.
#'@param n.cores        List of number cores for each node. 
#'@param verbose        Display some debug information.
#'@param out.file.path  File path where output will be written to.
open<-function(cluster.type = "PSOCK", user = NULL, nodes = NULL, n.cores = NULL, verbose = F, out.file.path = "") {
  
  if(cluster.type == "PSOCK") {
    if(is.null(nodes)) {
      stop("List of nodes should be defined when build a PSOCK cluster.")
    }
    
    if(is.null(nodes)) {
      stop("List of number of cores/node should be defined when build a PSOCK cluster.")
    }
  }
  
  if(length(nodes) != length(n.cores)) {
    stop("Number of nodes should be the same length as the list of number of cores/node.")
  }
  
  if(cluster.type == "MPI") {
    if(verbose) {
      message("Building the MPI cluster...")
    }
    library(Rmpi)
    cluster<-makeCluster(mpi.universe.size(), type="MPI")
    config<-list("mpi.universe"=mpi.universe.size()
               , "cluster"=cluster)
  } else if(cluster.type == "PSOCK") {
    print("Creating the PSOCK specification...")
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
      message("Building the PSOCK cluster...")
    }
    cluster<-parallel::makeCluster(type='PSOCK',
                                   master=primary,
                                   spec=spec,
                                   outfile=out.file.path)
  
    config<-list("master"=primary
               , "machine.addresses"=machine.addresses
               , "n.cores"=n.cores
               , "cluster"=cluster)
  } else {
    stop(paste0("The given cluster.type ",cluster.type," is not recognized."))
  }

  # Register the cluster
  if(verbose) {
    message("\nRegistering the workers in doPar backend...")
  }
  registerDoSNOW(cluster)
  if(verbose & cluster.type == "PSOCK") {
    message(paste("Specifications of the cluster:", foreach::getDoParWorkers(), "cores", "spread over", length(machine.addresses), "nodes.\n"))
  }
  invisible(config)
}

#'@name close
#'@description    Close the given PSOCK cluster.
#'@param cluster  Cluster object returned by open.
#'@param verbose  Display some debug information.
close<-function(cluster, verbose = F) {
  # stop cluster and remove clients
  if(verbose) {
    message("\nClosing the PS cluster...")
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
#'@param cluster  Cluster object returned by open.
start_mnp<-function(cluster, cluster.type = "PSOCK") {
  if(cluster.type == "PSOCK") {
    if(cluster$config$type == "raw") {
      # Build the cluster
      cl <- open(cluster.type = cluster.type, user = cluster$config$def$user, nodes = cluster$config$def$nodes, n.cores = cluster$config$def$n.cores, verbose = cluster$config$def$verbose)
      cluster$bin<-cl
      return (cluster)
    } else if(cluster$config$type == "psc") {
      message("Please make sure that the PSC is closed after all computations are done!")
    } else {
      stop("The given cluster object is invalid!")
    }
  } else if(cluster.type == "MPI") {
    cl <- open(cluster.type = "MPI", verbose = T)
    cluster$bin<-cl
    return (cluster)
  }
}

#'@name stop_mnp
#'@description              Stop the multi-node parallelism cluster.
#'@param cluster            Cluster object returned by open.
#'@param monitor.progress   Whether to display the progress bar.
#'@param pb                 Progress bar object.
stop_mnp<-function(cluster, monitor.progress, pb) {
  # Close the monitor progress
  if(monitor.progress) {
    close(pb)
  }
  
  if(cluster$config$type == "raw") {
    # Close the PS cluster
    close(cluster = cluster$bin$cluster, verbose = T)
  }
}

#'@name mnp
#'@description              Run the given function f for each element in the given list l in the given multi-node cluster.
#'@param l                  List of elements which the given function f will be applied on.
#'@param f                  Function to apply on each of the elements in the given list l.
#'@param combine            Function to combine all the results from calling the function on each of the elements.
#'@param cluster            List defining the configuration of the cluster.
#'@param cluster.keep.open  Whether to close the cluster when the task is finished. If true, stopping the cluster has to be done by yourself.
#'@param monitor.progress   Whether to show a progress bar.
#'@param ...              
mnp<-function(l, f, combine, cluster, cluster.type = "PSOCK", cluster.keep.open = F, monitor.progress = T, verbose = F, packages = NULL, export.vars = NULL, env = NULL, ...) {
  if(monitor.progress) {
    # Monitoring the progress
    pb <- txtProgressBar(min=1, max=length(l), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  } else {
    opts <- NULL
  }
  
  if("bin" %in% names(cluster)) {
    print("PS cluster has already been built. Skip building.")
    cl<-cluster
  } else {
    cl<-start_mnp(cluster = cluster, cluster.type = cluster.type)
  }

  if(verbose) {
    print("Cluster information:")
    print(cl)
  }
  
  if(!is.null(export.vars)) {
    parallel::clusterExport(cl = cl$bin$cluster, varlist = export.vars, envir = env)
  }
  
  out <- tryCatch({
    suppressPackageStartupMessages(out <- doRNG::"%dorng%"(foreach::foreach(x=l, .combine=combine, .options.snow = opts, .packages=packages), f(x)))
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
    if(!cluster.keep.open) {
      stop_mnp(cluster = cl, monitor.progress = monitor.progress, pb = pb)
    }
  })
  return (out)
}