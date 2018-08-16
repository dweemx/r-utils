# Build a Parallel Socket Cluster
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

close_PSC<-function(cluster, verbose = F) {
    # stop cluster and remove clients
    if(verbose) {
        message("\nStop cluster and remove clients.")
    }
    stopCluster(cluster)
    
    # insert serial backend, otherwise error in repetetive tasks
    if(verbose) {
        message("\nInsert serial backend.")
    }
    registerDoSEQ()

    # clean up a bit. (https://github.com/tobigithub/R-parallel/blob/gh-pages/R/code-setups/Install-doSNOW-parallel-DeLuxe.R)
    invisible(gc)
    remove(cluster); 
}

goPaR<-function(title, l, ex, combine) {
    if(monitorProgress) {
        # Monitoring the progress
        pb <- txtProgressBar(min=1, max=length(trajValNewPts), style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
    } else {
        opts <- NULL
    }
    
    library(doSNOW)
    library(doRNG)
    library(foreach)
    library(doParallel)
    
    if(cluster$config$type == "raw") {
        # Load utils to build the cluster
        source("https://raw.githubusercontent.com/mase5/r-utils/master/parallel_utils.R")
        # Build the cluster
        message("Building the PS cluster...")
        cl <- open_PSC(user = cluster$config$def$user, nodes = cluster$config$def$nodes, n.cores = cluster$config$def$n.cores, verbose = cluster$config$def$verbose)
        cl.c <- cl$config$def$cluster
    } else if(cluster$config$type == "psc") {
        message("Please make sure that the PSC is closed after all computations are done!")
        cl.c <- cluster$config$def$cluster
    } else {
        stop("The given cluster object is invalid!")
    }

    suppressPackageStartupMessages(out <- doRNG::"%dorng%"(foreach::foreach(x=l, .combine=combine, .options.snow = opts), ex))
    
    # Close the monitor progress
    if(monitorProgress) {
        close(pb)
    }
    
    if(cluster$config$type == "raw") {
        # Close the PS cluster
        message("Closing the PS cluster...")
        close_PSC(cluster = cl.c, verbose = T)
    }

    return (out)
}