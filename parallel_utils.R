# Build a Parallel Socket Cluster
open_PSC<-function(user, nodes, n.cores, verbose = F, out.file.path = "_PSC_log.txt") {
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