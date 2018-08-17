#' Title interpolation of the expression data along one trajectory
#' 
#' @author lindsaysmoore, mase5
#' @description This a parallel version of the interWeight function available in cellAlign R package (cfr.: https://github.com/shenorrLab/cellAlign)
#'
#' @param expDataBatch    raw single cell gene expression data (genes on rows, cells on columns)
#' @param trajCond        a vector of pseudo-time scores of the data-points whose length equals to the number of samples
#' @param winSz           window size of the interpolation
#' @param numPts          number of desired interpolated points
#' @param cluster         List defining the cluster configuration
#' @param monitorProgress Whether to display progress message
cellAlignInterWeight<-function (expDataBatch, trajCond, winSz, numPts, cluster, monitorProgress = T) {
  if (sum(is.na(trajCond)) > 0) {
    expDataBatch = expDataBatch[, -which(is.na(trajCond))]
    trajCond = trajCond[-which(is.na(trajCond))]
  }
  trajValNewPts = seq(from = min(trajCond), to = max(trajCond), 
                      length.out = numPts)
  source("https://raw.githubusercontent.com/mase5/r-utils/master/mnp.R")
  cat("\nCalculating new trajectory points...\n")
  tic()
  ValNewPts<-mnp(l = trajValNewPts, f = function(trajPt) {
    dist2Others = trajCond - trajPt
    weightedData = exp(-(dist2Others^2)/(winSz^2))
    weightedData = weightedData/sum(weightedData)
    return(as.matrix(expDataBatch) %*% weightedData)
  }, combine = cbind, cluster = cluster, cluster.keep.open = T, monitor.progress = monitorProgress, expDataBatch)
  tt<-toc()
  elapsed.time.1<-as.numeric(tt$toc-tt$tic)
  print(paste0("Elapsed time: ",elapsed.time.1))
  a = 1
  closestInt = sapply(trajCond, function(trajVal) {
    return(which.min(abs(trajValNewPts - trajVal)))
  })
  cat("\nCalculating error per gene...\n")
  tic()
  errPerGene<-mnp(l = 1:nrow(expDataBatch), f = function(rowInd) {
    return(abs(expDataBatch[rowInd, ] - ValNewPts[rowInd,
                                                  closestInt]))
  }, combine = rbind, cluster = cluster, cluster.keep.open = T, monitor.progress = T, expDataBatch, ValNewPts, closestInt)
  tt<-toc()
  elapsed.time.2<-as.numeric(tt$toc-tt$tic)
  print(paste0("Elapsed time: ",elapsed.time.2))
  
  cat("\nCalculating error for interpolated values...\n")
  tic()
  errInterpolated<-mnp(l = trajValNewPts, f = function(trajPt) {
    dist2Others = trajCond - trajPt
    weightedData = exp(-(dist2Others^2)/(winSz^2))
    weightedData = weightedData/sum(weightedData)
    return(as.matrix(errPerGene) %*% weightedData)
  }, combine = cbind, cluster = cluster, cluster.keep.open = F, monitor.progress = T, errPerGene)
  tt<-toc()
  elapsed.time.3<-as.numeric(tt$toc-tt$tic)
  print(paste0("Elapsed time: ",elapsed.time.3))
  total.elapsed.time<-elapsed.time.1+elapsed.time.2+elapsed.time.3
  print(paste0("Total elapsed time: ", total.elapsed.time))
  
  rownames(errInterpolated) = rownames(ValNewPts)
  return(list(interpolatedVals = ValNewPts, error = errInterpolated, traj = trajValNewPts))
}