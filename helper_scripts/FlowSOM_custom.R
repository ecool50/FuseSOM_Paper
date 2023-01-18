FlowSOM_custom <- function(input, pattern = ".fcs", 
                    compensate = FALSE, spillover = NULL, 
                    transform = FALSE, toTransform = NULL, 
                    transformFunction = flowCore::logicleTransform(), 
                    transformList = NULL, scale = FALSE, 
                    scaled.center = TRUE, scaled.scale = TRUE, silent = TRUE, 
                    colsToUse = NULL, nClus = 10, maxMeta = NULL, importance = NULL, 
                    seed = NULL, clusterAlg = 'hc', distance = 'euclidean', ...){
  # Method to run general FlowSOM workflow. 
  # Will scale the data and uses consensus meta-clustering by default.
  #
  # Args:
  #    input: dirName, fileName, array of fileNames, flowFrame or 
  #           array of flowFrames
  #    colsToUse: column names or indices to use for building the SOM
  #    maxMeta: maximum number of clusters for meta-clustering
  #
  # Returns:
  #    list with the FlowSOM object and an array with final clusterlabels
  
  # library(FlowSOM)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  
  t <- system.time(fsom <- ReadInput(input, pattern = pattern, 
                                     compensate = compensate, 
                                     spillover = spillover, 
                                     transform = transform, 
                                     toTransform = toTransform, 
                                     transformFunction = transformFunction, 
                                     transformList = transformList,
                                     scale = scale,
                                     scaled.center = scaled.center, 
                                     scaled.scale = scaled.scale, 
                                     silent = silent))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- FlowSOM::BuildSOM(fsom, colsToUse, silent = silent, 
                                    importance = importance))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- FlowSOM::BuildMST(fsom, silent = silent))
  if(!silent) message(t[3], "\n")
  if(is.null(maxMeta)){
    t <- system.time(cl <- as.factor(
      metaClustering_consensus_custom(fsom$map$codes, nClus, seed = seed,
                                      clusterAlg = clusterAlg, distance = distance)))
  } else {
    t <- system.time(cl <- as.factor(MetaClustering_custom(fsom$map$codes,
                                                    "metaClustering_consensus_custom", 
                                                    maxMeta,
                                                    seed = seed,
                                                    clusterAlg = clusterAlg, distance = distance)))
  }
  fsom$map$nMetaclusters <- length(levels(cl))
  fsom$metaclustering <- cl
  # fsom <- FlowSOM:::UpdateDerivedValues(fsom)
  fsom$info$parameters <- match.call()
  fsom$info$date <- as.character(Sys.time())
  fsom$info$version <- as.character(utils::packageVersion("FlowSOM"))
  if(!silent) message(t[3], "\n")
  return(fsom)
  
}



#' MetaClustering
#'
#' Cluster data with automatic number of cluster determination for 
#' several algorithms
#'
#' @param data   Matrix containing the data to cluster
#' @param method Clustering method to use
#' @param max    Maximum number of clusters to try out
#' @param seed   Seed to pass on to given clustering method
#' @param ...    Extra parameters to pass along
#' 
#' @return Numeric array indicating cluster for each datapoint
#' @seealso   \code{\link{metaClustering_consensus}}
#'
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate = TRUE,transform = TRUE,
#'                             scale = TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse = c(9, 12, 14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply metaclustering
#'    metacl <- MetaClustering(flowSOM.res$map$codes,
#'                             "metaClustering_consensus",
#'                             max = 10)
#'    
#'    # Get metaclustering per cell
#'    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[, 1]]    
#'
#' @export
MetaClustering_custom <- function(data, method, max = 20, seed = NULL, 
                                  clusterAlg = 'hc', distance = 'euclidean',...){
  res <- DetermineNumberOfClusters(data, max, method, seed = seed, ...)
  method <- get(method)
  method(data, k = res, seed = seed, clusterAlg = clusterAlg, distance = distance)
}

DetermineNumberOfClusters_custom <- function(data, max, method, plot = FALSE,
                                      smooth = 0.2, seed = NULL,clusterAlg = 'hc', distance = 'euclidean', ...){
  # Try out a clustering algorithm for several numbers of clusters and 
  # select optimal
  #
  # Args:
  #     data:     Matrix containing the data to cluster
  #     max:        Maximum number of clusters to try
  #     method: Clustering method to use
  #     plot:     Whether to plot the results for different k
  #     smooth: Smoothing option to find elbow: 
  #             0: no smoothing, 1: maximal smoothing
  #     seed:   Seed to pass on to given method
  #
  # Returns:
  #     Optimal number of clusters
  if(method ==    "metaClustering_consensus_custom"){
    results <- consensus_custom(data,max, seed, 
                                clusterAlg = clusterAlg, distance = distance, ...)
    res <- rep(0,max)
    res[1] <- SSE_custom(data,rep(1,nrow(data)))
    for(i in 2:max){
      c <- results[[i]]$consensusClass
      res[i] <- SSE_custom(data, c)
    }
  } else {
    method <- get(method)
    res <- rep(0,max)
    for(i in 1:max){
      c <- method(data, k = i,...)
      res[i] <- SSE(data, c)
    }
  }
  
  for(i in 2:(max - 1)){
    res[i] <- (1 - smooth) * res[i] + 
      (smooth / 2) * res[i - 1] + 
      (smooth / 2) * res[i + 1]
  }
  
  if(plot) plot(1:max, res, type = "b", xlab = "Number of Clusters", 
                ylab = "Within groups sum of squares")
  findElbow(res)
}

findElbow <- function(data){
  n <- length(data)    
  data <- as.data.frame(cbind(1:n,data))
  colnames(data) <- c("X","Y")
  
  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  optimal
}

#' MetaClustering
#' 
#' Cluster data using hierarchical consensus clustering with k clusters
#'
#' @param data Matrix containing the data to cluster
#' @param k    Number of clusters
#' @param seed Seed to pass to consensusClusterPlus
#' 
#' @return  Numeric array indicating cluster for each datapoint
#' @seealso \code{\link{MetaClustering}}
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate = TRUE,transform = TRUE,
#'                             scale = TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse = c(9, 12, 14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply consensus metaclustering
#'    metacl <- metaClustering_consensus(flowSOM.res$map$codes, k = 10)    
#'
#' @export
metaClustering_consensus_custom <- function(data, k = 7, seed = NULL, 
                                     clusterAlg = 'hc', distance = 'pearson'){
  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1,
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = clusterAlg, # "hc","km","kmdist","pam"
    distance = distance ,
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))
  # 
  results[[k]]$consensusClass
  # if(distance == 'euclidean'){
  #   d <- parallelDist(data)
  # }else{
  #   d <- DiffCorr::cor.dist(data)
  # }
  # 
  # np <- NNGraphParam(k=20, cluster.fun="louvain")
  # res <- clusterRows(d, np)
  
  # print(paste(dim(data)))
  # res <- ADPclustering(as.matrix(d), ClusterNo = k)$Cls
  # return(res)
  
  
  
}

consensus_custon <- function(data, max, seed = NULL,
                             clusterAlg = 'hc', distance = 'euclidean'){
  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = max, reps = 100, pItem = 0.9, pFeature = 1,
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = 'hc', # "hc","km","kmdist","pam"
    distance = distance,
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))
}

# metaClustering_hclust <- function(data, k = 7, seed = NULL){
#   set.seed(seed)
#   d <- stats::dist(data, method = "minkowski")
#   fit <- stats::hclust(d, method = "ward.D2")
#   stats::cutree(fit, k = k)
# }
# 
# metaClustering_kmeans <- function(data, k = 7, seed = NULL){
#   set.seed(seed)
#   stats::kmeans(data, centers = k)$cluster
# }
# 
# metaClustering_som <- function(data, k = 7, seed = NULL){
#   set.seed(seed)
#   s <- SOM(data, xdim = k, ydim = 1, rlen = 100)
#   s$unit.classif
# }

SSE_custom <- function(data,clustering){
  if(!is(clustering, "numeric"))
    clustering <- as.numeric(as.factor(clustering))
  c_wss <- 0
  for(j in seq_along(clustering)){
    if(sum(clustering == j) > 1){
      c_wss <- c_wss + (nrow(data[clustering == j, , drop = FALSE]) - 1) *
        sum(apply(data[clustering == j, , drop = FALSE], 2, stats::var))
    }
  }
  c_wss
}

#' F measure
#' 
#' Compute the F measure between two clustering results
#'
#' @param realClusters Array containing real cluster labels for each sample
#' @param predictedClusters Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE (default), print some information about 
#'                  precision and recall
#' 
#' @return  F measure score
#' @examples
#' # Generate some random data as an example
#' realClusters <- sample(1:5,100,replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' FMeasure(realClusters,predictedClusters)
#' @export
FMeasure <- function(realClusters, predictedClusters,silent = FALSE){
  if (sum(predictedClusters) == 0)
    return(0);
  a <- table(realClusters, predictedClusters);
  p <- t(apply(a, 1, function(x) x / colSums(a)))
  r <- apply(a, 2, function(x) x / rowSums(a))
  if(!silent) message("Precision: ",
                      sum(apply(p, 1, max) * (rowSums(a) / sum(a))),
                      "\nRecall: ", 
                      sum(apply(r, 1, max) * (rowSums(a) / sum(a))), 
                      "\n")
  f <- 2 * r * p / (r + p)
  f[is.na(f)] <- 0
  sum(apply(f, 1, max) * (rowSums(a) / sum(a)))
}
