Rphenograph_custom <- function(data, k=30, nthreads=64, distmet='euclidean'){
  library(igraph)
  library(ggplot2)
  
  if(is.data.frame(data))
    data <- as.matrix(data)
  
  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")
  
  if(k<1){
    stop("k must be a positive integer!")
  }else if (k > nrow(data)-2){
    stop("k must be smaller than the total number of points!")
  }
  
  message("Run Rphenograph starts:","\n", 
          "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
          "  -k is set to ", k)
  
  cat("  Finding nearest neighbors...")
  t1 <- system.time(neighborMatrix <- find_neighbors_new(data, k=k+1, nthreads = nthreads, metric = distmet)[,-1])
  cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
  t2 <- system.time(links <- Rphenograph:::jaccard_coeff(neighborMatrix))
  
  cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  t3 <- system.time(g <- graph.data.frame(relations, directed=FALSE))
  
  # Other community detection algorithms: 
  #    cluster_walktrap, cluster_spinglass, 
  #    cluster_leading_eigen, cluster_edge_betweenness, 
  #    cluster_fast_greedy, cluster_label_prop  
  cat("DONE ~",t3[3],"s\n", " Run louvain clustering on the graph ...")
  t4 <- system.time(community <- cluster_louvain(g))
  cat("DONE ~",t4[3],"s\n")
  
  message("Run Rphenograph DONE, totally takes ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
  cat("  Return a community class\n  -Modularity value:", modularity(community),"\n")
  cat("  -Number of clusters:", length(unique(membership(community))))
  
  return(list(g, membership(community)))
}


#' K Nearest Neighbour Search
#'
#' Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#' 
#' Use the nn2 function from the RANN package, utilizes the Approximate Near Neighbor (ANN) C++ library, 
#' which can give the exact near neighbours or (as the name suggests) approximate near neighbours 
#' to within a specified error bound. For more information on the ANN library please 
#' visit http://www.cs.umd.edu/~mount/ANN/.
#' 
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours
#' 
#' @return a n-by-k matrix of neighbor indices
#' 
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' neighbors <- find_neighbors(data, k=10)
#' 
#' @importFrom RANN nn2
#' @export
find_neighbors_new <- function(data, k, nthreads = 64, metric = 'euclidean'){
  nearest <-
    nnd_knn(k = k,
            data,
            metric = metric,
            verbose = TRUE,
            n_threads = nthreads,
            low_memory = F
    )
  return(nearest[[1]])
  
}
