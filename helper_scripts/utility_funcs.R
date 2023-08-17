# a function to run Kmeans
kmeansRun <- function(dat, markers, distmet = 'pearson'){
  celltypes <- dat$CellType
  kres <- amap::Kmeans(dat[, markers], 
                       length(unique(celltypes)),method = distmet, 
                       iter.max = 100)
  km.clusts <- kres$cluster
  ARI = aricode::ARI(km.clusts, dat$CellType)
  NMI = aricode::NMI(km.clusts, dat$CellType)
  FM = dendextend::FM_index(km.clusts, dat$CellType)[[1]]
  Fmeasure <-  FlowSOM::FMeasure(as.numeric(km.clusts), as.numeric(as.factor(dat$CellType)))
  
  return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM, FMeasure = Fmeasure,
              Metric = distmet))
  
}

# A function to Hclust
HclustRun <- function(dat, markers, distmet = 'pearson', link = 'average'){
  celltypes <- dat$CellType
  Hres <- amap::hcluster(dat[, markers],
                            method = distmet, nbproc = 8, link = link)
  hc.clusts <- cutree(Hres, length(unique(celltypes)))
  ARI = aricode::ARI(hc.clusts, celltypes)
  NMI = aricode::NMI(hc.clusts, celltypes)
  FM = dendextend::FM_index(hc.clusts, celltypes)
  Fmeasure = FlowSOM::FMeasure(as.numeric(hc.clusts), as.numeric(as.factor(celltypes)))
  
  return(list(ARI = ARI, NMI = NMI, FM = FM, FMeasure = Fmeasure,
              Metric = distmet))
  
}

# A function to run Phenograph
PhenoRun <- function(dat, distmet = 'euclidean', markers, nthreads = 8, k = 20){
  pheno_clusters <- ReductionWrappers::phenograph(dat[, markers], n_jobs = nthreads,
                                primary_metric = distmet, k = k)
  # pheno <- Rphenograph_custom(dat[, markers],
  #                           k = 20, distmet = distmet)
  # pheno_clusters <- as.character(pheno)
  # ARI = aricode::ARI(pheno_clusters, dat$CellType)
  # NMI = aricode::NMI(pheno_clusters, dat$CellType)
  # FM = dendextend::FM_index(pheno_clusters, dat$CellType)[[1]]
  # Fmeasure = FlowSOM::FMeasure(as.numeric(pheno_clusters), as.numeric(as.factor(dat$CellType)))
  # 
  # return(list(ARI = ARI, NMI = NMI, FM = FM, Fmeasure = Fmeasure))
  return(pheno_clusters)
}

# A function to run flowsom
FloSOMRun <- function(dat, distance = 'euclidean', clusterAlg = 'hc', markers, 
                      normalization = 'None', labels = F){
  gc()
  dat_new <- dat[, markers]
  if(normalization == 'None'){
    dat_new <- dat_new
  } else if(normalization == 'percentile'){
    dat_new <- percentile(dat_new)
  } else if(normalization == 'min_max'){
    dat_new <- min_max(dat_new)
  } else if(normalization == 'zscore'){
    dat_new <- scale(dat_new, center = T, scale = T)
  } else {
    dat_new <- arsinh(dat_new)
  }
  
  
  ff <- DFtoFF(as.data.frame(dat_new))
  flo_res <- FlowSOM_custom(ff,scale=F,nClus=length(unique(dat$CellType)), 
                     distance = distance, clusterAlg = clusterAlg)
  flo_res.clusters <- GetMetaclusters(flo_res)
  ARI = aricode::ARI(flo_res.clusters, dat$CellType)
  NMI = aricode::NMI(flo_res.clusters, dat$CellType)
  FM = dendextend::FM_index(flo_res.clusters, dat$CellType)[[1]]
  Fmeasure <-  FlowSOM::FMeasure(as.numeric(flo_res.clusters), as.numeric(as.factor(dat$CellType)))
  
  if(labels) {
    return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM, FMeasure = Fmeasure, Labels = flo_res.clusters))
  }

  return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM, FMeasure = Fmeasure))
  # return(flo_res.clusters)
}

subsetMatrix <- function(mat, markers, rep = 5, seed = 1994, nrows = 5000) {
  # subset the matrix
  set.seed(seed)
  sub.mat <- list()
  
  # remove rows with 0 stdev
  mat <- mat[!apply(mat[, markers], 1, sd, na.rm=TRUE) == 0, ]
  set.seed(seed)

  for(i in 1:rep) {
    rows <- sample(nrow(mat), nrows)
    sub.mat[[i]] <- mat[rows, ]
    
  }
  return (sub.mat)
}


percentile <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  q_1 <- percentiles[[1]]
  q_99 <- percentiles[[2]]
  
  x <- (x - q_1)/(q_99)
  return(x)
}

min_max <- function(x){
  x <- as.matrix(x)
  percentiles <- quantile(x, probs = c(0.01, 0.99))
  min_val <- percentiles[[1]]
  max_val <- percentiles[[2]]
  x <- (x - min_val)/(max_val - min_val)
  return(x)
}

arsinh <- function(x, cofactor=5){
  x <- asinh(x/cofactor)
  return(x)
}


# A function to generate the plots
GeneratePlots <- function(kmeans.pear, kmeans.eucl, hclust.pear, hclust.eucl,
                          pheno.pear, pheno.eucl, flo.pear, flo.eucl, datName = NULL, numIter = 5){
  # ARI
  res.ari <- data.frame(kmeans_corr = sapply(kmeans.pear,"[[",1), kmeans_eucl = sapply(kmeans.eucl,"[[",1), 
                         hclust_corr = sapply(hclust.pear,"[[",1), hclust_eucl = sapply(hclust.eucl,"[[",1), 
                         pheno_corr = sapply(pheno.pear,"[[",1), pheno.eucl = sapply(flo.eucl,"[[",1),
                         flow_corr = sapply(flo.pear,"[[",1), flow_eucl = sapply(flo.eucl,"[[",1))
  res.ari <- as.data.table(res.ari)
  
  res.ari.melted <- melt(res.ari)
  colnames(res.ari.melted) <- c('Method', 'ARI')
  
  # NMI
  res.nmi <- data.frame(kmeans_corr = sapply(kmeans.pear,"[[",2), kmeans_eucl = sapply(kmeans.eucl,"[[",2), 
                        hclust_corr = sapply(hclust.pear,"[[",2), hclust_eucl = sapply(hclust.eucl,"[[",2), 
                        pheno_corr = sapply(pheno.pear,"[[",2), pheno.eucl = sapply(flo.eucl,"[[",2),
                        flow_corr = sapply(flo.pear,"[[",2), flow_eucl = sapply(flo.eucl,"[[",2))
  res.nmi <- as.data.table(res.nmi)
  
  res.nmi.melted <- melt(res.nmi)
  colnames(res.nmi.melted) <- c('Method', 'NMI')
  
  # FM
  res.fm <- data.frame(kmeans_corr = sapply(kmeans.pear,"[[",3), kmeans_eucl = sapply(kmeans.eucl,"[[",3), 
                        hclust_corr = sapply(hclust.pear,"[[",3), hclust_eucl = sapply(hclust.eucl,"[[",3), 
                        pheno_corr = sapply(pheno.pear,"[[",3), pheno.eucl = sapply(flo.eucl,"[[",3),
                        flow_corr = sapply(flo.pear,"[[",3), flow_eucl = sapply(flo.eucl,"[[",3))
  res.fm <- as.data.table(res.fm)
  
  res.fm.melted <- melt(res.fm)
  colnames(res.fm.melted) <- c('Method', 'FM')
  
  # generate plots
  p.ari <- ggplot(res.ari.melted, aes(x=Method, y=ARI)) + 
    geom_boxplot(aes(color = Method)) + 
    ggtitle(paste0(datName, " Subsetted to 10K Cells"), subtitle = paste0('Number of Iteration = ', numIter)) +
    xlab("Method") + ylab("Adjusted Rand Index") + 
    theme(
      plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      plot.subtitle = element_text(color = "red", size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(color="Black", size=14, face="bold"),
      axis.title.y = element_text(color="Black", size=14, face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
  
  p.nmi <- ggplot(res.nmi.melted, aes(x=Method, y=NMI)) + 
    geom_boxplot(aes(color = Method)) + 
    ggtitle(paste0(datName, " Subsetted to 10K Cells"), subtitle = paste0('Number of Iteration = ', numIter)) +
    xlab("Method") + ylab("Normalized Mutual Information") + 
    theme(
      plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      plot.subtitle = element_text(color = "red", size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(color="Black", size=14, face="bold"),
      axis.title.y = element_text(color="Black", size=14, face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
  
  p.fm <- ggplot(res.fm.melted, aes(x=Method, y=FM)) + 
    geom_boxplot(aes(color = Method)) + 
    ggtitle(paste0(datName, " Subsetted to 10K Cells"), subtitle = paste0('Number of Iteration = ', numIter)) +
    xlab("Method") + ylab("Fowlkes–Mallows Index") + 
    theme(
      plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      plot.subtitle = element_text(color = "red", size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(color="Black", size=14, face="bold"),
      axis.title.y = element_text(color="Black", size=14, face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
  
  p.final <- cowplot::plot_grid(p.ari, p.nmi, p.fm)
  
  return(p.final)
  
  
}


generateHeatmap <- function(dat, markers, features_clust){
  features <- dat
  features_heatmap <- aggregate(.~as.character(features_clust),
                                features[,markers],
                                mean)
  rownames(features_heatmap) <- features_heatmap[,1]
  features_heatmap <- features_heatmap[,-1]

  features_heatmap <- sweep(features_heatmap,2, colMeans(features_heatmap), "-")
  features_heatmap <- sweep(features_heatmap,2, apply(features_heatmap,2,sd), "/")
  features_heatmap[features_heatmap>2] <- 2
  features_heatmap[features_heatmap< -2] <- -2
  
  annotation_row = data.frame(Clusters = rownames(features_heatmap))
  
  rn <- rownames(features_heatmap)
  features_heatmap <- as.matrix(features_heatmap)
  rownames(features_heatmap) <- rn
  rownames(annotation_row) <- rownames(features_heatmap)
  
  gaps_row <- which(!duplicated(substr(rownames(features_heatmap),1,2)))[-1]-1
  
  p <- ggplotify::as.ggplot(pheatmap(features_heatmap, gaps_row = gaps_row, 
                                     annotation_row = annotation_row, annotation_legend = FALSE, 
                                     cluster_rows = FALSE, cluster_cols = F))
  return(p)
}


RunPipeline <- function(dat, markers){
  subsets <- subsetMatrix(dat)
  # Kmeans with eulidean
  print(paste('=====================Now running Kmeans on Euclidean and Correlation====================='))
  kmeans.pear <- lapply(subsets, kmeansRun, markers = markers) %>% bind_rows()
  kmeans.pear$Method <- rep('kmeans_pear', nrow(kmeans.pear))
  
  # kmeans with correlation
  kmeans.eucl <- lapply(subsets, kmeansRun, markers = markers,
                                 distmet = 'euclidean')  %>% bind_rows()
  kmeans.eucl$Method <- rep('kmeans_eucl', nrow(kmeans.eucl))
  
  print(paste('=====================Now running Hclust on Euclidean and Correlation====================='))
  # hclust with pearson
  hclust.pear <- lapply(subsets, HclustRun, markers = markers)  %>% bind_rows()
  hclust.pear$Method <- rep('hclust_pear', nrow(hclust.pear))
  # hclust with euclidean
  hclust.eucl <- lapply(subsets, HclustRun, markers = markers,
                                 distmet = 'euclidean')  %>% bind_rows()
  hclust.eucl$Method <- rep('hclust_eucl', nrow(hclust.eucl))
  
  print(paste('=====================Now running FlowSOM on Euclidean and Correlation====================='))
  # flosom withe pearson
  flo.pear <- lapply(subsets, FloSOMRun, distance = 'pearson',
                              markers = markers)  %>% bind_rows()
  flo.pear$Method <- rep('FloSOM_pear', nrow(flo.pear))
  # flowsom with euclidean
  flo.eucl <- lapply(subsets, FloSOMRun, markers = markers)  %>% bind_rows()
  flo.eucl$Method <- rep('FloSOM_eucl', nrow(flo.eucl))
  
  print(paste('=====================Now running Phenograph on Euclidean and Correlation====================='))
  # phenograph with pearson
  pheno.pear <- lapply(subsets, PhenoRun, distmet = 'correlation',
                                markers = markers)  %>% bind_rows()
  pheno.pear$Method <- rep('Phengraph_pear', nrow(pheno.pear))
  # phenograph with euclidean
  pheno.eucl <- lapply(subsets, PhenoRun, markers = markers)  %>% bind_rows()
  pheno.eucl$Method <- rep('Phengraph_eucl', nrow(pheno.eucl))
  
  
  print(paste('=====================Now running FloSOM Custom on Euclidean and Correlation====================='))
  # phenograph with pearson
  eli.pear <- lapply(subsets, RunSOM, distmet = 'pearson',
                       markers = markers)  %>% bind_rows()
  eli.pear$Method <- rep('Elijah_pear', nrow(eli.pear))
  # phenograph with euclidean
  eli.eucl <- lapply(subsets, RunSOM, markers = markers,
                     distmet = 'euclidean')  %>% bind_rows()
  eli.eucl$Method <- rep('Elijah_eucl', nrow(eli.eucl))
  
  print(paste('================================Now Parsing Pipeline Results================================'))
  # combine the results
  res <- rbind(kmeans.eucl, kmeans.pear,hclust.eucl,
                        hclust.pear, flo.eucl, flo.pear, pheno.eucl, pheno.pear,
                        eli.eucl, eli.pear)
  return(res)
  
}



RunPipelineFull <- function(dat, markers){
  subsets <- subsetMatrixFull(dat)
  print(paste('=====================Now running FlowSOM ========================'))
  # flosom withe pearson
  flo <- lapply(subsets, FloSOMRun,
                     markers = markers, distance = 'pearson')  %>% bind_rows()
  flo$Method <- rep('FlowSOM', nrow(flo))
  gc()
  
  
  print(paste('=====================Now running FuseSOM================'))
  # phenograph with pearson
  fuse.res <- lapply(subsets, RunSOM,
                     markers = markers)  %>% bind_rows()
  fuse.res$Method <- rep('FuseSOM', nrow(fuse.res))
  
  
  print(paste('===============Now Parsing Pipeline Results======================'))
  # combine the results
  res <- rbind(flo, fuse.res)
  return(res)
}


#' Train the SOM, specifiy grid size and other optional parameters based on the
#' SOM Toolbox
#' @references
#'   http://www.cis.hut.fi/somtoolbox/package/docs2/som_topol_struct.html
#' @noRd
som_train <- function(x, xdim = 10, ydim = 10, rlen = 10, alpha = c(0.05,0.01), 
                      topo = "hexagonal") {
  # Create SOM grid and map data into the grid
  neurons <- 5 * sqrt(nrow(x))
  eigenvalues <- eigen(stats::cor(x))$values
  eigenratio <- eigenvalues[1] / eigenvalues[2]
  xdim <- xdim %||% sqrt(neurons / eigenratio)
  ydim <- ydim %||% neurons / xdim
  grid = kohonen::somgrid(xdim, ydim, topo = 'hexagonal')
  # x <- t(scale(t(x), scale = F))
  # init <- aweSOM::somInit(as.matrix(x), xdim, ydim, method = 'pca')
  kohonen::som(x, grid = grid, 
               rlen = 100, 
               alpha = c(0.01, 0.04), 
               dist.fcts = 'sumofsquares'
               )
}


GeneratePrototypes <- function(data, markers, normalization = 'None'){
  
  dat_new <- data[, markers]
  if(normalization == 'None'){
    dat_new <- dat_new
  } else if(normalization == 'percentile'){
    dat_new <- percentile(dat_new)
  } else if(normalization == 'min_max'){
    dat_new <- min_max(dat_new)
  } else if(normalization == 'zscore'){
    dat_new <- scale(dat_new, center = T, scale = T)
  } else {
    dat_new <- arsinh(dat_new)
  }
  print(paste('=====================Now generating SOM grid======================='))
  npcs <- SC3:::estkTW(data[, markers])
  # k = length(unique(data$CellType))
  k = 10
  print(paste('Optimal Grid Size is: ', npcs))
  # 
  if((npcs*npcs) < k){
    npcs <- 10
  } else{
    npcs <- npcs + 2
  }
  
  sg <- yasomi::somgrid(xdim=npcs,ydim=npcs,topo="hex")
  
  init.res <- yasomi::sominit.pca(as.matrix(dat_new), somgrid = sg)
  print(paste('=====================Now tunning and running the SOM model====================='))
  
  som_model <- yasomi::batchsom(as.matrix(dat_new), sg,
                                prototypes = init.res$prototypes,
                                verbose = F
  )
  return(som_model)
}

ClusterPrototypes <- function(som_model, distmet = 'pearson', numClusters){
  prototypes <- som_model$prototypes
  print(paste('=====================Now clustering prototypes====================='))
  
  # all the single distances
  if(distmet == 'pearson'){
    S <- cor2dist(stats::cor(t(prototypes), method = 'pearson'))
  } else if(distmet == 'cosine'){
    S <- cor2dist(coop::tcosine(prototypes))
  } else if(distmet == 'spearman'){
    S <- cor2dist(stats::cor(t(prototypes), method = 'spearman'))
    # the double distance metrics
  } else if(distmet == 'euclidean'){
    S <- as.matrix(dist(prototypes))
  } else if(distmet == 'pear_spear'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    spear <- stats::cor(t(prototypes), method = 'spearman')
    S <- as.matrix(fuse(cor2dist(pear),cor2dist(spear)))
  } else if(distmet == 'cos_spear'){
    cosi <- coop::tcosine(prototypes)
    spear <- stats::cor(t(prototypes), method = 'spearman')
    S <- as.matrix(fuse(cor2dist(cosi),cor2dist(spear)))
  } else if(distmet == 'pear_cos'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    cosi <- coop::tcosine(prototypes)
    S <- as.matrix(fuse(cor2dist(cosi),cor2dist(pear)))
  } else if(distmet == 'pear_eucl'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(pear), eucl))
  } else if(distmet == 'spear_eucl'){
    spear <- stats::cor(t(prototypes), method = 'spearman')
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(spear), eucl))
  } else if(distmet == 'cos_eucl'){
    cosi <- coop::tcosine(prototypes)
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(cosi), eucl))
    # the triple distance metrics
  } else if(distmet == 'pear_spear_cos'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    spear <- stats::cor(t(prototypes), method = 'spearman')
    cosi <- coop::tcosine(prototypes)
    S <- as.matrix(fuse(cor2dist(pear), cor2dist(spear), cor2dist(cosi)))
    
  } else if(distmet == 'pear_spear_eucl'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    spear <- stats::cor(t(prototypes), method = 'spearman')
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(pear), cor2dist(spear), eucl))
    
  } else if(distmet == 'spear_cos_eucl'){
    spear <- stats::cor(t(prototypes), method = 'spearman')
    cosi <- coop::tcosine(prototypes)
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(spear), cor2dist(cosi), eucl))
    
  } else if(distmet == 'pear_cos_eucl'){
    spear <- stats::cor(t(prototypes), method = 'spearman')
    cosi <- coop::tcosine(prototypes)
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(spear), cor2dist(cosi), eucl))
    # all four metrics
  } else if(distmet == 'all'){
    pear <- stats::cor(t(prototypes), method = 'pearson')
    cosi <- coop::tcosine(prototypes)
    spear <- stats::cor(t(prototypes), method = 'spearman')
    eucl <- as.matrix(dist(prototypes))
    S <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear), eucl))
    
  }
  res <- FCPS::HierarchicalClustering(S, ClusterNo = numClusters,
                                      Type = 'AverageL')$Cls
  print(paste('=====================Now mapping clusters to data ====================='))
  
  cluster_assignment <- res[som_model$classif]
  return(cluster_assignment)
}


ComputeMetrics <- function(trueLabels, clusterLabels){
  
  ARI = aricode::ARI(trueLabels, clusterLabels)
  NMI = aricode::NMI(trueLabels, clusterLabels)
  FM = dendextend::FM_index(trueLabels, clusterLabels)[[1]]
  Fmeasure <-  FlowSOM::FMeasure(as.numeric(as.factor(clusterLabels)), 
                                 as.numeric(as.factor(trueLabels)))
  
  return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM, FMeasure = Fmeasure))
}

RunSOM <- function(data, markers, normalization = 'None', labels=F){
  # print(paste('=====================Now Generating Prototypes  ========================'))
  # som_model <- GeneratePrototypes(data, markers = markers, 
  #                                 normalization = normalization)
  k = length(unique(data$CellType))
  # 
  # print(paste('=====================Now running on Pearson and Cosine================'))
  # # Cosine
  # both.clusters <- ClusterPrototypes(som_model, distmet = 'all', 
  #                                    numClusters = k)
  
  clusters <- FuseSOM::runFuseSOM(data, markers = markers, 
                                  numClusters = k, verbose = F)$clusters
  both.metrics <- ComputeMetrics(data$CellType, clusters)
  
  if(labels){
    return(list(Metrics=both.metrics, Labels=clusters))
  }
  
  return(both.metrics)
}

RunDistAnalysis <- function(data, markers){
  print(paste('=====================Now Generating Prototypes  ========================'))
  som_model <- GeneratePrototypes(data, markers = markers)
  k = length(unique(data$CellType))
  
  # pearson
  print(paste('=====================Now running on Pearson ========================'))
  pear.clusters <- ClusterPrototypes(som_model, distmet = 'pearson', 
                                     numClusters = k)
  pear.metrics <- ComputeMetrics(data$CellType, pear.clusters) %>% bind_rows()
  pear.metrics$Method <- 'Pearson'
  
  # cosine
  print(paste('=====================Now running on Cosine================'))
  cosine.clusters <- ClusterPrototypes(som_model, distmet = 'cosine', 
                                       numClusters = k)
  cosine.metrics <- ComputeMetrics(data$CellType, cosine.clusters) %>% bind_rows()
  cosine.metrics$Method <- 'Cosine'
  
  # spearman
  print(paste('=====================Now running on Spearman================'))
  spearman.clusters <- ClusterPrototypes(som_model, distmet = 'spearman', 
                                         numClusters = k)
  spearman.metrics <- ComputeMetrics(data$CellType, spearman.clusters) %>% bind_rows()
  spearman.metrics$Method <- 'Spearman'
  
  # euclidean
  print(paste('=====================Now running on Euclidean================'))
  
  eucl.clusters <- ClusterPrototypes(som_model, distmet = 'euclidean', 
                                     numClusters = k)
  eucl.metrics <- ComputeMetrics(data$CellType, eucl.clusters) %>% bind_rows()
  eucl.metrics$Method <- 'Euclidean'
  # pear_cos
  print(paste('=====================Now running on Pearson and Cosine================'))
  pear_cos.clusters <- ClusterPrototypes(som_model, distmet = 'pear_cos', 
                                     numClusters = k)
  pear_cos.metrics <- ComputeMetrics(data$CellType, pear_cos.clusters) %>% bind_rows()
  pear_cos.metrics$Method <- 'Pearson + Cosine'
  
  # pear_spear
  print(paste('=====================Now running on Pearson and Spearman================'))
  pear_spear.clusters <- ClusterPrototypes(som_model, distmet = 'pear_spear', 
                                           numClusters = k)
  pear_spear.metrics <- ComputeMetrics(data$CellType, pear_spear.clusters) %>% bind_rows()
  pear_spear.metrics$Method <- 'Pearson + Spearman'
  
  # cos_spear
  print(paste('=====================Now running on Cosine and Spearman================'))
  cos_spear.clusters <- ClusterPrototypes(som_model, distmet = 'cos_spear', 
                                          numClusters = k)
  cos_spear.metrics <- ComputeMetrics(data$CellType, cos_spear.clusters) %>% bind_rows()
  cos_spear.metrics$Method <- 'Cosine + Spearman'
  
  # pear_eucl
  print(paste(('=====================Now running on Pearson and Euclidean================')))
  pear_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'pear_eucl', 
                                          numClusters = k)
  pear_eucl.metrics <- ComputeMetrics(data$CellType, cos_spear.clusters) %>% bind_rows()
  pear_eucl.metrics$Method <- 'Pearson + Euclidean'
  
  # spear_eucl
  print(paste(('=====================Now running on Spearman and Euclidean================')))
  pear_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'spear_eucl', 
                                          numClusters = k)
  spear_eucl.metrics <- ComputeMetrics(data$CellType, pear_eucl.clusters) %>% bind_rows()
  spear_eucl.metrics$Method <- 'Spearson + Euclidean'
  
  # cos_eucl
  print(paste(('=====================Now running on Cosine and Euclidean================')))
  cos_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'cos_eucl', 
                                          numClusters = k)
  cos_eucl.metrics <- ComputeMetrics(data$CellType, cos_eucl.clusters) %>% bind_rows()
  cos_eucl.metrics$Method <- 'Cosine + Euclidean'
  
  # pear_spear_cos
  print(paste(('=====================Now running on Pearson, Spearman and Cosine================')))
  pear_spear_cos.clusters <- ClusterPrototypes(som_model, distmet = 'pear_spear_cos', 
                                         numClusters = k)
  pear_spear_cos.metrics <- ComputeMetrics(data$CellType, pear_spear_cos.clusters) %>% bind_rows()
  pear_spear_cos.metrics$Method <- 'Pearson + Spearman + Cosine'
  
  # pear_spear_eucl
  print(paste(('=====================Now running on Pearson, Spearman and Euclidean================')))
  pear_spear_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'pear_spear_eucl', 
                                               numClusters = k)
  pear_spear_eucl.metrics <- ComputeMetrics(data$CellType, pear_spear_eucl.clusters) %>% bind_rows()
  pear_spear_eucl.metrics$Method <- 'Pearson + Spearman + Euclidean'
  
  # spear_cos_eucl
  print(paste(('=====================Now running on Spearman, Cosine and Euclidean================')))
  spear_cos_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'spear_cos_eucl', 
                                                numClusters = k)
  spear_cos_eucl.metrics <- ComputeMetrics(data$CellType, spear_cos_eucl.clusters) %>% bind_rows()
  spear_cos_eucl.metrics$Method <- 'Spearman + Cosine + Euclidean'
  
  # pear_cos_eucl
  print(paste(('=====================Now running on Pearson, Cosine and Euclidean================')))
  pear_cos_eucl.clusters <- ClusterPrototypes(som_model, distmet = 'pear_cos_eucl', 
                                               numClusters = k)
  pear_cos_eucl.metrics <- ComputeMetrics(data$CellType, spear_cos_eucl.clusters) %>% bind_rows()
  pear_cos_eucl.metrics$Method <- 'Pearson + Cosine + Euclidean'
  
  # everything
  print(paste('=====================Now running on Pearson, Spearman and Cosine================'))
  all.clusters <- ClusterPrototypes(som_model, distmet = 'all', 
                                    numClusters = k)
  all.metrics <- ComputeMetrics(data$CellType, all.clusters) %>% bind_rows()
  all.metrics$Method <- 'All'
  
  print(paste('==================Now Combining Results========================='))
  # combine the results
  res <- rbind(pear.metrics, cosine.metrics, spearman.metrics, eucl.metrics, 
               pear_cos.metrics, pear_spear.metrics, pear_eucl.metrics,
               spear_eucl.metrics, cos_eucl.metrics, pear_spear_cos.metrics,
               pear_spear_eucl.metrics, spear_cos_eucl.metrics,
               spear_cos_eucl.metrics, pear_cos_eucl.metrics, all.metrics)
  return(res)
}

subsetMatrixFull <- function(mat, p = 0.8, rep = 5, seed = 1) {
  # subset the matrix
  tmp <- colnames(mat)
  mat <- mat[,order(tmp)]
  colnames(mat) <- tmp[order(tmp)]
  subsetIndex <- c()
  sub.mat <- list()
  
  set.seed(seed)
  for (i in 1:rep) {
    subsetIndex[i] <- caret::createDataPartition(colnames(mat), p = p)
    sub.mat[[i]] <- mat[,subsetIndex[[i]]]
  }
  return (sub.mat)
}

ComputeMinClusterSize <- function(som_model){
  prop <- max(table(som_model$classif))/sum(table(som_model$classif))
  nrows <- nrow(som_model$prototypes)
  n_min <- max(5,round(prop*nrows))
  
  return(n_min)
}

ComputeNumClusters <- function(x, minClusterSize, seed = 1994, alpha = 0.001){
  set.seed(seed)
  n <- nrow(x)
  p <- ncol(x)
  nd_type <- rep("", n-1)
  p_emp <- rep(0, n-1)
  
  numClusters <- 0
  
  # perform consensus clustering
  cat("\nNow Computing The Number Of Clusters\n")
  output <- fastcluster::hclust(dist(x, method = 'maximum'), method = "average")
  
  # do some processing
  hc_dat <- output
  idx_hc <- .idx_hc(output, n)
  cutoff <- .fwer_cutoff.matrix(idx_hc, alpha)
  pd_map <- .pd_map(output, n)
  nd_type <- rep("", n-1)
  
  # run significance testing on each node
  for (k in 1:(n-1)) {
    ## indices for subtree
    idx_vals <- idx_hc[k, ]
    idx_sub <- unlist(idx_hc[k, ])
    n_sub <- length(idx_sub)
    
    ## only calc p-values for branches w/ more than n_mini
    if (n_sub < minClusterSize) {
      nd_type[k] <- "n_small"
      next
    }
    
    if ((alpha < 1) && (k > 1) && (nd_type[pd_map[k]] != "sig")) {
      nd_type[k] <- "no_test"
      p_emp[k] <- 1
      next
    }
    
    # Generate initial assingments
    t <- c(idx_vals[[1]], idx_vals[[2]])
    x_comb <- x[t,t]
    assignments <- kmeans(x_comb, 2)$cluster
    # assignments <- cutree(fastcluster::hclust(dist(x_comb, method = 'maximum'), method = "ward"),2)
    
    
    # Compute the pvalues
    # compute the discriminant projections
    x_new <- fpc::discrcoord(x=x_comb, clvecd = assignments)$proj[, 1]
    res.pval <- dip.test(x_new)$p.value
    
    # print(res.pval)
    
    # update results
    if(alpha < 1){
      if(res.pval < alpha ){
        nd_type[k] <- "sig"
        numClusters = numClusters + 1
      } else{
        nd_type[k] <- "not_sig"
      }
      p_emp[k] <- res.pval
    }
  }
  return(numClusters)
}


EstimateK <- function(som_model, # the self organizing map model
                      method = c('Discriminant','Distance', 'Stability'),# the method to use
                      kseq = 2:20 # number of clusters to access
)
  {
  # extract the prototypes from the model
  prototypes <- som_model$prototypes
  
  k_discr <- NULL
  if('Discriminant' %in% method){
    
    print('==================Now Computing The Number Of Clusters Using Discriminant Analysis==================')
    nmin <- 20
    
    pear <- stats::cor(t(prototypes), method = 'pearson')
    cosi <- coop::tcosine(prototypes)
    spear <- stats::cor(t(prototypes), method = 'spearman')
    
    S <- as.matrix(fuse(cor2dist(pear),cor2dist(cosi),cor2dist(spear)))
    
    k_discr = ComputeNumClusters(S,nmin)
  }
  
  k_dist <- NULL
  if('Distance' %in% method){
    print('==================Now Computing The Number Of Clusters Using Distance Analysis==================')
    k_dist <- cDistance(prototypes, kseq = kseq)
    k_dist$k_Gap <- ComputeElbow(k_dist$Gaps)
    k_dist$k_Slope <- ComputeElbow(k_dist$Slopes)
    k_dist$k_Jump <- ComputeElbow(k_dist$Jumps)
    k_dist$k_WCD <- ComputeElbow(k_dist$WCD)
    k_dist$k_Sil <- ComputeElbow(k_dist$Silhouettes)
  }
  
  k_stab <- NULL
  if('Stability' %in% method){
    print('==================Now Computing The Number Of Clusters Using Stability Analysis==================')
    k_stab <- cStability(prototypes,kseq = kseq, nB = 10)
  }
  
  outlist <- list('Discriminant' = k_discr,
                  'Distance' = k_dist,
                  'Stability' = k_stab)
  
  return(outlist)

}

ComputeElbow <- function(vals){
  diffs <- diff(vals)
  diffs <- diffs[-1]
  optkb <- which.max(abs(diffs)) + 1
  return(optkb)
}

ClustersRE <- function(k, k_est){
  RE_Discr <- (k - k_est$Discriminant)/k
  RE_Gap <- (k - k_est$Distance$kGap)/k
  RE_Slope <- (k - k_est$Distance$kSlope)/k
  RE_Jump <- (k - k_est$Distance$kJump)/k
  RE_WCD <- (k - k_est$Distance$kWCD)/k
  RE_Sil <- (k - k_est$Distance$kSil)/k
  
  return(list(Discrimant = RE_Discr, Gap = RE_Gap, Slope = RE_Slope,
              Jump = RE_Jump, WCD = RE_WCD, Silhouette = RE_Sil))
  
}

RunFuseSOM <- function(data, assay = NULL, markers = NULL, numClusters){
  
  # if we have a dataframe or a matrix
  if(class(data) %in% c('data.frame', 'matrix')){
    message("You have provided a dataset of class data.frame or matrix")
    # if no markers are given, make sure all the columns are numeric
    if(is.null(markers)){
      num_numeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
      if(num_numeric != ncol(data)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      }
    }else{
      # extract the markers of interest
      data <- data[, markers]
    }
  }
  
  # if we have a single cell experiment object
  if(class(data) == "SingleCellExperiment"){
    message("You have provided a dataset of class SingleCellExperiment")
    # make sure an assay is provided
    if(is.null(assay)){
      stop("If a SingleCellExperiment, make sure the appropriate assay is provided as well")
    }
    data <- t(assay(data, assay))
    
    # again if no markers are given, make sure all the columns are numeric
    if(is.null(markers)){
      num_numeric  <- sum(apply(data, 2, function(x) is.numeric(x)))
      if(num_numeric != ncol(data)){
        stop("If markers of interest are not provided, make sure the data contains all numeric columns")
      }
    }else{
      # extract the markers of interest
      data <- data[, markers]
    }
    
    
  }
  
  # now we can run the FuseSOM algorithm
  message("Everything looks good. Now running the FuseSOM algorithm")
  data <- apply(data, 2, function(x) as.numeric(x))
  som_model <- GeneratePrototypes(data)
  clusters <- ClusterPrototypes(som_model, numClusters = numClusters)
  
  message("The FuseSOM algorithm has completed successfully")
  
  return(list(model = som_model, clusters = clusters))
  
}

run_prox <- function(dat, markers, prox_metrics){
  res.ari <- vector(mode = "numeric", length = 3)
  res.nmi <- vector(mode = "numeric", length = 3)
  res.fm <- vector(mode = "numeric", length = 3)
  res.fmeasure <- vector(mode = "numeric", length = 3)
  celltype <- dat$CellType
  dat <- dat[, markers]
  k = length(unique(celltype))
  for(i in 1:length(prox_metrics)){
    d <- parallelDist(as.matrix(dat), method = prox_metrics[[i]])
    labels <- cutree(fastcluster::hclust(d, method = 'average'),k)
    
    ARI <- aricode::ARI(labels, celltype)
    NMI <- aricode::NMI(labels, celltype)
    FM <- dendextend::FM_index(labels, celltype)
    FMeasure <- FlowSOM::FMeasure(as.numeric(labels), 
                                  as.numeric(as.factor(celltype)))
    
    res.ari[[i]] <- ARI
    res.nmi[[i]] <- NMI
    res.fm[[i]] <- FM
    res.fmeasure[[i]] <- FMeasure
  }
  
  return(list(ARI = res.ari, NMI = res.nmi, `FM-Index` = res.fm,
              FMeasure = res.fmeasure, Metric = prox_metrics))
  
}

run_corr <- function(dat, markers, corr_metrics){
  res.ari <- vector(mode = "numeric", length = 3)
  res.nmi <- vector(mode = "numeric", length = 3)
  res.fm <- vector(mode = "numeric", length = 3)
  res.fmeasure <- vector(mode = "numeric", length = 3)
  celltype <- dat$CellType
  dat <- as.matrix(dat[, markers])
  k = length(unique(celltype))
  for(i in 1:length(corr_metrics)){
    if(corr_metrics[[i]] == 'cosine'){
      d <- as.dist(cosineDist(dat))
    } else{
      d <- as.dist(cor.dist(dat, methods = corr_metrics[[i]]))
    }
    
    labels <- cutree(fastcluster::hclust(d, method = 'average'),k)
    
    ARI <- aricode::ARI(labels, celltype)
    NMI <- aricode::NMI(labels, celltype)
    FM <- dendextend::FM_index(labels, celltype)
    FMeasure <- FlowSOM::FMeasure(as.numeric(labels), 
                                  as.numeric(as.factor(celltype)))
    
    res.ari[[i]] <- ARI
    res.nmi[[i]] <- NMI
    res.fm[[i]] <- FM
    res.fmeasure[[i]] <- FMeasure
  }
  return(list(ARI = res.ari, NMI = res.nmi, `FM-Index` = res.fm,
              FMeasure = res.fmeasure, Metric = corr_metrics))
  
}

run_pheno <- function(dat, markers, metric){
  # res.ari <- vector(mode = "numeric", length = 2)
  # res.nmi <- vector(mode = "numeric", length = 2)
  # res.fm <- vector(mode = "numeric", length = 2)
  # res.fmeasure <- vector(mode = "numeric", length = 2)
  celltype <- dat$CellType
  dat <- as.matrix(dat[, markers])
  k = length(unique(celltype))
  
  labels <- as.vector(PhenoRun(dat, distmet = metric, markers = markers))
  
    
  ARI <- aricode::ARI(labels, celltype)
  NMI <- aricode::NMI(labels, celltype)
  FM <- dendextend::FM_index(labels, celltype)
  FMeasure <- FlowSOM::FMeasure(as.numeric(labels), 
                                  as.numeric(as.factor(celltype)))
    
    # res.ari[[i]] <- ARI
    # res.nmi[[i]] <- NMI
    # res.fm[[i]] <- FM
    # res.fmeasure[[i]] <- FMeasure
  return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM,
              FMeasure = FMeasure, Metric = metric))
  
  
}

run_flowsom <- function(dat, markers, metric){
  res.ari <- vector(mode = "numeric", length = 2)
  res.nmi <- vector(mode = "numeric", length = 2)
  res.fm <- vector(mode = "numeric", length = 2)
  res.fmeasure <- vector(mode = "numeric", length = 2)
  celltype <- dat$CellType
  # dat <- as.matrix(dat[, markers])
  k = length(unique(celltype))
  ff <- DFtoFF(as.data.frame(dat[, markers]))
  flo_res <- FlowSOM_custom(ff,scale=F,nClus=length(unique(dat$CellType)), 
                           distance = metric, clusterAlg = 'hc')
  labels <- GetMetaclusters(flo_res)

  ARI = aricode::ARI(labels, celltype)
  NMI = aricode::NMI(labels, celltype)
  FM = dendextend::FM_index(labels, celltype)[[1]]
  Fmeasure <-  FlowSOM::FMeasure(as.numeric(labels), as.numeric(as.factor(celltype)))
  
  return(list(ARI = ARI, NMI = NMI, `FM-Index` = FM,
              FMeasure = Fmeasure, Metric = metric))
  
  
}

manhattan_dist = function(x){ dist(x,method="manhattan")}
maximum_dist = function(x){ dist(x,method="maximum")}
cosine_dist = function(x){as.dist(cosineDist(x))}


runClustNum <- function(dat, markers){
  # generate the model
  som_model <- FuseSOM::generatePrototypes(data = dat[, markers])
  
  # get the true number of clusters
  k_true <- length(unique(dat$CellType))
  
  # estimate the number of clusters
  k_est <- FuseSOM::estimateNumCluster(som_model, kSeq = 2:(k_true+5))
  
  # compute the relative error
  message('Now computing the cluster estimation relative error')
  RE <- ClustersRE(k_true, k_est) %>% bind_rows()
  
  # compute the metrics based on the estimated number of clusters
  message('Now computing the cluster estimation performance metrics')
  metrics_disc <- ComputeMetrics(dat$CellType,
                                          FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Discriminant)) %>% 
    bind_rows()
  metrics_gap <- ComputeMetrics(dat$CellType, 
                                         FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Distance$kGap)) %>% 
    bind_rows()
  metrics_slope <- ComputeMetrics(dat$CellType,
                                           FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Distance$kSlope)) %>% 
    bind_rows()
  metrics_jump <- ComputeMetrics(dat$CellType,
                                          FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Distance$kJump)) %>% 
    bind_rows()
  metrics_wcd <- ComputeMetrics(dat$CellType,
                                         FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Distance$kWCD)) %>% 
    bind_rows()
  metrics_sil <- ComputeMetrics(dat$CellType,
                                         FuseSOM::clusterPrototypes(som_model, numClusters = k_est$Distance$kSil)) %>% 
    bind_rows()
  
  message('Done')
  res <- rbind(metrics_disc, metrics_gap, metrics_slope,
                        metrics_jump, metrics_wcd, metrics_sil)
  res$Method <- c('Discriminant', 'Gap', 'Slope', 'Jump', 'WCD',
                           'Silhouette')
  res <- cbind(res, t(RE))
  
  colnames(res) <- c("ARI","NMI","FM-Index","FMeasure","Method" ,"RE" )
  
  return(res)
  
}

split_by_pca = function(x) {
  x = scale(x)
  res = svd(x)
  v = res$v
  pc_mat = res$u[, 1 : 1, drop = F]
  res = x - pc_mat %*% (t(pc_mat) %*% x)
  return(res)
}

customWeighting <- function(method, data, markers, prototypes){
  # generate the prototypes
  
  # compute the distances on the generated prototypes
  pear <- cor(t(prototypes$prototypes), method = "pearson")
  cosi <- tcosine(prototypes$prototypes)
  spear <- cor(t(prototypes$prototypes), method = "spearman")
  eucl <- dist(prototypes$prototypes) %>%
    as.matrix()
  
  message(paste('Now running for method:',method))
  
  if(method == 'pca'){
    # Vectorize the distance matrices and combine
    combined_matrix <- rbind(as.vector(pear), 
                             as.vector(cosi),
                             as.vector(spear),
                             as.vector(eucl))
    
    # Perform PCA
    pca_result <- prcomp(combined_matrix)
    
    # The loadings of the first principal component will be our weights
    weights <- pca_result$sdev
    
    # Normalize the weights to sum up to 1
    weights <- weights / sum(weights)
    
    # Combine distance matrices using the weights
    # finalDist <- weights[1] * cor2dist(pear) + weights[2] * cor2dist(cosi) + 
    #   weights[3] * cor2dist(spear) + weights[4] * eucl
    
    finalDist <- as.matrix(fuse(cor2dist(pear), cor2dist(cosi),
                                cor2dist(spear), eucl, weights = weights))
  } else if(method == 'fuse'){
    finalDist <- as.matrix(fuse(cor2dist(pear), cor2dist(cosi),
                                cor2dist(spear), eucl))
  } else {
    # cluster the distances individually
    cluster_pear <- HierarchicalClustering(cor2dist(pear),
                                           ClusterNo = length(unique(data$CellType)),
                                           Type = "AverageL")$Cls
    cluster_cosi <- HierarchicalClustering(cor2dist(cosi),
                                           ClusterNo = length(unique(data$CellType)),
                                           Type = "AverageL")$Cls
    cluster_spear <- HierarchicalClustering(cor2dist(spear),
                                            ClusterNo = length(unique(data$CellType)),
                                            Type = "AverageL")$Cls
    cluster_eucl <- HierarchicalClustering(eucl,
                                           ClusterNo = length(unique(data$CellType)),
                                           Type = "AverageL")$Cls
    
    # compute the scores for each of the metho
    score1 <- cluster.stats(as.dist(cor2dist(pear)), cluster_pear)[[method]]
    score2 <- cluster.stats(as.dist(cor2dist(cosi)), cluster_cosi)[[method]]
    score3 <- cluster.stats(as.dist(cor2dist(spear)), cluster_spear)[[method]]
    score4 <- cluster.stats(as.dist(eucl), cluster_eucl)[[method]]
    
    # Normalize the scores to sum to 1
    scores_all <- c(score1, score2, score3, score4)
    weights <- scores_all/sum(scores_all)
    
    # Combine distance matrices using the weights
    # finalDist <- weights[1] * cor2dist(pear) + weights[2] * cor2dist(cosi) + 
    #   weights[3] * cor2dist(spear) + weights[4] * eucl
    
    finalDist <- as.matrix(fuse(cor2dist(pear), cor2dist(cosi),
                                cor2dist(spear), eucl, weights = weights))
    
  }
  
  # Then perform clustering on the combined matrix
  clusters <- HierarchicalClustering(finalDist, 
                                     ClusterNo = length(unique(data$CellType)), 
                                     Type = "AverageL")$Cls
  
  clusters <- clusters[prototypes$classif]
  
  res <- ComputeMetrics(data$CellType, clusters)
  
  return(res)
}


