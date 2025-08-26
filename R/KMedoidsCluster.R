#' K-Medoids clustering
#'
#' Perform K-Medoids clustering on the local CLR normalized data.
#' 
#' @param clr.norm The local CLR normalized data comes from the output of 
#'        \link{LocalCLRNorm} function.  
#' @param fast_mode A logical value indicating whether to perform fast mode. If TRUE, 
#'        CLARA clustering is performed; otherwise, PAM clustering is performed. 
#'        Fast mode (TRUE) is recommended when there are more than 10,000 cells. 
#'        Default: FALSE.
#' @param prop Proportion of cells used in CLARA clustering. This parameter will
#'        only be used when fast_mode is TRUE. Default: 0.1.
#' @param optional A logical value indicating whether extra clusters should be created. 
#'        If TRUE, optional clustering is performed, resulting in more clusters 
#'        than the number of hashes. Default: FALSE. 
#' @param extra_cluster The number of extra clusters. The total number of clusters 
#'        after clustering will equal the number of hashtags plus the number of 
#'        extra clusters. This parameter will only be used when the optional 
#'        parameter is TRUE. Default: 1.  
#' @details 
#' It is recommended to set `fast_mode` to TRUE when the hashing data contain more 
#' than 10,000 cells. In fast mode, CLARA clustering is performed, whereas in the 
#' default mode, PAM clustering is used. Since CLARA clustering uses only a 
#' proportion of cells, the `prop` parameter specifies the proportion of cells to 
#' be used for clustering.
#' 
#' Optional clustering is recommended when the number of clusters—equal to the 
#' number of hashtags—is not sufficient to achieve good demultiplexing results. 
#' This is particularly useful for low-quality data. If the number of clusters 
#' equal to the number of hashtags does not yield satisfactory demultiplexing, try 
#' gradually increasing the number of clusters by setting `extra_cluster` from 1 upward.
#' 
#' @return A list of clustering information. 
#' \describe{
#'   \item{medoids}{Medoids of clusters.}
#'   \item{clustering}{Cluster of each cell.}
#'   \item{\dots}{Other outputs from the underlying `pam` or `clara` function, such as `id.med`, `objective`, `isolation`, etc.}
#' }
#' 
#' @references 
#' Maechler M, Rousseeuw P, Struyf A, Hubert M, Hornik K (2025). 
#' cluster: Cluster Analysis Basics and Extensions. 
#' 
#' @examples 
#' \donttest{
#' hto.kmed.cl <- KmedCluster(hto.clr.mtx)
#' # Perform clustering in fast mode when there are more than 100,000 cells or when the default clustering is too slow
#' hto.kmed.cl <- KmedCluster(hto.clr.mtx, fast_mode = TRUE)
#' # PPerform optional clustering when more clusters than hashtags are required
#' hto.kmed.cl <- KmedCluster(hto.clr.mtx, optional = TRUE)
#' # Perform optional clustering by specifying the number of extra clusters
#' hto.kmed.cl <- KmedCluster(hto.clr.mtx, optional = TRUE, extra_cluster = 3)
#' hto.kmed.cl$medoids
#' head(hto.kmed.cl$clustering)
#' }
#' 
#' @importFrom cluster pam clara
#' @export
KmedCluster <- function(clr.norm, fast_mode = FALSE, prop = 0.1, optional = FALSE, extra_cluster = 1){
  if(fast_mode){
    hash.kmded.cluster <- clara(t(clr.norm), k=nrow(clr.norm), metric = "euclidean", samples = prop*ncol(clr.norm))
  }else{
    hash.kmded.cluster <- pam(t(clr.norm), nrow(clr.norm), metric = "euclidean")
  }
  if(isTRUE(optional)){
    if(fast_mode){
      hash.kmded.cluster <- clara(t(clr.norm), k=nrow(clr.norm)+extra_cluster, metric = "euclidean", samples = prop*ncol(clr.norm))
    }else{
      hash.kmded.cluster <- pam(t(clr.norm), nrow(clr.norm)+extra_cluster, metric = "euclidean")
    }
  }
  return(hash.kmded.cluster)
}


#' Euclidean distance to the cluster centroid 
#'
#' Calculate the Euclidean distance from each cell to its assigned cluster centroid.
#' @param clr.norm The local CLR normalized data comes from the output of 
#'        \link{LocalCLRNorm} function.  
#' @param kmed.cl K-Medoids clustering information from the output of
#'        \link{KmedCluster} function. 
#'        
#' @details
#' The Euclidean distance \eqn{d} between each cell \eqn{c} to its cluster \eqn{l}
#' centroid \eqn{M^{(l)}} is:
#' \eqn{d = \| h_c - M^{(l)} \|_2, \quad c \in l}         
#' 
#' @return A list of Euclidean distances. The length of the list equals the number 
#' of clusters, and each element corresponds to the Euclidean distances of the 
#' cells in that cluster. 
#' 
#' @examples 
#' \donttest{
#' hto.cl.dist <- KmedCluster(hto.clr.mtx, hto.kmed.cl)
#' length(hto.cl.dist)
#' head(hto.cl.dist[[1]])
#' }   
#' 
#' @export 
EuclideanClusterDist <- function(clr.norm, kmed.cl){
  eu.dist <- apply(kmed.cl$medoids, 1, FUN = function(x){ 
    sqrt(colSums((clr.norm - matrix(rep(x, ncol(clr.norm)), nrow = nrow(clr.norm), byrow = FALSE))^2))
  })
  cl.dist <- lapply(sort(unique(kmed.cl$clustering)), FUN = function(x){ 
    eu.dist[which(kmed.cl$clustering %in% x),x]
  })
  return(cl.dist)
}


#' Examine the number of clusters
#' 
#' Examine the number of clusters that should be used in the data. 
#' 
#' @param clr.norm The local CLR normalized data comes from the output of 
#'        \link{LocalCLRNorm} function.
#' @param kmed.cl1 K-Medoids clustering results with fewer clusters (smaller k) 
#'        from the output of the \link{KmedCluster} function. 
#' @param kmed.cl2 K-Medoids clustering results with more clusters (larger k) 
#'        from the output of the \link{KmedCluster} function. 
#' @param fc_cut Cut-off for determining the number of clusters. It is defined as 
#'        the log fold change of the CLR expression of a specific hashtag between 
#'        the top clusters for two k values. Default: 1. 
#'        
#' @details
#' For each hashtag \eqn{i}, the median local CLR value \eqn{h} across all clusters 
#' \eqn{l} is defined as
#' \deqn{med_{i,l} = median_{c \in l}\{ h_{i,c} \}}
#'
#' The top cluster with the largest median CLR is defined as
#' \deqn{l_{i,\max,k} = \arg\max_{l=1,\dots,k} \{ med_{i,l} \}}
#'
#' Two different numbers of clusters (\eqn{k_2 > k_1}) are then compared by 
#' their differences in CLR values for each hashtag between the top clusters,
#' \deqn{f_i = med_{i,l_{i,\max,k_2}} - med_{i,l_{i,\max,k_1}}, \quad k_2 > k_1}
#'
#' If the maximum of \eqn{f_i} across all hashtags,
#' i.e., \deqn{\max_{i=1,\dots,N} \{ f_i \}},
#' exceeds a given threshold (default: 1), then \eqn{k_2} is recommended.
#' 
#' @return A message is printed indicating the largest CLR difference, its value, 
#' and the number of clusters selected.
#' 
#' @examples 
#' \donttest{
#' hto.kmed.cl1 <- KmedCluster(hto.clr.mtx)
#' hto.kmed.cl2 <- KmedCluster(hto.clr.mtx, optional = TRUE, extra_cluster = 1)
#' ExamineClusterNumber(hto.clr.mtx, hto.kmed.cl1, hto.kmed.cl2)
#' }   
#' 
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' 
#' @export 
ExamineClusterNumber <- function(clr.norm, kmed.cl1, kmed.cl2, fc_cut = 1){
  cl.clr.fc <- apply(clr.norm, 1, FUN = function(x){
    cl1.cluster.df <- data.frame("exp" = x, "cluster" = kmed.cl1$clustering)
    cl2.cluster.df <- data.frame("exp" = x, "cluster" = kmed.cl2$clustering)
    cl1.cluster.df2 <- cl1.cluster.df %>% 
      group_by(cluster) %>%
      summarise("m_exp" = median(exp))  
    cl1.top.clr <- max(cl1.cluster.df2$m_exp)
    cl2.cluster.df2 <- cl2.cluster.df %>% 
      group_by(cluster) %>%
      summarise("m_exp" = median(exp)) 
    cl2.top.clr <- max(cl2.cluster.df2$m_exp)
    return(cl2.top.clr-cl1.top.clr)
  }) 
  if(max(cl.clr.fc) > fc_cut){
    print(paste0("The largest CLR difference is:", names(which.max(cl.clr.fc)), " ", round(max(cl.clr.fc), 3), ". Choose optional clustering method (kmed.cl2 with", length(unique(kmed.cl2$clustering)), " clusters)."))
  }else if(min(cl.clr.fc) < -fc_cut){
    print(paste0("The smallest CLR difference is:", names(which.min(cl.clr.fc)), " ", round(min(cl.clr.fc), 3), ". Choose default clustering method (kmed.cl1 with", length(unique(kmed.cl1$clustering)), " clusters)."))
  }else{
    print(paste0("Use the default of clustering.(kmed.cl1 with", length(unique(kmed.cl1$clustering)), " clusters)"))
  }  
}


