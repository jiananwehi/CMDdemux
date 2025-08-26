#' Calculate Mahalanobis distance
#'
#' Calculate the Mahalanobis distance of each cell to its cluster centroid.
#' 
#' @param clr.norm The local CLR normalized data comes from the output of the 
#'        \link{LocalCLRNorm} function.
#' @param core.drop.assign Classification of core and non-core cells from the 
#'        output of the \link{DefineNonCore} function.
#' @param kmed.cl K-Medoids clustering results from the output of the 
#'        \link{KmedCluster} function.
#' @param cluster.assign The cluster assignments and corresponding hashtags from 
#'        the output of \link{LabelClusterHTO}.      
#'        
#' @details
#' The covariance matrix \eqn{\Gamma} between pairs of hashtags (\eqn{i,i'}) is
#' \deqn{
#' \Gamma_{i,i'} = \frac{\sum_l \sum_{c \in l} (h_{i,c} - M_i^{(l)}) (h_{i',c} - M_{i'}^{(l)})}{\sum_l n_l},
#' }
#' Here, \eqn{n_l} denotes the number of cells in cluster \eqn{l}, and \eqn{h_{i,c}} 
#' denotes the local CLR value of hashtag \eqn{i} in cell \eqn{c}.
#' 
#' Since \eqn{\Gamma} is singular and not invertible, we use the Moore-Penrose 
#' generalized inverse, denoted \eqn{\Gamma^+}, to compute a pseudo-inverse of 
#' \eqn{\Gamma}.
#' 
#' The Mahalanobis distance \eqn{MD_{c,l}} for each cell \eqn{c} to the centroid 
#' \eqn{M^{(l)}} of the K-medoids cluster \eqn{l} is defined as
#' \deqn{
#' MD_{c,l} = (h_c - M^{(l)})^T \Gamma^+ (h_c - M^{(l)}), \quad
#' h_c = (h_{1,c}, \dots, h_{N,c})^T.
#' }
#' 
#' @return A matrix of Mahalanobis distances, where rows correspond to hashtags 
#' and columns correspond to cells.
#' 
#' @references
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S-PLUS. 
#' Fourth Edition. Springer.
#' 
#' Mahalanobis, P. C. (2018). On the generalized distance in statistics. 
#' Sankhyā: The Indian Journal of Statistics, Series A (2008-), 80, S1-S7.
#' 
#' @examples 
#' \donttest{
#' hto.md.mat <- CalculateMD(hto.cl.mtx, hto.noncore, hto.kmed.cl, hto.cluster.assign)
#' dim(hto.md.mat)
#' hto.md.mat[,1:10]
#' }
#' 
#' @importFrom MASS ginv
#' @importFrom stats mahalanobis
#' 
#' @export
CalculateMD <- function(clr.norm, core.drop.assign, kmed.cl, cluster.assign){
  cov.mat <- matrix(0, nrow = nrow(clr.norm), ncol = nrow(clr.norm))
  for(i in 1:nrow(cov.mat)){
    for(j in 1:ncol(cov.mat)){
      cluster.value <- lapply(unique(kmed.cl$clustering), FUN = function(x){
        sum((clr.norm[i,which(core.drop.assign %in% paste0("Cluster", x))] - kmed.cl$medoids[x,i])*(clr.norm[j,which(core.drop.assign %in% paste0("Cluster", x))] - kmed.cl$medoids[x,j]))
      })
      cov.mat[i,j] <- sum(unlist(cluster.value))
    }
  }
  cov.mat <- cov.mat/length(which(!core.drop.assign %in% "non-core"))
  rownames(cov.mat) <- rownames(clr.norm)
  colnames(cov.mat) <- rownames(clr.norm)
  pcov.mat <- ginv(cov.mat)
  md.list <- lapply(1:nrow(clr.norm), FUN = function(x){
    center.barcode <- rownames(kmed.cl$medoids)[cluster.assign[x]]
    center <- clr.norm[,center.barcode]
    md <- sqrt(mahalanobis(t(clr.norm), center = center, cov = pcov.mat, inverted = TRUE))
  })
  md.mat <- matrix(unlist(md.list), nrow = nrow(clr.norm), byrow = TRUE)
  colnames(md.mat) <- colnames(clr.norm)
  rownames(md.mat) <- rownames(clr.norm)
  return(md.mat)
}