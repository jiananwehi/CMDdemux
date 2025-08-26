#' Local CLR normalization
#'
#' Normalize the hashing data using a local (within-cell) CLR strategy.  
#'
#' @param hash.count The raw hashing count matrix has rows representing hashtags 
#'        and columns representing cells.
#' @param pseudo_count A pseudo count is added to each entry of the matrix to 
#'        avoid the problem of 0 in the denominator. Default value is 1. 
#'        
#' @details
#' For hashtag \eqn{i} and cell \eqn{c}, the HTO counts are normalized as:
#' \deqn{\pi_{i,c} = \frac{HTO_{i,c} + 1}{HTO_{+,c} + N}, \quad i = 1, \dots, N}
#' where \eqn{HTO_{+,c} = \sum_{i} HTO_{i,c}} and \eqn{N} is the total number of hashtags.
#' A pseudo count of 1 is added here.
#'
#' The geometric mean is calculated as:
#' \deqn{GM_{\pi_c} = \left( \prod_{i=1}^N \pi_{i,c} \right)^{\frac{1}{N}}}
#'
#' The local CLR value is then:
#' \deqn{h_{i,c} = \log \left( \frac{\pi_{i,c}}{GM_{\pi_c}} \right)}
#' 
#' @return A normalized hash count matrix with rows representing hashtags and 
#' columns representing cells.
#'  
#' @examples 
#' \donttest{
#' hto.clr.mtx <- LocalCLRNorm(hto.mtx)
#' dim(hto.clr.mtx)
#' hto.clr.mtx[, 1:10]
#' }
#' 
#' @importFrom Matrix colSums
#' 
#' @export
LocalCLRNorm <- function(hash.count, pseudo_count = 1){
  hash.count <- as.matrix(hash.count)
  if(is.null(rownames(hash.count))){
    rownames(hash.count) <- paste0("HTO", 1:nrow(hash.count))
  }
  if(is.null(colnames(hash.count))){
    colnames(hash.count) <- paste0("Cell", 1:ncol(hash.count))
  }
  lib.norm <- t(t((as.matrix(hash.count+pseudo_count)))/(colSums(hash.count+pseudo_count)))
  lib.log.norm <- log(lib.norm)
  geo_mean <- apply(lib.norm, 2, FUN = function(x){ 
    mean(log(x))
  })
  clr.norm.list <- lapply(1:length(geo_mean), FUN = function(x){ 
    lib.log.norm[,x] - geo_mean[x]
  })
  clr.norm <- matrix(unlist(clr.norm.list), nrow = nrow(lib.norm), ncol = ncol(lib.norm), byrow = FALSE)
  rownames(clr.norm) <- rownames(lib.norm)
  colnames(clr.norm) <- colnames(lib.norm)
  return(clr.norm)
}
