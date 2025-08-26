#' Define non-core cells
#'
#' Non-core cells are defined based on their Euclidean distances to cluster centroids.
#' 
#' @param cl.dist The Euclidean distance for each cell to its cluster centroid 
#'        from the output of \link{EuclideanClusterDist} function.
#' @param kmed.cl K-Medoids clustering results from the output of 
#'        \link{KmedCluster} function.
#' @param eu_cut_q A vector of cut-offs for defining non-core cells using the 
#'        quantile of Euclidean distances within each cluster. Must be the same 
#'        length as the number of clusters. This parameter can be explored by 
#'        visualizing the distribution of Euclidean distances using the 
#'        \link{EuclideanDistPlot} function. Default: 0.9.
#' @param eu_cut A vector of cut-offs for defining non-core cells using the raw 
#'        Euclidean distance within each cluster. Must be the same length as the 
#'        number of clusters. Use either \code{eu_cut} or \code{eu_cut_q} to 
#'        determine the cut-offs of non-core cells. Default: NULL.
#' @param optional Logical; whether optional clustering has been performed to 
#'        obtain more clusters than the number of hashtags. This parameter should 
#'        be the same as the \code{optional} parameter in the \link{KmedCluster} 
#'        function. Default: FALSE.
#' @param clr.norm Local CLR-normalized data from the output of \link{LocalCLRNorm} 
#'        function. This parameter only needs to be set when \code{optional = TRUE}. 
#'        Default: NULL.
#'        
#' @details 
#' Core cells are those with smaller Euclidean distances to their cluster centroids, 
#' whereas non-core cells have larger Euclidean distances. The cut-off for defining 
#' non-core cells can be based either on the quantile of Euclidean distances or 
#' on a raw Euclidean distance value. Cells with Euclidean distances larger than 
#' the cut-off in each cluster are defined as non-core cells. 
#' 
#' When optional clustering is performed, cells in additional clusters that do not 
#' correspond to any hashtag are also defined as non-core cells. Each hashtag should 
#' have its maximum CLR expression in one cluster. Extra clusters that do not show 
#' maximum CLR expression for any hashtag are therefore classified as non-core, 
#' and all cells within them are considered non-core cells.
#' 
#' @return A vector of cell labels. Non-core cells are labeled as "non-core", and 
#' core cells are labeled with their corresponding cluster. 
#' 
#' @examples 
#' \donttest{
#' hto.kmed.cl1 <- KmedCluster(hto.clr.mtx)
#' hto.noncore <- DefineNonCore(hto.cl.list, hto.kmed.cl1, c(0.9, 0.75, 0.94))
#' hto.kmed.cl2 <- KmedCluster(hto.clr.mtx, optional = TRUE, extra_cluster = 1)
#' # Define non-core cells in the data using optional clustering
#' hto.noncore <- DefineNonCore(hto.cl.list, hto.kmed.cl2, c(0.9, 0.95, 0.83, 0.65), optional = TRUE, clr.norm = hto.clr.norm)
#' head(hto.noncore)
#' }
#' 
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' 
#' @export 
DefineNonCore <- function(cl.dist, kmed.cl, eu_cut_q = 0.9, eu_cut = NULL, optional = FALSE, clr.norm = NULL){
  cl.new <- kmed.cl$clustering
  flag_q <- TRUE
  if(!is.null(eu_cut) & length(eu_cut) >= length(cl.dist)){
    flag_q <- FALSE
  }
  noncore.barcodes <- lapply(sort(unique(kmed.cl$clustering)), FUN = function(x){ 
    tryCatch(
      {
        if(!flag_q){
          if(eu_cut[x] <= max(cl.dist[[x]])){
            outlier <- eu_cut[x]
          }else{
            outlier <- quantile(cl.dist[[x]], probs = eu_cut_q[1]) 
          }
        }else{
          if(length(eu_cut_q) < length(cl.dist)){
            outlier <- quantile(cl.dist[[x]], probs = eu_cut_q[1]) 
          }else{
            outlier <- quantile(cl.dist[[x]], probs = eu_cut_q[x]) 
          }
        }
      }, error = function(e){
        print("Wrong input parameters for this function.")
        outlier <- quantile(cl.dist[[x]], probs = 0.9) 
      }
    )
    outlier.barcode <- names(cl.dist[[x]])[which(cl.dist[[x]] > outlier)]
    return(outlier.barcode)
  })
  if(isTRUE(optional)){
    clr.max <- apply(clr.norm, 1, FUN = function(x){
      data.df <- data.frame("exp" = x, "cl" = kmed.cl$clustering)
      data.df$cl <- factor(data.df$cl)
      cl.median <- data.df %>% 
        group_by(cl) %>%
        summarise("value" = median(exp))
      return(cl.median$cl[which.max(cl.median$value)])
    }) 
    all.cl <- unique(kmed.cl$clustering)
    outlier.cl <- all.cl[which(!all.cl %in% clr.max)]
    noncore.barcodes <- append(noncore.barcodes, list(names(kmed.cl$clustering)[which(kmed.cl$clustering %in% outlier.cl)]))
  }
  cl.new[as.vector(unlist(noncore.barcodes))] <- "non-core"
  cl.new[which(!cl.new %in% "non-core")] <- paste0('Cluster', cl.new[which(!cl.new %in% "non-core")])
  return(cl.new)
}


#' Label each cluster to its HTO
#'
#' Labels each cluster according to its sample of origin.
#' 
#' @param clr.norm The local CLR normalized data comes from the output of the 
#'        \link{LocalCLRNorm} function.
#' @param kmed.cl K-Medoids clustering results from the output of the 
#'        \link{KmedCluster} function.
#' @param core.drop.assign Classification of core and non-core cells from the 
#'        output of the \link{DefineNonCore} function.
#' @param label_method Specifies the method for labelling clusters. Options: 
#'        "medoids" (assigns clusters based on cluster medoids) or "expression" 
#'        (assigns clusters based on the maximum-expression hashtag). Default is 
#'        "medoids".
#' 
#' @details
#' In the medoids-based approach, a hashtag \eqn{i} is assigned to the cluster 
#' \eqn{l} whose medoid has the highest local CLR value \eqn{h}, i.e., 
#' \eqn{l = \arg\max_i \big(h_{i, M^{(l)}}\big)}.
#' 
#' In the expression-based approach, non-core cells are first excluded. A hashtag 
#' \eqn{i} is assigned to the cluster \eqn{l} with the highest median CLR value 
#' \eqn{h}, i.e., \eqn{l = \arg\max_i \big(\mathrm{med}_{i,l}\big)}.
#' 
#' When a cluster shows the maximum median local CLR value in two (or more) 
#' hashtags, this indicates conflicting hashtags, denoted \eqn{\mathcal{C}_i}. 
#' This means the cluster expresses multiple hashtags at high levels. To assign 
#' each conflicting hashtag \eqn{i} to a unique cluster, we calculate a normalized 
#' score \eqn{D}, defined as
#' \deqn{
#' D_i = \frac{\big| \mathrm{Med}(h_{i,l_1}) - \mathrm{Med}(h_{i,l_2}) \big|}
#'            { \tfrac{1}{2} \big( \mathrm{IQR}(h_{i,l_1}) + \mathrm{IQR}(h_{i,l_2}) \big)},
#' \quad i \in \mathcal{C}_i,
#' }
#' where \eqn{l_1} and \eqn{l_2} are the top two clusters with the highest 
#' median local CLR expression for hashtag \eqn{i}. The median local CLR value is 
#' normalized by the interquartile range (IQR) of each top cluster. The cluster 
#' is assigned to the hashtag with the maximum score \eqn{D}, i.e., 
#' \eqn{l = \arg\max_{i \in \mathcal{C}_i} D_i}. The remaining conflicting hashtags 
#' in \eqn{\mathcal{C}_i} are assigned to the unassigned clusters with the next 
#' highest median CLR expression.
#' 
#' @return A named integer vector, where values represent clusters and names 
#' correspond to hashtags. 
#' 
#' @examples 
#' \donttest{
#' # Label clusters using medoids information
#' hto.cluster.assign <- LabelClusterHTO(hto.clr.mtx, hto.kmed.cl, hto.noncore, "medoids")
#' # Label clusters using expression information
#' hto.cluster.assign <- LabelClusterHTO(hto.clr.mtx, hto.kmed.cl, hto.noncore, "expression")
#' hto.cluster.assign
#' }
#' 
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' 
#' @export 
LabelClusterHTO <- function(clr.norm, kmed.cl, core.drop.assign, label_method = "medoids"){
  medoids_fail <- FALSE
  if(label_method == "medoids"){
    cluster_hto <- apply(kmed.cl$medoids, 2, which.max)
    if(length(unique(cluster_hto)) < nrow(clr.norm)){
      medoids_fail <- TRUE
    }
  }
  if(label_method == "expression" | isTRUE(medoids_fail)){
    sub.clr.norm <- clr.norm[,which(!core.drop.assign %in% "non-core")]
    cluster_hto <- apply(sub.clr.norm, 1, FUN = function(x){
      df <- data.frame("data"= x, "cluster" = core.drop.assign[which(!core.drop.assign %in% "non-core")])
      df$cluster <- as.numeric(gsub("Cluster", "", df$cluster))
      sub.df <- df %>% group_by(cluster) %>% 
        summarise("clr" = median(data))
      return(sub.df$cluster[which.max(sub.df$clr)])
    })
    # Deal with conflicts between assignments (multiple HTOs assign to the same cluster)
    if(length(unique(cluster_hto)) < nrow(clr.norm)){
      duplicate_table <- table(cluster_hto)
      duplicate_clusters <- names(duplicate_table)[which(duplicate_table > 1)]
      update_cluster_hto <- cluster_hto
      update_cluster_hto[which(update_cluster_hto %in% duplicate_clusters)] <- NA
      for(i in duplicate_clusters){
        duplicate_htos <- names(cluster_hto)[which(cluster_hto == i)]
        duplicate_hto_medians <- apply(sub.clr.norm[duplicate_htos,], 1, FUN = function(x){
          df <- data.frame("data"= x, "cluster" = core.drop.assign[which(!core.drop.assign %in% "non-core")])
          df$cluster <- as.numeric(gsub("Cluster", "", df$cluster))
          sub.df <- df %>% group_by(cluster) %>% 
            summarise("clr" = median(data), "iqr" = iqr(data))
          return(sub.df)
        })
        # Define the duplicate cluster belongs to the HTO that has largest normalized difference of HTO expression (median(cl1)-median(cl2)/(1/2(IQR(cl1)+IQR(cl2)))) between to the top and second largest clusters
        duplicate_hto_dif <- lapply(duplicate_hto_medians, FUN = function(x){
          clr_sort_exp <- as.data.frame(x[order(x$clr, decreasing = TRUE),])
          max_second_dif <- abs(clr_sort_exp[1,"clr"]-clr_sort_exp[2, "clr"])/((clr_sort_exp[1,"iqr"]+clr_sort_exp[2,"iqr"])/2)
          return(max_second_dif)
        })
        define_duplicate_hto <- names(duplicate_hto_dif)[which.max(duplicate_hto_dif)]
        update_cluster_hto[define_duplicate_hto] <- as.numeric(i)
        # Other duplicate HTOs: belong to the non-assigned cluster with largest HTO expression
        other_duplicate_htos <-  duplicate_htos[!duplicate_htos == define_duplicate_hto]
        for(j in other_duplicate_htos){
          defined_clusters <- as.numeric(unique(c(update_cluster_hto[!is.na(update_cluster_hto)], duplicate_clusters)))
          clr_sort_df <- duplicate_hto_medians[[j]][order(duplicate_hto_medians[[j]]$clr, decreasing = TRUE),]
          non_defined_cluster <- clr_sort_df$cluster[which(!clr_sort_df$cluster %in% defined_clusters)][1]
          update_cluster_hto[j] <- non_defined_cluster
        }
      }
      cluster_hto <- update_cluster_hto
    }
  }
  return(cluster_hto)
}


#' Assign outlier droplets
#' 
#' Classify cells as singlets or outliers (negatives and doublets) based on Mahalanobis 
#' distance. 
#' 
#' @param md.mat Matrix of Mahalanobis distance from the output of the 
#'        \link{CalculateMD} function.
#' @param md_cut_q Numeric; cut-off to determine outlier cells (negatives and doublets)
#'        based on the quantile of Mahalanobis distance. This parameter is the same 
#'        as in the \link{MDDensPlot} function and is indicated by the red line 
#'        in its plot. Default is 0.95.
#' @param md_cut Numeric; cut-off to determine outlier cells (negatives and doublets)
#'        based on the Mahalanobis distance directly. You only need to specify one of
#'        \code{md_cut_q}, \code{md_cut}, or \code{md_auto} to determine the cut-offs.
#'        This parameter is the same as in the \link{MDDensPlot} function and is 
#'        indicated by the red line in its plot. Default is NULL.
#' @param md_auto Logical; whether to use an automatic method to determine cut-offs
#'        of outlier cells (negatives and doublets). If TRUE, the 0.975 quantile of 
#'        a chi-square distribution will be used as the cut-off. You only need to
#'        specify one of \code{md_cut_q}, \code{md_cut}, or \code{md_auto}. Default 
#'        is FALSE.
#'
#' @details
#' In automatic mode, it is assumed that the Mahalanobis distance follows a 
#' chi-square distribution with \eqn{N-1} degrees of freedom, where \eqn{N} is 
#' the number of hashtags. The cut-off is set at the 0.975 quantile of this distribution.
#' 
#' Cells with Mahalanobis distances larger than the cut-off, as determined by one of 
#' \code{md_cut_q}, \code{md_cut}, or \code{md_auto}, are classified as outlier cells 
#' (negatives and doublets). Others are singlets.
#' 
#' @return A character vector, with one element per cell, indicating whether the 
#' cell is a "Singlet" or an "Outlier" (negatives or doublets).
#' 
#' @examples 
#' \donttest{
#' hto.outlier.assign <- AssignOutlierDrop(hto.md.mat)
#' # Set the cut-off for outlier cells using the quantile of the Mahalanobis distance
#' hto.outlier.assign <- AssignOutlierDrop(hto.md.mat, md_cut_q = 0.52)
#' # Set the cut-off for outlier cells using the value of Mahalanobis distance directly
#' hto.outlier.assign <- AssignOutlierDrop(hto.md.mat, md_cut = 5)
#' # Set the cut-off for outlier cells using an automatic method
#' hto.outlier.assign <- AssignOutlierDrop(hto.md.mat, md_auto = TRUE)
#' head(hto.outlier.assign)
#' table(hto.outlier.assign)
#' }  
#' 
#' @export 
AssignOutlierDrop <- function(md.mat, md_cut_q = 0.95, md_cut = NULL, md_auto = FALSE){
  min.md <- apply(md.mat, 2, FUN = function(x){min(x)})
  if(md_auto){
    outlier.cut <- sqrt(qchisq(p = 0.975, df = (nrow(md.mat)-1)))
  }else{
    tryCatch(
      {
        if(!is.null(md_cut)){
          if(md_cut <= max(min.md)){
            outlier.cut <- md_cut
          }else{
            outlier.cut <- quantile(min.md, probs = md_cut_q) 
          }
        }else{
          outlier.cut <- quantile(min.md, probs = md_cut_q) 
        }
      }, error = function(e){
        print("Wrong input parameters for this function.")
        outlier.cut <- quantile(min.md, probs = 0.95) 
      }
    )
  }
  outlier.barcodes <- names(min.md)[which(min.md > outlier.cut)]
  global.assign <- rep("Singlet", ncol(md.mat))
  names(global.assign) <- colnames(md.mat)
  global.assign[outlier.barcodes] <- "Outlier"
  return(global.assign)
}


#' @importFrom dplyr group_by summarise arrange
#' @importFrom magrittr %>%
.ExamineUnlabelHash <- function(clr.norm, kmed.cl, unlabel_cut = 0.5){
  cl.clr.fc <- apply(clr.norm, 1, FUN = function(x){
    cluster.df <- data.frame("exp" = x, "cluster" = kmed.cl$clustering) 
    cluster.df2 <- cluster.df %>% 
      group_by(cluster) %>%
      summarise("m_exp" = median(exp)) %>%
      arrange(desc(m_exp))
    top2.clusters <- cluster.df2$cluster[c(1,2)]
    cv.top2 <- cluster.df2$m_exp[1]-cluster.df2$m_exp[2]
    return(cv.top2)
  })
  unlabel_hash <- names(cl.clr.fc)[which(cl.clr.fc < unlabel_cut)]
  if(length(unlabel_hash) == 0){
    return(NULL)
  }else{
    return(unlabel_hash)
  }
}


#' Demultiplex cells
#' 
#' Determine the sample identity of each cell based on hashtag information. This
#' function further classfies outlier cells into negatives and doublets.
#' 
#' @param md.mat Matrix of Mahalanobis distance from the output of the 
#'        \link{CalculateMD} function. 
#' @param hash.count The raw hashing count matrix has rows representing hashtags 
#'        and columns representing cells.
#' @param drop.outlier.assign A classification vector of outlier cells and singlets, 
#'        from the output of the \link{AssignOutlierDrop} function.  
#' @param use_gex_data Logical; whether to use gene expression data. It is recommended
#'        to use gene expression data if available. Default is FALSE. 
#' @param gex.count A gene expression count matrix with rows representing genes 
#'        and columns representing cells. Only required when \code{use_gex_data} 
#'        is TRUE. Default is NULL. 
#' @param num_modes Integer; the number of modes (including antimodes) in the 
#'        distribution of hashtag library size, used to determine the cut-off 
#'        between negatives and doublets. This parameter can be visualized by the 
#'        \link{OutlierHTOPlot} function. Default is 3.
#' @param cut_no Integer; the cut-off set at the \code{cut_no}th mode (including 
#'        antimodes) of the hashtag library size distribution. This parameter can 
#'        be visualized by the \link{OutlierHTOPlot} function. This number corresponds 
#'        to the number shown on the red line in the plot. Default is 2.
#' @param lib_cut Double; the cut-off value set directly on hashtag library sizes. 
#'        This parameter can be visualized by the \link{OutlierHTOPlot} function. 
#'        It is an alternative to \code{cut_no}, and should only be specified when 
#'        \code{num_modes} and \code{cut_no} do not provide an ideal cut-off. 
#'        Default is NULL.
#' @param optional Logical; indicating whether optional (extra) clusters have been 
#'        created. If TRUE, reclassification of cells in the extra cluster(s) will 
#'        be performed. This parameter should match the \code{optional} parameter 
#'        in the \link{KmedCluster} function. This parameter only works when 
#'        \code{use_gex_data} is TRUE and \code{gex.count} is provided. Default 
#'        is FALSE.
#' @param kmed.cl K-Medoids clustering results from the output of the 
#'        \link{KmedCluster} function. This parameter only needs to be provided
#'        when \code{optional} is TRUE. Default is NULL. 
#' @param cluster.assign An integer vector indicating the corresponding hashtag 
#'        for each cluster, with names of hashtags. This parameter comes from the
#'        output of the \link{LabelClusterHTO} function. It only needs to be 
#'        provided when \code{optional} is TRUE. Default is NULL. 
#' @param clr.norm Local CLR-normalized data from the output of the 
#'        \link{LocalCLRNorm} function. This parameter only needs to be provided 
#'        when unlabelled hashtags exist in the data. Default is NULL.
#' @param unlabel_cl_cut Double; cut-off to determine unlabelled hashtag(s). 
#'        Default is 0.5.
#' @param unlabel_raw_cut Double; cut-off of log raw hashtag counts to distinguish
#'        negatives from singlets on unlabelled hashtag(s). Default is 4.6.
#' @param unlabel_clr_cut Double; cut-off of local CLR values on unlabelled 
#'        hashtag(s) to distinguish negatives from singlets. Default is -0.5.
#'
#' @details
#' Theoretically, the log hashtag library sizes for outlier cells should be bimodally
#' distributed, with cells having smaller library sizes corresponding to negatives,
#' and those with larger library sizes corresponding to doublets. The 
#' \link{OutlierHTOPlot} function provides a visualization of the modes and antimodes 
#' of the distribution of hashtag library sizes for outlier cells. The parameter 
#' \code{num_modes} specifies the number of modes (including antimodes) applied to 
#' the distribution, and \code{cut_no} specifies the cut-off at the \code{cut_no}th 
#' mode (including antimodes) to separate negatives from doublets. Users can explore 
#' different settings by specifying varying numbers of modes in the 
#' \link{OutlierHTOPlot} function. The value of \code{cut_no} corresponds to the 
#' index of the red line(s) in the plot. Negatives are cells with hashtag library 
#' sizes smaller than the cut-off, while those with larger sizes are classified as 
#' doublets. 
#'
#' It is recommended that users provide gene expression data for demultiplexing 
#' when available, since hashtag information alone can be limited. Gene expression 
#' data can be supplied using the \code{use_gex_data} and \code{gex.count} parameters. 
#' When gene expression data is provided, the median mRNA library size across singlets 
#' is used as a threshold. Initial negatives with mRNA library sizes larger than this 
#' threshold are reclassified as singlets, while initial doublets with mRNA library 
#' sizes smaller than the threshold are reclassified as singlets. 
#' 
#' If optional clustering has been performed, cells in extra clusters that do not 
#' belong to any hashtag will be reclassified. Two thresholds are used for doublets: 
#' (1) the 0.9 quantile of mRNA library sizes across all cells, and (2) the median 
#' mRNA library size across all initial doublets outside the extra clusters. Cells 
#' in the extra cluster with mRNA library sizes greater than either threshold are 
#' reclassified as doublets. The threshold for negatives is the median mRNA library 
#' size across all initial singlets outside the extra clusters. Cells in the extra 
#' cluster with mRNA library sizes smaller than this threshold are reclassified as 
#' negatives. All other cells in the extra cluster are classified as singlets. 
#' 
#' In the case of unlabelled hashtag(s), the first step is to examine whether such 
#' hashtag(s) exist. If any hashtag shows a local CLR difference between the top 
#' two clusters (highest and second highest local CLR values) smaller than the 
#' threshold \code{unlabel_cl_cut}, it is defined as an unlabelled hashtag. The 
#' next step is to rescue singlets from unlabelled hashtag(s). Initial negatives 
#' and initial singlets belonging to unlabelled hashtags are pooled together as 
#' temporary negatives. For each cell \eqn{c} in this group, the difference \eqn{D} 
#' between the local CLR value \eqn{h} of the unlabelled hashtag \eqn{u} and the 
#' maximum local CLR value among all hashtags is calculated as:
#' \eqn{D_c = h_{u,c} - \max_{i=1,...,N}(h_{i,c})}. 
#' The log expression value of the unlabelled hashtag, \eqn{\log(HTO_{u,c})}, is 
#' also recorded. Cells in the temporary negatives with \eqn{D} larger than 
#' \code{unlabel_clr_cut} and \eqn{\log(HTO_{u,c})} larger than \code{unlabel_raw_cut} 
#' are reclassified as singlets belonging to the unlabelled hashtag.
#'
#' @return A data frame of demultiplexing results, where rows represent cells 
#' and row names correspond to cell barcodes. The columns are:
#' \describe{
#'   \item{demux_id}{The hashtag assigned to the cell, defined as the one with the 
#'   minimum Mahalanobis distance.}
#'   \item{global_class}{Classification of each cell as a singlet, doublet, or negative.}
#'   \item{demux_global_class}{Demultiplexing result at the global level. Singlets 
#'   are named by their corresponding hashtag.}
#'   \item{doublet_class}{Demultiplexing result for doublets. Doublets are named by 
#'   the two hashtags with the minimum and second minimum Mahalanobis distances 
#'   for the cell.}
#' }
#' 
#' @references
#' Ameijeiras-Alonso, J., Crujeiras, R. M., & Rodriguez-Casal, A. (2021). 
#' Multimode: an R package for mode assessment. Journal of Statistical Software, 97, 1-32.
#' 
#' @examples 
#' \donttest{
#' # Gene expression data not available:
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign)
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign,
#'   num_modes = 20, cut_no = 26)
#' # Gene expression data is available:
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign, 
#'   num_modes = 20, cut_no = 26, use_gex_data = TRUE, gex.count = hto.gex.count)
#' # Use HTO library size value directly to determine negatives and doublets:
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign,
#'   lib_cut= 5.6, use_gex_data = TRUE, gex.count = hto.gex.count)
#' # Data with optional clustering:
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign,
#'   num_modes = 7, cut_no = 4, use_gex_data = TRUE, gex.count = hto.gex.count, 
#'   optional = TRUE, kmed.cl = hto.kmed.cl, cluster.assign = hto.cluster.assign)
#' # Data with unlabelled hashtag(s):
#' hto.cmddemux.assign <- CMDdemuxClass(hto.md.mat, hto.mtx, hto.outlier.assign, 
#'   num_modes = 18, cut_no = 16, use_gex_data = TRUE, gex.count = hto.gex.count, 
#'   optional = TRUE, kmed.cl = hto.kmed.cl, cluster.assign = hto.cluster.assign,
#'   clr.norm = hto.clr.mtx, unlabel_cl_cut = 0.5, unlabel_raw_cut = 4.6, 
#'   unlabel_clr_cut = -0.5)
#' head(hto.cmddemux.assign)
#' } 
#' 
#' @importFrom multimode locmodes
#' @importFrom Matrix colSums
#' 
#' @export 
CMDdemuxClass <- function(
    md.mat, 
    hash.count, 
    drop.outlier.assign,
    use_gex_data = FALSE,
    gex.count = NULL, 
    num_modes = 3, 
    cut_no = 2, 
    lib_cut = NULL, 
    optional = FALSE, 
    kmed.cl = NULL, 
    cluster.assign = NULL, 
    clr.norm = NULL, 
    unlabel_cl_cut = 0.5, 
    unlabel_raw_cut = 4.6, 
    unlabel_clr_cut = -0.5
){
  demux.id <- apply(md.mat, 2, FUN = function(x){rownames(md.mat)[which.min(x)]})
  hto.lib <- log(colSums(hash.count)+1)
  outlier.lib <- hto.lib[names(drop.outlier.assign)[which(drop.outlier.assign %in% "Outlier")]]
  tryCatch(
    {
      if(is.null(lib_cut)){
        suppressWarnings(cut_off <- locmodes(outlier.lib,num_modes,display=FALSE)$locations[cut_no])
      }else{
        cut_off <- lib_cut[1]
      }
    }, error = function(e){
      print("Wrong input parameters for this function.")
      suppressWarnings(cut_off <- locmodes(outlier.lib,3,display=FALSE)$locations[4])
    }
  )
  negative.barcodes <- names(outlier.lib)[which(outlier.lib <= cut_off)]
  doublet.barcodes <- names(outlier.lib)[which(outlier.lib > cut_off)]
  global.class <- drop.outlier.assign
  global.class[negative.barcodes] <- "Negative"
  global.class[doublet.barcodes] <- "Doublet"
  demux.global.class <- demux.id
  demux.global.class[negative.barcodes] <- "Negative"
  demux.global.class[doublet.barcodes] <- "Doublet"
  if(!is.null(gex.count) & isTRUE(use_gex_data)){
    gex.count <- gex.count[,colnames(hash.count)[which(colnames(hash.count) %in% colnames(gex.count))]]
    gex.lib <- log(colSums(gex.count)+1)
    singlet.barcodes <- names(global.class)[which(global.class %in% "Singlet")]
    singlet.median.lib <- median(gex.lib[singlet.barcodes], na.rm = TRUE)
    global.class[intersect(negative.barcodes, names(gex.lib))] <- ifelse(gex.lib[intersect(negative.barcodes, names(gex.lib))] < singlet.median.lib, "Negative", "Singlet")
    change.barcodes <- names(global.class[negative.barcodes])[which(global.class[negative.barcodes] %in% "Singlet")]
    demux.global.class[change.barcodes] <- demux.id[change.barcodes]
    global.class[intersect(doublet.barcodes, names(gex.lib))] <- ifelse(gex.lib[intersect(doublet.barcodes, names(gex.lib))] < singlet.median.lib, "Singlet", "Doublet")
    change.barcodes2 <- names(global.class[doublet.barcodes])[which(global.class[doublet.barcodes] %in% "Singlet")]
    demux.global.class[change.barcodes2] <- demux.id[change.barcodes2]
  }
  
  # Reclassify droplets in outlier cluster
  if(optional & !is.null(gex.count) & use_gex_data){
    all.clusters <- unique(kmed.cl$clustering)
    outlier.clusters <- all.clusters[which(!all.clusters %in% cluster.assign)]
    outlier.drops <- intersect(names(kmed.cl$clustering)[which(kmed.cl$clustering %in% outlier.clusters)], names(drop.outlier.assign)[which(drop.outlier.assign %in% "Outlier")])
    outlier.db.cut1 <- quantile(gex.lib, 0.9)
    outlier.db.barcodes1 <- outlier.drops[which(gex.lib[outlier.drops] > outlier.db.cut1)]
    nonoutlier.db.bacodes <- names(demux.global.class[which(demux.global.class %in% "Doublet" & !names(demux.global.class) %in% outlier.drops)])
    outlier.db.cut2 <- median(gex.lib[nonoutlier.db.bacodes], na.rm = TRUE)
    outlier.db.barcodes2 <- outlier.drops[which(gex.lib[outlier.drops] > outlier.db.cut2)]
    outlier.db.barcodes <- unique(c(outlier.db.barcodes1, outlier.db.barcodes2))
    demux.global.class[outlier.db.barcodes] <- "Doublet"
    global.class[outlier.db.barcodes] <- "Doublet"
    outlier.nondb.barcodes <- outlier.drops[!outlier.drops %in% outlier.db.barcodes]
    nonoutlier.singlet.barcodes <- singlet.barcodes[!singlet.barcodes %in% outlier.drops]
    nonoutlier.singlet.lib <- lapply(unique(demux.id[nonoutlier.singlet.barcodes]), FUN = function(x){median(gex.lib[nonoutlier.singlet.barcodes[which(demux.id[nonoutlier.singlet.barcodes] %in% x)]])})
    names(nonoutlier.singlet.lib) <- unique(demux.id[nonoutlier.singlet.barcodes])
    demux.global.class[outlier.nondb.barcodes] <- ifelse(gex.lib[outlier.nondb.barcodes] > nonoutlier.singlet.lib[demux.id[outlier.nondb.barcodes]], demux.id[outlier.nondb.barcodes], "Negative")
    demux.global.class[outlier.nondb.barcodes][which(is.na(demux.global.class[outlier.nondb.barcodes]))] <- "Negative"
    global.class[outlier.nondb.barcodes] <- ifelse(gex.lib[outlier.nondb.barcodes] > nonoutlier.singlet.lib[demux.id[outlier.nondb.barcodes]], "Singlet", "Negative")
    global.class[outlier.nondb.barcodes][which(is.na(global.class[outlier.nondb.barcodes]))] <- "Negative"
  }
  
  # Find unlabelled hash and rescue droplets with unlabelled hash
  if(!is.null(kmed.cl) & !is.null(clr.norm)){
    unlabelled.hash <- .ExamineUnlabelHash(clr.norm, kmed.cl, unlabel_cl_cut)
    if(!is.null(unlabelled.hash)){
      unlabelled.barcodes <- names(demux.global.class)[which(demux.global.class %in% unlabelled.hash)]
      global.class[unlabelled.barcodes] <- "Negative"
      demux.global.class[unlabelled.barcodes] <- "Negative"
      tmp.negative.barcodes <- names(global.class)[which(global.class %in% "Negative")]
      for(i in 1:length(unlabelled.hash)){
        unlabelled.hash.idx <- which(rownames(clr.norm) == unlabelled.hash[i])
        unlabelled.fc <- apply(clr.norm[,tmp.negative.barcodes], 2, FUN = function(x){
          unlabelled.clr <- x[unlabelled.hash.idx]
          top <- max(x)
          return(unlabelled.clr-top)
        })
        unlabelled.raw.hash <- log(hash.count[unlabelled.hash[i],tmp.negative.barcodes]+1)
        unlabelled.rescue.barcodes <- names(unlabelled.fc)[which(unlabelled.fc > unlabel_clr_cut & unlabelled.raw.hash > unlabel_raw_cut)]
        global.class[unlabelled.rescue.barcodes] <- "Singlet"
        demux.global.class[unlabelled.rescue.barcodes] <- unlabelled.hash[i]
      }
    }
  }
  
  doublet.class <- demux.global.class
  doublet.barcodes <- colnames(md.mat)[which(demux.global.class == "Doublet")]
  doublet.second.id <- apply(md.mat[,doublet.barcodes], 2, FUN = function(x){rownames(md.mat)[order(x)[2]]})
  doublet.sort.id <- lapply(1:length(doublet.barcodes), FUN = function(x){
    ifelse(which(rownames(hash.count) %in% demux.id[doublet.barcodes[x]]) <  which(rownames(hash.count) %in% doublet.second.id[x]), paste0(demux.id[doublet.barcodes][x], "_", doublet.second.id[x]), paste0(doublet.second.id[x], "_", demux.id[doublet.barcodes][x]))
  })
  doublet.class[doublet.barcodes] <- unlist(doublet.sort.id)
  class.df <- data.frame("demux_id" = demux.id, "global_class" = global.class, "demux_global_class" = demux.global.class, "doublet_class" = doublet.class)
  return(class.df)
}







