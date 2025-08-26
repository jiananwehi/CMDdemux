#' Density plot of Euclidean distances
#'
#' Draw a density plot of the Euclidean distance distribution for cells in each 
#' cluster.
#' 
#' @param cl.dist The Euclidean distance for each cell to its cluster centroid 
#'        from the output of \link{EuclideanClusterDist} function.
#' @param eu_cut_q A vector of cut-offs for defining non-core cells using the 
#'        quantile of Euclidean distances within each cluster. Must be the same 
#'        length as the number of clusters. The red line in the plot represents 
#'        the defined cut-off. Default: 0.9.
#' @param eu_cut A vector of cut-offs for defining non-core cells based on the 
#'        within-cluster Euclidean distances. Must be the same length as the number 
#'        of clusters. This only needs to be specified when the \code{eu_cut_q} 
#'        parameter cannot provide suitable cut-offs. The red line in the plot 
#'        represents the defined cut-off. Default: NULL.
#'        
#' @return A combined layout of density plots, where each plot shows the distribution 
#' of Euclidean distances for a cluster.
#' 
#' @examples 
#' \donttest{
#' EuclideanDistPlot(hto.cl.dist, c(0.9, 0.9, 0.75, 0.92, 0.95, 0.94, 0.95, 0.92))
#' }
#' 
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' 
#' @export
EuclideanDistPlot <- function(cl.dist, eu_cut_q = 0.9, eu_cut = NULL){
  plot.list <- list()
  flag_q <- TRUE
  if(!is.null(eu_cut) & length(eu_cut) >= length(cl.dist)){
    flag_q <- FALSE
  }
  for(i in 1:length(cl.dist)){
    plot.df <- data.frame("x" = cl.dist[[i]])
    tryCatch(
      {
        if(!flag_q){
          if(eu_cut[i] <= max(cl.dist[[i]])){
            cut_off <- eu_cut[i]
          }else{
            cut_off <- quantile(cl.dist[[i]], probs = eu_cut_q[1]) 
          }
        }else{
          if(length(eu_cut_q) < length(cl.dist)){
            cut_off <- quantile(cl.dist[[i]], probs = eu_cut_q[1]) 
          }else{
            cut_off <- quantile(cl.dist[[i]], probs = eu_cut_q[i]) 
          }
        }
      },error = function(e){
        print("Wrong input parameters for this function.")
        cut_off <<- quantile(cl.dist[[i]], probs = 0.9) 
      }
    )
    p <- ggplot(plot.df, aes(x)) + 
      geom_density(color = "#0073C2FF", fill = "#0073C2FF", alpha = 0.4) + 
      geom_vline(xintercept = cut_off, color = "red", linewidth = 1) + 
      labs(title = paste0("cluster", i), x = "Euclidean distance", y = "Density") + 
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face='bold')) +
      ggtitle(paste0("Cluster", i))
    plot.list[[i]] <- p 
  }
  plots <- wrap_plots(plot.list, ncol = 3)
  return(plots)
}


#' @import ggplot2
.CLRClusterPlot <- function(clr.norm, kmed.cl, drop.assign, i, j){
  plot_color_palette <- c("#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", "#FF9896", "#9467BD", "#8C564B", "#E377C2", "#BCBD22", "#AEC7E8", "#17BECF", "#FFBB78", "#98DF8A", "#C5B0D5","#C49C94",  "#F7B6D2",  "#DBDB8D",  "#9EDAE5")
  title.x <- paste("h[", i, "]", sep = "") 
  title.y <- paste("h[", j, "]", sep = "")
  plot.df <- data.frame("x" = clr.norm[i,], "y" = clr.norm[j,], "color" = drop.assign)
  if("non-core" %in% drop.assign){
    all.levels <- c(paste0("Cluster", sort(unique(kmed.cl$clustering))), "non-core")
    all.colors <- c(plot_color_palette[1:length(unique(kmed.cl$clustering))], "grey")
    names(all.colors) <- all.levels
    plot.df$color <- factor(plot.df$color, levels = all.levels[all.levels %in% unique(drop.assign)])
    plot.colors <- all.colors[levels(plot.df$color)]
  }else{
    all.levels <- c(paste0("Cluster", sort(unique(kmed.cl$clustering))))
    all.colors <- plot_color_palette[1:length(unique(kmed.cl$clustering))]
    names(all.colors) <- all.levels
    plot.df$color <- factor(plot.df$color, levels = all.levels[all.levels %in% unique(drop.assign)])
    plot.colors <- all.colors[levels(plot.df$color)]
  }
  sub.df <- data.frame("x" = clr.norm[i,rownames(kmed.cl$medoids)], "y" = clr.norm[j,rownames(kmed.cl$medoids)])
  p <- ggplot(plot.df, aes(x, y, color = color)) + 
    geom_point() +
    scale_color_manual(values = plot.colors) +
    geom_point(sub.df, mapping = aes(x, y), color = "black", shape = 4, stroke = 2) +
    theme_classic() +
    theme(legend.title = element_blank()) + 
    labs(x = parse(text = title.x), y = parse(text = title.y))
  return(p)
}


#' Scatter plot of pairwise local CLR values
#'
#' Scatter plots of pairwise local CLR values, where core cells are colored by 
#' cluster, non-core cells are shown in grey, and cluster medoids are indicated.
#' 
#' @param clr.norm The local CLR normalized data comes from the output of the 
#'        \link{LocalCLRNorm} function.
#' @param kmed.cl K-Medoids clustering results from the output of \link{KmedCluster} 
#'        function.
#' @param core.drop.assign Classification of core and non-core cells from the 
#'        output of the \link{DefineNonCore} function.
#' @param return_cluster_only Logical. If TRUE, one plot is randomly selected per 
#'        cluster. If FALSE, plots of all pairwise local CLR values (h) are returned. 
#'        Default: TRUE.
#' @param seed Set a random seed. Default: 10000. 
#' 
#' @return A combined layout of scatter plots of pairwise local CLR values (h), 
#' where core cells are colored by cluster and non-core cells are shown in grey. 
#' Cluster medoids are indicated by cross markers.
#' 
#' @examples 
#' \donttest{
#' CLRPlot(hto.clr.mtx, hto.kmed.cl, hto.noncore)
#' }
#' 
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' 
#' @export
CLRPlot <- function(clr.norm, kmed.cl, core.drop.assign, return_cluster_only = TRUE, seed = 10000){
  plot.list <- list()
  n <- 1
  if(return_cluster_only){
    cluster_hto <- apply(kmed.cl$medoids, 2, which.max)
    i_j_pair <- c()
    for(i in match(names(sort(cluster_hto)), rownames(clr.norm))){
      set.seed(seed)
      j <- sample(c(1:nrow(clr.norm))[-i], 1)
      plot.list[[n]] <- .CLRClusterPlot(clr.norm, kmed.cl, core.drop.assign, i, j)
      n <- n + 1
    }
  }else{
    for(i in 1:(nrow(clr.norm)-1)){
      for(j in (i+1):nrow(clr.norm)){
        plot.list[[n]] <- .CLRClusterPlot(clr.norm, kmed.cl, core.drop.assign, i, j)
        n <- n + 1
      }
    }
  }
  plots <- wrap_plots(plot.list, ncol = 3)
  return(plots)
}

#' Density plot of Mahalanobis distances
#'
#' Draw a density plot of the minimum Mahalanobis distance distribution across 
#' all cells.
#' 
#' @param md.mat Matrix of Mahalanobis distance from the output of the 
#'        \link{CalculateMD} function.
#' @param md_cut_q Numeric; cut-off to determine outlier cells (negatives and doublets)
#'        based on the quantile of Mahalanobis distance. This parameter is indicated 
#'        by the red line in the plot. Default is 0.95.
#' @param md_cut Numeric; cut-off to determine outlier cells (negatives and doublets)
#'        based on the Mahalanobis distance directly. You only need to specify one of
#'        \code{md_cut_q}, \code{md_cut}, or \code{md_auto} to determine the cut-offs.
#'        This parameter is indicated by the red line in the plot. Default is NULL.
#' @param md_auto Logical; whether to use an automatic method to determine cut-offs
#'        of outlier cells (negatives and doublets). If TRUE, the 0.975 quantile of 
#'        a chi-square distribution will be used as the cut-off. You only need to
#'        specify one of \code{md_cut_q}, \code{md_cut}, or \code{md_auto}. Default 
#'        is FALSE.
#' 
#' @return A density plot of Mahalanobis distances. 
#'       
#' @examples 
#' \donttest{
#' # Visualize the cut-off for outlier cells based on the quantile of the Mahalanobis distance
#' MDDensPlot(hto.md.mat, md_cut_q = 0.94)
#' # Visualize the cut-off for outlier cells based on the value of Mahalanobis distance directly
#' MDDensPlot(hto.md.mat, md_cut = 5)
#' # Visualize the cut-off for outlier cells based on an automatic method
#' MDDensPlot(hto.md.mat, md_auto = TRUE)
#' }
#' 
#' @import ggplot2
#' 
#' @export
MDDensPlot <- function(md.mat, md_cut_q = 0.95, md_cut = NULL, md_auto = FALSE){
  min.md <- apply(md.mat, 2, FUN = function(x){min(x)})
  if(md_auto){
    cut_off <- sqrt(qchisq(p = 0.975, df = (nrow(md.mat)-1)))
  }else{
    tryCatch(
      {
        if(!is.null(md_cut)){
          if(md_cut <= max(min.md)){
            cut_off <- md_cut
          }else{
            cut_off <- quantile(min.md, probs = md_cut_q) 
          }
        }else{
          cut_off <- quantile(min.md, probs = md_cut_q) 
        }
      }, error = function(e){
        print("Wrong input parameters for this function.")
        cut_off <<- quantile(min.md, probs = 0.95) 
      }
    )
  }
  plot.df <- data.frame("x" = min.md)
  p <- ggplot(plot.df, aes(x)) + 
    geom_density(color = "#0073C2FF", fill = "#0073C2FF", alpha = 0.4) + 
    geom_vline(xintercept = cut_off, color = "red", linewidth = 1) +
    labs(x = "Mahalanobis distance", y = "Density") + 
    theme_classic() 
  return(p)
}


#' Density plot of hashtag library sizes of outlier cells
#'
#' Draw a density plot of hashtag library sizes for outlier cells (negatives and 
#' doublets). 
#' 
#' @param hash.count The raw hashing count matrix has rows representing hashtags 
#'        and columns representing cells.
#' @param drop.outlier.assign A classification vector of outlier cells and singlets, 
#'        from the output of the \link{AssignOutlierDrop} function.  
#' @param num_modes Integer; the number of modes (including antimodes) in the 
#'        distribution of hashtag library size, used to determine the cut-off 
#'        between negatives and doublets. The modes and antimodes are displayed 
#'        with red lines in the plot. Default is 3.
#' @param lib_cut Double; the cut-off value applied directly to hashtag library 
#'        sizes. This parameter is indicated by a red line in the plot. It is an 
#'        alternative to \code{num_modes} and should only be specified when 
#'        \code{num_modes} does not provide an ideal cut-off. The default is NULL.
#'        
#' @return A density plot of hashtag library sizes for outlier cells (negatives 
#' and doublets).
#' 
#' @examples 
#' \donttest{
#' OutlierHTOPlot(hto.mtx, hto.outlier.assign)
#' OutlierHTOPlot(hto.mtx, hto.outlier.assign, num_modes = 10)
#' # Specify the cut-off directly using hashtag library size
#' OutlierHTOPlot(hto.mtx, hto.outlier.assign, lib_cut= 5.6)
#' }
#' 
#' @import ggplot2
#' @importFrom Matrix colSums
#' 
#' @export
OutlierHTOPlot <- function(hash.count, drop.outlier.assign, num_modes = 3, lib_cut = NULL){
  hto.lib <- log(colSums(hash.count)+1)
  outlier.lib <- hto.lib[names(drop.outlier.assign)[which(drop.outlier.assign %in% "Outlier")]]
  tryCatch(
    {
      if(is.null(lib_cut)){
        suppressWarnings(cut_off <- locmodes(outlier.lib,num_modes,display=FALSE)$locations)
      }else{
        cut_off <- lib_cut
      }
    }, error = function(e){
      print("Wrong input parameters for this function.")
      suppressWarnings(cut_off <<- locmodes(outlier.lib,3,display=FALSE)$locations)
    }
  )
  plot.df <- data.frame("lib" = outlier.lib)
  p <- ggplot(plot.df, aes(lib)) + 
    geom_density(color = "#0073C2FF", fill = "#0073C2FF", alpha = 0.4) + 
    geom_vline(xintercept = cut_off, color = "red", linewidth = 1) +
    labs(x = "log HTO library", y = "Density") + 
    theme_classic() 
  if(is.null(lib_cut)){
    p <- p + annotate("text", x = cut_off, y = rep(max(density(plot.df$lib)$y), length(cut_off)), label = c(1:length(cut_off)), size = 5)
  }
  return(p)
}


#' Check the number of clusters with a boxplot
#'
#' Use a boxplot to assess whether the chosen number of clusters is appropriate. 
#' Ideally, one box should be significantly higher than the others in each panel.
#' 
#' @param clr.norm The local CLR normalized data comes from the output of the 
#'        \link{LocalCLRNorm} function.
#' @param kmed.cl K-Medoids clustering results from the output of \link{KmedCluster} 
#'        function.
#'
#' @return A combined layout of boxplots, where each panel represents the 
#' expression of a specific hashtag across clusters.
#' 
#' @examples 
#' \donttest{
#' CheckCLRBoxPlot(hto.clr.mtx, hto.kmed.cl)
#' }
#' 
#' @import ggplot2
#' 
#' @export
CheckCLRBoxPlot <- function(clr.norm, kmed.cl){
  plot.df <- NULL
  plot_color_palette <- c("#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", "#FF9896", "#9467BD", "#8C564B", "#E377C2", "#BCBD22", "#AEC7E8", "#17BECF", "#FFBB78", "#98DF8A", "#C5B0D5","#C49C94",  "#F7B6D2",  "#DBDB8D",  "#9EDAE5")
  for(i in 1:nrow(clr.norm)){
    data.df <- data.frame("exp" = clr.norm[i,], "cluster" = kmed.cl$clustering, "hash" = rownames(clr.norm)[i])
    plot.df <- rbind(plot.df, data.df)
  }
  plot.df$hash <- factor(plot.df$hash, levels = rownames(clr.norm))
  plot.df$cluster <- factor(plot.df$cluster)
  p <- ggplot(plot.df, aes(cluster, exp, color = cluster)) + 
    facet_wrap(.~hash, scales = "free_y") +
    geom_boxplot(linewidth = 1.5) + 
    scale_color_manual(values=plot_color_palette) +
    theme_classic() +
    labs(color = "", y = "Local CLR")
  return(p)
}


#' Check suspicious demultiplexed cells 
#'
#' Use pairwise information from hashtag library sizes, mRNA library sizes, 
#' the second Mahalanobis distance (2nd MD), and the ratio between the minimum 
#' Mahalanobis distance (1st MD) and the 2nd MD to assess whether droplets are 
#' correctly classified.
#' 
#' @param hash.count The raw hashing count matrix has rows representing hashtags 
#'        and columns representing cells.
#' @param gex.count A gene expression count matrix with rows representing genes 
#'        and columns representing cells. Only required when gene expression data
#'        is available. Default is NULL.
#' @param md.mat Matrix of Mahalanobis distance from the output of the 
#'        \link{CalculateMD} function.
#' @param cmddemux.assign The demultiplexing results are returned by the 
#'        \link{CMDdemuxClass} function.
#' @param check_barcodes A character vector of cell barcodes to be examined. 
#'        Default: NULL. 
#' @param check_barcodes_col The color used to highlight the specified cells 
#'        in the plot. Default: "blue".
#' @param check_barcode_size The size of the points representing the specified 
#'        cells in the plot. Default: 2.
#'        
#' @return A combined layout of pairwise information from hashtag library sizes, 
#' mRNA library sizes, the 2nd MD, and the ratio between the 1st MD and the 2nd 
#' MD. Cells are color-coded by classification (singlets, negatives, doublets), 
#' with selected barcodes highlighted.
#' 
#' @examples 
#' \donttest{
#' hto.check.barcodes <- colnames(hto.mtx)[1:30]
#' CheckAssign2DPlot(hto.mtx, hto.md.mat, hto.cmddemux.assign, gex.count = hto.gex.count, check_barcodes = hto.check.barcodes)
#' }
#' 
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom Matrix colSums
#' 
#' @export
CheckAssign2DPlot <- function(
    hash.count, 
    gex.count = NULL,
    md.mat, 
    cmddemux.assign, 
    check_barcodes = NULL, 
    check_barcodes_col = "blue", 
    check_barcode_size = 2
){
  common.barcodes <- intersect(colnames(hash.count), colnames(gex.count))
  check_barcodes <- check_barcodes[check_barcodes %in% common.barcodes]
  hash.count <- hash.count[,common.barcodes]
  md.mat <- md.mat[,common.barcodes]
  cmddemux.assign <- cmddemux.assign[common.barcodes,]
  first.md <- apply(md.mat, 2, FUN = function(x){min(x)})
  second.md <- apply(md.mat, 2, FUN = function(x){sort(x)[2]})
  md.ratio <- first.md/second.md
  if(!is.null(gex.count)){
    gex.count <- gex.count[,common.barcodes]
    data.df <- data.frame("hto" = log(colSums(hash.count)+1), "gex" = log(colSums(gex.count)+1), "md2" = second.md, "ratio" = md.ratio, "CMDdemux" = cmddemux.assign$global_class)
    data.col.names <- c("Log HTO library size", "Log mRNA library size", "Second Mahalanobis distance", "Ratio of 1st MD/2nd MD", "CDMdemux")
  }else{
    data.df <- data.frame("hto" = log(colSums(hash.count)+1), "md2" = second.md, "ratio" = md.ratio, "CMDdemux" = cmddemux.assign$global_class)
    data.col.names <- c("Log HTO library size", "Second Mahalanobis distance", "Ratio of 1st MD/2nd MD", "CDMdemux")
  }
  
  plot.list <- list()
  n <- 1
  for(i in 1:(ncol(data.df)-2)){
    for(j in (i+1):(ncol(data.df)-1)){
      if(!(data.col.names[i] == "Second Mahalanobis distance" & data.col.names[j] == "Ratio of 1st MD/2nd MD")){
        plot.df <- data.frame("x" = data.df[,i], "y" = data.df[,j], "CMDdemux" = data.df$CMDdemux)
        rownames(plot.df) <- rownames(data.df)
        plot.df$CMDdemux <- factor(plot.df$CMDdemux, levels = c("Singlet", "Doublet", "Negative"))
        p <- ggplot(plot.df, aes(x, y, color = CMDdemux)) + 
          geom_point(size = 1) + 
          scale_color_manual(values=c("black", "red", "grey")) + 
          theme_classic() +
          labs(x = data.col.names[i], y = data.col.names[j])
        if(!is.null(check_barcodes) & length(check_barcodes) > 0){
          sub.plot.df <- plot.df[check_barcodes,]
          p <- p + geom_point(data = sub.plot.df, mapping = aes(x, y), size = check_barcode_size, shape = 15, color = check_barcodes_col)
        }
        plot.list[[n]] <- p
        n <- n + 1
      }
    }
  }
  combined_plot <- wrap_plots(plot.list) + plot_layout(guides = "collect")
  return(combined_plot)
}


#' Check suspicious demultiplexed cells in a 3D plot
#'
#' Assess the accuracy of droplet classification by visualizing hashtag library 
#' sizes, mRNA library sizes, and the ratio of the minimum to the second 
#' Mahalanobis distance.
#' 
#' @param hash.count The raw hashing count matrix has rows representing hashtags 
#'        and columns representing cells.
#' @param gex.count A gene expression count matrix with rows representing genes 
#'        and columns representing cells. Only required when gene expression data
#'        is available. 
#' @param md.mat Matrix of Mahalanobis distance from the output of the 
#'        \link{CalculateMD} function.
#' @param cmddemux.assign The demultiplexing results are returned by the 
#'        \link{CMDdemuxClass} function.
#' @param check_barcodes A character vector of cell barcodes to be examined. 
#'        Default: NULL. 
#' @param check_barcodes_col The color used to highlight the specified cells 
#'        in the plot. Default: "blue".
#' @param point_size The size of points in the plot. Default: 2. 
#' @param check_barcode_size The size of the points representing the specified 
#'        cells in the plot. Default: 5.
#'        
#' @return A 3D plot with library sizes, mRNA library sizes, and the ratio of the minimum to the second 
#' Mahalanobis distance. Cells are color-coded by classification (singlets, 
#' negatives, doublets), with selected barcodes highlighted.
#' 
#' @examples 
#' \donttest{
#' hto.check.barcodes <- colnames(hto.mtx)[1:30]
#' CheckAssign3DPlot(hto.mtx, hto.gex.count, hto.md.mat, hto.cmddemux.assign, hto.check.barcodes)
#' }
#' 
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom Matrix colSums
#' 
#' @export
CheckAssign3DPlot <- function(
    hash.count, 
    gex.count, 
    md.mat, 
    cmddemux.assign, 
    check_barcodes = NULL, 
    check_barcodes_col = "blue", 
    point_size = 2, 
    check_barcode_size = 5){
  common.barcodes <- intersect(colnames(hash.count), colnames(gex.count))
  check_barcodes <- check_barcodes[check_barcodes %in% common.barcodes]
  hash.count <- hash.count[,common.barcodes]
  gex.count <- gex.count[,common.barcodes]
  md.mat <- md.mat[,common.barcodes]
  cmddemux.assign <- cmddemux.assign[common.barcodes,]
  first.md <- apply(md.mat, 2, FUN = function(x){min(x)})
  second.md <- apply(md.mat, 2, FUN = function(x){sort(x)[2]})
  md.ratio <- first.md/second.md
  plot.df <- data.frame("hto" = log(colSums(hash.count)+1), "gex" = log(colSums(gex.count)+1), "ratio" = md.ratio, "CMDdemux" = cmddemux.assign$global_class)
  plot.df$CMDdemux <- factor(plot.df$CMDdemux, levels = c("Singlet", "Doublet", "Negative"))
  sub.plot.df <- plot.df[check_barcodes,]

  p <- plot_ly(plot.df, x = ~ratio, y = ~hto, z = ~gex, color = ~CMDdemux, colors = c("black", "red", "grey"),  type = 'scatter3d', mode = 'markers', size = point_size)%>% 
    add_trace(data = sub.plot.df, x = ~ratio, y = ~hto, z = ~gex, type = 'scatter3d', mode = 'markers', marker = list(color = check_barcodes_col, size = check_barcode_size, symbol = 'circle'), showlegend = FALSE) %>% 
    layout(scene = list(xaxis = list(title = "Ratio of 1st MD/2nd MD"), yaxis = list(title = 'Log HTO library'), zaxis = list(title = "Log mRNA library")))
  return(p)
}














