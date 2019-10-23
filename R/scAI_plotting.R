
#' ggplot theme in scAI
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
scAI_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}



#' Visualize the inferred biologically relevant factors
#' We plot the heatmap of the three learned low-rank matrices using hierarchical clustering.
#' @param object scAI object
#' @param color.by the name of the variable in object.pData; defining cell groups (not necessary)
#' @param colors.use defined colors of the cell groups
#' @param do.sampling whether perform sampling of loci when generating heatmap of the loci-factor matrix
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats setNames
#' @importFrom grid grid.grabExpr grid.newpage pushViewport grid.draw unit gpar viewport popViewport
lmHeatmap <- function(object, color.by, colors.use = NULL,do.sampling = T ){

  H <- as.matrix(object@fit$H)
  H <- sweep(H,2,colSums(H),FUN = `/`)

  label <- object@pData[[color.by]]
  df<- data.frame(group = label); rownames(df) <- colnames(H)

  if (is.null(colors.use)) {
    colors.use <- brewer.pal(length(unique(label)), "Set1")
  }
  cell.cols.assigned <- setNames(colors.use, unique(as.character(df$group)))
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),annotation_name_side = "left",simple_anno_size = grid::unit(0.2, "cm"))
  colormap = structure(rev(brewer.pal(9,"RdBu")))

  ht1 = Heatmap(H,name = "H",
                clustering_method_columns = "average",
                clustering_distance_columns = "euclidean",
                col = colormap,
                cluster_rows = FALSE, show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left", row_names_rot = 0,row_names_gp = gpar(fontsize = 10),
                width = unit(6, "cm"), height = unit(4, "cm"),
                top_annotation = col_annotation,
                column_title = "Cell loading matrix",
                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                heatmap_legend_param = list(title = "H", at = c(0, 0.5, 1),legend_width = unit(0.0001, "cm"),legend_height = unit(2, "cm"),labels_gp = gpar(font = 6))
  )


  # heatmap for W1
  W1 <- as.matrix(object@fit$W[[1]])
  W1 <- sweep(W1,1,rowSums(W1),FUN = `/`)
  W1[is.na(W1)] <- 0
  colormap = structure(rev(brewer.pal(11,"RdBu")))
  ht2 = Heatmap(W1,name = "W1",
                clustering_method_rows = "average",
                col = colormap,
                cluster_columns = FALSE, show_column_names = T, show_row_names = F, column_names_gp = gpar(fontsize = 10),
                width = unit(4, "cm"), height = unit(8, "cm"),
                column_title = "Gene loading matrix (scRNA-seq)",
                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                row_title = "Genes", row_title_rot = 90,row_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "W1", at = c(0, 0.5, 1),legend_width = unit(0.0001, "cm"),legend_height = unit(2, "cm"),labels_gp = gpar(font = 6))
  )

  # heatmap for W1
  W2 <- as.matrix(object@fit$W[[2]])
  W2 <- sweep(W2,1,rowSums(W2),FUN = `/`)
  W2[is.na(W2)] <- 0
  if (nrow(W2) > 5000 & do.sampling) {
    loci.use <- sample(1:nrow(W2), 5000, replace=F)
    W2 <- W2[sort(loci.use),]
  }

  colormap = structure(rev(brewer.pal(9,"Spectral")))
  ht3 = Heatmap(W2,name = "W2",
                clustering_method_rows = "average",
                col = colormap,
                cluster_columns = FALSE, show_column_names = T, show_row_names = F, column_names_gp = gpar(fontsize = 10),
                width = unit(4, "cm"), height = unit(8, "cm"),
                column_title = "Locus loading matrix (scATAC-seq)",
                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                row_title = "Loci", row_title_rot = 90,row_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "W2", at = c(0, 0.5, 1),legend_width = unit(0.0001, "cm"),legend_height = unit(2, "cm"),labels_gp = gpar(font = 6))
  )
  gb_ht1 = grid.grabExpr(draw(ht1))
  gb_ht2 = grid.grabExpr(draw(ht2))
  gb_ht3 = grid.grabExpr(draw(ht3))
  grid.newpage()
  pushViewport(viewport(x = 0.2,y = 1, width = 0.5, height = 0.3, just = c("left", "top")))
  grid.draw(gb_ht1)
  popViewport()

  pushViewport(viewport(x = 0.1, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
  grid.draw(gb_ht2)
  popViewport()

  pushViewport(viewport(x = 0.5, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
  grid.draw(gb_ht3)
  popViewport()
}



#' visualize cells in 2D-dimensional space
#'
#' @param object scAI object
#' @param cell_coords 2D embedding coordinates of cells
#' @param color.by the name of the variable in pData, defining cell groups, cells are colored based on the labels
#' @param labels.order defining the factor level of cell groups
#' @param colors.use defining the color for each cell group
#' @param brewer.use use RColorBrewer palette instead of default ggplot2 color
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param label.size font size of the legend
#' @param cell.size size of the dots
#' @param font.size font size
#' @param do.label label the cluster in 2D space
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#'
#' @return ggplot2 object with 2D plot
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot geom_point aes scale_color_manual facet_wrap element_text theme guides element_blank element_rect geom_line
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% summarize
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats median
cellVisualization <- function(object, cell_coords, color.by, labels.order = NULL, colors.use = NULL, brewer.use = T,
                              xlabel = "tSNE1", ylabel = "tSNE2", title = NULL,
                              label.size = 4, cell.size = 0.3, font.size = 10, do.label = F, show.legend = T, show.axes = T) {


  labels <- object@pData[[color.by]]

  if (is.null(labels.order) == FALSE) {
    labels <- factor(labels, levels = labels.order)
  } else if (class(labels) != "factor") {
    labels <- as.factor(labels)
  }

  df <- data.frame(x = cell_coords[, 1], y = cell_coords[, 2], group = labels)

  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(colour = labels), size = cell.size) + scAI_theme_opts() +
    theme(text = element_text(size = 10)) + labs(title = title, x = xlabel, y = ylabel) +
    guides(colour = guide_legend(override.aes = list(size = label.size)))
  numCluster = length(unique((labels)))
  if (is.null(colors.use)) {
    if (brewer.use) {
      if (numCluster < 9) {
        colors <- RColorBrewer::brewer.pal(numCluster, "Set1")
      } else {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(numCluster)
      }
      names(colors) <- levels(labels)
      gg <- gg + scale_color_manual(values = colors)
    }
  } else {
    gg <- gg + scale_color_manual(values = colors.use)
  }

  if (do.label) {
    centers <- df %>% dplyr::group_by(group) %>% dplyr::summarize(x = median(x = x), y = median(x = y))
    gg <- gg + ggrepel::geom_text_repel(data = centers, mapping = aes(x, y, label = group), size = label.size)
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}



#' Ranking the features (genes/loci) and show the top markers in each factor
#'
#' @param object scAI object
#' @param assay define an assay to show, e.g., assay = "RNA"
#' @param feature.show a vector of the features that are labeled on the plot
#' @param feature.show.names instead of the default name in feature.show, one can show the manual feature name such as the enriched motif
#' @param top.p showing the features in top ranking
#' @param features.diff a table includes the differential features, returned from identifyfactorMakrers.R
#' @param ylabel ylabel shown on the y-axis
#'
#' @return
#' @export
#'
#' @examples
featureRankingPlot <- function(object, assay, feature.show = NULL, feature.show.names = NULL, top.p = 0.5, features.diff = NULL, ylabel = "Weight") {
  W <- object@fit$W[[assay]]
  features <- rownames(W)
  K = ncol(W)
  W <- sweep(W,1,rowSums(W),FUN = `/`)
  W[is.na(W)] <- 0

  Wg <- vector("list", K)
  for (i in 1:K) {
    W_order <- sort(W[,i],decreasing=F, index.return = T)
    features_ordered <- features[W_order$ix]
    if (!is.null(features.diff)) {
      features.diffi <- as.character(features.diff$features[features.diff$factors == i])
    }else {
      features.diffi <- as.character(features)
    }

    if (!is.null(feature.show)) {
      features.diffi <- intersect(features.diffi, feature.show)
    }
    idx <- match(features.diffi, features_ordered)
    data_show <- matrix(0, nrow(W), 1); data_show[idx] <- 1
    if (!is.null(top.p) & top.p < 1) {
      idx_bottom <- seq_len(floor((1-top.p)*nrow(W))); data_show[idx_bottom] <- 0
    }

    Wg[[i]] <- cbind(Weight =  as.numeric(W_order$x), factor = i, Ranking = seq_len(nrow(W)), Show = as.numeric(data_show), Genes = features_ordered)
  }
  data <- Wg[[1]]
  for (i in 2:K) {
    data <- rbind(data, Wg[[i]])
  }

  df <- as.data.frame(data, stringsAsFactors=FALSE)
  colnames(df) <- c("Weight", "factor", "Ranking", "Show","Genes")
  df$factor <- paste('Factor',df$factor, sep = " ")
  df$Weight <- as.numeric(as.character(df$Weight))
  df$Ranking <- as.numeric(as.character(df$Ranking))
  df$Show <- as.numeric(as.character(df$Show))

  if (!is.null(feature.show.names)) {
    idx <- which(df$Genes %in% feature.show)
    df$Genes[idx] <- feature.show.names
  }

  data_topFeature = df[df$Show == 1,]

  gg <- ggplot(df, aes(Ranking, Weight)) +
    geom_line(colour = "grey80",size = 1) + facet_wrap(~ factor, ncol = 5, scales = "free")+
    scAI_theme_opts()+
    theme(text = element_text(size = 10), axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
    theme(strip.background = element_rect(fill="grey80")) +
    ylab(ylabel) +
    geom_point(size = 3, shape = 1, data = data_topFeature) +
    ggrepel::geom_text_repel(aes(label = Genes), data = data_topFeature, segment.color = "grey50", segment.alpha = 1,
                             direction = "y",nudge_x = -150, hjust = 1,size = 3,segment.size = 0.3) # hjust = 1 for right-align
  gg
}




#' VscAI visualize the genes, loci and factors that separate cell states on two dimensions alongside the cells
#'
#' @param object scAI object
#' @param gene.use embedded genes
#' @param loci.use embedded loci
#' @param loci.use.names alternative names of embedded loci, e.g, the corresponding motif
#' @param color.by the name of the variable in pData, defining cell groups, cells are colored based on the labels
#' @param labels.order defining the factor level
#' @param colors.use defining the color for each cell group
#' @param brewer.use use RColorBrewer palette instead of default ggplot2 color
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param label.size font size of the legend
#' @param cell.size size of the dots
#' @param font.size size of font
#' @param do.label label the cluster in 2D space
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#'
#' @return ggplot2 object with 2D plot
#' @export
#'
#' @examples
#' @importFrom ggplot2 guide_legend guides labs element_text theme xlab ylab scale_fill_manual scale_color_manual scale_shape_manual scale_size_manual

VscAIplot <- function(object, gene.use, loci.use, loci.use.names, color.by,
                      labels.order = NULL, colors.use = NULL, brewer.use = T, xlabel = "VscAI1",
                      ylabel = "VscAI2", title = NULL, label.size = 3, cell.size = 0.3, font.size = 10,
                      do.label = T, show.legend = T, show.axes = T) {

  cell_coords <- object@embed$VscAI$cells
  factor_coords <- object@embed$VscAI$factors
  gene_coords <- object@embed$VscAI$genes
  loci_coords <- object@embed$VscAI$loci

  labels <- object@pData[[color.by]]

  if (is.null(labels.order) == FALSE) {
    labels <- factor(labels, levels = labels.order)
  } else if (class(labels) != "factor") {
    labels <- as.factor(labels)
  }

  df.cell <- data.frame(x = cell_coords[, 1], y = cell_coords[, 2], group = labels)

  gg <- ggplot(data = df.cell, aes(x, y)) +
    geom_point(aes(colour = labels), size = cell.size) +
    scAI_theme_opts() + theme(text = element_text(size = 10)) +
    labs(title = title) + xlab(xlabel) + ylab(ylabel) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    guides(fill = guide_legend(title = "Cell groups")) + scale_fill_manual("Cell groups")

  numCluster = length(unique((labels)))
  if (is.null(colors.use)) {
    if (brewer.use) {
      if (numCluster < 9) {
        colors <- RColorBrewer::brewer.pal(numCluster, "Set1")
      } else {
        colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(numCluster)
      }
      names(colors) <- levels(labels)
      gg <- gg + scale_color_manual(values = colors)
    }
  } else {
    gg <- gg + scale_color_manual(values = colors.use)
  }


  # embedding factors
  if (do.label) {
    df.factor <- data.frame(factor_coords, label.name = paste0("F", seq_len(length(factor_coords[, 1]))), Embedding = "Factors")
    df.features <- df.factor
  }

  # embedding genes
  if (!is.null(gene.use)) {
    df.genes <- data.frame(gene_coords[gene.use, ], label.name = gene.use,
                           Embedding = "Genes")
    df.features <- rbind(df.features, df.genes)
  }

  # embedding loci
  if (!is.null(loci.use)) {
    df.loci <- data.frame(loci_coords[loci.use, ], label.name = loci.use.names,
                          Embedding = "Loci")
    df.features <- rbind(df.features, df.loci)
  }


  gg <- gg + geom_point(data = df.features, aes(x, y, shape = Embedding, size = Embedding)) +
    scale_shape_manual(values = c(1, 16, 5)) +
    scale_size_manual(values = c(3, 2, 2)) +
    ggrepel::geom_text_repel(data = df.features, aes(label = label.name), size = label.size,
                             segment.color = "grey50", segment.size = 0.3, box.padding = grid::unit(0.35, "lines"), point.padding = grid::unit(0.2, "lines"))


  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}



#' visualize cells on the 2D space with gene expression or chromatian accessibility overlayed
#'
#' @param object scAI object
#' @param assay define an assay to show, e.g., assay = "RNA"
#' @param feature.use a vector of features
#' @param method dimensional reduction method, e.g., VscAI, tsne, umap
#' @param nCol number of columns of the plot
#' @param xlabel label shown on x-axis
#' @param ylabel label shown on y-axis
#' @param cell.size the size of points (cells)
#' @param show.legend whether show individual legend
#' @param show.legend.combined  whether just show one legend
#' @param show.axes whether show the axes
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 guide_colorbar scale_colour_gradientn
featureVisualization <- function(object, assay, feature.use, method = "VscAI", nCol = NULL,
                                 xlabel = "VscAI1", ylabel = "VscAI2", cell.size = 0.3,
                                 show.legend = T, show.legend.combined = F, show.axes = T) {

  data <- object@norm.data[[assay]]

  feature.use <- intersect(feature.use, rownames(data))
  data.use <- data[feature.use, ]

  if (is.null(nCol)) {
    if (length(feature.use) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(feature.use), 3)
    }
  }
  if (method == "VscAI") {
    cell_coords <- object@embed$VscAI$cells
  } else if (method == "tsne") {
    cell_coords <- object@embed$tsne
    xlabel = "tSNE1"
    ylabel = "tSNE2"
  } else if (method == "umap") {
    cell_coords <- object@embed$umap
    xlabel = "UMAP1"
    ylabel = "UMAP2"
  }

  colormap <- colorRampPalette(c("#FFFFEF", "#FFFF00", "#FF0000", "#0A0000"))(64)
  colormap[1] <- "#E5E5E5"

  df <- data.frame(x = cell_coords[, 1], y = cell_coords[, 2])
  numFeature = length(feature.use)
  gg <- vector("list", numFeature)
  for (i in seq_len(numFeature)) {
    feature.name <- feature.use[i]
    df$feature.data <- data.use[i, ]
    g <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(colour = feature.data), size = cell.size) +
      scale_colour_gradientn(colours = colormap, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey") +
      labs(title = feature.name) + scAI_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) + labs(x = xlabel, y = ylabel)

    if (!show.legend) {
      g <- g + theme(legend.position = "none")
    }

    if (show.legend.combined & i == numFeature) {
      g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.5, "in"), legend.title = NULL)
    }

    if (!show.axes) {
      g <- g + theme_void()
    }
    gg[[i]] <- g
  }
  gg.combined <- cowplot::plot_grid(plotlist = gg, ncol = nCol)

  gg.combined
}


#' visualize cells on the 2D space with features overlayed
#'
#' @param object scAI object
#' @param feature.use a vector of features
#' @param feature.scores a matrix containing the feature scores
#' @param method dimensional reduction method, e.g., VscAI, tsne, umap
#' @param colormap RColorbrewer palette to use
#' @param color.direction Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
#' @param nCol number of columns of the plot
#' @param xlabel label shown on x-axis
#' @param ylabel label shown on y-axis
#' @param cell.size the size of points (cells)
#' @param show.legend whether show individual legend
#' @param show.legend.combined  whether just show one legend
#' @param show.axes whether show the axes
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 guide_colorbar scale_color_distiller
featureScoreVisualization <- function(object, feature.use = NULL, feature.scores, method = "VscAI",
                                      colormap = "RdPu", color.direction = 1,
                                      nCol = NULL, xlabel = "VscAI1", ylabel = "VscAI2",
                                      show.axes = T,  cell.size = 0.3,
                                      show.legend = T, show.legend.combined = F) {

  data.use <- as.matrix(feature.scores[ ,feature.use])

  if (is.null(nCol)) {
    if (length(feature.use) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(feature.use), 3)
    }
  }

  if (method == "VscAI") {
    cell_coords <- object@embed$VscAI$cells
  } else if (method == "tsne") {
    cell_coords <- object@embed$tsne
    xlabel = "tSNE1"
    ylabel = "tSNE2"
  } else if (method == "umap") {
    cell_coords <- object@embed$umap
    xlabel = "UMAP1"
    ylabel = "UMAP2"
  }

  df <- data.frame(x = cell_coords[, 1], y = cell_coords[, 2])
  numFeature = length(feature.use)
  gg <- vector("list", numFeature)
  for (i in seq_len(numFeature)) {
    feature.name <- feature.use[i]
    df$feature.data <- data.use[ ,i]

    g <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(colour = feature.data), size = cell.size) +
      scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "lightgrey") +
      labs(title = feature.name) + scAI_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in")) + labs(x = xlabel, y = ylabel)

    if (!show.legend) {
      g <- g + theme(legend.position = "none")
    }

    if (show.legend.combined & i == numFeature) {
      g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.5, "in"), legend.title = NULL)
    }

    if (!show.axes) {
      g <- g + theme_void()
    }
    gg[[i]] <- g
  }
  gg.combined <- cowplot::plot_grid(plotlist = gg, ncol = nCol)

  gg.combined
}



#' generate a heatmap for the expression of differential features across different cell groups
#'
#' @param object scAI object
#' @param assay define an assay to show, e.g., assay = "RNA"
#' @param feature.use a vector of features to show
#' @param group.by the name of the variable in pData, defining cell groups. cells are grouped together
#' @param rescaling whether rescale each feature across all the cells
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom circlize colorRamp2
featureHeatmap <- function(object, assay, feature.use,  group.by, rescaling = T) {

  data <- as.matrix(object@norm.data[[assay]])
  groups = object@pData[[group.by]]

  feature.use <- intersect(feature.use, rownames(data))
  data.use <- data[feature.use,]

  if(rescaling) {
    data.use = t(scale(t(data.use), center = T))
  }

  cell.order <- order(groups)
  data.use <- data.use[,cell.order]
  numCluster <- length(unique(groups))

  colorGate = structure(brewer.pal(numCluster, "Set1"), names = as.character(unique(groups)))

  col_annotation = HeatmapAnnotation(group = sort(groups),col = list(group = colorGate),
                                     annotation_name_side = "left",simple_anno_size = unit(0.2, "cm"))
  Heatmap(data.use,name = "zscore",
          col = colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"),space = "LAB"),
          cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
          show_row_names = TRUE, row_names_side = "left", row_names_rot = 0,row_names_gp = gpar(fontsize = 8),
          width = unit(6, "cm"),
          bottom_annotation = col_annotation,
          heatmap_legend_param = list(title = NULL, legend_width = unit(0.0001, "cm"),labels_gp = gpar(font = 6))
  )

}
