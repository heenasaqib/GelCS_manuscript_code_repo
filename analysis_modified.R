#' added two more parameteres to the cellchat netAnalysis_signalingRole_heatmap function
#' @param color.rownames assign color to rownames
#' @param return.data return data is TRUE
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' 
netAnalysis_signalingRole_heatmap_modified <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                              color.use = NULL, color.heatmap = "BuGn", color.rownames = NULL, return.data = FALSE,
                                              title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){

  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }

  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA


  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }

  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)

  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows, cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size, col = color.rownames),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 12),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  if (return.data) {
    return(list(mat, ht1))
  } else {
    return(ht1)
  }
}

#' Heatmap showing the contribution of signals (signaling pathways or ligand-receptor pairs)
#' to cell groups in terms of differential outgoing or incoming signaling between two datasets
#'
#' In this heatmap, colobar represents the relative signaling strength of differential signaling
#' pathway across cell groups (NB: values are row-scaled). Postive bar means stronger pathway
#' signaling in the second dataset
#'
#' The top colored bar plot shows the total differential signaling strength of a cell group by
#' summarizing all signaling pathways displayed in the heatmap. Postive bar means stronger overall
#' cell signaling in the second dataset
#'
#' The right grey bar plot shows the total signaling strength of a signaling pathway by
#' summarizing all cell groups displayed in the heatmap. Postive grey bar means stronger overall
#' signaling pathway signaling in the second dataset.
#'
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_barplot rowAnnotation
#' @importFrom stats setNames
#'
#' @return
#' @export
#'
netAnalysis_diff_signalingRole_heatmap <- function(object, signaling = NULL, signaling.exclude = NULL, comparison = c(1,2), pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                              color.use = NULL, color.heatmap = "BuGn", color.rownames = NULL, return.data = FALSE, 
                                              title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  if (is.list(object)) {
    object <- mergeCellChat(object, add.names = names(object))
  }
  if (!is.list(object@net[[1]])) {
    stop("This function cannot be applied to a single cellchat object from one dataset!")
  }
  if (length(names(object@net)) != 2) {
    stop("This function cannot be applied to cellchat object from more than two datasets")
  }

  pattern <- match.arg(pattern)
  dataset.name <- names(object@net)

  message(paste0("Visualizing differential pathway ", pattern, " signaling changes from ", dataset.name[comparison[1]], " to ", dataset.name[comparison[2]]))

  if (is.null(signaling)) {
    signaling <- union(object@netP[[comparison[1]]]$pathways, object@netP[[comparison[2]]]$pathways)
  }
  if (!is.null(signaling.exclude)) {
    signaling <- setdiff(signaling, signaling.exclude)
  }

  cell.levels <- levels(object@idents$joint)
  mat.all.merged <- list()
  mat.ori.all.merged <- list()
  for (ii in 1:length(comparison)) {
    if (length(slot(object, slot.name)[[comparison[ii]]]$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores for each dataset seperately! ")
    }
    if (sum(c("outdeg", "indeg") %in% names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]])) !=2) {
      stop(paste0("`outdeg, indeg` should be one of ", paste(names(slot(object, slot.name)[[comparison[ii]]]$centr[[1]]),collapse=", "), '\n', "`outdeg_unweighted` is only supported for version >= 1.1.2"))
    }

    centr <- slot(object, slot.name)[[comparison[ii]]]$centr
    outgoing <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    incoming <- matrix(0, nrow = length(cell.levels), ncol = length(centr))
    dimnames(outgoing) <- list(cell.levels, names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[,i] <- centr[[i]]$outdeg
      incoming[,i] <- centr[[i]]$indeg
    }
    mat.out <- t(outgoing)
    mat.in <- t(incoming)

    mat.all <- array(0, dim = c(length(signaling),ncol(mat.out),2))
    mat.ori.all <- array(0, dim = c(length(signaling),ncol(mat.out),2))
    mat.t <-list(mat.out, mat.in)
    for (i in 1:length(comparison)) {
      mat = mat.t[[i]]
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
      dimnames(mat) <- list(signaling, colnames(mat1))

      mat.ori <- mat
      mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
      mat[is.na(mat)] <- 0
      mat.ori[is.na(mat.ori)] <- 0
      mat.all[,,i] = mat
      mat.ori.all[,,i] = mat.ori
    }
    dimnames(mat.all) <- list(dimnames(mat)[[1]], dimnames(mat)[[2]], c("outgoing", "incoming"))
    dimnames(mat.ori.all) <- list(dimnames(mat.ori)[[1]], dimnames(mat.ori)[[2]], c("outgoing", "incoming"))
    mat.all.merged[[ii]] <- mat.all
    mat.ori.all.merged[[ii]] <- mat.ori.all
  }

  mat.diff <- mat.all.merged[[2]] -  mat.all.merged[[1]] # subtracting between the two datasets
  mat.ori.diff <- mat.ori.all.merged[[2]] -  mat.ori.all.merged[[1]]

  outgoing.diff <- mat.diff[ , , 1]
  incoming.diff <- mat.diff[ , , 2]
  all.diff <- outgoing.diff + incoming.diff

  outgoing.ori.diff <- mat.ori.diff[ , , 1]
  incoming.ori.diff <- mat.ori.diff[ , , 2]
  ori.all.diff <- outgoing.ori.diff + incoming.ori.diff

  if (pattern == "outgoing") {
    mat <- outgoing.diff
    mat <- sweep(mat, 1L, apply(abs(mat), 1, max), '/', check.margin = FALSE)
    mat.ori <- outgoing.ori.diff
    mat[mat == 0] <- NA
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- incoming.diff
    mat <- sweep(mat, 1L, apply(abs(mat), 1, max), '/', check.margin = FALSE)
    mat.ori <- incoming.ori.diff
    mat[mat == 0] <- NA
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- all.diff
    mat <- sweep(mat, 1L, apply(abs(mat), 1, max), '/', check.margin = FALSE)
    mat.ori <- ori.all.diff
    mat[mat == 0] <- NA
    legend.name <- "Overall"
  }

  if (is.null(title)) {
    title <- paste0(legend.name, " signaling changes ", " (", dataset.name[comparison[1]], " vs. ", dataset.name[comparison[2]], ")")
  }

  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }

  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)

  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative differential strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows, cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size, col = color.rownames), column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  if (return.data) {
    return(list(mat, ht1))
  } else {
    return(ht1)
  }
}

#' Rank signaling networks based on the information flow or the number of interactions
#'
#' This function can also be used to rank signaling from certain cell groups to other cell groups
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param measure "weight" or "count". "weight": comparing the total interaction weights (strength); "count": comparing the number of interactions;
#' @param mode "single","comparison"
#' @param comparison a numerical vector giving the datasets for comparison; a single value means ranking for only one dataset and two values means ranking comparison for two datasets
#' @param color.use defining the color for each cell group
#' @param stacked whether plot the stacked bar plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a vector giving the signaling pathway to show
#' @param pairLR a vector giving the names of L-R pairs to show (e.g, pairLR = c("IL1A_IL1R1_IL1RAP","IL1B_IL1R1_IL1RAP"))
#' @param signaling.type a char giving the types of signaling from the three categories c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
#' @param do.stat whether do a paired Wilcoxon test to determine whether there is significant difference between two datasets. Default = FALSE
#' @param cutoff.pvalue the cutoff of pvalue when doing Wilcoxon test; Default = 0.05
#' @param tol a tolerance when considering the relative contribution being equal between two datasets. contribution.relative between 1-tol and 1+tol will be considered as equal contribution
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @param do.flip whether flip the x-y axis
#' @param x.angle,y.angle,x.hjust,y.hjust parameters for rotating and spacing axis labels
#' @param axis.gap whetehr making gaps in y-axes
#' @param ylim,segments,tick_width,rel_heights parameters in the function gg.gap when making gaps in y-axes
#' e.g., ylim = c(0, 35), segments = list(c(11, 14),c(16, 28)), tick_width = c(5,2,5), rel_heights = c(0.8,0,0.1,0,0.1)
#' https://tobiasbusch.xyz/an-r-package-for-everything-ep2-gaps
#' @param show.raw whether show the raw information flow. Default = FALSE, showing the scaled information flow to provide compariable data scale; When stacked = TRUE, use raw information flow by default.
#' @param return.data whether return the data.frame consisting of the calculated information flow of each signaling pathway or L-R pair
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param font.size font size

#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankNet_modified <- function(object, slot.name = "netP", measure = c("weight","count"),
                    mode = c("comparison", "single"), comparison = c(1,2), color.use = NULL,
                    stacked = FALSE, sources.use = NULL, targets.use = NULL,  signaling = NULL,
                    pairLR = NULL, signaling.type = NULL, do.stat = FALSE, cutoff.pvalue = 0.05,
                    tol = 0.05, thresh = 0.05, show.raw = FALSE, return.data = FALSE,
                    x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8,
                    do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1, y.hjust = 1,
                    significant.only = FALSE, axis.gap = FALSE, ylim = NULL, segments = NULL,
                    tick_width = NULL, rel_heights = c(0.9,0,0.1)) {
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  } else if (measure == "count") {
    ylabel = "Number of interactions"
  }
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)
    prob = object1$prob
    prob[object1$pval > thresh] <- 0
    if (measure == "count") {
      prob <- 1*(prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        } else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        } else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[ ,idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }

    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      pSum <- -1/log(pSum)
      pSum[is.na(pSum)] <- 0
      idx1 <- which(is.infinite(pSum) | pSum < 0)
      values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
      position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    } else if (measure == "count") {
      pSum <- pSum.original
    }

    pair.name <- names(pSum)

    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }

    if (!is.null(signaling.type)) {
      LR <- subset(object@DB$interaction, annotation %in% signaling.type)
      if (slot.name == "netP") {
        signaling <- unique(LR$pathway_name)
      } else if (slot.name == "net") {
        pairLR <- LR$interaction_name
      }
    }

    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab(ylabel) + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }

  } else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- object.list$prob
      prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1*(prob > 0)
      }
      prob.list[[i]] <- prob
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[ ,idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- -1/log(pSum.original[[i]])
        pSum[[i]][is.na(pSum[[i]])] <- 0
        idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
        pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      } else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    if (measure == "weight") {
      values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix
      for (i in 1:length(comparison)) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
        } else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
        }
      }
    }



    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0, group = object.names[comparison[i]], row.names = pair.name.all)
      df[[i]][pair.name[[i]],3] <- pSum[[i]]
      df[[i]][pair.name[[i]],2] <- pSum.original[[i]]
    }


    # contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution/abs(df[[1]]$contribution), digits=1))
    # #  contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution.scaled/abs(df[[1]]$contribution.scaled), digits=1))
    # contribution.relative2 <- as.numeric(format(df[[length(comparison)-1]]$contribution/abs(df[[1]]$contribution), digits=1))
    # contribution.relative[is.na(contribution.relative)] <- 0
    # for (i in 1:length(comparison)) {
    #   df[[i]]$contribution.relative <- contribution.relative
    #   df[[i]]$contribution.relative2 <- contribution.relative2
    # }
    # df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    # idx <- with(df[[1]], order(-contribution.relative,  -contribution.relative2, contribution, -contribution.data2))
    #
    contribution.relative <- list()
    for (i in 1:(length(comparison)-1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison)-i+1]]$contribution/df[[1]]$contribution, digits=1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2,contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }



    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL

    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)

    if (is.null(color.use)) {
      color.use =  ggPalette(length(comparison))
    }

    # https://stackoverflow.com/questions/49448497/coord-flip-changes-ordering-of-bars-within-groups-in-grouped-bar-plot
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)

    # perform statistical analysis
    # if (do.stat) {
    #   pvalues <- c()
    #   for (i in 1:length(pair.name.all)) {
    #     df.prob <- data.frame()
    #     for (j in 1:length(comparison)) {
    #       if (pair.name.all[i] %in% pair.name[[j]]) {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(prob.list[[j]][ , , pair.name.all[i]]), group = comparison[j]))
    #       } else {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(matrix(0, nrow = nrow(prob.list[[j]]), ncol = nrow(prob.list[[j]]))), group = comparison[j]))
    #       }
    #
    #     }
    #     df.prob$group <- factor(df.prob$group, levels = comparison)
    #     if (length(comparison) == 2) {
    #       pvalues[i] <- wilcox.test(prob ~ group, data = df.prob)$p.value
    #     } else {
    #       pvalues[i] <- kruskal.test(prob ~ group, data = df.prob)$p.value
    #     }
    #   }
    #   df$pvalues <- pvalues
    # }
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][ , , pair.name.all[i]])
          } else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) >3 & sum(is.na(prob.values)) == 0) {
          pvalues <- wilcox.test(prob.values[ ,1], prob.values[ ,2], paired = TRUE)$p.value
        } else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }

    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), ]
      }
    }

    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2], ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
        if (significant.only == TRUE) {
          pathway.select <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                            (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)
          df <- df[pathway.select, ] # subset df to only include the relevant pathway
          colors.text <- colors.text[colors.text!="black"]
        }
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
      }
    } else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }

    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    if (stacked) {
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position ="fill") # +
      # xlab("") + ylab("Relative information flow") #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      #  scale_y_discrete(breaks=c("0","0.5","1")) +
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      } else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }

      gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)
    } else {
      if (show.raw) {
        gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      } else {
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      }

      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg,
                       ylim = ylim,
                       segments = segments,
                       tick_width = tick_width,
                       rel_heights = rel_heights)
      }
    }
    gg <- gg +  CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }

    } else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text)))

    }

    gg <- gg + theme(axis.text=element_text(size=font.size),
                     axis.title=element_text(size=font.size*1.5),
                     legend.text=element_text(size=font.size),
                     )
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, hjust=x.hjust),
                     axis.text.y = element_text(angle = y.angle, hjust=y.hjust))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }
  }

  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}


get_signalingRole_data <- function(object, signaling = NULL, slot.name = "netP", 
x.measure = "outdeg", y.measure = "indeg") {
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  if (sum(c(x.measure, y.measure) %in% names(slot(object, slot.name)$centr[[1]])) !=2) {
    stop(paste0("`x.measure, y.measure` should be one of ", paste(names(slot(object, slot.name)$centr[[1]]),collapse=", "), '\n', "`outdeg_unweighted` is only supported for version >= 1.1.2"))
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]][[x.measure]]
    incoming[,i] <- centr[[i]][[y.measure]]
  }
  if (is.null(signaling)) {
    message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
  } else {
    message("Signaling role analysis on the cell-cell communication network from user's input")
    signaling <- signaling[signaling %in% object@netP$pathways]
    if (length(signaling) == 0) {
      stop('There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`')
    }
    outgoing <- outgoing[ , signaling, drop = FALSE]
    incoming <- incoming[ , signaling, drop = FALSE]
  }
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  num.link <- aggregateNet(object, signaling = signaling, return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link)-diag(num.link)
  df <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells),
                   Count = num.link)
  df$labels <- factor(df$labels, levels = names(incoming.cells))
  return(df)
}
