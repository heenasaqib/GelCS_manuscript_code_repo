#' Subset CellChat object using a portion of cells
#'
#' @param object  A CellChat object (either an object from a single dataset or a merged objects from multiple datasets)
#' @param cells.use a char vector giving the cell barcodes to subset. If cells.use = NULL, USER must define `idents.use`
#' @param idents.use a subset of cell groups used for analysis
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param invert whether invert the idents.use
#' @param thresh threshold of the p-value for determining significant interaction. A parameter as an input of the function `computeCommunProbPathway`
#' @importFrom methods slot new
#'
#' @return
#' @export
#'
subsetCellChat_mod <- function(object, cells.use = NULL, idents.use = NULL, group.by = NULL, invert = FALSE, thresh = 0.05) {
  if (!is.null(idents.use)) {
    if (is.null(group.by)) {
      labels <- object@idents
      if (object@options$mode == "merged") {
        message("Use the joint cell labels from the merged CellChat object")
        labels <- object@idents$joint
      }
    } else {
      labels <- object@meta[[group.by]]
    }
    if (!is.factor(labels)) {
      labels <- factor(labels)
    }
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(labels)]

    if (invert) {
      level.use <- level.use[!(level.use %in% idents.use)]
    } else {
      level.use <- level.use[level.use %in% idents.use]
    }
    cells.use.index <- which(as.character(labels) %in% level.use)
    cells.use <- names(labels)[cells.use.index]
  } else if (!is.null(cells.use)) {
    labels <- object@idents
    if (object@options$mode == "merged") {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
    cells.use.index <- which(names(labels) %in% cells.use)
  } else {
    stop("USER should define either `cells.use` or `idents.use`!")
  }
  cat("The subset of cell groups used for CellChat analysis are ", level.use, '\n')

  if (nrow(object@data) > 0) {
    data.subset <- object@data[, cells.use.index]
  } else {
    data.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  if (nrow(object@data.project) > 0) {
    data.project.subset <- object@data.project[, cells.use.index]
  } else {
    data.project.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  data.signaling.subset <- object@data.signaling[, cells.use.index]

  meta.subset <- object@meta[cells.use.index, , drop = FALSE]


  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)
    net.subset <- vector("list", length = length(object@net))
    netP.subset <- vector("list", length = length(object@netP))
    idents.subset <- vector("list", length = length(idents))
    names(net.subset) <- names(object@net)
    names(netP.subset) <- names(object@netP)
    names(idents.subset) <- names(object@idents[1:(length(object@idents)-1)])
    images.subset <- vector("list", length = length(idents))
    names(images.subset) <- names(object@idents[1:(length(object@idents)-1)])

    for (i in 1:length(idents)) {
      cat("Update slots object@images, object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      images <- object@images[[i]]
      for (images.j in names(images)) {
        values <- images[[images.j]]
        if (images.j %in% c("coordinates")) {
          values.new <- values[cells.use.index, ]
          images[[images.j]] <- values.new
        }
        if (images.j %in% c("distance")) {
          values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
          images[[images.j]] <- values.new
        }
      }
      images.subset[[i]] <- images

      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- values[group.existing.index, group.existing.index, ]
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- values[group.existing.index, group.existing.index]
          net[[net.j]] <- values.new
        }
       # net[[net.j]] <- values.new
      }
      net.subset[[i]] <- net

      # cat("Update slot object@netP...", '\n')
      # netP <- object@netP[[i]]
      # for (netP.j in names(netP)) {
      #   values <- netP[[netP.j]]
      #   if (netP.j %in% c("pathways")) {
      #     values.new <- values
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("prob")) {
      #     values.new <- values[group.existing.index, group.existing.index, ]
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("centr")) {
      #     for (k in 1:length(values)) {
      #       values.new <- lapply(values, function(x) {
      #         values.new2 <- lapply(x, function(x) {
      #           values.new1 <- x[group.existing.index]
      #           names(values.new1) <- group.existing
      #           return(values.new1)
      #         })
      #         names(values.new2) <- names(x)
      #         return(values.new2)
      #       })
      #       names(values.new) <- names(values)
      #     }
      #   }
      #   netP[[netP.j]] <- values.new
      # }
      netP = computeCommunProbPathway(net = net.subset[[i]], pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
      netP$centr = netAnalysis_computeCentrality(net =  net.subset[[i]]$prob)
      netP.subset[[i]] <- netP
      idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells.use]
      idents.subset[[i]] <- factor(idents.subset[[i]], levels = levels(idents[[i]])[levels(idents[[i]]) %in% level.use])
    }
    idents.subset$joint <- factor(object@idents$joint[cells.use.index], levels = level.use)

  } else {
    cat("Update slots object@images, object@net, object@netP in a single dataset...", '\n')

    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)

    images <- object@images
    for (images.j in names(images)) {
      values <- images[[images.j]]
      if (images.j %in% c("coordinates")) {
        values.new <- values[cells.use.index, ]
        images[[images.j]] <- values.new
      }
      if (images.j %in% c("distance")) {
        values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
        images[[images.j]] <- values.new
      }
    }
    images.subset <- images


    # cat("Update slot object@net...", '\n')
    net <- object@net
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
        net[[net.j]] <- values.new
      }
    }
    net.subset <- net

    # cat("Update slot object@netP...", '\n')
    # netP <- object@netP
    # for (netP.j in names(netP)) {
    #   values <- netP[[netP.j]]
    #   if (netP.j %in% c("pathways")) {
    #     values.new <- values
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("prob")) {
    #     values.new <- values[group.existing.index, group.existing.index, ]
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("centr")) {
    #     for (k in 1:length(values)) {
    #       values.new <- lapply(values, function(x) {
    #         values.new2 <- lapply(x, function(x) {
    #           values.new1 <- x[group.existing.index]
    #           names(values.new1) <- group.existing
    #           return(values.new1)
    #         })
    #         names(values.new2) <- names(x)
    #         return(values.new2)
    #       })
    #       names(values.new) <- names(values)
    #     }
    #   }
    #   netP[[netP.j]] <- values.new
    # }
    netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, thresh = thresh)
    netP$centr = netAnalysis_computeCentrality(net = net.subset$prob)
    netP.subset <- netP
    idents.subset <- object@idents[cells.use.index]
    idents.subset <- factor(idents.subset, levels = level.use)
  }


  object.subset <- methods::new(
    Class = "CellChat",
    data = data.subset,
    data.signaling = data.signaling.subset,
    data.project = data.project.subset,
    images = images.subset,
    net = net.subset,
    netP = netP.subset,
    meta = meta.subset,
    idents = idents.subset,
    var.features = object@var.features,
    LR = object@LR,
    DB = object@DB,
    options = object@options
  )
  return(object.subset)
}