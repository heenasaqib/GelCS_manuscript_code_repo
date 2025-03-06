netVisual_aggregate_check <- function(object, signaling, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord","spatial"),
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                slice.use = NULL, alpha.image = 0.15, point.size = 1.5,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  }  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1,2), sum)
    if (vertex.weight == "incoming"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$indeg
    } else if (vertex.weight == "outgoing"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$outdeg
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    meta.t <- object@meta
    meta.t$labels <- labels
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, meta = meta.t, slice.use = slice.use, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)

  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }

  return(gg)

}

searchPair <- function(signaling = c(), pairLR.use, key = c("pathway_name","ligand"), matching.exact = FALSE, pair.only = TRUE) {
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(
    X = 1:length(signaling),
    FUN = function(x) {
      if (!matching.exact) {
        index <- grep(signaling[x], pairLR.use[[key]])
      } else {
        index <- which(pairLR.use[[key]] %in% signaling[x])
      }
      if (length(index) > 0) {
        if (pair.only) {
          pairLR <- dplyr::select(pairLR.use[index, ], interaction_name, pathway_name, ligand, receptor)
        } else {
          pairLR <- pairLR.use[index, ]
        }
        return(pairLR)
      } else {
        stop(cat(paste("Cannot find ", signaling[x], ".", "Please input a correct name!"),'\n'))
      }
    }
  )
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4*i-3, 4*i-2, 4*i-1, 4*i)]), ncol=4, byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[,1]
  } else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i*ncol(pairLR.use)-(ncol(pairLR.use)-1)):(i*ncol(pairLR.use))]), ncol=ncol(pairLR.use), byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[,1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}
