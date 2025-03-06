library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(reticulate)
library(Matrix)
library(dplyr)
library(plyr)
library(Seurat)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(viridis)
library(BiocManager)
library(gage)
library(fgsea)
library(data.table)
library(ggforce)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


#----------------------------------------------------------------------
# Project 2: CELLCHAT
#----------------------------------------------------------------------
### Optional: Make cellchat object for Gel100 and Gel50, focusing on msc+mø
# cell.type.use <- c("M1-like", "M2-like 1", "M2-like 2", "Fibroblasts", "Stromal cells")
# alldata <- subset(alldata, idents = cell.type.use)
# levels(alldata) <- cell.type.use

##### Part I: Data Input & Processing and Initialization of CellChat Object #####
alldata_cellchat <- list(data = alldata@assays$RNA@data,
                         meta = data.frame(condition = alldata@meta.data$orig.ident,
                                           cell.type = alldata@active.ident))

data.input <- alldata_cellchat$data # Normalized count matrix
meta <- alldata_cellchat$meta # meta data of cells in rows

# Make cellchat object for Gel100 and Gel50, not focusing on any cells
cell.use <- rownames(meta)[meta$condition == "Gel100"] # Gel100 Gel50_CS50
data.input.100 <- data.input[, cell.use]
meta.100 <- meta[cell.use, ]

cell.use <- rownames(meta)[meta$condition == "Gel50_CS50"] # Gel100 Gel50_CS50
data.input.50 <- data.input[, cell.use]
meta.50 <- meta[cell.use, ]

# Create CellChat Object
cellchat.100 <- createCellChat(object = data.input.100, meta = meta.100, group.by = "cell.type")
cellchat.50  <- createCellChat(object = data.input.50, meta = meta.50, group.by = "cell.type")
cellchat.100 <- setIdent(cellchat.100, ident.use = "cell.type")
cellchat.50  <- setIdent(cellchat.50, ident.use = "cell.type")

# Set the Ligand-Receptor Interaction Database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

which(CellChatDB[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1887,]
which(CellChatDB[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1900,]
cellchat.100@DB <- CellChatDB
cellchat.50@DB <- CellChatDB

# Preprocessing the Expression Data for Cell-Cell Communication Analysis
cellchat.100 <- subsetData(cellchat.100)
cellchat.50  <- subsetData(cellchat.50)
future::plan("multisession", workers = 4)
cellchat.100 <- identifyOverExpressedGenes(cellchat.100)
cellchat.100 <- identifyOverExpressedInteractions(cellchat.100)
cellchat.100 <- projectData(cellchat.100, PPI.mouse) # Project gene expression data onto PPI
cellchat.50  <- identifyOverExpressedGenes(cellchat.50)
cellchat.50  <- identifyOverExpressedInteractions(cellchat.50)
cellchat.50  <- projectData(cellchat.50, PPI.mouse)

##### Part II: Inference of cell-cell communication network #####
cellchat.100 <- computeCommunProb(cellchat.100, raw.use = FALSE, population.size =  TRUE)
cellchat.50  <- computeCommunProb(cellchat.50, raw.use = FALSE, population.size =  TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.100 <- filterCommunication(cellchat.100, min.cells = 10)
cellchat.50  <- filterCommunication(cellchat.50, min.cells = 10)


# Extract the inferred cellular communication network as a data frame
df.net.100 <- subsetCommunication(cellchat.100)
df.net.50  <- subsetCommunication(cellchat.50)

# Infer the cell-cell communication at a signaling pathway level
cellchat.100 <- computeCommunProbPathway(cellchat.100)
cellchat.50  <- computeCommunProbPathway(cellchat.50)

# Calculate the aggregated cell-cell communication network
cellchat.100 <- aggregateNet(cellchat.100)
cellchat.50 <- aggregateNet(cellchat.50)
cellchat.100 <- netAnalysis_computeCentrality(cellchat.100, slot.name = "netP")
cellchat.50 <- netAnalysis_computeCentrality(cellchat.50, slot.name = "netP")

object.list <- list(GEL100 = cellchat.100, GEL50 = cellchat.50)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat.100, file = "cellchat_Gel100.rds")
saveRDS(cellchat.50, file = "cellchat_Gel50.rds")
# saveRDS(cellchat, file = "cellchat_Gel100vsGEL50_Combined.rds")

