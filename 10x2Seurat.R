# This script is used for preprocessing scRNA dataset for CellChat purposes
# Follow "Inference and analysis of cell-cell communication using CellChat" tutorial uptil part II. 

library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(reticulate)
library(Matrix)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(viridis)
library(BiocManager)
library(gage)
library(data.table)
library(ggforce)
library(gplots)
options(stringsAsFactors = FALSE)


# Set WD -----
setwd("~/Desktop/School Work/Stanford/Yang Lab/2023 Spring rotation_Heena/CS scRNAseq data share/CS0_50_RAW data")

# Open and create seurat object for CS0 = dat -----
dat.sham<-Read10X(
  data.dir = "CS0/filtered_feature_bc_matrix",
  gene.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) 

#Mapping the Ensembl Gene ID back to the first symbol
x<- select(x = org.Mm.eg.db, 
           keys = rownames(dat.sham), 
           column = "SYMBOL", 
           keytype = "ENSEMBL",
           multiVals = "first",asNA=F)

head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL) # issues with excess matches so collapse into a single col

for (name in rownames(dat.sham)) {
  rownames(dat.sham)[match(name, rownames(dat.sham))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(dat.sham)) #check


dat <- CreateSeuratObject(
  dat.sham,
  project = "Gel100",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^Mt") ## if human, "MT"
rownames(dat)
range(dat[["percent.mt"]])## !
VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Open and create seurat object for CS50= dat1 ------
dat.sham1<-Read10X(
  data.dir = "CS50/filtered_feature_bc_matrix", 
  gene.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) 

#Mapping the Ensembl Gene ID back to the first symbol
x<- select(x = org.Mm.eg.db, 
           keys = rownames(dat.sham1), 
           column = "SYMBOL", 
           keytype = "ENSEMBL",
           multiVals = "first",asNA=F)

head(x)
x$FINAL <- ifelse(is.na(x$SYMBOL), x$ENSEMBL, x$SYMBOL) # issues with excess matches so collapse into a single col

for (name in rownames(dat.sham1)) {
  rownames(dat.sham1)[match(name, rownames(dat.sham1))]<-x$FINAL[match(name, x$ENSEMBL)]
}
head(rownames(dat.sham1)) #check


dat1 <- CreateSeuratObject(
  dat.sham1,
  project = "Gel50_CS50",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)

dat1[["percent.mt"]] <- PercentageFeatureSet(dat1, pattern = "^Mt") ## if human, "MT"
rownames(dat1)
range(dat1[["percent.mt"]])## !
VlnPlot(dat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



# Combine CS0 and CS50 = all data ------
alldata <- merge(dat, y=dat1, add.cell.ids=c("Gel100","Gel50_CS50"), project="alldata")
alldata
head(colnames(alldata))
table(alldata$orig.ident) 

###______This is redundant: percent.mt has already been calculated in previous section_____###
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^Mt") ## if human, "MT"
rownames(alldata)
range(alldata[["percent.mt"]])## !
###__I compared columns from previous section to this one, it shows that they are the same__###

# Quality control - since percent.mt < 5 for all, this param not really filtering any out --------
VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
alldata <- subset(alldata, subset = nFeature_RNA > 500 & nFeature_RNA<7000 & nCount_RNA>500 & percent.mt < 5)

dat <- subset(dat, subset = nFeature_RNA > 500 & nFeature_RNA<7000 & nCount_RNA>500 & percent.mt < 5)
dat1 <- subset(dat1, subset = nFeature_RNA > 500 & nFeature_RNA<7000 & nCount_RNA>500 & percent.mt < 5)

VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(alldata$orig.ident)

# Normalization ---------
alldata <- NormalizeData(alldata)
rm(dat, dat.sham, dat.sham1, dat1)

# Variable gene selection -------
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(alldata), 20) # displays 20 most highly variable genes
top20

# Run PCA -----
# scale
alldata <- ScaleData(alldata, features = rownames(alldata))

# dim reduction
alldata <- RunPCA(alldata, features = VariableFeatures(object = alldata))

# Clusters and UMAP ------------
alldata <- FindNeighbors(alldata, dims = 1:20) ## iterativly decide the dims tbh
alldata <- FindClusters(alldata, resolution = 0.3) ## semi - supervise resolution, experiment! youll return through this cycle alot
alldata <- RunUMAP(alldata, dims = 1:20, n.neighbors = 20,min.dist = .01,spread = 6) ## twiddle the dials to get umaps of face validity
#dat <- RunUMAP(dat, dims = 1:15, n.neighbors = 4,min.dist = .001,spread = 5) ## twiddle the dials to get umaps of face validity


## plot 1 : dimension reduction plots
## plot 2 : overlap two samples UMAP or separate UMAP of two samples 
Sub_Gel100 <-subset(x = alldata, subset = orig.ident == "Gel100")
Sub_Gel50_CS50 <-subset(x = alldata, subset = orig.ident == "Gel50_CS50")



# Clusters marker identification -----
# based on r=0.3
alldata <- RenameIdents(object = alldata, 
                        `0` = "M1-like",
                        `1` = "M2-like 1",
                        `2` = "Neutrophil 1",
                        `3` = "NK Cells",
                        `4` = "B Cells",
                        `5` = "M2-like 2",
                        `6` = "Dendritic 1",
                        `7` = "Stromal 2",
                        `8` = "T Cells",
                        `9` = "Proliferating Cells",
                        `10` = "γδ T Cells",
                        `11` = "Dendritic 2",
                        `12` = "Endothelial Cells",
                        `13` = "Neutrophil 2",
                        `14` = "Stromal 1")

gg <- DimPlot(alldata, reduction = "umap", pt.size = 0.2, split.by = "orig.ident")
ggsave("../../figure_output/all/manuscript/umap_condition.jpeg", plot = gg, width = 12, height = 6, units = c("in"), dpi = 1200)
saveRDS(alldata, file = "gelCS_alldata.rds")

# # Optional: Combine stromal cells and fibroblasts to make a msc cluster
# alldata <- RenameIdents(alldata,
#                         `Stromal cells` = "MSC",
#                         `Fibroblasts`   = "MSC")
                        
# DimPlot(alldata, reduction = "umap", pt.size = 0.2)
