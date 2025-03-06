# Introduction to this script ----
# The goal of this script is to perform GSEA downstream scRNA-seq analysis
# We will have different modules for inter Mø comparison combining Gel100 and Gel50 samples
# We will focus on two methods of obtaining DEGs:
# 1). From Seurat Object, perform FindMarkers() or FindAllMarkers() to obtain a list of DEGs
# 2). From 10X `raw`, perform data preparation including normalization and QC to obtain DEGs (hold)
# November 2023

# Load Libraries ----
library(Seurat)
library(tidyverse)
library(DT)
library(clusterProfiler)
#library(biomaRt)
library(org.Mm.eg.db)
library(GSEABase)
library(stringr)
library(enrichplot)
library(msigdbr) # access to msigdb collections directly within R
library(ggplot2)

options("DT.TOJSON_ARGS" = list(na = "string"))
options(stringsAsFactors = F)
source("enrichment_analysis_func.R")
# Set WD ----
setwd("~/Desktop/School Work/Stanford/Yang Lab/2023 Spring rotation_Heena/")

# Load Dataset ----
alldata <- readRDS("CS scRNAseq data share/CS0_50_RAW data/gelCS_alldata.rds")

# Create Cluster Subset ----
cluster0 <- subset(x = alldata, idents = "M1-like")
cluster1 <- subset(x = alldata, idents = "M2-like 1")
cluster5 <- subset(x = alldata, idents = "M2-like 2")
cluster015 <- subset(x = alldata, idents = c("M1-like", "M2-like 1", "M2-like 2"))

clustermsc <- subset(x = alldata, idents = c("Stromal 1","Stromal 2")) # Stromal 1 -> msc; Stromal 2 -> fibroblast

##### GO Analysis ##### ----
# change ident to orig.ident for comparison within the groups
# Default test.use == wilcox; see manual for FindAllMarkers() to see other options
DEGs_møs <- FindAllMarkers(cluster015, min.pct = 0.05)

DEGs_msc <- FindAllMarkers(clustermsc, min.pct = 0.05)

# Filter out low significance and FC genes
p_val_thresh <- 0.05 
logFC_thresh <- 0.5 #default = 0.1; increasing this param reduces weak signals

DEGs_møs_ori <- DEGs_møs
DEGs_møs <- DEGs_møs %>%
            filter(p_val_adj < p_val_thresh & abs(avg_log2FC) > logFC_thresh) %>%
            mutate(regulation = ifelse(avg_log2FC<0, "down", "up")) %>%
            mutate(rank = abs(avg_log2FC)*(-log(p_val_adj)))
            # add <rank> parameter for filtering 'best' DEGs

# Unify keytypes of `gene` - currently, most are in 'SYMBOL', few are in 'ENSEMBL'
# Convert all 'ENSEMBL' to 'SYMBOL'
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # for mouse genes
ensembl_ids <- DEGs_møs$gene[grep("^ENSMUSG|^ENSG", DEGs_møs$gene)]

gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = ensembl_ids,
                   mart = ensembl)

# Rename the columns for clarity
colnames(gene_info) <- c("ensembl_id", "gene_symbol")

# Replace the ENSEMBL IDs in df$genes with the corresponding gene symbols
DEGs_møs$gene_SYMBOL <- sapply(DEGs_møs$gene, function(x) {
    symbol <- gene_info$gene_symbol[match(x, gene_info$ensembl_id)]
    if (!is.na(symbol)) return(symbol) else return(x)
})

# Convert SYMBOL gene annotations to ENTREZID for GO/GSEA analysis
ENTREZID <- bitr(DEGs_møs$gene_SYMBOL,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Mm.eg.db)

# Assign ENTREZID to each gene in DEGs_møs
temp <- inner_join(DEGs_møs, ENTREZID, by=c("gene_SYMBOL"="SYMBOL"))

# Subset for subpopulations for macrophages
DEGs_M1 <- subset(temp, cluster == "M1-like")
DEGs_M2l1 <- subset(temp, cluster == "M2-like 1")
DEGs_M2l2 <- subset(temp, cluster == "M2-like 2")

DEGs_fib <- subset(temp, cluster == "Fibroblasts")
DEGs_stm <- subset(temp, cluster == "Stromal cells")

rownames(DEGs_M1) <- DEGs_M1$gene
rownames(DEGs_M2l1) <- DEGs_M2l1$gene
rownames(DEGs_M2l2) <- DEGs_M2l2$gene

rownames(DEGs_fib) <- DEGs_fib$gene
rownames(DEGs_stm) <- DEGs_stm$gene

write.csv(DEGs_M1, "figure_output/all/manuscript/DEGs_M1.csv", row.names = TRUE)
write.csv(DEGs_M2l1, "figure_output/all/manuscript/DEGs_M2like1.csv", row.names = TRUE)
write.csv(DEGs_M2l2, "figure_output/all/manuscript/DEGs_M2like2.csv", row.names = TRUE)
write.csv(DEGs_fib, "DEGs_fib.csv", row.names = TRUE)
write.csv(DEGs_stm, "DEGs_stm.csv", row.names = TRUE)

# Heatmap showing DEG within Mø
DEGs_møs %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(cluster015, features = top10$gene) + NoLegend()

# Visualize datatable---
datatable(DEGs_M2l1,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 2: DEGs in M2-like 1 Mø',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 20,
                         lengthMenu = c("10", "25", "50", "100"),
                         columnDefs = list(list(targets = 9,  # Adjust column index as needed
                         render = JS(
                            "function(data, type, row, meta) {
                                if (type === 'sort' || type === 'type') {
                                    return (data === 'Inf' ? 1e6 : data);
                                    // Replace Infinity with a large number for sorting
                                }
                                return data; // For display, return the data as is
                            }"
                         )))),
          editable = TRUE) %>%
  formatRound(columns=c(2:4), digits=3) %>% #formats column style - rounding
  formatSignif(columns=c(1,5), digits=3)    #formats column style - scientific notation

# double check: any intersection of upregulated genes across each cell type
temp1 <- DEGs_M1$gene_SYMBOL[DEGs_M1$regulation=="up"]
temp2 <- DEGs_M2l1$gene_SYMBOL[DEGs_M2l1$regulation=="up"]
temp3 <- DEGs_M2l2$gene_SYMBOL[DEGs_M2l2$regulation=="up"]

cat("Intersected Upregulated Genes among M1 and M2l1 Møs Are: \n",  intersect(temp1, temp2),'\n')
cat("Intersected Upregulated Genes among M2l1 and M2l2 Møs Are: \n", intersect(temp2, temp3), '\n')
cat("Intersected Upregulated Genes among M1 and M2l2 Møs Are: \n", intersect(temp1, temp3), '\n')


# GO Analysis for M1-like phenotype
ego_M1.down.BP  <- ego(DEGs_M1, reg.direction = "down", ont = "BP", keyType = "ENTREZID")
ego_M1.down.CC  <- ego(DEGs_M1, reg.direction = "down", ont = "CC", keyType = "ENTREZID")
ego_M1.down.MF  <- ego(DEGs_M1, reg.direction = "down", ont = "MF", keyType = "ENTREZID")
ego_M1.down.all <- ego(DEGs_M1, reg.direction = "down", ont = "all", keyType = "ENTREZID")

ego_M1.up.BP    <- ego(DEGs_M1, reg.direction = "up", ont = "BP", keyType = "ENTREZID")
ego_M1.up.CC    <- ego(DEGs_M1, reg.direction = "up", ont = "CC", keyType = "ENTREZID")
ego_M1.up.MF    <- ego(DEGs_M1, reg.direction = "up", ont = "MF", keyType = "ENTREZID")
ego_M1.up.all   <- ego(DEGs_M1, reg.direction = "up", ont = "all", keyType = "ENTREZID")

# Plot ALL EGO, split by "ONTOLOGY"
gg1 <- dotplot(ego_M1.down.all$ego.obj, split="ONTOLOGY", title = "M1 Downregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")
gg2 <- dotplot(ego_M1.up.all$ego.obj, split="ONTOLOGY", title = "M1 Upregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")
gg1

# Plot individual EGO
gg1 <- dotplot(ego_M1.down.BP$ego.obj, showCategory = 15,
       title = "M1 Downregulated Gene for Biological Processes")
gg2 <- dotplot(ego_M1.up.BP$ego.obj, showCategory = 15,
         title = "M1 Upregulated Gene for Biological Processes")
gg1+gg2

ggsave("figure_output/all/manuscript/GO_mø/GO_M1_down_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M1_up_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

# Plot EGO Functional Map
gg1 <- goplot(ego_M1.down.BP$ego.obj, title = "M1 Downregulated Gene for Biological Processes")
gg2 <- goplot(ego_M1.up.BP$ego.obj, title = "M1 Upregulated Gene for Biological Processes")

ggsave("figure_output/all/manuscript/GO_mø/GO_M1_down_MAP_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M1_up_MAP_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

length(ego_M1.up.MF$ego.obj@gene)

write.csv(ego_M1.down.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M1_down_BP.csv", row.names = TRUE)
write.csv(ego_M1.up.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M1_up_BP.csv", row.names = TRUE)

temp <- simplify(ego_M2l2.up.BP$ego.obj, cutoff=0.6)
write.csv(temp, "figure_output/all/manuscript/GO_mø/GO_raw/GO_M2l2_up_BP_simplified.csv", row.names = TRUE)


##### GO Analysis for M2-like 1 phenotype #####
ego_M2l1.down.BP  <- ego(DEGs_M2l1, reg.direction = "down", ont = "BP")
ego_M2l1.down.CC  <- ego(DEGs_M2l1, reg.direction = "down", ont = "CC")
ego_M2l1.down.MF  <- ego(DEGs_M2l1, reg.direction = "down", ont = "MF")
ego_M2l1.down.all <- ego(DEGs_M2l1, reg.direction = "down", ont = "all")

ego_M2l1.up.BP    <- ego(DEGs_M2l1, reg.direction = "up", ont = "BP")
ego_M2l1.up.CC    <- ego(DEGs_M2l1, reg.direction = "up", ont = "CC")
ego_M2l1.up.MF    <- ego(DEGs_M2l1, reg.direction = "up", ont = "MF")
ego_M2l1.up.all   <- ego(DEGs_M2l1, reg.direction = "up", ont = "all")

gg1 <- dotplot(ego_M2l1.down.all$ego.obj, split="ONTOLOGY", title = "M2-like 1 Downregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")
gg2 <- dotplot(ego_M2l1.up.all$ego.obj, split="ONTOLOGY", title = "M2-like 1 Upregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_down_all.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_up_all.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

write.csv(ego_M2l1.down.all$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l1_down_all.csv", row.names = TRUE)
write.csv(ego_M2l1.up.all$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l1_up_all.csv", row.names = TRUE)
write.csv(ego_M2l1.down.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l1_down_BP.csv", row.names = TRUE)
write.csv(ego_M2l1.up.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l1_up_BP.csv", row.names = TRUE)

# Plot individual EGO
gg1 <- dotplot(ego_M2l1.down.BP$ego.obj, showCategory = 15,
       title = "M2-like 1 Downregulated Gene for Biological Processes")
gg2 <- dotplot(ego_M2l1.up.BP$ego.obj, showCategory = 15,
         title = "M2-like 1 Upregulated Gene for Biological Processes")
gg1+gg2

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_down_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_up_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

# Plot EGO Functional Map
gg1 <- goplot(ego_M2l1.down.BP$ego.obj, title = "M2-like 1 Downregulated Gene for Biological Processes")
gg2 <- goplot(ego_M2l1.up.BP$ego.obj, title = "M2-like 1 Upregulated Gene for Biological Processes")

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_down_MAP_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l1_up_MAP_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")


##### GO Analysis for M2-like 2 phenotype
ego_M2l2.down.BP  <- ego(DEGs_M2l2, reg.direction = "down", ont = "BP")
ego_M2l2.down.CC  <- ego(DEGs_M2l2, reg.direction = "down", ont = "CC")
ego_M2l2.down.MF  <- ego(DEGs_M2l2, reg.direction = "down", ont = "MF")
ego_M2l2.down.all <- ego(DEGs_M2l2, reg.direction = "down", ont = "all")

ego_M2l2.up.BP    <- ego(DEGs_M2l2, reg.direction = "up", ont = "BP")
ego_M2l2.up.CC    <- ego(DEGs_M2l2, reg.direction = "up", ont = "CC")
ego_M2l2.up.MF    <- ego(DEGs_M2l2, reg.direction = "up", ont = "MF")
ego_M2l2.up.all   <- ego(DEGs_M2l2, reg.direction = "up", ont = "all")

gg1 <- dotplot(ego_M2l2.down.all$ego.obj, split="ONTOLOGY", title = "M2-like 2 Downregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")
gg2 <- dotplot(ego_M2l2.up.all$ego.obj, split="ONTOLOGY", title = "M2-like 2 Upregulated Gene EGO") +
       facet_grid(ONTOLOGY~., scale="free")

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_down_all.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_up_all.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

write.csv(ego_M2l2.down.all$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l2_down_all.csv", row.names = TRUE)
write.csv(ego_M2l2.up.all$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l2_up_all.csv", row.names = TRUE)
write.csv(ego_M2l2.down.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l2_down_BP.csv", row.names = TRUE)
write.csv(ego_M2l2.up.BP$ego.obj, "figure_output/all/manuscript/GO_mø/GO_M2l2_up_BP.csv", row.names = TRUE)

# Plot individual EGO
gg1 <- dotplot(ego_M2l2.down.BP$ego.obj, showCategory = 15,
       title = "M2-like 2 Downregulated Gene for Biological Processes")
gg2 <- dotplot(ego_M2l2.up.BP$ego.obj, showCategory = 15,
         title = "M2-like 2 Upregulated Gene for Biological Processes")
gg1+gg2

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_down_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_up_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")

# Plot EGO Functional Map
gg1 <- goplot(ego_M2l2.down.BP$ego.obj, title = "M2-like 2 Downregulated Gene for Biological Processes")
gg2 <- goplot(ego_M2l2.up.BP$ego.obj, title = "M2-like 2 Upregulated Gene for Biological Processes")

ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_down_MAP_BP.jpeg", plot = gg1,  width = 10, height = 12, unit = "in")
ggsave("figure_output/all/manuscript/GO_mø/GO_M2l2_up_MAP_BP.jpeg", plot = gg2,  width = 10, height = 12, unit = "in")



##### test out: reading csv and plotting dot plot #####
gg2 <- dotplot(ego_M2l2.up.BP$ego.obj, showCategory = 15,
       title = "M2-like 2 Upregulated Gene for Biological Processes")

temp <- read.csv("figure_output/all/manuscript/GO_mø/GO_raw/GO_M2l2_up_BP.csv")
data <- as.data.frame(temp)
data <- data %>% top_n(15, Count) %>% arrange(desc(Count))
S1 <- ggplot(data, aes(x= GeneRatio, y=Description, size=Count, color=p.adjust)) +
      geom_point(alpha = 0.8) + theme_classic() +
      scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab",
      limit = c(0.000000000000000007, 0.002)) + scale_size(range = c(2, 8))
S1 



##### GO Analysis for MSC #####
ego_fibroblast.down.BP  <- ego(DEGs_fib, reg.direction = "down", ont = "BP")
ego_fibroblast.up.BP    <- ego(DEGs_fib, reg.direction = "up", ont = "BP")

ego_stromal.down.BP  <- ego(DEGs_stm, reg.direction = "down", ont = "BP")
ego_stromal.up.BP    <- ego(DEGs_stm, reg.direction = "up", ont = "BP")

ego_fibroblast.down.BP$gg.obj
ego_fibroblast.up.BP$gg.obj

ego_stromal.down.BP$gg.obj
ego_stromal.up.BP$gg.obj

write.csv(ego_fibroblast.down.BP$ego.obj, "GO_Fibroblast_down_BP.csv", row.names = TRUE)
write.csv(ego_fibroblast.up.BP$ego.obj, "GO_Fibroblast_up_BP.csv", row.names = TRUE)

# GOST plot - Example
gost.res <- gost(genes.use, organism = "mmusculus", correction_method = "fdr")
mygostplot <- gostplot(gost.res, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0071345", "GO:0002833", "GO:1902533", "REAC:R-MMU-168249"),
  filename = "Manhattan Plot for M1-like Cells.jpeg",
  width = NA,
  height = NA)

source("enrichment_analysis.R")
# KEGG Analysis ------
# KEGG Analysis for M1 phenotype
ekegg_M1.down <- ekegg(DEGs_M1, reg.direction = "down")
ekegg_M1.up   <- ekegg(DEGs_M1, reg.direction = "up")

# KEGG Analysis for M2-like 1 phenotype
ekegg_M2l1.down <- ekegg(DEGs_M2l1, reg.direction = "down")
ekegg_M2l1.up   <- ekegg(DEGs_M2l1, reg.direction = "up")

# KEGG Analysis for M2-like 2 phenotype
ekegg_M2l2.down <- ekegg(DEGs_M2l2, reg.direction = "down")
ekegg_M2l2.up   <- ekegg(DEGs_M2l2, reg.direction = "up")

# KEGG Visualization
gg_down <- ekegg_M2l1.down$gg.obj + ggtitle("M2l1 Downregulated Gene KEGG Enrichment")
gg_up   <- ekegg_M2l1.up$gg.obj + ggtitle("M2l1 Upregulated Gene KEGG Enrichment")

gg_down
gg_up

# GSEA Analysis ------

# Load background gene
geneset.m5 <- read.gmt("background gene/m5.all.v2023.1.Mm.symbols.gmt") # ontology
geneset.m8 <- read.gmt("background gene/m8.all.v2023.2.Mm.symbols.gmt") # cell type signature
geneset.c7 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C7") %>% # mouse immunogenic gene from C7
              dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
geneset.c2 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C2") %>% # mouse immunogenic gene from C7
              dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 


geneset.m5$term = str_remove(geneset$term,"HALLMARK_") # remove any portion of gene name

# GSEA on M1 Phenotype
geneList <- DEGs_M1$avg_log2FC
names(geneList) <- DEGs_M1$gene_SYMBOL
geneList <- sort(geneList, decreasing = TRUE)

egmt.m5 <- GSEA(geneList, TERM2GENE=geneset.m5,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.m8 <- GSEA(geneList, TERM2GENE=geneset.m8,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c7 <- GSEA(geneList, TERM2GENE=geneset.c7,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c2 <- GSEA(geneList, TERM2GENE=geneset.c2,verbose=F,pvalueCutoff = 0.5, eps = 0)

write.csv(egmt.m8, "GSEA_M1_M8.csv", row.names = TRUE)
write.csv(egmt.c7, "GSEA_M1_C7.csv", row.names = TRUE)
write.csv(egmt.c2, "GSEA_M1_C2.csv", row.names = TRUE)

dotplot(egmt.m5,split=".sign", showCategory = 10) + facet_grid(~.sign)
dotplot(egmt.m8,split=".sign", showCategory = 8) + facet_grid(~.sign)
dotplot(egmt.c7,split=".sign", showCategory = 8) + facet_grid(~.sign)

gseaplot2(egmt.m5, geneSetID = 1, title = egmt.m5$Description[1])
gseaplot2(egmt.m8, geneSetID = 5, title = egmt.m8$Description[5])

# Ridge plot to visualize GSEA correlation among features
ridgeplot(egmt.m8, showCategory = 10, fill = "p.adjust", core_enrichment = TRUE)

# GSEA on M2-like 1 Phenotype
geneList <- DEGs_M2l1$avg_log2FC
names(geneList) <- DEGs_M2l1$gene_SYMBOL
geneList <- sort(geneList, decreasing = TRUE)

egmt.m5 <- GSEA(geneList, TERM2GENE=geneset.m5,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.m8 <- GSEA(geneList, TERM2GENE=geneset.m8,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c7 <- GSEA(geneList, TERM2GENE=geneset.c7,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c2 <- GSEA(geneList, TERM2GENE=geneset.c2,verbose=F,pvalueCutoff = 0.5, eps = 0)

write.csv(egmt.m8, "GSEA_M2l1_M8.csv", row.names = TRUE)
write.csv(egmt.c7, "GSEA_M2l1_C7.csv", row.names = TRUE)
write.csv(egmt.c2, "GSEA_M2l1_C2.csv", row.names = TRUE)

dotplot(egmt.m5,split=".sign", showCategory = 10) + facet_grid(~.sign)
dotplot(egmt.m8,split=".sign", showCategory = 8) + facet_grid(~.sign)
dotplot(egmt.c7,split=".sign", showCategory = 8) + facet_grid(~.sign)

gseaplot2(egmt.m5, geneSetID = 1, title = egmt.m5$Description[1])
gseaplot2(egmt.m8, geneSetID = 5, title = egmt.m8$Description[5])

# Ridge plot to visualize GSEA correlation among features
ridgeplot(egmt.m8, showCategory = 10, fill = "p.adjust", core_enrichment = TRUE)

# GSEA on M2-like 2 Phenotype
geneList <- DEGs_M2l2$avg_log2FC
names(geneList) <- DEGs_M2l2$gene_SYMBOL
geneList <- sort(geneList, decreasing = TRUE)

egmt.m5 <- GSEA(geneList, TERM2GENE=geneset.m5,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.m8 <- GSEA(geneList, TERM2GENE=geneset.m8,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c7 <- GSEA(geneList, TERM2GENE=geneset.c7,verbose=F,pvalueCutoff = 0.5, eps = 0)
egmt.c2 <- GSEA(geneList, TERM2GENE=geneset.c2,verbose=F,pvalueCutoff = 0.5, eps = 0)

write.csv(egmt.m8, "GSEA_M2l2_M8.csv", row.names = TRUE)
write.csv(egmt.c7, "GSEA_M2l2_C7.csv", row.names = TRUE)
write.csv(egmt.c2, "GSEA_M2l2_C2.csv", row.names = TRUE)

dotplot(egmt.m5,split=".sign", showCategory = 10) + facet_grid(~.sign)
dotplot(egmt.m8,split=".sign", showCategory = 8) + facet_grid(~.sign)
dotplot(egmt.c7,split=".sign", showCategory = 8) + facet_grid(~.sign)
dotplot(egmt.c2,split=".sign", showCategory = 8) + facet_grid(~.sign)

gseaplot2(egmt.m5, geneSetID = 1, title = egmt.m5$Description[1])
gseaplot2(egmt.m8, geneSetID = 5, title = egmt.m8$Description[5])

# Ridge plot to visualize GSEA correlation among features
ridgeplot(egmt.m8, showCategory = 10, fill = "p.adjust", core_enrichment = TRUE)

# plot GSEA table
df.egmt.c2 <- as.data.frame(egmt.c2)[1:10]
df.egmt.c2 <- df.egmt.c2 %>%
              filter(abs(NES)>2 & p.adjust < p_val_thresh)

ggplot(df.egmt.c2, aes(x = reorder(Description, NES), y = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for horizontal bars
    theme_minimal() +
    labs(x = "Gene Set", y = "Normalized Enrichment Score (NES)", title = "M2-like 2 GSEA Results C2")


df.egmt.c7 <- as.data.frame(egmt.c7)
df.egmt.c7 <- df.egmt.c7 %>%
              filter(abs(NES)>2.5 & p.adjust < p_val_thresh)

ggplot(df.egmt.c7, aes(x = reorder(Description, NES), y = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for horizontal bars
    theme_minimal() +
    labs(x = "Gene Set", y = "Normalized Enrichment Score (NES)", title = "M2-like 2 GSEA Results C7")


df.egmt.m8 <- as.data.frame(egmt.m8)
df.egmt.m8 <- df.egmt.m8 %>%
              filter(abs(NES)>1 & p.adjust < p_val_thresh)

ggplot(df.egmt.m8, aes(x = reorder(Description, NES), y = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for horizontal bars
    theme_minimal() +
    labs(x = "Gene Set", y = "Normalized Enrichment Score (NES)", title = "M2-like 2 GSEA Results M8")


# DEPCRECATED CODE -----
# genes.use <- DEGs_M1 %>% 
#              filter(regulation == 'down') %>% # Adjust this as needed
#              pull(gene_SYMBOL)
# ego_M1_BP <- enrichGO(gene         = genes.use,
#                    OrgDb        = org.Mm.eg.db, # Replace with appropriate organism database
#                    keyType      = "SYMBOL", # Ensure this matches your gene ID type
#                    ont          = "BP", # Can be "BP", "CC", or "MF"
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.05,
#                    readable = TRUE)
# ego_M1_CC <- enrichGO(gene         = genes.use,
#                    OrgDb        = org.Mm.eg.db, # Replace with appropriate organism database
#                    keyType      = "SYMBOL", # Ensure this matches your gene ID type
#                    ont          = "CC", # Can be "BP", "CC", or "MF"
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.05,
#                    readable = TRUE)
# ego_M1_MF <- enrichGO(gene         = genes.use,
#                    OrgDb        = org.Mm.eg.db, # Replace with appropriate organism database
#                    keyType      = "SYMBOL", # Ensure this matches your gene ID type
#                    ont          = "MF", # Can be "BP", "CC", or "MF"
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.05,
#                    readable = TRUE)
# head(summary(ego_M1_BP))
# dotplot(ego_M1_BP, showCategory=20)
# dotplot(ego_M1_CC, showCategory=20)
# dotplot(ego_M1_MF, showCategory=20)




# #you can also generate a table of your gost results
# publish_gosttable(
#   gost.res,
#   highlight_terms = c("GO:0071345", "GO:0002833"),
#   use_colors = TRUE,
#   show_columns = c("source", "term_name", "term_size", "intersection_size"),
#   filename = "GOST Table for M1-like Cells.jpeg",
#   ggplot=TRUE)

# ranked_genes <- DEGs_M1 %>% arrange(desc(rank)) %>% select(gene, rank)

