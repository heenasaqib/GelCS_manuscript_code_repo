# Introduction to this script ----
# The goal of this script is to continue GO analysis from the Enrichment Analysis Macrophage.R script
# We aim to customize the GO dotplot for macrophages enrichedGO() results

# Load Libraries ----
library(Seurat)
library(tidyverse)
library(DT)
# library(fgsea)
library(clusterProfiler)
# library(gprofiler2)
library(biomaRt)
library(org.Mm.eg.db)
library(GSEABase)
library(stringr)
library(enrichplot)
library(msigdbr) # access to msigdb collections directly within R
library(ggplot2)
library(DOSE)

options("DT.TOJSON_ARGS" = list(na = "string"))
options(stringsAsFactors = F)
source("enrichment_analysis_func.R")
# Set WD ----
setwd("~/Desktop/School Work/Stanford/Yang Lab/2023 Spring rotation_Heena/")

# Read files ----
## Read DEG csv files ##
DEGs_M1 <- read.csv("figure_output/all/manuscript/GO_mø/DEGs_raw/DEGs_M1.csv",row.names = 1)
DEGs_M2l1 <- read.csv("figure_output/all/manuscript/GO_mø/DEGs_raw/DEGs_M2like1.csv",row.names = 1)
DEGs_M2l2 <- read.csv("figure_output/all/manuscript/GO_mø/DEGs_raw/DEGs_M2like2.csv",row.names = 1)

# Analysis on Intersected Genes ----
#### double check: any intersection of upregulated genes across each cell type ####
temp1 <- DEGs_M1$gene_SYMBOL[DEGs_M1$regulation=="up"]
temp2 <- DEGs_M2l1$gene_SYMBOL[DEGs_M2l1$regulation=="up"]
temp3 <- DEGs_M2l2$gene_SYMBOL[DEGs_M2l2$regulation=="up"]

cat("Intersected Upregulated Genes among M1 and M2l1 Møs Are: \n",  intersect(temp1, temp2),'\n')
cat("Intersected Upregulated Genes among M2l1 and M2l2 Møs Are: \n", intersect(temp2, temp3), '\n')
cat("Intersected Upregulated Genes among M1 and M2l2 Møs Are: \n", intersect(temp1, temp3), '\n')

#### double check: see whether the intersected genes collectively lead to any enriched GO ####
pAdjustMethod <- "BH"
qvalueCutoff <- 0.05

# enrichedGO of Intersected Upregulated Genes among M2l1 and M2l2 Møs
ego <- enrichGO(gene    = intersect(temp2, temp3),
                  OrgDb        = org.Mm.eg.db, # Replace with appropriate organism database
                  keyType      = 'SYMBOL', # Ensure this matches your gene ID type 'ENTREZID'
                  ont          = "BP", # Can be "BP", "CC", or "MF"
                  pAdjustMethod = pAdjustMethod,
                  qvalueCutoff = qvalueCutoff,
                  readable = TRUE)

dotplot(ego, showCategory = 15)
# use simplify() to remove redundancy of enriched GO terms
ego.simplified <- simplify(ego, cutoff = 0.6)
dotplot(ego.simplified)
write.csv(ego, "figure_output/all/manuscript/GO_mø/GO_M2l1_M2l2_intersect.csv", row.names = TRUE)
write.csv(ego.simplified, "figure_output/all/manuscript/GO_mø/GO_M2l1_M2l2_intersect_simplified.csv", row.names = TRUE)

# enrichedGO of Intersected Upregulated Genes among M1 and M2l2 Møs
ego2 <- enrichGO(gene    = intersect(temp1, temp3),
                  OrgDb        = org.Mm.eg.db, # Replace with appropriate organism database
                  keyType      = 'SYMBOL', # Ensure this matches your gene ID type 'ENTREZID'
                  ont          = "BP", # Can be "BP", "CC", or "MF"
                  pAdjustMethod = pAdjustMethod,
                  qvalueCutoff = qvalueCutoff,
                  readable = TRUE)

# use simplify() to remove redundancy of enriched GO terms
ego2.simplified <- simplify(ego2, cutoff = 0.6)
write.csv(ego2, "figure_output/all/manuscript/GO_mø/GO_M1_M2l2_intersect.csv", row.names = TRUE)
write.csv(ego2.simplified, "figure_output/all/manuscript/GO_mø/GO_M1_M2l2_intersect_simplified.csv", row.names = TRUE)

# Histogram too visualize count in M2l1 and M2l2 intersected gene GO
# Note: didn't do the same for M1 and M2l2 intersected GO because result can be
#       manually interpretted
ego.simplified <- as.data.frame(ego.simplified)
plot <- ggplot(ego.simplified, aes(x = Count)) + geom_histogram(binwidth=1, fill="#69b3a2", color = "black") + 
        hrbrthemes::theme_ipsum()
plot
# Remark: From the plot above, I have chosen to focus on pathways with 2 or more counts

datatable(ego.simplified, caption = 'Table 1: Intersecting DEG GO for M2l1 and M2l2 Mø')
ego.simplified.filtered <- ego.simplified %>% filter(Count>=2, p.adjust <=0.01)
# Remark: Now we have around 10 GOs representing intersected genes between M2l1 and M2l2

# Extract GO ID from ego.simplified.filtered and ego2.simplified[1]
GO_ID.intersect <- c(ego.simplified.filtered$ID, ego2.simplified$ID[1])
ego.M1.simplified <- simplify(ego_M1.up.BP$ego.obj, cutoff=0.6)
ego.M2l1.simplified <- simplify(ego_M2l1.up.BP$ego.obj, cutoff=0.6)
ego.M2l2.simplified <- simplify(ego_M2l2.up.BP$ego.obj, cutoff=0.6)

GO_ID.intersect %in% ego_M1.up.BP$ego.obj$ID
GO_ID.intersect %in% ego_M2l1.up.BP$ego.obj$ID
GO_ID.intersect %in% ego.M1.simplified$ID
GO_ID.intersect %in% ego.M2l1.simplified$ID

(GO_ID.intersect %in% ego_M2l1.up.BP$ego.obj$ID) | (GO_ID.intersect %in% ego_M1.up.BP$ego.obj$ID)
(GO_ID.intersect %in% ego.M1.simplified$ID) | (GO_ID.intersect %in% ego.M2l1.simplified$ID)
GO_ID.intersect

# Remark: I have decided to move forward with none simplified version of mø ego
#         and have manually selected overlapping GO to represent
GO_ID.intersect <- GO_ID.intersect[c(T,T,F,F,F,T,T,T,T,F,T,T)]

# Create a dataframe of compiliing GO pathway of interest from all mø ego dataframe
# Need features: ID, Description, GeneRatio per macrophage, p.adjus per macrophage
df.M1 <- as.data.frame(ego_M1.up.BP$ego.obj)
df.M2l1 <- as.data.frame(ego_M2l1.up.BP$ego.obj)
df.M2l2 <- as.data.frame(ego_M2l2.up.BP$ego.obj)

GO_ID <- ego.simplified[ego.simplified$ID %in% GO_ID.intersect, 1:2]
df.M1.GO <- df.M1[df.M1$ID %in% GO_ID.intersect & df.M1$p.adjust<0.01, c(1:3, 6)]
df.M2l1.GO <- df.M2l1[df.M2l1$ID %in% GO_ID.intersect & df.M2l1$p.adjust<0.01, c(1:3, 6)]
df.M2l2.GO <- df.M2l2[df.M2l2$ID %in% GO_ID.intersect & df.M2l2$p.adjust<0.01, c(1:3, 6)]

df.M1.GO$Phenotype <- "M1 like"
df.M2l1.GO$Phenotype <- "M2 like 1"
df.M2l2.GO$Phenotype <- "M2 like 2"

df <- rbind(df.M1.GO, df.M2l1.GO, df.M2l2.GO) %>% arrange(desc(GeneRatio))
df$GeneRatio <- parse_ratio(df$GeneRatio)
# Plot enrichment dotplot using ggplot()
gg <- ggplot(df, aes(x= Phenotype, y=Description, size=GeneRatio, color=p.adjust, group = Phenotype)) +
      geom_point(alpha = 0.8) + theme_dose(20) + ylab(NULL) + xlab(NULL) +
      scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab") + 
      scale_size(range = c(2, 8)) 
gg


#### Read and analyze from manually filtered GO file
df.mø.GO <- read.csv("figure_output/all/manuscript/GO_mø/GO_raw/Macrophage_GO_updated.csv")
df.mø.GO$GeneRatio <- parse_ratio(df.mø.GO$GeneRatio)

# Plot individual phenotype/population with facet_drig()
idx <- order(df.mø.GO$GeneRatio, decreasing = TRUE)
df.mø.GO$Description <- factor(df.mø.GO$Description,
                         levels=rev(unique(df.mø.GO$Description[idx])))

gg <- ggplot(df.mø.GO, aes(x= GeneRatio, y=Description, size=GeneRatio, color=p.adjust)) +
      geom_point(alpha = 0.8) + theme_dose(12) + ylab(NULL) + xlab(NULL) +
      scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab") + 
      scale_size(range = c(2, 8)) + facet_grid(Population~., scale='free')
gg

# Plot phenotype/population together in same grid
# temp <- df.mø.GO %>% arrange(is.unique, Population) %>% mutate(idx = 1:n())

# temp <- df.mø.GO %>% group_by(Population, is.unique) %>%
#                      mutate(group_id = cur_group_id()) %>% ungroup() 
                     

# temp1 <- df.mø.GO %>% filter(is.unique == "TRUE") %>% group_by(Population) %>%
#                      mutate(group_id = cur_group_id()) %>% ungroup()
# temp2 <- df.mø.GO %>% filter(is.unique == "FALSE") %>% group_by(Population) %>%
#                      mutate(group_id = cur_group_id()+3) %>% ungroup()
# df.mø.GO <- as.data.frame(rbind(temp1, temp2))
# df.mø.GO$idx <- 1:nrow(df)

# # df.mø.GO <- df.mø.GO[order(df.mø.GO$idx),]

# df.mø.GO$Description <- factor(df.mø.GO$Description,
#                          levels=rev(df.mø.GO$Description[df$idx]))

# unique(df.mø.GO$Description[idx])
df.mø.GO$is.unique <- !(df.mø.GO$Description %in% 
                        df.mø.GO$Description[duplicated(df.mø.GO$Description)])
idx <- order(-df.mø.GO$is.unique, df.mø.GO$Population, -df.mø.GO$GeneRatio, df.mø.GO$Population)

df.mø.GO$Description <- factor(df.mø.GO$Description,
                         levels=rev(unique(df.mø.GO$Description[idx])))


gg <- ggplot(df.mø.GO, aes(x= Population, y=Description, size=GeneRatio, color=p.adjust, group = Population)) +
      geom_point(alpha = 0.8) + theme_dose(16) + ylab(NULL) + xlab(NULL) +
      scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab") + 
      scale_size(range = c(2, 8)) 
gg
# ggsave("figure_output/all/manuscript/GO_mø/GO_Mø_up_updated.jpeg", plot = gg,  width = 10, height = 10, unit = "in")

