# This script is used for mø-msc crosstalk
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(reticulate)
library(Matrix)
library(dplyr)
library(Seurat)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(viridis)
library(BiocManager)
library(gage)
library(data.table)
library(ggforce)
library(CellChat)
library(patchwork)
library(gplots)
library(ggpubr)
library(httpgd)
library(languageserver)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)
source("analysis_modified.R")

# Setup Directory
setwd("~/Desktop/School Work/Stanford/Yang Lab/2023 Spring rotation_Heena")

data_dir <- "CS scRNAseq data share/CS0_50_RAW data/"

# cellchat <- readRDS("cellchat_Gel100vsGEL50_Combined.rds")
cellchat.100 <- readRDS("CS scRNAseq data share/CS0_50_RAW data/cellchat_Gel100.rds")
cellchat.50 <- readRDS("CS scRNAseq data share/CS0_50_RAW data/cellchat_Gel50.rds")

object.list <- list(GEL100 = cellchat.100, GEL50 = cellchat.50)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Macrophage - MSC crosstalk----
msc <- c("Stromal cells", "Fibroblasts")
macrophages <- c("M1-like", "M2-like 1", "M2-like 2")

tol <- 0.05
cutoff.pvalue <- 0.05
color.use <-  ggPalette(2)
color.use <- rev(color.use)

# Compare the macrophage information flow to msc+fibroblast+mø
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = macrophages, targets.use = c(msc, macrophages))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = macrophages, targets.use = c(msc, macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="Macrophages as Senders and MSCs+MØ as Reciever Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure 
ggsave("figure_output/mø_msc/approach1/from_mø/relativePathwayInfoFlow_møAsSenders_msc+MøAsReciever.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.1 <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                  (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)

df <- df[pathway.select.1, ]

pathway.sig.1 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
pathway.sig <- pathway.sig.1

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]


png("figure_output/mø_msc/approach1/from_mø/signalingRole_specificPathways.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_mø/signalingRole_specificPathways_outgoing.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_mø/signalingRole_specificPathways_incoming.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_mø/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach1/from_mø/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig.1))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach1/from_mø/diff_sign_pathways_score_by_group.csv", row.names = FALSE)







# Compare the msc+fibroblast information flow to msc+fibroblast+mø-----
gg1 <- rankNet_check(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = msc, targets.use = c(msc, macrophages))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = msc, targets.use = c(msc, macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="MSC as Senders and MSCs+MØ as Reciever Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure
ggsave("figure_output/mø_msc/approach1/from_msc/relativePathwayInfoFlow_mscAsSenders_msc+MøAsReciever.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.2 <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                  (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)
df <- df[pathway.select.2, ]
pathway.sig.2 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
pathway.sig <- pathway.sig.2

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]

png("figure_output/mø_msc/approach1/from_msc/signalingRole_specificPathways.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_msc/signalingRole_specificPathways_outgoing.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_msc/signalingRole_specificPathways_incoming.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/from_msc/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach1/from_msc/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach1/from_msc/diff_sign_pathways_score_by_group.csv", row.names = FALSE)







# Compare the msc+mø internal information flow------
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = c(msc, macrophages), targets.use = c(msc, macrophages))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = c(msc, macrophages), targets.use = c(msc, macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="MSCs+MØ Internal Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure
ggsave("figure_output/mø_msc/approach1/internal/relativePathwayInfoFlow_msc+MøInternal.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.3 <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                  (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)
df <- df[pathway.select.3, ]
pathway.sig.3 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
pathway.sig <- pathway.sig.3

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]

png("figure_output/mø_msc/approach1/internal/signalingRole_specificPathways.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/internal/signalingRole_specificPathways_outgoing.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/internal/signalingRole_specificPathways_incoming.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach1/internal/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach1/internal/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach1/internal/diff_sign_pathways_score_by_group.csv", row.names = FALSE)






# Extract information on genes involved in each pathway --- pathway.sig.3 = union(pathway.sig.2, pathway.sig.1)
genes_per_pathway <- data.frame(matrix(ncol = 2, nrow = length(pathway.sig.3)))
colnames(genes_per_pathway) <- c('pathway_name','genes_involved')
genes_per_pathway$pathway_name <- pathway.sig.3
genes_per_pathway <- genes_per_pathway %>%
                     group_by(pathway_name) %>%
                     mutate(genes_involved = toString(extractEnrichedLR(cellchat, signaling = pathway_name, geneLR.return = TRUE, enriched.only = FALSE)$geneLR))

write.csv(genes_per_pathway, "figure_output/mø_msc/approach1/diff_sign_pathways_genes.csv", row.names = FALSE)







source("classify_modified.R")
source("analysis_modified.R")
# Approach 2+3: Subset cellchat to only mø-msc
df <-subsetCommunication(cellchat, sources.use = c(msc, macrophages), targets.use = c(msc, macrophages), slot.name = 'netP')
write.csv(df[[1]], "figure_output/mø_msc/subsetCommunication_GEL100.csv", row.names = FALSE)
write.csv(df[[2]], "figure_output/mø_msc/subsetCommunication_GEL50.csv", row.names = FALSE)


# Approach 2:
cellchat100.mmsc <- subsetCellChat_mod(cellchat.100, idents.use = c(macrophages, msc))
cellchat50.mmsc <- subsetCellChat_mod(cellchat.50, idents.use = c(macrophages, msc))

cellchat100.mmsc <- netAnalysis_computeCentrality(cellchat100.mmsc)
cellchat50.mmsc <- netAnalysis_computeCentrality(cellchat50.mmsc)

# Approach 3:
cellchat100.mmsc <- readRDS("CS scRNAseq data share/CS0_50_RAW data/cellchat_Gel100_mmsc5.rds")
cellchat50.mmsc <- readRDS("CS scRNAseq data share/CS0_50_RAW data/cellchat_Gel50_mmsc5.rds")

cellchat100.mmsc <- updateCellChat(cellchat100.mmsc)
cellchat50.mmsc  <- updateCellChat(cellchat50.mmsc)

object.list <- list(GEL100 = cellchat100.mmsc, GEL50 = cellchat50.mmsc)
cellchat.mmsc <- mergeCellChat(object.list, add.names = names(object.list))

# netVisual plots for MSC_Mø
gg1 <- netVisual_heatmap(cellchat.mmsc, font.size = 15, font.size.title = 15)
gg2 <- netVisual_heatmap(cellchat.mmsc, measure = "weight", font.size = 15, font.size.title = 15)
figure <- gg1 + gg2
figure

weight.max.count  <- getMaxWeight(object.list, attribute = c("idents", "count"))
weight.max.weight <- getMaxWeight(object.list, attribute = c("idents", "weight"))
par(mfrow = c(2, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.count[2], edge.width.max = 3, vertex.weight = 3,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = weight.max.weight[2], edge.width.max = 3, vertex.weight = 3,
                   title.name = paste0("Interaction Weights/Strength - ", names(object.list)[i]))
}

# Macrophage as Senders to MSC_Mø---
# Macrophage - MSC crosstalk----
# msc <- c("MSC")
msc <- c("Stromal cells", "Fibroblasts")
macrophages <- c("M1-like", "M2-like 1", "M2-like 2")

tol <- 0.05
cutoff.pvalue <- 0.05
color.use <-  ggPalette(2)
color.use <- rev(color.use)

gg1 <- rankNet_modified(cellchat.mmsc , mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = macrophages, targets.use = c(msc, macrophages))
gg2 <- rankNet_modified(cellchat.mmsc , mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = macrophages, targets.use = c(msc, macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="Macrophages as Senders and MSC+Møs as Reciever Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure 
ggsave("figure_output/mø_msc/approach4/from_mø/relativePathwayInfoFlow_møAsSenders_msc+møAsReciever.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.1 <- ((df$contribution.relative.1 < 1 - tol) & (df$pvalues < cutoff.pvalue)) |
                  ((df$contribution.relative.1 > 1 + tol) & (df$pvalues < cutoff.pvalue))

df <- df[pathway.select.1, ]

pathway.sig.3.1 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
# pathway.sig.3.1[length(pathway.sig.3.1)+1] <- "JAM"
pathway.sig <- pathway.sig.3.1

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]


png("figure_output/mø_msc/approach4/from_mø/signalingRole_specificPathways.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_mø/signalingRole_specificPathways_outgoing.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_mø/signalingRole_specificPathways_incoming.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_mø/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach4/from_mø/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig.1))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach4/from_mø/diff_sign_pathways_score_by_group.csv", row.names = FALSE)




# MSC as Senders to MSC_Mø---
gg1 <- rankNet_modified(cellchat.mmsc, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = msc, targets.use = c(macrophages))
gg2 <- rankNet_modified(cellchat.mmsc, mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = msc, targets.use = c(macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="MSC as Senders and Møs as Reciever Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure
ggsave("figure_output/mø_msc/approach4/from_msc/relativePathwayInfoFlow_mscAsSenders_møAsReciever.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.2 <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                  (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)
df <- df[pathway.select.2, ]
pathway.sig.3.2 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
pathway.sig <- pathway.sig.3.2

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]

png("figure_output/mø_msc/approach4/from_msc/signalingRole_specificPathways.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 10, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_msc/signalingRole_specificPathways_outgoing.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_msc/signalingRole_specificPathways_incoming.png",width=8,height=8,units="in",res=1200)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 10, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/from_msc/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 10, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach4/from_msc/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach4/from_msc/diff_sign_pathways_score_by_group.csv", row.names = FALSE)




# Compare the msc+mø internal information flow------
gg1 <- rankNet_modified(cellchat.mmsc, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE, sources.use = c(msc, macrophages), targets.use = c(msc, macrophages))
gg2 <- rankNet_modified(cellchat.mmsc, mode = "comparison", stacked = FALSE, do.stat = TRUE, sources.use = c(msc, macrophages), targets.use = c(msc, macrophages))
figure <- (gg1$gg.obj + gg2) + plot_annotation(title="MSCs+MØ Internal Information Flow") & theme(plot.title = element_text(hjust = 0.5, size = 20))
figure
ggsave("figure_output/mø_msc/approach4/internal/relativePathwayInfoFlow_msc+MøInternal.jpeg", plot = figure, width = 22, height = 13, unit = "in")

df <- gg1$signaling.contribution
pathway.select.3 <- (df$contribution.relative < 1 - tol) & (df$pvalues < cutoff.pvalue) |
                  (df$contribution.relative > 1 + tol) & (df$pvalues < cutoff.pvalue)
df <- df[pathway.select.3, ]
pathway.sig.3.3 <- unique(df$name) # each pathway occurs twice for evaluation of each condition
pathway.sig <- pathway.sig.3.3

colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2],
                       ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
colors.text <- colors.text[1:length(pathway.sig)]

png("figure_output/mø_msc/approach4/internal/signalingRole_specificPathways.png",width=8,height=9,units="in",res=1600)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text,
                                         width = 5, height = 16, color.heatmap = "OrRd")

ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "all",
                                         signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text,
                                         width = 5, height = 16, color.heatmap = "OrRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/internal/signalingRole_specificPathways_outgoing.png",width=8,height=9,units="in",res=1600)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 16)
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "outgoing", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 16)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/internal/signalingRole_specificPathways_incoming.png",width=8,height=9,units="in",res=1600)
ht1 <- netAnalysis_signalingRole_heatmap_modified(object.list[[1]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[1], color.rownames = colors.text, width = 5, height = 16, color.heatmap = "PuRd")
ht2 <- netAnalysis_signalingRole_heatmap_modified(object.list[[2]], pattern = "incoming", signaling = pathway.sig, title = names(object.list)[2], color.rownames = colors.text, width = 5, height = 16, color.heatmap = "PuRd")
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png("figure_output/mø_msc/approach4/internal/differential_signalingRole_specificPathways.png",width=12,height=15,units="in",res=1200)
ht1 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("outgoing"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 16, return.data = FALSE)
ht2 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("incoming"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 16, return.data = FALSE)
ht3 <- netAnalysis_diff_signalingRole_heatmap(cellchat.mmsc, signaling = pathway.sig, comparison = c(1,2), pattern = c("all"), color.heatmap = "RdBu", color.rownames = colors.text, width = 5, height = 16, return.data = FALSE)
ComplexHeatmap::draw(ht1 + ht2 + ht3, ht_gap = unit(2, "cm"))
dev.off()

# Author Heena: Calculate the score of each pathway in each material composition
# and determine which material composition has higher sore
df <- df %>% group_by(name) %>% mutate(total.contribution = sum(contribution)) %>% 
             ungroup() %>% mutate(score = contribution / total.contribution * contribution.scaled) %>%
             group_by(name) %>% mutate(is.max = ifelse(score == max(score), TRUE, FALSE))

# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
# Extract subset of dataframe for pathways (`name`) with material composition `group`
# that has higher score (`is.max` = TRUE)
df.pathways.comp.max <- df[df$is.max == TRUE,] %>%
                   dplyr::select(-c(total.contribution, is.max))
df.pathways.comp.max <- as.data.frame(df.pathways.comp.max)

write.csv(df.pathways.comp.max, "figure_output/mø_msc/approach4/internal/diff_sign_pathways_with_score.csv", row.names = FALSE)

df.pathways.sig <- head(arrange(df.pathways.comp.max, dplyr::desc(score)), n = length(pathway.sig))

# Create a new `init` dataframe extracting info of only pathway, score of each group
init <- data.frame(matrix(ncol = 3, nrow = length(pathway.sig)))
colnames(init) <- c('pathway_name','score_GEL100', 'score_GEL50')
init$pathway_name <- pathway.sig
init$score_GEL100 <- df$score[df$group == "GEL100"]
init$score_GEL50 <- df$score[df$group == "GEL50"]
init <- init %>% mutate(dominance = ifelse(score_GEL100 > score_GEL50, "GEL100", "GEL50"))
write.csv(init, "figure_output/mø_msc/approach4/internal/diff_sign_pathways_score_by_group.csv", row.names = FALSE)


# Extract information on genes involved in each pathway --- pathway.sig.3 = union(pathway.sig.2, pathway.sig.1)
genes_per_pathway <- data.frame(matrix(ncol = 2, nrow = length(pathway.sig.3.3)))
colnames(genes_per_pathway) <- c('pathway_name','genes_involved')
genes_per_pathway$pathway_name <- pathway.sig.3.3
genes_per_pathway <- genes_per_pathway %>%
                     group_by(pathway_name) %>%
                     mutate(genes_involved = toString(extractEnrichedLR(cellchat.mmsc, signaling = pathway.sig.3.3, geneLR.return = TRUE, enriched.only = FALSE)$geneLR))

write.csv(genes_per_pathway, "figure_output/mø_msc/approach4/diff_sign_pathways_genes.csv", row.names = FALSE)



# Hierarchy Plot Example
gg <- list()
pathways.show <- c("OSM")
vertex.receiver <- seq(4,5)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  gg[[i]] <- netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, vertex.weight = 0.5,
  layout = "hierarchy", edge.weight.max = weight.max[1], arrow.size = 1, signaling.name = paste(pathways.show, names(object.list)[i]))
}





source("classify_modified.R")
source("analysis_modified.R")
source("visualization_modified.R")




ggarrange(gg[[1]],gg[[2]], ncol = 1, nrow = 2, labels = c("A","B"))

ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
  color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

for (i in 1:length(object.list)) {
  ht[[i]] <- netAnalysis_signalingRole_network(object.list[[i]], signaling = pathways.show, 
  color.heatmap = "Reds")
}


## circle plot
pathways.show <- c("OSM") 
par(mfrow=c(1,1))
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle")
