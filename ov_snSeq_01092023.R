# ------ 1: preprocess the count matrix from 10x Cellranger and MTD pipelines ------ #

### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)

# change the current plan to access parallelization
plan("multisession", workers = 4)

setwd("~/ovary_SC")
# import 10X results
cor_ct <- Read10X(data.dir = "/bigrock/FurmanLab/Fei/ovary_SC/ov-cor-ct/outs/filtered_feature_bc_matrix/")
cor_dox <- Read10X(data.dir = "/bigrock/FurmanLab/Fei/ovary_SC/ov-cor-dox/outs/filtered_feature_bc_matrix/")
med_ct <- Read10X(data.dir = "/bigrock/FurmanLab/Fei/ovary_SC/ov-med-ct/outs/filtered_feature_bc_matrix/")
med_dox <- Read10X(data.dir = "/bigrock/FurmanLab/Fei/ovary_SC/ov-med-dox/outs/filtered_feature_bc_matrix/")

# make a sample list
Sample.list <- list("cor_ct"=cor_ct,
                    "cor_dox"=cor_dox,
                    "med_ct"=med_ct,
                    "med_dox"=med_dox)
# transfer to dataframe for reprocessing
Sample.list <- lapply(Sample.list, as.data.frame)

# create seuratobject for QC
Sample.list.qc<-list()
for (i in names(Sample.list)){
  Sample.list.qc[[i]] <- CreateSeuratObject(Sample.list[[i]], project=i,min.cells = 3, min.feature = 250)
}
# function for Pre-process Seurat object: QC, PCA and UMAP
SCT_m<-function(l){
  l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
  l <- subset(l, subset= nFeature_RNA>250 & percent.mt < 10)
  l <- SCTransform(l, vst.flavor = "v2", method = "glmGamPoi",verbose = FALSE)
  l <- RunPCA(l, verbose = FALSE)
  l <- RunUMAP(l,dims = 1:30, verbose = FALSE)
  l <- FindNeighbors(l, dims = 1:30, verbose = FALSE)
  l <- FindClusters(l, verbose = FALSE)
}

# SCTransform for each dataset independently
Sample.list.qc<-lapply(Sample.list.qc,SCT_m)

# Remove doublet
library(DoubletFinder)
for (i in names(Sample.list.qc)){
  # pK Identification (no ground-truth)
  sweep.res.list_sample <- paramSweep_v3(Sample.list.qc[[i]], PCs = 1:30, sct = T)
  sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  pK<-as.numeric(as.character(bcmvn_sample$pK))[bcmvn_sample$BCmetric==max(bcmvn_sample$BCmetric)]
  ## Homotypic Doublet Proportion Estimate
  annotations <- Sample.list.qc[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Sample.list.qc[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  re.pAnn<-names(Sample.list.qc[[i]]@meta.data)[(length(Sample.list.qc[[i]]@meta.data)-1)]
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = re.pAnn, sct = T)
}

# filter samples based on results of QC & Doublet
for (i in names(Sample.list.qc)){
  sample.meta.data<-Sample.list.qc[[i]]@meta.data
  singlet<-row.names(sample.meta.data)[sample.meta.data[length(sample.meta.data)]=="Singlet"]
  Sample.list[[i]]<-Sample.list[[i]][singlet]
}
# create seuratobject for integration
for (i in names(Sample.list)){
  Sample.list[[i]] <- CreateSeuratObject(Sample.list[[i]], project=i)
  Sample.list[[i]] <- SCTransform(Sample.list[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}
features <- SelectIntegrationFeatures(object.list = Sample.list, nfeatures = 3000)
# run PCA before RPCA integration
Sample.list <- PrepSCTIntegration(object.list = Sample.list, anchor.features = features)
Sample.list <- lapply(X = Sample.list, FUN = RunPCA, features = features)

# Perform RPCA integration
plan("multisession", workers = 16)
options(future.globals.maxSize= +Inf) # increase maxSize of future function

Anchors <- FindIntegrationAnchors(object.list = Sample.list, normalization.method = "SCT",
                                  anchor.features = features, reduction = "rpca")
Sample.combined <- IntegrateData(anchorset = Anchors, normalization.method = "SCT")

# Dimensional reduction
Sample.combined <- RunPCA(Sample.combined, npcs = 20 ,verbose = FALSE)
# ElbowPlot(Sample.combined,ndims = 30)
Sample.combined <- RunUMAP(Sample.combined, reduction = "pca", dims = 1:20)
# Cluster the cells
Sample.combined <- FindNeighbors(Sample.combined, dims = 1:20)
Sample.combined <- FindClusters(Sample.combined, resolution = 0.3)
DimPlot(Sample.combined, reduction = "umap",split.by = "orig.ident",label =T)
ggsave("ov_clusters.pdf",width=12,height=6)

# Prepare object to run differential expression on SCT assay with multiple models
Sample.combined <- PrepSCTFindMarkers(Sample.combined)

# add more metadata
Sample.combined$tissue <- case_when(Sample.combined$orig.ident ==
                                             "cor_ct" ~ "cortex",
                                    Sample.combined$orig.ident ==
                                             "cor_dox" ~ "cortex",
                                    Sample.combined$orig.ident ==
                                      "med_ct" ~ "medulla",
                                    Sample.combined$orig.ident ==
                                      "med_dox" ~ "medulla")

# find markers for every cluster compared to all remaining cells
markers.all <- FindAllMarkers(Sample.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="tissue")
write.csv(markers.all,"markers_clusters.csv")

# plot top markers in heatmap
library(dplyr)
markers.all %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(Sample.combined, features=top10$gene)
ggsave("ov_markers_top10.pdf",width=8,height=16)

# find DEGs across conditions
Idents(Sample.combined) <- "orig.ident"
cor_dox.vs.ct <- FindMarkers(Sample.combined, assay = "SCT", ident.1 = "cor_dox", ident.2 = "cor_ct",
                             min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
write.csv(cor_dox.vs.ct,"DEGs_cor_dox_vs_ct.csv")
med_dox.vs.ct <- FindMarkers(Sample.combined, assay = "SCT", ident.1 = "med_dox", ident.2 = "med_ct",
                             min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
write.csv(med_dox.vs.ct,"DEGs_med_dox_vs_ct.csv")

library(SingleR)
HPCA.ref <- celldex::HumanPrimaryCellAtlasData()

# set the active assay back to “RNA,” and re-do the normalization
DefaultAssay(Sample.combined) <- "RNA"
Sample.combined <- NormalizeData(Sample.combined)
# convert Seurat object to single cell experiment (SCE) for convenience
sce <- as.SingleCellExperiment(DietSeurat(Sample.combined))
HPCA.main <- SingleR(test = sce,assay.type.test = 1,ref = HPCA.ref,labels = HPCA.ref$label.main)
# see the summary of general cell type annotations
table(HPCA.main$pruned.labels)
# add the annotations to the Seurat object metadata
Sample.combined@meta.data$HPCA.main <- HPCA.main$pruned.labels

# most frequent predicted cell type for each cluster
Sample.combined$predicted.cellType <- as.character(Sample.combined$seurat_clusters)
for (c in unique(Sample.combined$seurat_clusters)){
  cell_type <- Sample.combined@meta.data %>%
    filter(seurat_clusters==c) %>%
    select(HPCA.main)
  cell_type_counts<-table(cell_type)
  most_frequent_cell_type<-names(which.max(cell_type_counts))
  Sample.combined$predicted.cellType[Sample.combined$predicted.cellType==c]<-most_frequent_cell_type
}
Idents(Sample.combined) <- "predicted.cellType"
DimPlot(Sample.combined, reduction = "umap",split.by = "orig.ident",label =T)
ggsave("ov_cellTypes.pdf",width=12,height=6)

# quick check cell annotation results for a cluster
table(Sample.combined@meta.data %>%
  filter(seurat_clusters==8) %>%
  select(HPCA.main))

# refine cell annotation results by scRNAseq ovary literatures collected by researcher
Genes_list_human_ovary <- readxl::read_excel("Genes_list_human_ovary.xlsx")
cell4type <- unique(Genes_list_human_ovary$cluster)
cluster4type <- unique(markers.all$cluster)
# correlation test
correlation.df<-data.frame()
for (c in cell4type){
  cell_FC_ref <- Genes_list_human_ovary %>% filter(cluster==c) %>% select(c(gene,avg_logFC,cluster))
  for (cl in cluster4type){
    cell_FC_target <- markers.all %>% filter(cluster==cl) %>% select(c(gene,avg_log2FC,cluster))
    FC4cor <- inner_join(cell_FC_ref,cell_FC_target,by="gene")
    correlation <- cor(FC4cor$avg_log2FC,FC4cor$avg_logFC)
    cor.df <- data.frame(cluster=cl,cell_type=c,correlation=correlation,"overlaps_ref_genes"=paste0(nrow(FC4cor),"/",nrow(cell_FC_ref)))
    correlation.df <- rbind(correlation.df,cor.df)
  }
}
write.csv(correlation.df,"correlation_celltype_annotation.csv")

# cell type enrichment test
for (c in cell4type){
  cell_FC_ref_pos <- Genes_list_human_ovary %>% filter(cluster==c & avg_logFC>0) %>% select(gene)
  celltype_features <- list(as.character(cell_FC_ref_pos$gene))
  ov_module <- AddModuleScore(
    object = ov_module,
    features = celltype_features,
    name = c)
}

celltype.df <- ov_module@meta.data %>% select(c("seurat_clusters","endothelial.cells1","stromal.cells1",
                                 "perivascular.cells1","immune.cells1","smooth.muscle.cells1","theca...stroma1"))
box_ls<-list()
for (cs in names(celltype.df)[-1]){
  box_plot_cs <- ggplot(celltype.df, aes_string(x = "seurat_clusters", y = cs)) + 
    geom_boxplot() + 
    xlab("clusters") + 
    ylab(paste0("cell scores of ",cs)) + 
    ggtitle(paste0("scores of ",cs))
  box_ls[[cs]] <- box_plot_cs
}
library(gridExtra)
grid_plot<-grid.arrange(grobs = box_ls, ncol = 2)
ggsave(paste0("Box Plot of celltype scores.pdf"),plot = grid_plot,width=12,height=11)

# annotate cell type by combining results from HPCA and manually curate the gene list from ovary publications
Sample.combined$celltype_refine <- case_when(Sample.combined$seurat_clusters ==
                                          "0" ~ "Stromal_cells",
                                        Sample.combined$seurat_clusters ==
                                          "1" ~ "Stromal_cells",
                                        Sample.combined$seurat_clusters ==
                                          "2" ~ "Smooth_muscle_cells",
                                        Sample.combined$seurat_clusters ==
                                          "3" ~ "Epithelial_cells",
                                        Sample.combined$seurat_clusters ==
                                          "4" ~ "Perivascular_cells",
                                        Sample.combined$seurat_clusters ==
                                          "5" ~ "Perivascular_cells",
                                        Sample.combined$seurat_clusters ==
                                          "6" ~ "Perivascular_cells",
                                        Sample.combined$seurat_clusters ==
                                          "7" ~ "Endothelial_cells",
                                        Sample.combined$seurat_clusters ==
                                          "8" ~ "MSC",
                                        Sample.combined$seurat_clusters ==
                                          "9" ~ "Tissue_stem_cells",
                                        Sample.combined$seurat_clusters ==
                                          "10" ~ "Epithelial_cells")

# set the levels in order
Sample.combined$celltype_refine <- factor(Sample.combined$celltype_refine,
                                          levels=sort(unique(Sample.combined$celltype_refine)))

# plot UMAP based on refined cell types
Idents(Sample.combined) <- "celltype_refine"
DimPlot(Sample.combined, reduction = "umap",split.by = "orig.ident",label =T,repel =T)
ggsave("ov_cellTypes_refined.pdf",width=12,height=6)

# ------ find DEG for all cell types after Dox treatment ------ #
# with default filter min.pct=0.1 and logfc=0.25
Sample.combined$celltype_refine.treatment <- paste(Sample.combined$celltype_refine, Sample.combined$orig.ident,
                                            sep = "_")
Idents(Sample.combined)<-"celltype_refine.treatment"
plan("multisession", workers = 4)
cell4DEG<-unique(Sample.combined$celltype_refine)
for (d in cell4DEG){
  DEG.cellType <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_cor_dox"), ident.2 =paste0(d,"_cor_ct"))
  write.csv(DEG.cellType, paste0(d,"_cor_DoxvsCt.csv"))
}
for (d in cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]){
  DEG.cellType <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_med_dox"), ident.2 =paste0(d,"_med_ct"))
  write.csv(DEG.cellType, paste0(d,"_med_DoxvsCt.csv"))
}

# ------ Cell type ratio change stack barplot ------ #
sub.prop.all<-data.frame()
for (l in unique(Sample.combined$orig.ident)){
  sub.treat<-Sample.combined@meta.data[Sample.combined$orig.ident==l,]
  sub.prop<-data.frame(table(sub.treat$celltype_refine)/sum(table(sub.treat$celltype_refine)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}
write.csv(sub.prop.all,"sub.prop.all.csv")

library(tidyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("cor_ct","cor_dox","med_ct","med_dox"))
sub.prop.all$labs<-round(sub.prop.all$Freq,3) # prepare cell proportion label
sub.prop.all$labs[sub.prop.all$labs<0.05]<-""
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="Ovary") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_all.pdf",width=6,height=5)


# normalize and scale RNA assays
DefaultAssay(Sample.combined) <- "RNA"
Sample.combined <- NormalizeData(Sample.combined)
Sample.combined <- ScaleData(Sample.combined, features = all.genes)

# load the senescence gene sets
Senescence_genesets <- read.csv("Senescence_genesets.csv")

sc_module <- Sample.combined

for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  geneset_features <- list(sc.set)
  # score each senescent gene set
  sc_module <- AddModuleScore(
    object = sc_module,
    features = geneset_features,
    name = names(Senescence_genesets[sc]))
}

## cor_Dox_vs_Control
library(tibble)
# sc_score <- data.frame()
# sc_score_c <- data.frame(matrix(ncol = 0, nrow = 7))
# row.names(sc_score_c) <- names(Senescence_genesets)
# sc_score_p <- data.frame()
# sc_score_p_c <- data.frame(matrix(ncol = 0, nrow = 7))
# row.names(sc_score_p_c) <- names(Senescence_genesets)
sc_score.ls<-list()
sc_score_p.ls<-list()
for (cell in cell4DEG){
  sc_score <- data.frame()
  sc_score_p <- data.frame()
  for (sc in 1:ncol(Senescence_genesets)){
    names(Senescence_genesets[sc])
    ctl <- sc_module@meta.data %>% 
      filter(celltype_refine.treatment==paste0(cell,"_cor_ct")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
    dox <- sc_module@meta.data %>% 
      filter(celltype_refine.treatment==paste0(cell,"_cor_dox")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
    compare_result <- t.test(dox,ctl)
    avg_dox <- mean(dox[,1])
    avg_ctl <- mean(ctl[,1])
    #log2fc <- log2(avg_dox/avg_ctl)
    p.value <- compare_result$p.value
    # sc_score_df <- data.frame(SC_gene_set=names(Senescence_genesets[sc]),log2fc,p.value,
    #                           avg_dox, avg_ctl,
    #                           comparison=cell)
    sc_score_df <- data.frame(names(Senescence_genesets[sc]),
                              avg_dox, avg_ctl)
    names(sc_score_df) <- c("SC_gene_set",paste0(cell,"_Dox"),paste0(cell,"_Ctl"))
    sc_score <- rbind(sc_score,sc_score_df)
    sc_score_df_p <- data.frame(names(Senescence_genesets[sc]),
                                p.value, p.value)
    names(sc_score_df_p) <- c("SC_gene_set",paste0(cell,"_Dox_p_value"),paste0(cell,"_Ctl_p_value"))
    sc_score_p <- rbind(sc_score_p,sc_score_df_p)
  }
  sc_score.ls[[cell]] <- sc_score
  sc_score_p.ls[[cell]] <- sc_score_p
  # sc_score_cr <- column_to_rownames(sc_score,"SC_gene_set")
  # sc_score_c <- cbind(sc_score_c,sc_score_cr)
  # sc_score_p_cr <- column_to_rownames(sc_score_p,"SC_gene_set")
  # sc_score_p_c <- cbind(sc_score_p_c,sc_score_p_cr)
}
# for Overall comparison
sc_score <- data.frame()
sc_score_p <- data.frame()
for (sc in 1:ncol(Senescence_genesets)){
  ctl <- sc_module@meta.data %>% 
    filter(orig.ident=="cor_ct") %>% 
    select(paste0(names(Senescence_genesets[sc]),"1"))
  dox <- sc_module@meta.data %>% 
    filter(orig.ident=="cor_dox") %>% 
    select(paste0(names(Senescence_genesets[sc]),"1"))
  compare_result <- t.test(dox,ctl)
  avg_dox <- mean(dox[,1])
  avg_ctl <- mean(ctl[,1])
  #log2fc <- log2(avg_dox/avg_ctl)
  p.value <- compare_result$p.value
  sc_score_df <- data.frame(names(Senescence_genesets[sc]),
                            avg_dox, avg_ctl)
  names(sc_score_df) <- c("SC_gene_set","Overall_Dox","Overall_Ctl")
  sc_score <- rbind(sc_score,sc_score_df)
  sc_score_df_p <- data.frame(names(Senescence_genesets[sc]),
                              p.value, p.value)
  names(sc_score_df_p) <- c("SC_gene_set","Overall_Dox","Overall_Ctl")
  sc_score_p <- rbind(sc_score_p,sc_score_df_p)
}
sc_score.ls[["Overall"]] <- sc_score
sc_score_p.ls[["Overall"]] <- sc_score_p
library(tidyverse)
sc_score <- sc_score.ls %>% reduce(left_join, by = "SC_gene_set")
sc_score_p <- sc_score_p.ls %>% reduce(left_join, by = "SC_gene_set")

sc_score[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score[,1])
sc_score_p[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score_p[,1])
library(tibble)
# score4hp <- sc_score %>% select(c("SC_gene_set","log2fc","comparison"))
# score4hp <- pivot_wider(score4hp, names_from = comparison, values_from = log2fc)
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
# score4hp_p <- sc_score %>% select(c("SC_gene_set","p.value","comparison"))
# score4hp_p <- pivot_wider(score4hp_p, names_from = comparison, values_from = p.value)
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_cor.pdf", width=9, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          at=c(-0.4,-0.2,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p[i, j] > 0.05) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(score4hp[i, j],3), x, y,gp=gpar(fontsize=6.5))
                }})
ht<-draw(ht)
dev.off()

# draw only Dox samples
col_fun = colorRamp2(c(-0.1, 0, 0.4), c("blue", "white", "red"))
pdf(file="SnC-score_dox_ov_cor.pdf", width=7, height = 5)
score4hp <- score4hp %>% as.data.frame() %>% select(grep("_Dox",colnames(score4hp),value=T))
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
              grid.text(round(score4hp[i, j],3), x, y,gp=gpar(fontsize=6))})
ht<-draw(ht)
dev.off()


# for medula
sc_score.ls<-list()
sc_score_p.ls<-list()
for (cell in cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]){
  sc_score <- data.frame()
  sc_score_p <- data.frame()
  for (sc in 1:ncol(Senescence_genesets)){
    names(Senescence_genesets[sc])
    ctl <- sc_module@meta.data %>% 
      filter(celltype_refine.treatment==paste0(cell,"_med_ct")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
    dox <- sc_module@meta.data %>% 
      filter(celltype_refine.treatment==paste0(cell,"_med_dox")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
    compare_result <- t.test(dox,ctl)
    avg_dox <- mean(dox[,1])
    avg_ctl <- mean(ctl[,1])
    #log2fc <- log2(avg_dox/avg_ctl)
    p.value <- compare_result$p.value
    # sc_score_df <- data.frame(SC_gene_set=names(Senescence_genesets[sc]),log2fc,p.value,
    #                           avg_dox, avg_ctl,
    #                           comparison=cell)
    sc_score_df <- data.frame(names(Senescence_genesets[sc]),
                              avg_dox, avg_ctl)
    names(sc_score_df) <- c("SC_gene_set",paste0(cell,"_Dox"),paste0(cell,"_Ctl"))
    sc_score <- rbind(sc_score,sc_score_df)
    sc_score_df_p <- data.frame(names(Senescence_genesets[sc]),
                                p.value, p.value)
    names(sc_score_df_p) <- c("SC_gene_set",paste0(cell,"_Dox_p_value"),paste0(cell,"_Ctl_p_value"))
    sc_score_p <- rbind(sc_score_p,sc_score_df_p)
  }
  sc_score.ls[[cell]] <- sc_score
  sc_score_p.ls[[cell]] <- sc_score_p
  # sc_score_cr <- column_to_rownames(sc_score,"SC_gene_set")
  # sc_score_c <- cbind(sc_score_c,sc_score_cr)
  # sc_score_p_cr <- column_to_rownames(sc_score_p,"SC_gene_set")
  # sc_score_p_c <- cbind(sc_score_p_c,sc_score_p_cr)
}
# for Overall comparison
sc_score <- data.frame()
sc_score_p <- data.frame()
for (sc in 1:ncol(Senescence_genesets)){
  ctl <- sc_module@meta.data %>% 
    filter(orig.ident=="med_ct") %>% 
    select(paste0(names(Senescence_genesets[sc]),"1"))
  dox <- sc_module@meta.data %>% 
    filter(orig.ident=="med_dox") %>% 
    select(paste0(names(Senescence_genesets[sc]),"1"))
  compare_result <- t.test(dox,ctl)
  avg_dox <- mean(dox[,1])
  avg_ctl <- mean(ctl[,1])
  #log2fc <- log2(avg_dox/avg_ctl)
  p.value <- compare_result$p.value
  sc_score_df <- data.frame(names(Senescence_genesets[sc]),
                            avg_dox, avg_ctl)
  names(sc_score_df) <- c("SC_gene_set","Overall_Dox","Overall_Ctl")
  sc_score <- rbind(sc_score,sc_score_df)
  sc_score_df_p <- data.frame(names(Senescence_genesets[sc]),
                              p.value, p.value)
  names(sc_score_df_p) <- c("SC_gene_set","Overall_Dox","Overall_Ctl")
  sc_score_p <- rbind(sc_score_p,sc_score_df_p)
}
sc_score.ls[["Overall"]] <- sc_score
sc_score_p.ls[["Overall"]] <- sc_score_p
library(tidyverse)
sc_score <- sc_score.ls %>% reduce(left_join, by = "SC_gene_set")
sc_score_p <- sc_score_p.ls %>% reduce(left_join, by = "SC_gene_set")

sc_score[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score[,1])
sc_score_p[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score_p[,1])
library(tibble)
# score4hp <- sc_score %>% select(c("SC_gene_set","log2fc","comparison"))
# score4hp <- pivot_wider(score4hp, names_from = comparison, values_from = log2fc)
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
# score4hp_p <- sc_score %>% select(c("SC_gene_set","p.value","comparison"))
# score4hp_p <- pivot_wider(score4hp_p, names_from = comparison, values_from = p.value)
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.4), c("blue", "white", "red"))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_med.pdf", width=9, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p[i, j] > 0.05) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(score4hp[i, j],3), x, y,gp=gpar(fontsize=8))
                }})
ht<-draw(ht)
dev.off()

# draw only Dox samples
pdf(file="SnC-score_dox_ov_med.pdf", width=7, height = 5)
score4hp <- score4hp %>% as.data.frame() %>% select(grep("_Dox",colnames(score4hp),value=T))
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                grid.text(round(score4hp[i, j],3), x, y,gp=gpar(fontsize=6))})
ht<-draw(ht)
dev.off()


# senescent gene lists tests
## SenMayo
geneset_features <- list(c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2",
                           "BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26",
                           "CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF1",
                           "CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12",
                           "CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR",
                           "EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG",
                           "HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3",
                           "IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18",
                           "IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2",
                           "ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13",
                           "MMP14","MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF",
                           "PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4",
                           "SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2",
                           "TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA",
                           "VEGFC","VGF","WNT16","WNT2"))
### SenMayo scores
ov_module <- AddModuleScore(
  object = Sample.combined,
  features = geneset_features,
  name = 'SenMayo')
ggplot(filter(ov_module@meta.data,tissue=="cortex"), aes(x=SenMayo1, fill=orig.ident)) +
  geom_density(alpha=0.4) + theme_bw()
ggsave("ov_cor_SenMayo.pdf", width = 6, height = 5)
ggplot(filter(ov_module@meta.data,tissue=="medulla"), aes(x=SenMayo1, fill=orig.ident)) +
  geom_density(alpha=0.4) + theme_bw()
ggsave("ov_med_SenMayo.pdf", width = 6, height = 5)

Sample.combined@scale.data

### fgsea for pathway enrichment by reranking gene list
# library(devtools)
# install_github("ctlab/fgsea")
library(fgsea)
# cortex
setwd("/opt/home/buckcenter.org/fwu/ovary_SC/cell_type_refine/Cortex")
file_list <- list.files(pattern="*.csv")
df_list <- list()
for (file in file_list) {
  df <- read.csv(file)
  cell4gsea <- gsub("_cor_DoxvsCt.csv","",file)
  df <- df[order(df$avg_log2FC, decreasing = T),]
  df <- df %>% filter(p_val_adj<0.05) %>% select(X,avg_log2FC)
  ranks <- deframe(df)
  df_list[[cell4gsea]] <- ranks
}
setwd("/opt/home/buckcenter.org/fwu/ovary_SC/")

# Load the pathways into a named list
pathways.SnC <- gmtPathways("Cellular_senescence_genesets.gmt")

# fgseaRes <- fgsea(pathways=pathways.SnC, stats=df_list[[1]])
fgseaRes <- lapply(df_list,function(x) fgsea(pathways=pathways.SnC, stats=x))

# Define a function to modify column names
modify_colnames <- function(df, name) {
  colnames(df)[-1] <- paste0(name, "_", colnames(df)[-1])
  return(df)
}

# Use lapply to apply the function to each data frame in the list
fgseaRes <- lapply(names(fgseaRes), function(x) modify_colnames(fgseaRes[[x]], x))

fgseaRes.df <- fgseaRes %>% reduce(left_join, by = "pathway")
fgseaRes.df <- column_to_rownames(fgseaRes.df,"pathway")
fgseaRes.NES <- fgseaRes.df %>% select(grep("_NES", colnames(fgseaRes.df), value = TRUE))
colnames(fgseaRes.NES) <- gsub("_NES","",colnames(fgseaRes.NES))
fgseaRes.p <- fgseaRes.df %>% select(grep("_padj", colnames(fgseaRes.df), value = TRUE))
colnames(fgseaRes.p) <- gsub("_padj","",colnames(fgseaRes.p))

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

fgseaRes.NES <- as.matrix(fgseaRes.NES)
fgseaRes.p <- as.matrix(fgseaRes.p)
pdf(file="SnC_NES_ov_cor.pdf", width=7, height = 5)
ht <- Heatmap(fgseaRes.NES,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-NES",
                                          at=c(-2,-1,0,1,2),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(fgseaRes.NES)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(fgseaRes.p[i, j] > 0.05 | is.na(fgseaRes.p[i, j])==T) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(fgseaRes.NES[i, j],2), x, y,gp=gpar(fontsize=6.5))
                }})
ht<-draw(ht)
dev.off()

# medulla
# cortex
setwd("/opt/home/buckcenter.org/fwu/ovary_SC/cell_type_refine/Medulla")
file_list <- list.files(pattern="*.csv")
df_list <- list()
for (file in file_list) {
  df <- read.csv(file)
  cell4gsea <- gsub("_med_DoxvsCt.csv","",file)
  df <- df[order(df$avg_log2FC, decreasing = T),]
  df <- df %>% filter(p_val_adj<0.05) %>% select(X,avg_log2FC)
  ranks <- deframe(df)
  df_list[[cell4gsea]] <- ranks
}
setwd("/opt/home/buckcenter.org/fwu/ovary_SC/")

# Load the pathways into a named list
pathways.SnC <- gmtPathways("Cellular_senescence_genesets.gmt")

# fgseaRes <- fgsea(pathways=pathways.SnC, stats=df_list[[1]])
fgseaRes <- lapply(df_list,function(x) fgsea(pathways=pathways.SnC, stats=x))

# Define a function to modify column names
modify_colnames <- function(df, name) {
  colnames(df)[-1] <- paste0(name, "_", colnames(df)[-1])
  return(df)
}

# Use lapply to apply the function to each data frame in the list
fgseaRes <- lapply(names(fgseaRes), function(x) modify_colnames(fgseaRes[[x]], x))

fgseaRes.df <- fgseaRes %>% reduce(full_join, by = "pathway")
fgseaRes.df <- column_to_rownames(fgseaRes.df,"pathway")
fgseaRes.NES <- fgseaRes.df %>% select(grep("_NES", colnames(fgseaRes.df), value = TRUE))
colnames(fgseaRes.NES) <- gsub("_NES","",colnames(fgseaRes.NES))
fgseaRes.p <- fgseaRes.df %>% select(grep("_padj", colnames(fgseaRes.df), value = TRUE))
colnames(fgseaRes.p) <- gsub("_padj","",colnames(fgseaRes.p))

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

fgseaRes.NES <- as.matrix(fgseaRes.NES)
fgseaRes.p <- as.matrix(fgseaRes.p)
pdf(file="SnC_NES_ov_med.pdf", width=7, height = 5)
ht <- Heatmap(fgseaRes.NES,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-NES",
                                          at=c(-2,-1,0,1,2),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(fgseaRes.NES)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(fgseaRes.p[i, j] > 0.05 | is.na(fgseaRes.p[i, j])==T) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(fgseaRes.NES[i, j],2), x, y,gp=gpar(fontsize=6.5))
                }})
ht<-draw(ht)
dev.off()

# check CD68 marker
Idents(Sample.combined) <- "seurat_clusters"
DefaultAssay(Sample.combined) <- "SCT"
FeaturePlot(Sample.combined, features=c("CD68"), pt.size = 0.2)
VlnPlot(Sample.combined, features = "CD68")

# ----------- to subcluster seurat_clusters #2 --------------- #
ov_sub2 <- FindSubCluster(Sample.combined,2, graph.name = "integrated_snn")
FeaturePlot(ov_sub2, features=c("CD68"), pt.size = 0.2)
Idents(ov_sub2) <- "sub.cluster"
DimPlot(ov_sub2)
DimPlot(ov_sub2, cells = names(ov_sub2$seurat_clusters[ov_sub2$seurat_clusters=="2"]))
DefaultAssay(ov_sub2)
# find markers for every cluster compared to all remaining cells
sub2.markers.all <- FindAllMarkers(ov_sub2, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="tissue")
write.csv(sub2.markers.all,"sub2_markers_clusters.csv")

# 2_0 and 2_2 would be Immune_cells CD68 marker positive
ov_sub2$celltype_refine2 <- ifelse(ov_sub2$sub.cluster == "2_0" | ov_sub2$sub.cluster == "2_2",
                                   "Immune_cells", as.character(ov_sub2$celltype_refine))
# set the levels in order
ov_sub2$celltype_refine2 <- factor(ov_sub2$celltype_refine2,
                                          levels=sort(unique(ov_sub2$celltype_refine2)))

Idents(ov_sub2) <- "celltype_refine2"
DimPlot(ov_sub2)
DimPlot(ov_sub2, reduction = "umap",split.by = "orig.ident",label =T,repel =T)
ggsave("ov_cellTypes_refined2.pdf",width=12,height=6)

# ------ Cell type ratio change stack barplot ------ #
sub.prop.all<-data.frame()
for (l in unique(ov_sub2$orig.ident)){
  sub.treat<-ov_sub2@meta.data[ov_sub2$orig.ident==l,]
  sub.prop<-data.frame(table(sub.treat$celltype_refine2)/sum(table(sub.treat$celltype_refine2)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}
library(tidyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("cor_ct","cor_dox","med_ct","med_dox"))
sub.prop.all$labs<-round(sub.prop.all$Freq,3) # prepare cell proportion label
sub.prop.all$labs[sub.prop.all$labs<0.05]<-""
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="Ovary") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_all2.pdf",width=6,height=5)

# ------ find DEG for all cell types after Dox treatment ------ #
# with default filter min.pct=0.1 and logfc=0.25
ov_sub2$celltype_refine2.treatment <- paste(ov_sub2$celltype_refine2, ov_sub2$orig.ident,
                                                   sep = "_")
Idents(ov_sub2)<-"celltype_refine2.treatment"
plan("multisession", workers = 8)
cell4DEG<-unique(ov_sub2$celltype_refine2)
for (d in cell4DEG){
  DEG.cellType <- FindMarkers(ov_sub2, assay = "SCT", ident.1 =paste0(d,"_cor_dox"), ident.2 =paste0(d,"_cor_ct"))
  write.csv(DEG.cellType, paste0(d,"_cor_DoxvsCt.csv"))
}
for (d in cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]){
  DEG.cellType <- FindMarkers(ov_sub2, assay = "SCT", ident.1 =paste0(d,"_med_dox"), ident.2 =paste0(d,"_med_ct"))
  write.csv(DEG.cellType, paste0(d,"_med_DoxvsCt.csv"))
}


# refine cell annotation results by scRNAseq ovary literatures collected by Cosmo Hahn
Genes_list_human_ovary.23 <- readxl::read_excel("Lengyel et al 2023 ovary DEGs.xlsx")
cell4type <- unique(Genes_list_human_ovary.23$`Annotated clusters`)
cluster4type <- unique(sub2.markers.all$cluster)
# correlation test
correlation.df<-data.frame()
for (c in cell4type){
  cell_FC_ref <- Genes_list_human_ovary.23 %>% filter(`Annotated clusters`==c) %>% select(c(Gene,avg_log2FC,`Annotated clusters`))
  cell_FC_ref <- as.data.frame(cell_FC_ref)
  for (cl in cluster4type){
    cell_FC_target <- sub2.markers.all %>% filter(cluster==cl) %>% select(c(gene,avg_log2FC,cluster))
    FC4cor <- inner_join(cell_FC_ref,cell_FC_target,by=c("Gene" = "gene"))
    correlation <- cor(FC4cor$avg_log2FC.x,FC4cor$avg_log2FC.y)
    cor.df <- data.frame(cluster=cl,cell_type=c,correlation=correlation,"overlaps_ref_genes"=paste0(nrow(FC4cor),"/",nrow(cell_FC_ref)))
    correlation.df <- rbind(correlation.df,cor.df)
  }
}
write.csv(correlation.df,"correlation_celltype_annotation23.csv")

# score the senescent cell again by using the refined2 cell types
DefaultAssay(ov_sub2)<- "RNA"
sc_module <- ov_sub2
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  geneset_features <- list(sc.set)
  # score each senescent gene set
  sc_module <- AddModuleScore(
    object = sc_module,
    features = geneset_features,
    name = names(Senescence_genesets[sc]))
}

FeaturePlot(ov_sub2, features = c("CD68"),split.by = "orig.ident",pt.size = 0.2)
ggsave("CD68.pdf",width=12,height=6)

FeaturePlot(ov_sub2, features = c("CD68"))
ggsave("CD68_merged.pdf",width=7,height=6)
            
table(ov_sub2$celltype_refine2.treatment)


# # volcano plot for overall dox vs ct in cortex
# library(EnhancedVolcano)
# EnhancedVolcano(cor_dox.vs.ct,
#                 lab = rownames(cor_dox.vs.ct),
#                 x = 'avg_log2FC',
#                 y = 'p_val_adj',
#                 #xlim =c(-1.5,2),
#                 #FCcutoff = 0.25,
#                 labSize = 5, 
#                 drawConnectors = T, arrowheads=F, min.segment.length=0.3,
#                 title = "Dox vs Control",
#                 subtitle = bquote(italic("Ovary Cortex")))
# ggsave("uGvs1G_unstimulated_Volcano.pdf",width=10,height=10)

ov_sub2[["SCT"]]@data[row.names(ov_sub2[["SCT"]]@data)=="CD68",]
?FetchData
CD68.p21 <- subset(ov_sub2, subset= orig.ident=="med_ct")
CD68 <- FetchData(CD68.p21,"CD68")
p21 <- FetchData(CD68.p21,"CDKN1A")
p16 <- FetchData(CD68.p21,"CDKN2A")
cor(CD68,p21,method="spearman")
cor.test(CD68[,1],p21[,1],method="spearman")

SC.bulk.test <- merge(CD68,p21)
SC.bulk.test <- merge(SC.bulk.test, p16)

library("ggpubr")
ggscatter(SC.bulk.test, x = "CD68", y = "CDKN1A", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "SC_means", ylab = "bulk_mean", title = "SC and Bulk")


