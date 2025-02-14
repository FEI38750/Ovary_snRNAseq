# following the ov_snSeq_01092023.R
# 10 days Dox treatment/culture

### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)

# change the current plan to access parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize= +Inf) # increase maxSize of future function

setwd("~/ovary_SC/batch3")
# import 10X results
# 3nd batch
ov_sample_list <- c("RTL-868-Cortex-CTRL-1","RTL-886-Cortex-CTRL-1","RTL-868-Medulla-CTRL-1","RTL-886-Medulla-CTRL-1","RTL-868-Cortex-DOXO-1","RTL-886-Cortex-DOXO-1",
    "RTL-868-Medulla-DOXO-1","RTL-886-Medulla-DOXO-1","RTL-872-Cortex-CTRL-1","RTL-872-Cortex-DOXO-1")
Sample.list.cmb <- list()
for (s in ov_sample_list){
  ov_sample <- Read10X(data.dir = paste0("/ovary_SC_3rd/",s,"/outs/filtered_feature_bc_matrix"))
  Sample.list.cmb[[s]] <- ov_sample
}

# create seuratobject for QC
Sample.list.qc.cmb<-list()
for (i in names(Sample.list.cmb)){
  Sample.list.qc.cmb[[i]] <- CreateSeuratObject(Sample.list.cmb[[i]], project=i,min.cells = 3, min.feature = 250)
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
Sample.list.qc.cmb<-lapply(Sample.list.qc.cmb,SCT_m)

# Remove doublet
library(DoubletFinder)
for (i in names(Sample.list.qc.cmb)){
  # pK Identification (no ground-truth)
  sweep.res.list_sample <- paramSweep_v3(Sample.list.qc.cmb[[i]], PCs = 1:30, sct = T)
  sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  pK<-as.numeric(as.character(bcmvn_sample$pK))[bcmvn_sample$BCmetric==max(bcmvn_sample$BCmetric)]
  ## Homotypic Doublet Proportion Estimate
  annotations <- Sample.list.qc.cmb[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Sample.list.qc.cmb[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  Sample.list.qc.cmb[[i]] <- doubletFinder_v3(Sample.list.qc.cmb[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  re.pAnn<-names(Sample.list.qc.cmb[[i]]@meta.data)[(length(Sample.list.qc.cmb[[i]]@meta.data)-1)]
  Sample.list.qc.cmb[[i]] <- doubletFinder_v3(Sample.list.qc.cmb[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = re.pAnn, sct = T)
}

# filter samples based on results of QC & Doublet
for (i in names(Sample.list.qc.cmb)){
  sample.meta.data<-Sample.list.qc.cmb[[i]]@meta.data
  singlet<-row.names(sample.meta.data)[sample.meta.data[length(sample.meta.data)]=="Singlet"]
  # Get the indices of the desired columns
  column_indices <- which(colnames(Sample.list.cmb[[i]]) %in% singlet)
  # Select the columns from the sparse matrix
  Sample.list.cmb[[i]]<-Sample.list.cmb[[i]][,column_indices]
}
# create seuratobject for integration
for (i in names(Sample.list.cmb)){
  Sample.list.cmb[[i]] <- CreateSeuratObject(Sample.list.cmb[[i]], project=i)
  Sample.list.cmb[[i]] <- SCTransform(Sample.list.cmb[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}

#Sample.list.cmb[["fist_batch"]] <- NULL
features.cmb <- SelectIntegrationFeatures(object.list = Sample.list.cmb, nfeatures = 3000)
# run PCA before RPCA integration
Sample.list.cmb <- PrepSCTIntegration(object.list = Sample.list.cmb, anchor.features = features.cmb)
Sample.list.cmb <- lapply(X = Sample.list.cmb, FUN = RunPCA, features = features.cmb)

# Perform RPCA integration
Anchors.cmb <- FindIntegrationAnchors(object.list = Sample.list.cmb, normalization.method = "SCT",
                                  anchor.features = features.cmb, reduction = "rpca")
Sample.cmb <- IntegrateData(anchorset = Anchors.cmb, normalization.method = "SCT")


# Dimensional reduction
Sample.cmb <- RunPCA(Sample.cmb, npcs = 20 ,verbose = FALSE)
ElbowPlot(Sample.cmb,ndims = 20)
Sample.cmb <- RunUMAP(Sample.cmb, reduction = "pca", dims = 1:20)
# Cluster the cells
Sample.cmb <- FindNeighbors(Sample.cmb, dims = 1:20)
Sample.cmb <- FindClusters(Sample.cmb, resolution = 0.3)
DimPlot(Sample.cmb, reduction = "umap",split.by = "orig.ident",label =T,ncol=4)
ggsave("ov_clusters_combined.pdf",width=14,height=8)

# Prepare object to run differential expression on SCT assay with multiple models
Sample.cmb <- PrepSCTFindMarkers(Sample.cmb)

# add metadata for tissue, batches, donor, and treatment
library(dplyr)
Sample.cmb$tissue <- case_when(Sample.cmb$orig.ident %in%
                                 c("RTL-868-Cortex-CTRL-1","RTL-886-Cortex-CTRL-1","RTL-872-Cortex-CTRL-1",
                                   "RTL-868-Cortex-DOXO-1","RTL-886-Cortex-DOXO-1","RTL-872-Cortex-DOXO-1") ~ "cortex",
                               .default = "medulla")

Sample.cmb$donor <- case_when(Sample.cmb$orig.ident %in%
                                  c("RTL-868-Cortex-CTRL-1","RTL-868-Medulla-CTRL-1","RTL-868-Cortex-DOXO-1","RTL-868-Medulla-DOXO-1") ~ "donor868",
                              Sample.cmb$orig.ident %in%
                                c("RTL-886-Cortex-CTRL-1","RTL-886-Medulla-CTRL-1","RTL-886-Cortex-DOXO-1","RTL-886-Medulla-DOXO-1") ~ "donor886",
                              Sample.cmb$orig.ident %in%
                                c("RTL-872-Cortex-CTRL-1","RTL-872-Cortex-DOXO-1") ~ "donor872")

Sample.cmb$treatment <- case_when(Sample.cmb$orig.ident %in%
                                 c("RTL-868-Cortex-DOXO-1","RTL-886-Cortex-DOXO-1","RTL-872-Cortex-DOXO-1",
                                   "RTL-868-Medulla-DOXO-1","RTL-886-Medulla-DOXO-1") ~ "dox",
                               .default = "ctl")

# find markers for every cluster compared to all remaining cells
Idents(Sample.cmb) <- "seurat_clusters"
markers.all <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(markers.all,"markers_all_clusters.csv")



# plot top markers in heatmap
library(dplyr)
markers.all %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(Sample.cmb, features=top10$gene)
ggsave("ov_markers_top10.pdf",width=8,height=16)



# find overall DEGs dox vs ctl for tissues
Sample.cmb$tissue_treat <- paste0(Sample.cmb$tissue,"_",Sample.cmb$treatment)
Idents(Sample.cmb) <- "tissue_treat"
# cortex overall DEGs
cor_dox.vs.ct <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "cortex_dox", ident.2 = "cortex_ctl",
                             min.pct = 0.1, logfc.threshold = 0.2, test.use = "MAST",latent.vars="donor")
write.csv(cor_dox.vs.ct,"Overall_DEGs_cor_dox_vs_ct.csv")
# medulla overall DEGs
med_dox.vs.ct <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "medulla_dox", ident.2 = "medulla_ctl",
                             min.pct = 0.1, logfc.threshold = 0.2, test.use = "MAST",latent.vars="donor")
write.csv(med_dox.vs.ct,"Overall_DEGs_med_dox_vs_ct.csv")

# UMAP before annotation
Idents(Sample.cmb) <- "seurat_clusters"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T)
ggsave("ov_seurat_clusters.pdf",width=12,height=6)

library(SingleR)
HPCA.ref <- celldex::HumanPrimaryCellAtlasData()

# set the active assay back to “RNA,” and re-do the normalization
DefaultAssay(Sample.cmb) <- "RNA"
Sample.cmb <- NormalizeData(Sample.cmb)
# convert Seurat object to single cell experiment (SCE) for convenience
sce <- as.SingleCellExperiment(Sample.cmb,assay = "RNA")
HPCA.main <- SingleR(test = sce,assay.type.test = 1,ref = HPCA.ref,labels = HPCA.ref$label.main)
# see the summary of general cell type annotations
table(HPCA.main$pruned.labels)
# add the annotations to the Seurat object metadata
Sample.cmb@meta.data$HPCA.main <- HPCA.main$pruned.labels

# most frequent predicted cell type for each cluster
Sample.cmb$predicted.cellType <- as.character(Sample.cmb$seurat_clusters)
for (c in unique(Sample.cmb$seurat_clusters)){
  cell_type <- Sample.cmb@meta.data %>%
    filter(seurat_clusters==c) %>%
    select(HPCA.main)
  cell_type_counts<-table(cell_type)
  most_frequent_cell_type<-names(which.max(cell_type_counts))
  Sample.cmb$predicted.cellType[Sample.cmb$predicted.cellType==c]<-most_frequent_cell_type
}
Idents(Sample.cmb) <- "predicted.cellType"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T)
ggsave("ov_cellTypes_HPCA.pdf",width=12,height=6)

DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue",label =T)
ggsave("ov_cellTypes_UMAP_tissue.pdf",width=12,height=6)

Idents(Sample.cmb) <- "tissue"
DimPlot(Sample.cmb, reduction = "umap",split.by = "treatment",label =F)
ggsave("ov_cellTypes_UMAP_tissue_treatment.pdf",width=10,height=5)

# count cell number for each tissue and condition
table(Sample.cmb$tissue_treat)



# quick check cell annotation results for a cluster
table(Sample.cmb@meta.data %>%
        filter(seurat_clusters==3) %>%
        select(HPCA.main)) %>% sort(decreasing = T)

# refine cell annotation results by scRNAseq ovary literatures collected by researcher
Genes_list_human_ovary <- readxl::read_excel("../Genes list - human ovary updated with PMO.xls")
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
ov_module <- Sample.cmb
for (c in cell4type){
  cell_FC_ref_pos <- Genes_list_human_ovary %>% filter(cluster==c & avg_logFC>0) %>% select(gene)
  celltype_features <- list(as.character(cell_FC_ref_pos$gene))
  ov_module <- AddModuleScore(
    object = ov_module,
    features = celltype_features,
    assay = "RNA",
    name = c)
}

celltype.df <- ov_module@meta.data %>% select("seurat_clusters",paste0(cell4type,"1"))
colnames(celltype.df) <- gsub(" ","_",colnames(celltype.df))
colnames(celltype.df) <- gsub("&_","",colnames(celltype.df))
colnames(celltype.df) <- gsub("/","_or_",colnames(celltype.df))
colnames(celltype.df) <- gsub("1","",colnames(celltype.df))

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
ggsave(paste0("Box Plot of celltype scores.pdf"),plot = grid_plot,width=12,height=14)


# annotate cell type by combining results from HPCA and manually curate the gene list from ovary publications
Sample.cmb$celltype_refine <- case_when(Sample.cmb$seurat_clusters ==
                                               "0" ~ "Stroma_1",
                                             Sample.cmb$seurat_clusters ==
                                               "1" ~ "Smooth_muscle_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "2" ~ "Stroma_1",
                                             Sample.cmb$seurat_clusters ==
                                               "3" ~ "unknown",
                                             Sample.cmb$seurat_clusters ==
                                               "4" ~ "Stroma_2",
                                             Sample.cmb$seurat_clusters ==
                                               "5" ~ "Epithelial_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "6" ~ "Stroma_1",
                                             Sample.cmb$seurat_clusters ==
                                               "7" ~ "Stroma_2",
                                             Sample.cmb$seurat_clusters ==
                                               "8" ~ "Perivascular_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "9" ~ "Stroma_2",
                                             Sample.cmb$seurat_clusters ==
                                               "10" ~ "Endothelial_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "11" ~ "Epithelial_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "12" ~ "Immune_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "13" ~ "Endothelial_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "14" ~ "Immune_cells",
                                             Sample.cmb$seurat_clusters ==
                                               "15" ~ "Smooth_muscle_cells")

# set the levels in order
Sample.cmb$celltype_refine <- factor(Sample.cmb$celltype_refine,
                                          levels=sort(unique(Sample.cmb$celltype_refine)))
# plot UMAP based on refined cell types
Idents(Sample.cmb) <- "celltype_refine"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
ggsave("ov_cellTypes_refined.pdf",width=12,height=6)

# find markers for each refined cell type compared to all remaining cell types
Idents(Sample.cmb) <- "celltype_refine"
CellType.markers.all <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(CellType.markers.all,"CellType_markers_all_clusters.csv")

# ------ Cell type ratio change stack barplot ------ #
sub.prop.all<-data.frame()
for (l in unique(Sample.cmb$tissue_treat)){
  sub.treat<-Sample.cmb@meta.data[Sample.cmb$tissue_treat==l,]
  sub.prop<-data.frame(table(sub.treat$celltype_refine)/sum(table(sub.treat$celltype_refine)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}
write.csv(sub.prop.all,"cell.prop.all.csv")

library(tidyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("cortex_ctl","cortex_dox","medulla_ctl","medulla_dox"))
sub.prop.all$labs<-round(sub.prop.all$Freq,3) # prepare cell proportion label
sub.prop.all$labs[sub.prop.all$labs<0.05]<-""
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="Ovary") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("Cell_prop_stacked_barplot_all.pdf",width=6,height=5)

# Calculate summary statistics
Cell_Results4plot_summary <- Sample.cmb@meta.data %>%
  group_by(tissue_treat,celltype_refine) %>%
  summarise(cell_num = length(celltype_refine)) %>%
  group_by(tissue_treat) %>%
  summarise(celltype_refine = celltype_refine,
            cell_num = cell_num,
            sum = sum(cell_num),
            Freq_mean = cell_num/sum(cell_num))

Cell_Results_donor <- Sample.cmb@meta.data %>%
  group_by(tissue_treat,celltype_refine,donor) %>%
  summarise(cell_num_donor = length(donor)) %>%
  group_by(tissue_treat,donor) %>%
  summarise(celltype_refine,
            cell_num_donor,
            sum_donor = sum(cell_num_donor),
            Freq_donor = cell_num_donor/sum_donor) %>%
  group_by(celltype_refine, tissue_treat) %>%
  summarise(donor,cell_num_donor,sum_donor,Freq_donor,
            se = sd(Freq_donor) / sqrt(n()))

Cell_Results4plot_summary <- left_join(Cell_Results4plot_summary,Cell_Results_donor,
                                       by=c("tissue_treat","celltype_refine"))

# Create the bar plot with T-shaped error bars
bar_plot <- ggplot(Cell_Results4plot_summary, aes(x = celltype_refine, y = Freq_mean, fill = tissue_treat)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Freq_mean - se, ymax = Freq_mean + se), position = position_dodge(0.9)) +
  ggtitle("Ovary Cell Type Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Add the individual "Freq" values as dots
dot_plot <- ggplot(Cell_Results4plot_summary, aes(x = celltype_refine, y = Freq_donor, fill = tissue_treat)) +
  geom_jitter(size = 1, width = 0.2 / length(unique(Cell_Results4plot_summary$tissue_treat))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Combine the bar plot with error bars and the dot plot
combined_plot <- bar_plot + 
  geom_point(data = Cell_Results4plot_summary, aes(x = celltype_refine, y = Freq_donor, fill = tissue_treat),
             position = position_dodge(0.9), size = 0.5, alpha=0.6)
print(combined_plot)
ggsave("Cell_type_frequency_bar.pdf",width=8,height=6)

# cell type ratio change: log2FC comparision unstim; remove average celltype ratio <1%
# Cortex
sub.prop.cortex_ctl<-subset(Cell_Results4plot_summary, tissue_treat=="cortex_ctl")
sub.prop.cortex_dox<-subset(Cell_Results4plot_summary, tissue_treat=="cortex_dox")
sub.prop.cortex.comm<-intersect(sub.prop.cortex_ctl$celltype_refine,sub.prop.cortex_dox$celltype_refine)
sub.prop.cortex_ctl <- subset(sub.prop.cortex_ctl, celltype_refine %in% sub.prop.cortex.comm)
sub.prop.cortex_dox <- subset(sub.prop.cortex_dox, celltype_refine %in% sub.prop.cortex.comm)
sub.prop.cor_DoxvsCtl <- data.frame(Cell_type=sub.prop.cortex_ctl$celltype_refine,
                                    group="Dox vs Ctl Cortex",
                                    Fold_change=sub.prop.cortex_dox$Freq_mean/sub.prop.cortex_ctl$Freq_mean)
sub.prop.cor_DoxvsCtl$log2FC<-log2(sub.prop.cor_DoxvsCtl$Fold_change)

sub.prop.cor_DoxvsCtl<-sub.prop.cor_DoxvsCtl%>% distinct()

# add t-test p-value
sub.prop.cortex_ctl.full <- subset(sub.prop.cortex_ctl, tissue_treat=="cortex_ctl")
sub.prop.cortex_dox.full <- subset(sub.prop.cortex_dox, tissue_treat=="cortex_dox")

t.test.num <- numeric()
for (c in sub.prop.cortex.comm){
  sub.prop.cortex_ctl.temp <- sub.prop.cortex_ctl.full %>% filter(celltype_refine == c) %>% select(Freq_donor)
  sub.prop.cortex_dox.temp <- sub.prop.cortex_dox.full %>% filter(celltype_refine == c) %>% select(Freq_donor)
  t.test.temp <- t.test(sub.prop.cortex_ctl.temp$Freq_donor, sub.prop.cortex_dox.temp$Freq_donor,paired =T)
  t.test.num <- c(t.test.num,t.test.temp$p.value)
}
sub.prop.cor_DoxvsCtl$p_value <- t.test.num
# bar plot
ggplot(sub.prop.cor_DoxvsCtl, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Dox vs Ctl Ovary Cortex") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  # geom_text(aes(label = ifelse(Cell_type == "CD14.Mono", "*", "")),
  #           position = position_dodge(0.9), vjust = 1.5,size=5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_Cortex_Dox_vs_Ctl.pdf",width=6,height=4)

# Medulla
sub.prop.medulla_ctl<-subset(Cell_Results4plot_summary, tissue_treat=="medulla_ctl" & Freq_mean>0.01)
sub.prop.medulla_dox<-subset(Cell_Results4plot_summary, tissue_treat=="medulla_dox" & Freq_mean>0.01)
sub.prop.medulla.comm<-intersect(sub.prop.medulla_ctl$celltype_refine,sub.prop.medulla_dox$celltype_refine)
sub.prop.medulla_ctl <- subset(sub.prop.medulla_ctl, celltype_refine %in% sub.prop.medulla.comm)
sub.prop.medulla_dox <- subset(sub.prop.medulla_dox, celltype_refine %in% sub.prop.medulla.comm)
sub.prop.med_DoxvsCtl <- data.frame(Cell_type=sub.prop.medulla_ctl$celltype_refine,
                                    group="Dox vs Ctl medulla",
                                    Fold_change=sub.prop.medulla_dox$Freq_mean/sub.prop.medulla_ctl$Freq_mean)
sub.prop.med_DoxvsCtl$log2FC<-log2(sub.prop.med_DoxvsCtl$Fold_change)

sub.prop.med_DoxvsCtl<-sub.prop.med_DoxvsCtl%>% distinct()

# add t-test p-value
sub.prop.medulla_ctl.full <- subset(sub.prop.medulla_ctl, tissue_treat=="medulla_ctl")
sub.prop.medulla_dox.full <- subset(sub.prop.medulla_dox, tissue_treat=="medulla_dox")

t.test.num <- numeric()
for (c in sub.prop.medulla.comm){
  sub.prop.medulla_ctl.temp <- sub.prop.medulla_ctl.full %>% filter(celltype_refine == c) %>% select(Freq_donor)
  sub.prop.medulla_dox.temp <- sub.prop.medulla_dox.full %>% filter(celltype_refine == c) %>% select(Freq_donor)
  t.test.temp <- t.test(sub.prop.medulla_ctl.temp$Freq_donor, sub.prop.medulla_dox.temp$Freq_donor,paired =T)
  t.test.num <- c(t.test.num,t.test.temp$p.value)
}
sub.prop.med_DoxvsCtl$p_value <- t.test.num
# bar plot
ggplot(sub.prop.med_DoxvsCtl, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Dox vs Ctl Ovary Medulla") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  # geom_text(aes(label = ifelse(Cell_type == "CD14.Mono", "*", "")),
  #           position = position_dodge(0.9), vjust = 1.5,size=5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_medulla_Dox_vs_Ctl.pdf",width=6,height=4)


# ------ find DEG for all cell types after Dox treatment ------ #
# with default filter min.pct=0.1 and logfc=0.25
Sample.cmb$celltype_treatment <- paste(Sample.cmb$celltype_refine, Sample.cmb$tissue_treat,
                                       sep = "_")
Idents(Sample.cmb)<-"celltype_treatment"
plan("multisession", workers = 4)
cell4DEG<-unique(Sample.cmb$celltype_refine)
for (d in cell4DEG){
  DEG.cellType <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 =paste0(d,"_cortex_dox"), ident.2 =paste0(d,"_cortex_ctl"),
                              min.pct = 0.1, logfc.threshold = 0.2,test.use = "MAST",latent.vars="donor")
  write.csv(DEG.cellType, paste0(d,"_cor_DoxvsCt.csv"))
}
for (d in cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]){
  DEG.cellType <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 =paste0(d,"_medulla_dox"), ident.2 =paste0(d,"_medulla_ctl"),
                              min.pct = 0.1, logfc.threshold = 0.2,test.use = "MAST",latent.vars="donor")
  write.csv(DEG.cellType, paste0(d,"_med_DoxvsCt.csv"))
}


# plot UMAP visualization for each participant for each condition
Idents(Sample.cmb) <- "donor"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
# subset by donor
# donor868
Idents(Sample.cmb) <- "donor"
Sample.cmb1 <- subset(Sample.cmb, idents="donor868")
Idents(Sample.cmb1) <- "celltype_refine"
DimPlot(Sample.cmb1, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
ggsave("ov_cellTypes_refined_donor868.pdf",width=12,height=6)
# donor872
Sample.cmb2 <- subset(Sample.cmb, idents="donor872")
Idents(Sample.cmb2) <- "celltype_refine"
DimPlot(Sample.cmb2, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
ggsave("ov_cellTypes_refined_donor872.pdf",width=12,height=6)
# donor886
Sample.cmb3 <- subset(Sample.cmb, idents="donor886")
Idents(Sample.cmb3) <- "celltype_refine"
DimPlot(Sample.cmb3, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
ggsave("ov_cellTypes_refined_donor886.pdf",width=12,height=6)


# ---------------- #

# score the senecent cells
library(AUCell)
# load the senescence gene sets
Senescence_genesets <- read.csv("/ovary_SC/Senescence_genesets.csv")

sc_module <- Sample.cmb
# scoring the senescence
geneset_features <- list()
#sc_scores <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  geneset_features[colnames(Senescence_genesets)[sc]] <- list(sc.set)
  # score each senescent gene set
  cells_AUC <- AUCell_run(
    Sample.cmb@assays[["RNA"]]@counts,geneset_features)
  #sc_scores[[names(Senescence_genesets)[sc]]] <- cells_AUC
}


# integrate to metadata
cell_AUC_results <- t(cells_AUC@assays@data@listData[["AUC"]])
sc_module@meta.data <- merge(sc_module@meta.data,cell_AUC_results, by="row.names")
row.names(sc_module@meta.data) <- sc_module@meta.data$Row.names
sc_module@meta.data <- sc_module@meta.data[,-1]

# for Overall comparison: cortex
sc_module@meta.data$tissue_treat <- paste0(sc_module@meta.data$tissue,"_",sc_module@meta.data$treatment)
sc_score <- data.frame()
sc_score_p <- data.frame()
for (sc in 1:ncol(Senescence_genesets)){
  ctl <- sc_module@meta.data %>% 
    filter(tissue_treat=="cortex_ctl") %>% 
    select(names(Senescence_genesets[sc]))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="cortex_dox") %>% 
    select(names(Senescence_genesets[sc]))
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
library(tibble)
sc_score <- column_to_rownames(sc_score,"SC_gene_set")
sc_score_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))

sc_score <- as.matrix(sc_score)
pdf(file="SnC-score_ov_cor_sc_overall.pdf", width=6, height = 5)
ht <- Heatmap(sc_score,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox and Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(sc_score)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(sc_score_p[i, j] > 0.05) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(sc_score[i, j],3), x, y,gp=gpar(fontsize=6.5))
                }})
ht<-draw(ht)
dev.off()

# for Overall comparison: medulla
sc_score <- data.frame()
sc_score_p <- data.frame()
for (sc in 1:ncol(Senescence_genesets)){
  ctl <- sc_module@meta.data %>% 
    filter(tissue_treat=="medulla_ctl") %>% 
    select(names(Senescence_genesets[sc]))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="medulla_dox") %>% 
    select(names(Senescence_genesets[sc]))
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
library(tibble)
sc_score <- column_to_rownames(sc_score,"SC_gene_set")
sc_score_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))

sc_score <- as.matrix(sc_score)
pdf(file="SnC-score_ov_med_sc_overall.pdf", width=6, height = 5)
ht <- Heatmap(sc_score,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox and Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(sc_score)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(sc_score_p[i, j] > 0.05) {
                  grid.text("·", x, y)
                }else{
                  grid.text(round(sc_score[i, j],3), x, y,gp=gpar(fontsize=6.5))
                }})
ht<-draw(ht)
dev.off()

# add the genes to the overall senescence heatmap: cortex
overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- cor_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- cor_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Down","cell_senescence_signatures",
                "REACTOME_CELLULAR_SENESCENCE","REACTOME_SASP","GOBP_CELLULAR_SENESCENCE","CellAge_Up","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"overlapped_SC_DEGs_cortex.csv")

library(ComplexHeatmap)
library(circlize)
#col_fun = colorRamp2(c(0, 1.5), c("white", "red"))

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

# select top20 SnC genes
overlapped_SC_DEGs_top20 <- names(sort(colSums(overlapped_SC_DEGs.df),decreasing =T)[1:20])

pdf(file="SnC-score_ov_cor_sc_overall_genes_color2_top20.pdf", width=9, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df[,overlapped_SC_DEGs_top20],
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex: Doxo vs ctrl. Senescence score-associated genes",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              # row_names_gp = gpar(fontsize = 9.5),
              # column_names_gp = gpar(fontsize = 8)
              )
ht<-draw(ht)
dev.off()

# add the genes to the overall senescence heatmap: Medulla
overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- med_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- med_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","REACTOME_CELLULAR_SENESCENCE","CellAge_Down","cell_senescence_signatures",
                "SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Up",
                "GOBP_CELLULAR_SENESCENCE","REACTOME_SASP","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"overlapped_SC_DEGs_medulla.csv")

library(ComplexHeatmap)
library(circlize)
#col_fun = colorRamp2(c(0, 1.5), c("white", "red"))

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

# select top20 SnC genes
overlapped_SC_DEGs_top20 <- names(sort(colSums(overlapped_SC_DEGs.df),decreasing =T)[1:20])

pdf(file="SnC-score_ov_med_sc_overall_genes_color2_top20.pdf", width=9, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df[,overlapped_SC_DEGs_top20],
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla: Doxo vs ctrl. Senescence score-associated genes",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              # row_names_gp = gpar(fontsize = 9.5),
              # column_names_gp = gpar(fontsize = 8)
              )
ht<-draw(ht)
dev.off()


# in cortex, check epithelial and stromal_1
# epithelial
library(tibble)
epithlial_cor_dox.vs.ct <- read.csv("Epithelial_cells_cor_DoxvsCt.csv")
epithlial_cor_dox.vs.ct <- column_to_rownames(epithlial_cor_dox.vs.ct,"X")

overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- epithlial_cor_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- epithlial_cor_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Down","cell_senescence_signatures",
                "REACTOME_CELLULAR_SENESCENCE","REACTOME_SASP","GOBP_CELLULAR_SENESCENCE","CellAge_Up","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"epithlial_cor_overlapped_SC_DEGs_cortex.csv")

library(ComplexHeatmap)
library(circlize)
#col_fun = colorRamp2(c(0, 1.5), c("white", "red"))

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
library(viridis)
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

# select top20 SnC genes
overlapped_SC_DEGs_top20 <- names(sort(colSums(overlapped_SC_DEGs.df),decreasing =T)[1:20])

pdf(file="SnC-score_ov_cor_epithlial_sc_overall_genes_color2_top20.pdf", width=9, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df[,overlapped_SC_DEGs_top20],
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Epithelial cells: Doxo vs ctrl",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              # row_names_gp = gpar(fontsize = 9.5),
              # column_names_gp = gpar(fontsize = 8)
              )
ht<-draw(ht)
dev.off()

# stromal_1
# Cortex
stromal_1_cor_dox.vs.ct <- read.csv("Stromal_1_cor_DoxvsCt.csv")
stromal_1_cor_dox.vs.ct <- column_to_rownames(stromal_1_cor_dox.vs.ct,"X")

overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- stromal_1_cor_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- stromal_1_cor_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Down","cell_senescence_signatures",
                "REACTOME_CELLULAR_SENESCENCE","REACTOME_SASP","GOBP_CELLULAR_SENESCENCE","CellAge_Up","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"Stromal_1_cor_overlapped_SC_DEGs_cortex.csv")

library(ComplexHeatmap)
library(circlize)
#col_fun = colorRamp2(c(0, 1.5), c("white", "red"))

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

# select top20 SnC genes
overlapped_SC_DEGs_top20 <- names(sort(colSums(overlapped_SC_DEGs.df),decreasing =T)[1:20])

pdf(file="SnC-score_ov_cor_Stromal_1_sc_overall_genes_color2_top20.pdf", width=9, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df[,overlapped_SC_DEGs_top20],
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Stroma 1: Doxo vs ctrl",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              # row_names_gp = gpar(fontsize = 9.5),
              # column_names_gp = gpar(fontsize = 8)
              )
ht<-draw(ht)
dev.off()

# Medulla
# stromal_1
stromal_1_med_dox.vs.ct <- read.csv("Stromal_1_med_DoxvsCt.csv")
stromal_1_med_dox.vs.ct <- column_to_rownames(stromal_1_med_dox.vs.ct,"X")

overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- stromal_1_med_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- stromal_1_med_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Down","cell_senescence_signatures",
                "REACTOME_CELLULAR_SENESCENCE","REACTOME_SASP","GOBP_CELLULAR_SENESCENCE","CellAge_Up","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"Stromal_1_med_overlapped_SC_DEGs_medulla.csv")

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

pdf(file="SnC-score_ov_med_Stromal_1_sc_overall_genes_color2.pdf", width=9, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Stromal_1 Dox and Control",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              row_names_gp = gpar(fontsize = 9.5),
              column_names_gp = gpar(fontsize = 9.5))
ht<-draw(ht)
dev.off()





# add cell type information for the senescent gene comparison
#sc_module <- AddMetaData(Sample.cmb,sc_module@meta.data[11:22])

# add the genes to the overall senescence heatmap: cortex
overlapped_SC_DEGs.ls <- list()
for (sc in 1:ncol(Senescence_genesets)){
  sc.set <- Senescence_genesets[,sc]
  sc.set <- sc.set[!is.na(sc.set) & sc.set != ""]
  SC_DEGs <- cor_dox.vs.ct %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% rownames()
  overlapped_SC_DEGs <- intersect(sc.set,SC_DEGs)
  overlapped_SC_DEGs <- cor_dox.vs.ct[overlapped_SC_DEGs,] %>% select(avg_log2FC)
  colnames(overlapped_SC_DEGs) <- names(Senescence_genesets[sc])
  overlapped_SC_DEGs <- rownames_to_column(overlapped_SC_DEGs,var="SC_genes")
  overlapped_SC_DEGs.ls[[names(Senescence_genesets[sc])]] <- overlapped_SC_DEGs
}

library(purrr)
overlapped_SC_DEGs.df <- reduce(overlapped_SC_DEGs.ls,full_join,by="SC_genes")
overlapped_SC_DEGs.df <- as.matrix(overlapped_SC_DEGs.df)
overlapped_SC_DEGs.df <- t(overlapped_SC_DEGs.df)
colnames(overlapped_SC_DEGs.df) <- overlapped_SC_DEGs.df[1,]
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[-1,]
overlapped_SC_DEGs.df <- matrix(as.numeric(overlapped_SC_DEGs.df),nrow=nrow(overlapped_SC_DEGs.df),ncol=ncol(overlapped_SC_DEGs.df),
                                dimnames=dimnames(overlapped_SC_DEGs.df))
overlapped_SC_DEGs.df[is.na(overlapped_SC_DEGs.df)] <- 0
# re-arrange the order according to the SC-score result
row_orders <- c("SASP_Atlas.IR.Fibroblasts.","SASP_Ovary","SASP_Atlas.IR.Epithelial.","Aging_Markers_Bikem","CellAge_Down","cell_senescence_signatures",
                "REACTOME_CELLULAR_SENESCENCE","REACTOME_SASP","GOBP_CELLULAR_SENESCENCE","CellAge_Up","SenMayo")
overlapped_SC_DEGs.df <- overlapped_SC_DEGs.df[row_orders,]
write.csv(overlapped_SC_DEGs.df,"overlapped_SC_DEGs_cortex.csv")

library(ComplexHeatmap)
library(circlize)
#col_fun = colorRamp2(c(0, 1.5), c("white", "red"))

# Define specific breakpoints
breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(length(breaks)) # color theme 1
colors <- viridis(length(breaks)) # color theme 2
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

pdf(file="SnC-score_ov_cor_sc_overall_genes_color1.pdf", width=18, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox and Control",
              cluster_columns = T,
              cluster_rows = F,
              heatmap_legend_param = list(title="Log2FC",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(overlapped_SC_DEGs.df)),
              row_names_gp = gpar(fontsize = 9.5),
              column_names_gp = gpar(fontsize = 8))
ht<-draw(ht)
dev.off()

sc_module@meta.data <- sc_module@meta.data %>% select(-celltype_refine)
sc_module@meta.data <- merge(sc_module@meta.data,as.data.frame(Sample.cmb$celltype_refine),by="row.names")
sc_module@meta.data <- column_to_rownames(sc_module@meta.data, "Row.names")
colnames(sc_module@meta.data)[colnames(sc_module@meta.data)=="Sample.cmb$celltype_refine"] <- "celltype_refine"

sc_module$celltype_treatment <- paste(sc_module$celltype_refine, sc_module$tissue_treat,
                                       sep = "_")

# cortex
library(tibble)
sc_score.ls<-list()
sc_score_p.ls<-list()
for (cell in levels(cell4DEG)){
  sc_score <- data.frame()
  sc_score_p <- data.frame()
  for (sc in 1:ncol(Senescence_genesets)){
    names(Senescence_genesets[sc])
    ctl <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_cortex_ctl")) %>% 
      select(names(Senescence_genesets[sc]))
    dox <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_cortex_dox")) %>% 
      select(names(Senescence_genesets[sc]))
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
    filter(tissue_treat=="cortex_ctl") %>% 
    select(names(Senescence_genesets[sc]))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="cortex_dox") %>% 
    select(names(Senescence_genesets[sc]))
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
library(purrr)
sc_score <- sc_score.ls %>% purrr::reduce(left_join, by = "SC_gene_set")
sc_score_p <- sc_score_p.ls %>% purrr::reduce(left_join, by = "SC_gene_set")

sc_score[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score[,1])
sc_score_p[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score_p[,1])
library(tibble)
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))
column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_cor_cmb.pdf", width=10, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox and Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
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

# log2FC heatmap
cell4DEG<-unique(Sample.cmb$celltype_refine)
cell4DEG_fc <- c(levels(cell4DEG),"Overall")
#cell4DEG_fc <- c(as.character(cell4DEG),"Overall")
score4hp_fc <- list()
for (fc in cell4DEG_fc){
  score4hp_fcc <- log2(score4hp[,paste0(fc,"_Dox")]/score4hp[,paste0(fc,"_Ctl")])
  score4hp_fc[[fc]] <- score4hp_fcc
}
score4hp_fc_cmp <- purrr::reduce(score4hp_fc,cbind)
colnames(score4hp_fc_cmp) <- names(score4hp_fc)
score4hp_p_fc <- score4hp_p[seq(1,length(cell4DEG_fc)*2,2)]
write.csv(score4hp_fc_cmp,"score4hp_fc_cmp_cor_10D.csv")
write.csv(score4hp_p_fc,"score4hp_p_fc_cor_10D.csv")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
#score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_cor_cmb_fc.pdf", width=7, height = 5)
ht <- Heatmap(score4hp_fc_cmp,
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="Log2FC of \nSnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p_fc[i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()
# Immune  and stromal cells
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf(file="SnC-score_ov_cor_Immu_Strom_fc.pdf", width=5, height = 5)
ht <- Heatmap(score4hp_fc_cmp[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)],
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="Log2FC of \nSnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p_fc[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)][i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()

cor_imm_str <- score4hp_fc_cmp[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)]
cor_imm_str.p <- score4hp_p_fc[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)]

# for medula
sc_score.ls<-list()
sc_score_p.ls<-list()
cell4DEG_medulla <- levels(cell4DEG)[!levels(cell4DEG) %in% c("Epithelial_cells","Tissue_stem_cells")]
for (cell in cell4DEG_medulla){
  sc_score <- data.frame()
  sc_score_p <- data.frame()
  for (sc in 1:ncol(Senescence_genesets)){
    names(Senescence_genesets[sc])
    ctl <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_medulla_ctl")) %>% 
      select(names(Senescence_genesets[sc]))
    dox <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_medulla_dox")) %>% 
      select(names(Senescence_genesets[sc]))
    compare_result <- t.test(dox,ctl)
    avg_dox <- mean(dox[,1])
    avg_ctl <- mean(ctl[,1])
    p.value <- compare_result$p.value
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
}
# for Overall comparison
sc_score <- data.frame()
sc_score_p <- data.frame()
for (sc in 1:ncol(Senescence_genesets)){
  ctl <- sc_module@meta.data %>% 
    filter(tissue_treat=="medulla_ctl") %>% 
    select(names(Senescence_genesets[sc]))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="medulla_dox") %>% 
    select(names(Senescence_genesets[sc]))
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
sc_score <- sc_score.ls %>% purrr::reduce(left_join, by = "SC_gene_set")
sc_score_p <- sc_score_p.ls %>% purrr::reduce(left_join, by = "SC_gene_set")

sc_score[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score[,1])
sc_score_p[,1]<-gsub("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","REACTOME_SASP",sc_score_p[,1])
library(tibble)
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))
column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_med_cmb.pdf", width=9, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox and Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="SnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
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

# log2FC heatmap
cell4DEG_fc <- c(levels(cell4DEG)[!levels(cell4DEG) %in% c("Epithelial_cells","Tissue_stem_cells")],"Overall")
score4hp_fc <- list()
for (fc in cell4DEG_fc){
  score4hp_fcc <- log2(score4hp[,paste0(fc,"_Dox")]/score4hp[,paste0(fc,"_Ctl")])
  score4hp_fc[[fc]] <- score4hp_fcc
}
score4hp_fc_cmp <- purrr::reduce(score4hp_fc,cbind)
colnames(score4hp_fc_cmp) <- names(score4hp_fc)
score4hp_p_fc <- score4hp_p[seq(1,length(cell4DEG_fc)*2,2)]
write.csv(score4hp_fc_cmp,"score4hp_fc_cmp_med_10D.csv")
write.csv(score4hp_p_fc,"score4hp_p_fc_med_10D.csv")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
#score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_med_cmb_fc.pdf", width=7, height = 5)
ht <- Heatmap(score4hp_fc_cmp,
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="Log2FC of \nSnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p_fc[i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()
# Immune  and stromal cells
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf(file="SnC-score_ov_med_Immu_Strom_fc.pdf", width=5, height = 5)
ht <- Heatmap(score4hp_fc_cmp[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)],
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_columns = F,
              heatmap_legend_param = list(title="Log2FC of \nSnC-score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(score4hp_p_fc[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)][i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()

med_imm_str <- score4hp_fc_cmp[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)]
med_imm_str.p <- score4hp_p_fc[,grep("Immune_cells|Stromal_cells",cell4DEG_fc)]

imm_str_cor_med <- cbind(cor_imm_str,med_imm_str)
imm_str_cor_med.p <- cbind(cor_imm_str.p,med_imm_str.p)

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
column_split <- c(rep("Cortex", 2), rep("Medulla", 2))
tissue.anno <- HeatmapAnnotation(empty = anno_empty(border = FALSE),
                                 foo = anno_block(gp = gpar(fill = 2:3),labels = c("Cortex","Medulla")),
                                 annotation_name_side = "left")

pdf(file="Immune_strom_heatmap_p.pdf", width=6.5, height = 5.5)
ht <- Heatmap(imm_str_cor_med,
              #col = col_fun,
              column_split = column_split,
              top_annotation = tissue.anno,
              column_title = NULL,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Log2FC of \nSnC-score",
                                          #at=c(-3,-2,-1,0,1,2,3),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(imm_str_cor_med)),
              #row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(imm_str_cor_med.p[i, j] > 0.05) {
                  grid.text("·", x, y)
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


## IPA pathway comparison
# load IPA pathway comparison z-score
library(readxl)
PW_Z <- read_excel("ov_comparison_IPA.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
colnames(PW_Z) <- gsub("_DoxvsCt.*","",colnames(PW_Z)) # truncate the column names
# load IPA pathway comparison adj.p
PW_P <- read_excel("ov_comparison_IPA_adjP.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")
colnames(PW_P) <- gsub("_DoxvsCt.*","",colnames(PW_P)) # truncate the column names

Interested.PW <- c("Pulmonary Fibrosis Idiopathic Signaling Pathway", "Pulmonary Healing Signaling Pathway", 
"Wound Healing Signaling Pathway", "Role Of Osteoclasts In Rheumatoid Arthritis Signaling Pathway", 
"Senescence Pathway", "Hepatic Fibrosis Signaling Pathway", "NF-κB Activation by Viruses", 
"Leukocyte Extravasation Signaling","Tumor Microenvironment Pathway","Ovarian Cancer Signaling", 
"Macrophage Alternative Activation Signaling Pathway", "IL-6 Signaling", "VEGF Signaling", 
"VEGF Family Ligand-Receptor Interactions", "Chemokine Signaling", "TGF-β Signaling")

# select top 50 pathways by z-score
top50 <- head(PW_Z, n = 50)

# Iterate through the Interested.PW vector
for (name in Interested.PW) {
  # Check if the name is in the row names of top50
  if (!(name %in% rownames(top50))) {
    # If not, find the corresponding row in PW_Z and append it to top50
    row_to_append <- PW_Z[rownames(PW_Z) == name, ]
    top50 <- rbind(top50, row_to_append)
  }
}

PW_Z <- top50

# convert entire datafrome to numberic
PW_Z <- dplyr::mutate_all(PW_Z, function(x) as.numeric(as.character(x)))
PW_Z[is.na(PW_Z)] <- 0

# reorder the p value matrix by the z-score matrix
PW_P <- PW_P[rownames(PW_Z),colnames(PW_Z)]

PW_Z <- as.matrix(PW_Z)
PW_P <- as.matrix(PW_P)

# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
rownames(PW_Z) <- stri_replace_all_regex(rownames(PW_Z),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
rownames(PW_P) <- stri_replace_all_regex(rownames(PW_P),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)

# draw heatmap
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
column_split <- c(rep("Cortex", 9), rep("Medulla", 7))
tissue.anno <- HeatmapAnnotation(empty = anno_empty(border = FALSE),
                                 foo = anno_block(gp = gpar(fill = 2:3),labels = c("Cortex","Medulla")),
                                 annotation_name_side = "left")

pdf(file="IPA_heatmap_p.pdf", width=10, height = 12)
ht <- Heatmap(PW_Z,
              col = col_fun,
              column_split = column_split,
              top_annotation = tissue.anno,
              column_title = NULL,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-3,-2,-1,0,1,2,3),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(PW_P[i, j] < 1.3) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht, column_title = "Pathway Comparison",column_title_gp = gpar(fontsize = 16))
dev.off()


## IPA diseases comparison
# load IPA pathway comparison z-score
library(readxl)
PW_D_Z <- read_excel("ov_comparison_IPA_diseases.xls",skip=1)
PW_D_Z <- column_to_rownames(PW_D_Z,"Diseases and Bio Functions")
colnames(PW_D_Z) <- gsub("_DoxvsCt.*","",colnames(PW_D_Z)) # truncate the column names
# load IPA pathway comparison adj.p
PW_D_P <- read_excel("ov_comparison_IPA_diseases_adjP.xls",skip=1)
PW_D_P <- column_to_rownames(PW_D_P,"Diseases and Bio Functions")
colnames(PW_D_P) <- gsub("_DoxvsCt.*","",colnames(PW_D_P)) # truncate the column names
# convert entire datafrome to numberic
PW_D_Z <- dplyr::mutate_all(PW_D_Z, function(x) as.numeric(as.character(x)))
PW_D_Z[is.na(PW_D_Z)] <- 0
# reorder the p value matrix by the z-score matrix
PW_D_P <- PW_D_P[rownames(PW_D_Z),colnames(PW_D_Z)]

PW_D_Z <- as.matrix(PW_D_Z)
PW_D_P <- as.matrix(PW_D_P)

# draw heatmap
library(circlize)
col_fun_D = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
column_split <- c(rep("Cortex", 9), rep("Medulla", 7))
tissue.annoD <- HeatmapAnnotation(empty = anno_empty(border = FALSE),
                                 foo = anno_block(gp = gpar(fill = 2:3),labels = c("Cortex","Medulla")),
                                 annotation_name_side = "left",
                                 height = unit(1.6, "cm"))

pdf(file="IPA_heatmap_diseases.pdf", width=8, height = 5)
htD <- Heatmap(PW_D_Z,
              col = col_fun_D,
              column_split = column_split,
              top_annotation = tissue.annoD,
              column_title = NULL,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-2,-1,0,1,2),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_D_Z)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(PW_D_P[i, j] < 1.3) {
                  grid.text("·", x, y)
                }})
htD<-draw(htD, column_title = "Ovarian Diseases and Bio Functions Comparison",column_title_gp = gpar(fontsize = 16))
dev.off()


# select recognized senescent genes from DEGs
# CellAge
CellAge_up <- Senescence_genesets[,"CellAge_Up"][nzchar(Senescence_genesets[,"CellAge_Up"])]
CellAge_down <- Senescence_genesets[,"CellAge_Down"][nzchar(Senescence_genesets[,"CellAge_Down"])]

#### draw dot plot for markers by Log2FC & p-value values
# cellAge up in cortex
features4plot <- CellAge_up
celltype4dot <- as.character(unique(sc_module$predicted.id))
celltype4dot<-append("Overall",celltype4dot)
lsmarker.dot<-data.frame()
for (d in celltype4dot){
  markers4dot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot<-rbind(lsmarker.dot,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot$p_val_adj[lsmarker.dot$p_val_adj==0] <- sort(unique(lsmarker.dot$p_val_adj))[2]

gs<-lsmarker.dot$X %>% unique()
lsmarker.dot$Cell_type<-factor(lsmarker.dot$Cell_type, levels = unique(lsmarker.dot$Cell_type)) # order the x axis

library(dplyr)
library(ggplot2)
library(stringr)
lsmarker.dot %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
  ggplot(aes(x=X, y = Cell_type, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('CellAge Genes Up') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("../Dot_Markers_CellAge_gene_up_cortex.pdf", width = 10, height = 5,limitsize = FALSE)

# cellAge down in cortex
features4plot <- CellAge_down
celltype4dot <- as.character(unique(sc_module$predicted.id))
celltype4dot<-append("Overall",celltype4dot)
lsmarker.dot<-data.frame()
for (d in celltype4dot){
  markers4dot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot<-rbind(lsmarker.dot,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot$p_val_adj[lsmarker.dot$p_val_adj==0] <- sort(unique(lsmarker.dot$p_val_adj))[2]

gs<-lsmarker.dot$X %>% unique()
lsmarker.dot$Cell_type<-factor(lsmarker.dot$Cell_type, levels = unique(lsmarker.dot$Cell_type)) # order the x axis

lsmarker.dot %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
  ggplot(aes(x=X, y = Cell_type, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('CellAge Genes Down') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("../Dot_Markers_CellAge_gene_down_cortex.pdf", width = 10, height = 5,limitsize = FALSE)

# cellAge up in medulla
features4plot <- CellAge_up
celltype4dot <- sc_module@meta.data %>% filter(tissue=="medulla") %>% select(predicted.id) %>% unique()
celltype4dot <- as.character(celltype4dot[,1])
celltype4dot<-append("Overall",celltype4dot)
lsmarker.dot<-data.frame()
for (d in celltype4dot[!celltype4dot %in% c("Tissue_stem_cells","Epithelial_cells")]){
  markers4dot <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot<-rbind(lsmarker.dot,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot$p_val_adj[lsmarker.dot$p_val_adj==0] <- sort(unique(lsmarker.dot$p_val_adj))[2]

gs<-lsmarker.dot$X %>% unique()
lsmarker.dot$Cell_type<-factor(lsmarker.dot$Cell_type, levels = unique(lsmarker.dot$Cell_type)) # order the x axis

library(dplyr)
library(ggplot2)
library(stringr)
lsmarker.dot %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
  ggplot(aes(x=X, y = Cell_type, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('CellAge Genes Up') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("../Dot_Markers_CellAge_gene_up_medulla.pdf", width = 6, height = 5,limitsize = FALSE)

# cellAge down in medulla
features4plot <- CellAge_down
celltype4dot <- sc_module@meta.data %>% filter(tissue=="medulla") %>% select(predicted.id) %>% unique()
celltype4dot <- as.character(celltype4dot[,1])
celltype4dot<-append("Overall",celltype4dot)
lsmarker.dot<-data.frame()
for (d in celltype4dot[!celltype4dot %in% c("Tissue_stem_cells","Epithelial_cells")]){
  markers4dot <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot<-rbind(lsmarker.dot,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot$p_val_adj[lsmarker.dot$p_val_adj==0] <- sort(unique(lsmarker.dot$p_val_adj))[2]

gs<-lsmarker.dot$X %>% unique()
lsmarker.dot$Cell_type<-factor(lsmarker.dot$Cell_type, levels = unique(lsmarker.dot$Cell_type)) # order the x axis

library(dplyr)
library(ggplot2)
library(stringr)
lsmarker.dot %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
  ggplot(aes(x=X, y = Cell_type, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('CellAge Genes Down') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("../Dot_Markers_CellAge_gene_down_medulla.pdf", width = 6, height = 5,limitsize = FALSE)



p16_pos <- WhichCells(Sample.cmb,expression = CDKN2A > 1)

Sample.cmb$p16 <- "p16_neg"
Sample.cmb@meta.data[p16_pos,"p16"] <- "p16_pos"

Sample.cmb$celltype_treatment_p16 <- paste0(Sample.cmb$celltype_treatment,"_",Sample.cmb$p16)
Idents(Sample.cmb)<-"celltype_treatment_p16"
DefaultAssay(Sample.cmb)<-"SCT"
p16pos_stromal_cor <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "Stromal_cells_cortex_ctl_p16_pos", ident.2 = "Stromal_cells_cortex_ctl_p16_neg",
               min.pct = 0, logfc.threshold = 0, test.use = "MAST",latent.vars="donor")
write.csv(p16pos_stromal_cor,"~/ovary_SC/p16pos_stromal_cor.csv")

# view the p16+ cells on UMAP
# select p16+ cells
Idents(sc_module)<-"tissue_treat"
cortex_ctl <- subset(sc_module,idents="cortex_ctl")
Idents(cortex_ctl)<-"predicted.id"
cortex_ctl_stromal <- subset(cortex_ctl,idents="Stromal_cells")
FeaturePlot(cortex_ctl_stromal, features=c("CDKN2A"))
ggsave("cortex_ctl_stromal_p16.pdf", width = 6, height = 5,limitsize = FALSE)


## Doxo Response is Cell and Region Specific
# immune cells
immune_cortex <- read.csv("Immune_cells_cor_DoxvsCt.csv")
immune_med <- read.csv("Immune_cells_med_DoxvsCt.csv")

immune_cortex_up <- immune_cortex %>% filter(avg_log2FC>0 & p_val_adj<0.05)
immune_med_up <- immune_med %>% filter(avg_log2FC>0 & p_val_adj<0.05)
intersect(immune_cortex_up$X,immune_med_up$X)

immune_cortex_down <- immune_cortex %>% filter(avg_log2FC<0 & p_val_adj<0.05)
immune_med_down <- immune_med %>% filter(avg_log2FC<0 & p_val_adj<0.05)
intersect(immune_cortex_down$X,immune_med_down$X)
# Stromal cells
Stromal_cortex <- read.csv("Stromal_cells_cor_DoxvsCt.csv")
Stromal_med <- read.csv("Stromal_cells_med_DoxvsCt.csv")

Stromal_cortex_up <- Stromal_cortex %>% filter(avg_log2FC>0 & p_val_adj<0.05)
Stromal_med_up <- Stromal_med %>% filter(avg_log2FC>0 & p_val_adj<0.05)
intersect(Stromal_cortex_up$X,Stromal_med_up$X)

Stromal_cortex_down <- Stromal_cortex %>% filter(avg_log2FC<0 & p_val_adj<0.05)
Stromal_med_down <- Stromal_med %>% filter(avg_log2FC<0 & p_val_adj<0.05)
intersect(Stromal_cortex_down$X,Stromal_med_down$X)

### senescent pathway analysis - IPA method ###
# cellAge up in cortex
CellAge <- read.csv("CellAge.csv")
SASP_EP <- read.csv("SASP_Atlas_IR_Epithelial.csv")
SASP_FI <- read.csv("SASP_Atlas_IR_Fibroblast.csv")

CellAge_u <- CellAge %>% filter(Effect>0) %>% select(Symbol)
CellAge_d <- CellAge %>% filter(Effect<0) %>% select(Symbol)

geneset_u <- CellAge_u
geneset_d <- CellAge_d

cell4DEG<-unique(Sample.cmb$celltype_refine)
cell4DEG_fc <- c(cell4DEG,"Overall")
library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))

    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
CellAge.result <- zscore.vc
CellAge.P <- pvalue.vc


SASP_EP.u <- SASP_EP %>% filter(Log2FC>0) %>% select(Symbol)
SASP_EP.d <- SASP_EP %>% filter(Log2FC<0) %>% select(Symbol)

geneset_u <- SASP_EP.u
geneset_d <- SASP_EP.d

library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))
    
    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
SASP_EP.result <- zscore.vc
SASP_EP.P <- pvalue.vc


SASP_FI.u <- SASP_FI %>% filter(Log2FC>0) %>% select(Symbol)
SASP_FI.d <- SASP_FI %>% filter(Log2FC<0) %>% select(Symbol)

geneset_u <- SASP_FI.u
geneset_d <- SASP_FI.d

library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))
    
    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
SASP_FI.result <- zscore.vc
SASP_FI.P <- pvalue.vc

z_score.df <- data.frame(Cell_Age=CellAge.result,SASP_Epithelial=SASP_EP.result,SASP_Fibroblast=SASP_FI.result)
p_value.df <- data.frame(Cell_Age=CellAge.P,SASP_Epithelial=SASP_EP.P,SASP_Fibroblast=SASP_FI.P)

col_fun = colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
#score4hp <- as.matrix(score4hp)
pdf(file="Z-score_ov_cor.pdf", width=5, height = 4)
ht <- Heatmap(z_score.df,
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox vs Control",
              cluster_rows = F,
              cluster_columns = F,
              heatmap_legend_param = list(title="Activation Score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(p_value.df[i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()

#--- senescent score for medulla ---#
geneset_u <- CellAge_u
geneset_d <- CellAge_d

cell4DEG_fc <- c(cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")],"Overall")
library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))
    
    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
CellAge.result <- zscore.vc
CellAge.P <- pvalue.vc


geneset_u <- SASP_EP.u
geneset_d <- SASP_EP.d

library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))
    
    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
SASP_EP.result <- zscore.vc
SASP_EP.P <- pvalue.vc


geneset_u <- SASP_FI.u
geneset_d <- SASP_FI.d

library(GeneOverlap)
zscore.vc <- numeric()
pvalue.vc <- numeric()
for (d in cell4DEG_fc){
  deg4score <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  if (nrow(deg4score) != 0){
    deg_u <- deg4score %>% filter(avg_log2FC>0 & p_val_adj<0.05) %>% select(X)
    deg_d <- deg4score %>% filter(avg_log2FC<0 & p_val_adj<0.05) %>% select(X)
    uu <- intersect(deg_u[,1],geneset_u[,1])
    dd <- intersect(deg_d[,1],geneset_d[,1])
    ud <- intersect(deg_u[,1],geneset_d[,1])
    du <- intersect(deg_d[,1],geneset_u[,1])
    aa <- intersect(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]))
    z_score <- ((length(uu)+length(dd))-(length(ud)+length(du)))/sqrt(length(aa))
    
    go.obj <- newGeneOverlap(c(deg_d[,1],deg_u[,1]),c(geneset_d[,1],geneset_u[,1]),spec ="hg19.gene")
    go.obj <- testGeneOverlap(go.obj)
    p_value <- go.obj@pval
    
    names(z_score) <- d
    names(p_value) <- d
    
    zscore.vc <- c(zscore.vc,z_score)
    pvalue.vc <- c(pvalue.vc,p_value)
  }
}
SASP_FI.result <- zscore.vc
SASP_FI.P <- pvalue.vc

z_score.df <- data.frame(Cell_Age=CellAge.result,SASP_Epithelial=SASP_EP.result,SASP_Fibroblast=SASP_FI.result)
p_value.df <- data.frame(Cell_Age=CellAge.P,SASP_Epithelial=SASP_EP.P,SASP_Fibroblast=SASP_FI.P)

col_fun = colorRamp2(c(-2, 0, 4), c("blue", "white", "red"))
#score4hp <- as.matrix(score4hp)
pdf(file="Z-score_ov_med.pdf", width=5, height = 4)
ht <- Heatmap(z_score.df,
              col = col_fun,
              #column_split = column_split,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox vs Control",
              cluster_rows = F,
              cluster_columns = F,
              heatmap_legend_param = list(title="Activation Score",
                                          #at=c(-0.1,0,0.2,0.4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(score4hp)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(p_value.df[i, j] > 0.05) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht)
dev.off()


# top20 genes of each cell type
library(EnhancedVolcano)
vol.plot.ls <- list()
cell4vol <- c("Overall",as.character(cell4DEG))
for (d in cell4vol){
  markers4plot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  # sort genes
  topgenes_u <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$X
  topgenes_d <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% arrange(avg_log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$X
  topgenes <- c(topgenes_u,topgenes_d)
  # Create a label column, label only selected genes
  markers4plot$label <- ifelse(markers4plot$X %in% topgenes, markers4plot$X, "")
  vol.plot <- EnhancedVolcano(markers4plot,
                              lab = markers4plot$label,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              #xlim =c(-1.5,2),
                              FCcutoff = 0.25,
                              pCutoff = 0.05,
                              labSize = 5, 
                              drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                              title = paste0(d," Dox vs Ctl"),
                              subtitle = bquote(italic("10 days")))
  vol.plot.ls[[d]] <- vol.plot
}
library(patchwork)
combined_plot <- wrap_plots(vol.plot.ls, ncol = floor(length(vol.plot.ls)/2))
print(combined_plot)
ggsave("Volcano_cor_DoxvsCt1.pdf",width=40,height=20)

vol.plot.ls <- list()
cell4vol <- c("Overall",as.character(cell4DEG[!cell4DEG %in% c("Epithelial_cells")]))
for (d in cell4vol){
  markers4plot <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  # sort genes
  topgenes_u <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$X
  topgenes_d <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% arrange(avg_log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$X
  topgenes <- c(topgenes_u,topgenes_d)
  # Create a label column, label only selected genes
  markers4plot$label <- ifelse(markers4plot$X %in% topgenes, markers4plot$X, "")
  vol.plot <- EnhancedVolcano(markers4plot,
                              lab = markers4plot$label,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              #xlim =c(-1.5,2),
                              FCcutoff = 0.25,
                              pCutoff = 0.05,
                              labSize = 5, 
                              drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                              title = paste0(d," Dox vs Ctl"),
                              subtitle = bquote(italic("10 days")))
  vol.plot.ls[[d]] <- vol.plot
}
library(patchwork)
combined_plot <- wrap_plots(vol.plot.ls, ncol = 4)
print(combined_plot)
ggsave("Volcano_med_DoxvsCt.pdf",width=40,height=20)

# top15 pathways of each cell type
# load IPA pathway comparison z-score
library(readxl)
PW_Z <- read_excel("ov_10D_IPA_compar.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
colnames(PW_Z) <- gsub("_DoxvsCt.*|_dox_vs_ct.*","",colnames(PW_Z)) # truncate the column names
# load IPA pathway comparison adj.p
PW_P <- read_excel("ov_10D_IPA_compar_adjP.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")
colnames(PW_P) <- gsub("_DoxvsCt.*|_dox_vs_ct.*","",colnames(PW_P)) # truncate the column names

PW_Z_ls <- list()
for (c in names(PW_Z)){
  PW_P_c <- PW_P %>% dplyr::select(c) %>% filter(get(c)>1.3)
  if (nrow(PW_P_c)>1){
    PW_Z_c <- PW_Z[rownames(PW_P_c),] %>% dplyr::select(c)
    PW_Z_c[,1] <- as.numeric(PW_Z_c[,1])
    PW_Z_c <- PW_Z_c %>% arrange(desc(get(c))) %>% filter(!is.na(get(c)))
    PW_Z_c_u <- head(PW_Z_c,15)
    PW_Z_c_u <- PW_Z_c_u %>% filter(get(c)>0)
    PW_Z_c_d <- tail(PW_Z_c,15)
    PW_Z_c_d <- PW_Z_c_d %>% filter(get(c)<0)
    # PW_Z_u_ls[[c]] <- PW_Z_c_u
    # PW_Z_d_ls[[c]] <- PW_Z_c_d
    
    PW_P_c_u <- PW_P_c[row.names(PW_Z_c_u),,drop=F]
    PW_P_c_d <- PW_P_c[row.names(PW_Z_c_d),,drop=F]
    
    PW_Z_c_ud <- rbind(PW_Z_c_u,PW_Z_c_d)
    PW_P_c_ud <-rbind(PW_P_c_u,PW_P_c_d)
    PW_c_ud <- cbind(PW_Z_c_ud,PW_P_c_ud)
    colnames(PW_c_ud) <- c("Z_Score","Adj_P")
    PW_Z_ls[[c]] <- PW_c_ud
    write.csv(PW_c_ud,paste0("top_",c,".csv"))
    # bar plot
    PW_c_ud$Pathway <- factor(row.names(PW_c_ud),levels=rev(row.names(PW_c_ud)))
    library(ggplot2)
    library(viridis)
    # Calculate the minimum and maximum values of Adj_P
    min_val <- min(PW_c_ud$Adj_P, na.rm = TRUE)
    max_val <- max(PW_c_ud$Adj_P, na.rm = TRUE)
    # Create an evenly spaced sequence of break points
    breaks <- seq(from = min_val, to = max_val, length.out = 5)
    library(scales)
    ggplot(PW_c_ud, aes(x=Z_Score,y=Pathway,fill=Adj_P))+
      geom_bar(stat="identity") +
      ggtitle(paste0("Top pathways in ",c)) +
      labs(x="Activation z-scores",y="Pathways",fill="-log10(adj.P)") +
      theme_bw() +
      scale_fill_viridis_c(breaks = breaks,
                           labels=label_number(accuracy=0.1))
    ggsave(paste0("top_",c,".pdf"),width = 8,height=7)
  }
}


# comparison of 6 and 10 days number of DEGs
comp_days <- data.frame(days=rep(c("6 days","10 days"),1,each=4),
           tissue=rep(c("Cortex","Medulla"),2,each=2),
           directions=rep(c("down","up"),4),
           nDEGs=c(172,242,445,354,668,276,68,197))
comp_days$days <- factor(comp_days$days, levels=c("6 days","10 days"))

ggplot(comp_days, aes(x = days, y = nDEGs, fill = directions)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = nDEGs), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25) +
  facet_wrap(~tissue, scales = "free_x") +
  labs(x = "Days", y = "Number of DEGs", fill = "Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(comp_days, aes(x = tissue, y = nDEGs, fill = directions)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = nDEGs), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25) +
  facet_wrap(~days, scales = "free_x") +
  labs(x = "Tissues", y = "Number of DEGs", fill = "Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# prep for comparision bar plot
# 10 days
cor10.u <- cor_dox.vs.ct %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% row.names()
cor10.d <- cor_dox.vs.ct %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% row.names()
med10.u <- med_dox.vs.ct %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% row.names()
med10.d <- med_dox.vs.ct %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% row.names()
# 6 days
cor_dox.vs.ct.6d <- read.csv("../Re_6days/Overall_cor_DoxvsCt.csv")
cor_dox.vs.ct.6d <- column_to_rownames(cor_dox.vs.ct.6d,"X")
med_dox.vs.ct.6d <- read.csv("../Re_6days/Overall_med_DoxvsCt.csv")
med_dox.vs.ct.6d <- column_to_rownames(med_dox.vs.ct.6d,"X")
cor6.u <- cor_dox.vs.ct.6d %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% row.names()
cor6.d <- cor_dox.vs.ct.6d %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% row.names()
med6.u <- med_dox.vs.ct.6d %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% row.names()
med6.d <- med_dox.vs.ct.6d %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% row.names()

comp_days <- data.frame(days=rep(c("6 days","10 days"),1,each=4),
                        tissue=rep(c("Cortex","Medulla"),2,each=2),
                        directions=rep(c("down","up"),4),
                        nDEGs=c(length(cor6.d),length(cor6.u),length(med6.d),length(med6.u),
                                length(cor10.d),length(cor10.u),length(med10.d),length(med10.u)))
comp_days$days <- factor(comp_days$days, levels=c("6 days","10 days"))
# Modify the data frame to reorder 'directions'
comp_days$directions <- factor(comp_days$directions, levels = c("up", "down"))
library("ggsci")
ggplot(comp_days, aes(x = tissue, y = nDEGs, fill = directions)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = nDEGs), 
            position = position_dodge(width = 0.9), 
            vjust = -0.25, size=3.5) +
  facet_wrap(~days, scales = "free_x") +
  labs(x = "Tissues", y = "Number of DEGs", fill = "Direction") +
  theme_bw() + scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Larger font for x-axis labels
        axis.title.x = element_text(size = 14), # Larger font for x-axis title
        axis.text.y = element_text(size = 12), # Adjust y-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis title size
        strip.text = element_text(size = 12), # Adjust facet label size
        strip.background = element_rect(fill = "grey90", color = "grey20", size = 0.5), # Optimized for publication
        legend.title = element_text(size = 12), # Adjust legend title size
        legend.text = element_text(size = 10), # Adjust legend text size
        plot.title = element_text(size = 16, hjust = 0.5)) # Adjust plot title size and alignment
ggsave("comp_days.pdf", width = 6, height=5)

# prep for comparison upset plot
library(UpSetR)
comp.ls <- list(Cortex_up_10days=cor10.u,Cortex_down_10days=cor10.d,
                Medulla_up_10days=med10.u,Medulla_down_10days=med10.d,
                Cortex_up_6days=cor6.u,Cortex_down_6days=cor6.d,
                Medulla_up_6days=cor6.u,Medulla_down_6days=cor6.d)
upset(fromList(comp.ls), nsets=8, order.by = "freq")

upset(fromList(comp.ls),nsets=8,
      set.metadata = list(data = comp_days, plots = list(list(type = "heat", 
                                                              column = "tissue", assign = 10, colors = c(Cortex = "green", Medulla = "navy")), 
                                                         list(type = "heat", column = "days",assign = 10))))
meatadata <- comp_days
meatadata$names <- paste0(meatadata$tissue,"_",meatadata$directions,"_",meatadata$days)
meatadata$names <- gsub(" ","",meatadata$names)
# the first column must be the sample names
meatadata <- meatadata %>%
  select(names, everything())
meatadata$tissue <- as.character(meatadata$tissue)
meatadata$days <- as.character(meatadata$days)

# select colors
pal_npg("nrc")(10)
scales::show_col(pal_npg("nrc")(10))
# upset plot
pdf(file="upsetplot.pdf",width=10, height=6)
upset(fromList(comp.ls),nsets=8, nintersects = NA,
      sets=rev(c("Cortex_up_6days","Cortex_down_6days","Medulla_up_6days","Medulla_down_6days",
            "Cortex_up_10days","Cortex_down_10days","Medulla_up_10days","Medulla_down_10days")),
      keep.order = TRUE,
      mb.ratio = c(0.7, 0.3),
      mainbar.y.label = "nDEGs Intersections", sets.x.label = "nDEGs",
      set.metadata = list(data = meatadata, plots = list(
        list(type = "heat",column = "tissue", assign = 6, colors = c(Cortex = "#F39B7FFF", Medulla = "#4DBBD5FF")),
        list(type = "heat", column = "days", assign = 5,colors = c(`6 days` = "#91D1C2FF", `10 days` = "#B09C85FF")),
        list(type = "text", column = "directions", assign = 8,colors = c(up = "#DC0000FF", down = "#3C5488FF")),
        list(type = "matrix_rows", 
             column = "days", colors = c(`6 days` = "#91D1C2FF", `10 days` = "#B09C85FF"), 
             alpha = 0.5)
        )),
        queries = list(list(query = intersects,params = list("Cortex_up_6days", "Medulla_up_6days"), color = "orange", active = T),
                       list(query = intersects,params = list("Cortex_down_6days", "Medulla_down_6days"), color = "orange", active = T),
                       list(query = intersects,params = list("Cortex_up_10days", "Medulla_up_10days"), color = "orange", active = T),
                       list(query = intersects,params = list("Cortex_down_10days", "Medulla_down_10days"), color = "orange", active = T)
             )
      )
dev.off()

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
sets <- names(movies[3:19])
avgRottenTomatoesScore <- round(runif(17, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets, avgRottenTomatoesScore))
names(metadata) <- c("sets", "avgRottenTomatoesScore")
metadata$avgRottenTomatoesScore <- as.numeric(as.character(metadata$avgRottenTomatoesScore))
Cities <- sample(c("Boston", "NYC", "LA"), 17, replace = T)
metadata <- cbind(metadata, Cities)
metadata$Cities <- as.character(metadata$Cities)
metadata[which(metadata$sets %in% c("Drama", "Comedy", "Action", "Thriller", 
                                    "Romance")), ]
upset(movies, set.metadata = list(data = metadata, plots = list(list(type = "heat", 
                                                                     column = "Cities", assign = 10, colors = c(Boston = "green", NYC = "navy", 
                                                                                                                LA = "purple")))))

                                                 
library(UpSetR)

# Assuming comp.ls is a list of sets and comp_days is a data frame with 'tissue' and 'days' columns
# Make sure the names in the colors vector match the unique values in comp_days$tissue

upset(
  fromList(comp.ls), 
  nsets = 8,
  set.metadata = list(
    data = comp_days, 
    plots = list(
      list(
        type = "heat", 
        column = "tissue", 
        assign = 8, 
        colors = c(Cortex = "green", Medulla = "navy")
      ), 
      list(
        type = "heat", 
        column = "days", 
        assign = 8
      )
    )
  )
)




# make heatmap for top20 genes
cell4hp <- c("Overall","Epithelial_cells","Stromal_1")
# get top genes list of cortex
topGenes.cor.ls <- list()
for (d in cell4hp){
  markers4plot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  # sort genes
  topgenes_u <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$X
  topgenes_d <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% arrange(avg_log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$X
  topgenes <- c(topgenes_u,topgenes_d)
  topGenes.cor.ls[[d]] <- topgenes
}

# get top genes list of medulla
topGenes.med.ls <- list()
for (d in c("Overall","Stromal_1")){
  markers4plot <- read.csv(paste0(d,"_med_DoxvsCt.csv"))
  # sort genes
  topgenes_u <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$X
  topgenes_d <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% arrange(avg_log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$X
  topgenes <- c(topgenes_u,topgenes_d)
  topGenes.med.ls[[d]] <- topgenes
}

# function to make data frame of overall top genes' expression
top20HP <- function(conditions,topgenes){
  cells4hp <- names(Sample.cmb$tissue_treat)[Sample.cmb$tissue_treat==conditions[1]]
  condition1 <- FetchData(object = Sample.cmb, vars = topgenes, cells = cells4hp, layer = "data")
  condition1 <- colMeans(condition1)
  cells4hp <- names(Sample.cmb$tissue_treat)[Sample.cmb$tissue_treat==conditions[2]]
  condition2 <- FetchData(object = Sample.cmb, vars = topgenes, cells = cells4hp, layer = "data")
  condition2 <- colMeans(condition2)
  top20HP.df <- data.frame(condition1,condition2)
  colnames(top20HP.df) <- conditions
  return(top20HP.df)
}

library(ComplexHeatmap)
library(circlize)
library(viridis)
breaks <- seq(0, 3, by = 0.5)
col_fun <- colorRamp2(breaks, rev(brewer.pal(n = 7, name = "RdYlBu")))

# top 20 in cortex_dox
conditions <- c("cortex_dox","cortex_ctl")
top20HP.df <- top20HP(conditions=conditions,topgenes=topGenes.cor.ls[["Overall"]])
pdf(file=paste0(conditions[1],"_top20_gene.pdf"), width=3, height = 7)
Heatmap(top20HP.df,
        col = col_fun,
        border=T,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = F,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 11),
        heatmap_legend_param = list(title="Gene \nExpression")
)
dev.off()

# top 20 in medulla_dox
conditions <- c("medulla_dox","medulla_ctl")
top20HP.df <- top20HP(conditions=conditions,topgenes=topGenes.med.ls[["Overall"]])
pdf(file=paste0(conditions[1],"_top20_gene.pdf"), width=3, height = 7)
Heatmap(top20HP.df,
        col = col_fun,
        border=T,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = F,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 11),
        heatmap_legend_param = list(title="Gene \nExpression")
)
dev.off()


# function to make data frame of top genes' expression in each cell types
top20HP.ct <- function(conditions,topgenes,cell_type){
  cells4hp <- Sample.cmb@meta.data %>% filter(celltype_refine==cell_type & tissue_treat==conditions[1]) %>% row.names()
  condition1 <- FetchData(object = Sample.cmb, vars = topgenes, cells = cells4hp, layer = "data")
  condition1 <- colMeans(condition1)
  cells4hp <- Sample.cmb@meta.data %>% filter(celltype_refine==cell_type & tissue_treat==conditions[2]) %>% row.names()
  condition2 <- FetchData(object = Sample.cmb, vars = topgenes, cells = cells4hp, layer = "data")
  condition2 <- colMeans(condition2)
  top20HP.df <- data.frame(condition1,condition2)
  colnames(top20HP.df) <- conditions
  return(top20HP.df)
}

# top genes in cortex in cell types
conditions <- c("cortex_dox","cortex_ctl")
for (d in c("Epithelial_cells","Stroma_1")){
  #markers4plot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  top20HP.df <- top20HP.ct(conditions=conditions,topgenes=topGenes.cor.ls[[d]],cell_type=d)
  
  pdf(file=paste0(conditions[1],"_",d,"_top20_gene.pdf"), width=3, height = 7)
  Heatmap(top20HP.df,
          col = col_fun,
          border=T,
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_rows = F,
          cluster_columns = F,
          row_names_gp = gpar(fontsize = 11),
          heatmap_legend_param = list(title="Gene \nExpression")
  )
  dev.off()
}

# top genes in medulla in cell types
conditions <- c("medulla_dox","medulla_ctl")
  #markers4plot <- read.csv(paste0(d,"_cor_DoxvsCt.csv"))
  top20HP.df <- top20HP.ct(conditions=conditions,topgenes=topGenes.cor.ls[["Stromal_1"]],cell_type="Stroma_1")
  
  pdf(file=paste0(conditions[1],"_",d,"_top20_gene.pdf"), width=3, height = 7)
  Heatmap(top20HP.df,
          col = col_fun,
          border=T,
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_rows = F,
          cluster_columns = F,
          row_names_gp = gpar(fontsize = 11),
          heatmap_legend_param = list(title="Gene \nExpression")
  )
  dev.off()

  
# RRHO
# get the DEG full list
Idents(Sample.cmb) <- "tissue_treat"
# cortex overall DEGs
cor_dox.vs.ct.all <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "cortex_dox", ident.2 = "cortex_ctl",
                               min.pct = 0, logfc.threshold = 0, test.use = "MAST",latent.vars="donor")
write.csv(cor_dox.vs.ct.all,"Overall_DEGs_cor_dox_vs_ct_all.csv")
# medulla overall DEGs
med_dox.vs.ct.all <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "medulla_dox", ident.2 = "medulla_ctl",
                               min.pct = 0, logfc.threshold = 0, test.use = "MAST",latent.vars="donor")
write.csv(med_dox.vs.ct.all,"Overall_DEGs_med_dox_vs_ct_all.csv")


library(tibble)
library(RRHO2)
# Create "gene" lists:
  dox_vs_ctl_cor.all.rank <- cor_dox.vs.ct.all %>% 
    select(avg_log2FC,p_val_adj) %>%
    mutate(p_val_adj=ifelse(p_val_adj==0,1e-305,p_val_adj)) %>%
    rownames_to_column("gene_symbol")
  dox_vs_ctl_cor.all.rank$weight <- dox_vs_ctl_cor.all.rank$avg_log2FC*-log10((dox_vs_ctl_cor.all.rank$p_val_adj))
  dox_vs_ctl_cor.all.rank$rank <- order(dox_vs_ctl_cor.all.rank$weight,decreasing=T)
  dox_vs_ctl_cor.all.rank <- dox_vs_ctl_cor.all.rank %>% select(gene_symbol,weight)
  
  dox_vs_ctl_med.all.rank <- med_dox.vs.ct.all %>% 
    select(avg_log2FC,p_val_adj) %>%
    mutate(p_val_adj=ifelse(p_val_adj==0,1e-305,p_val_adj)) %>%
    rownames_to_column("gene_symbol")
  dox_vs_ctl_med.all.rank$weight <- dox_vs_ctl_med.all.rank$avg_log2FC*-log10((dox_vs_ctl_med.all.rank$p_val_adj))
  dox_vs_ctl_med.all.rank$rank <- order(dox_vs_ctl_med.all.rank$weight,decreasing=T)
  dox_vs_ctl_med.all.rank <- dox_vs_ctl_med.all.rank %>% select(gene_symbol,weight)

  # Compute overlap and significance
  RRHO_obj <-  RRHO2_initialize(dox_vs_ctl_cor.all.rank, dox_vs_ctl_med.all.rank, labels = c("Dox_vs_Ctl in Cortex", "Dox_vs_Ctl in medulla"), log10.ind=TRUE)
  
  pdf(file="RRHO2_dox_vs_ctl_CorvsMed.pdf", width=7, height = 6)
  RRHO2_heatmap(RRHO_obj)
  dev.off()
  
  pdf(file="RRHO2_vennDiagra_DD_dox_vs_ctl_CorvsMed.pdf", width=6, height = 5)
  RRHO2_vennDiagram(RRHO_obj,type="dd")
  dev.off()
  
  pdf(file="RRHO2_vennDiagra_UU_dox_vs_ctl_CorvsMed.pdf", width=6, height = 5)
  RRHO2_vennDiagram(RRHO_obj,type="uu")
  dev.off()
