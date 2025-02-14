# compare the senescence score changes between 6 and 10 days.
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

setwd("~/ovary_SC/6days")

# import 10X results
ov_sample_list <- c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox",
                    "RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-ct-1","RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2","RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1",
                    "RTL-857-Medulla-DOXO-1","RTL-857-Medulla-CTRL-1")
for (s in ov_sample_list){print(s)}

Sample.list.cmb <- list()
for (s in ov_sample_list){
  if (s %in% c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox")){
    ov_sample <- Read10X(data.dir = paste0("/opt/home/buckcenter.org/fwu/ovary_SC/",s,"/outs/filtered_feature_bc_matrix"))
  } else if (s %in% c("RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-ct-1",
                 "RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2","RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1")){
    ov_sample <- Read10X(data.dir = paste0("/data/array2/fwu/ovary_SC_2nd/",s,"/outs/filtered_feature_bc_matrix"))
  } else {
    ov_sample <- Read10X(data.dir = paste0("/data/array2/fwu/ovary_SC_3rd/",s,"/outs/filtered_feature_bc_matrix"))
  }
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

# transfer previous manual annotations.
ov_sub2 <- readRDS("ov_sub2.rds")
Idents(ov_sub2)<-"sub.cluster"
# Identify anchors between the reference (annotated) and new Seurat objects:
anchors_list <- list()
for (i in seq_along(Sample.list.cmb)) {
  anchors_list[[i]] <- FindTransferAnchors(reference = ov_sub2, query = Sample.list.cmb[[i]], dims = 1:30,
                                           reference.assay="integrated",
                                           normalization.method = "SCT")
}
# Transfer the cell type annotations from the reference object to the new objects:
for (i in seq_along(Sample.list.cmb)) {
  predictions <- TransferData(anchorset = anchors_list[[i]], refdata = ov_sub2$celltype_refine2)
  Sample.list.cmb[[i]] <- AddMetaData(Sample.list.cmb[[i]],metadata = predictions)
}

#Sample.list.cmb[["fist_batch"]] <- NULL
features.cmb <- SelectIntegrationFeatures(object.list = Sample.list.cmb, nfeatures = 3000)
# run PCA before RPCA integration
Sample.list.cmb <- PrepSCTIntegration(object.list = Sample.list.cmb, anchor.features = features.cmb)
Sample.list.cmb <- lapply(X = Sample.list.cmb, FUN = RunPCA, features = features.cmb)

# find the anchor by using reference
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
                                 c("ov-cor-ct","ov-cor-dox","RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1") ~ "cortex",
                               .default = "medulla")

Sample.cmb$batches <- case_when(Sample.cmb$orig.ident %in%
                                  c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox") ~ "batch1",
                                Sample.cmb$orig.ident %in%
                                  c("RTL-857-Medulla-CTRL-1","RTL-857-Medulla-DOXO-1") ~ "batch3",
                                .default = "batch2")

Sample.cmb$donor <- case_when(Sample.cmb$orig.ident %in%
                                c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox") ~ "donor1",
                              Sample.cmb$orig.ident %in%
                                c("RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-ct-1","RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2") ~ "donor2",
                              Sample.cmb$orig.ident %in%
                                c("RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1","RTL-857-Medulla-CTRL-1","RTL-857-Medulla-DOXO-1") ~ "donor3")

Sample.cmb$treatment <- case_when(Sample.cmb$orig.ident %in%
                                    c("ov-cor-dox","ov-med-dox","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2","RTL_857_ov-cor-dox-1","RTL-857-Medulla-DOXO-1") ~ "dox",
                                  .default = "ctl")

Sample.cmb$tissue_treat <- paste0(Sample.cmb$tissue,"_",Sample.cmb$treatment)

DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T,ncol=4)
ggsave("ov_clusters_combined.pdf",width=14,height=8)

# find markers for every cluster compared to all remaining cells
markers.all <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(markers.all,"markers_all_clusters.csv")

Idents(Sample.cmb) <- "predicted.id"
markers.all.cellType <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(markers.all.cellType,"markers_all_cellTypes.csv")

# find overall DEGs dox vs ctl for tissues
Idents(Sample.cmb) <- "tissue_treat"
# cortex overall DEGs
cor_dox.vs.ct <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "cortex_dox", ident.2 = "cortex_ctl",
                             min.pct = 0.1, logfc.threshold = 0.2, test.use = "MAST",latent.vars="donor")
write.csv(cor_dox.vs.ct,"Overall_DEGs_cor_dox_vs_ct.csv")
# medulla overall DEGs
med_dox.vs.ct <- FindMarkers(Sample.cmb, assay = "SCT", ident.1 = "medulla_dox", ident.2 = "medulla_ctl",
                             min.pct = 0.1, logfc.threshold = 0.2, test.use = "MAST",latent.vars="donor")
write.csv(med_dox.vs.ct,"Overall_DEGs_med_dox_vs_ct.csv")

### --- cell type annotation: refine --- ###
# # quick check cell annotation results for a cluster
# table(Sample.cmb@meta.data %>%
#         filter(predicted.id=="MSC") %>%
#         select(seurat_clusters)) %>% sort(decreasing = T)
# table(Sample.cmb@meta.data %>%
#         filter(seurat_clusters==14) %>%
#         select(predicted.id)) %>% sort(decreasing = T)
cell_type_counts <- list()
for (c in unique(Sample.cmb$seurat_clusters)){
  cell_type <- Sample.cmb@meta.data %>%
    filter(seurat_clusters==c) %>%
    select(predicted.id)
  # Creating a sample frequency table
  freq_table <- table(cell_type)
  # Convert to a data frame
  freq_table_df <- as.data.frame(freq_table)
  # Rename the columns for clarity
  names(freq_table_df) <- c("Cell_type", "Frequency")
  freq_table_df$Cluster <- c
  # Sort in descending order
  sorted_df <- freq_table_df[order(-freq_table_df$Frequency), ]
  cell_type_counts[[c]] <- sorted_df
}
cell_type_counts <- purrr::reduce(cell_type_counts,rbind)
write.csv(cell_type_counts,"cell_type_counts_transferred_from_batch1.csv")
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
  Q1 <- quantile(celltype.df[,cs], 0.25)
  Q3 <- quantile(celltype.df[,cs], 0.75)
  IQR <- Q3 - Q1
  box_plot_cs <- ggplot(celltype.df, aes_string(x = "seurat_clusters", y = cs)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(Q1 - 1.5 * IQR, Q3 + 1.5 * IQR) + 
    xlab("clusters") + 
    ylab(paste0("cell scores of ", cs)) + 
    ggtitle(paste0("scores of ", cs))
  box_ls[[cs]] <- box_plot_cs
}
library(gridExtra)
grid_plot<-grid.arrange(grobs = box_ls, ncol = 2)
ggsave(paste0("Box Plot of celltype scores.pdf"),plot = grid_plot,width=12,height=14)

# annotate cell type by combining results from HPCA and manually curate the gene list from ovary publications
Sample.cmb$celltype_refine <- case_when(Sample.cmb$seurat_clusters ==
                                          "0" ~ "Stroma_1",
                                        Sample.cmb$seurat_clusters ==
                                          "1" ~ "Stroma_1",
                                        Sample.cmb$seurat_clusters ==
                                          "2" ~ "Stroma_2",
                                        Sample.cmb$seurat_clusters ==
                                          "3" ~ "Stroma_1",
                                        Sample.cmb$seurat_clusters ==
                                          "4" ~ "Perivascular_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "5" ~ "Stroma_1",
                                        Sample.cmb$seurat_clusters ==
                                          "6" ~ "Stroma_2",
                                        Sample.cmb$seurat_clusters ==
                                          "7" ~ "Stroma_1",
                                        Sample.cmb$seurat_clusters ==
                                          "8" ~ "Endothelial_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "9" ~ "Epithelial_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "10" ~ "Smooth_muscle_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "11" ~ "Perivascular_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "12" ~ "Unknown",
                                        Sample.cmb$seurat_clusters ==
                                          "13" ~ "Smooth_muscle_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "14" ~ "Smooth_muscle_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "15" ~ "Immune_cells",
                                        Sample.cmb$seurat_clusters ==
                                          "16" ~ "Endothelial_cells")
# set the levels in order
Sample.cmb$celltype_refine <- factor(Sample.cmb$celltype_refine,
                                     levels=sort(unique(Sample.cmb$celltype_refine)))
# plot UMAP based on refined cell types
Idents(Sample.cmb) <- "celltype_refine"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue_treat",label =T,repel =T)
ggsave("ov_cellTypes_refined.pdf",width=12,height=6)

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
  geom_errorbar(aes(ymin = Freq_mean - se, ymax = Freq_mean + se), position = position_dodge(0.9), alpha=0.6) +
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
             position = position_dodge(0.9), size = 0.6, alpha=0.6)
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
#plan("multisession", workers = 4)
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

# score the senecent cells
library(AUCell)
# load the senescence gene sets
Senescence_genesets <- read.csv("/opt/home/buckcenter.org/fwu/ovary_SC/Senescence_genesets_10122023.csv")

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
    Sample.cmb@assays[["RNA"]]$counts,geneset_features)
  #sc_scores[[names(Senescence_genesets)[sc]]] <- cells_AUC
}

#cells_AUC@assays@data@listData[["AUC"]]
#colnames(cells_AUC@assays@data@listData[["AUC"]])

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
col_fun = colorRamp2(c(0, 0.15), c("white", "red"))

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
col_fun = colorRamp2(c(0, 0.15), c("white", "red"))

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
col_fun = colorRamp2(c(0, 1), c("white", "red"))

pdf(file="SnC-score_ov_cor_sc_overall_genes.pdf", width=12, height = 5)
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
col_fun = colorRamp2(c(0, 1), c("white", "red"))

pdf(file="SnC-score_ov_med_sc_overall_genes.pdf", width=10, height = 5)
ht <- Heatmap(overlapped_SC_DEGs.df,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox and Control",
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

# sc_module$celltype_refine <- gsub("Stromal_cells","Stroma_1",sc_module$celltype_refine)
# sc_module$celltype_refine <- gsub("Theca_Stroma","Stroma_2",sc_module$celltype_refine)
# sc_module$celltype_treatment <- paste(sc_module$celltype_refine, sc_module$tissue_treat,
#                                        sep = "_")
# Idents(sc_module)<-"celltype_treatment"
# add cell type information for the senescent gene comparison
# cortex
library(tibble)
sc_score.ls<-list()
sc_score_p.ls<-list()
for (cell in cell4DEG){
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
column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9))

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

colnames(score4hp) <- gsub("Stromal_cells","Stroma_1",colnames(score4hp))
colnames(score4hp) <- gsub("Theca_Stroma","Stroma_2",colnames(score4hp))

# log2FC heatmap
cell4DEG<-unique(Sample.cmb$celltype_refine)
cell4DEG_fc <- c(as.character(cell4DEG),"Overall")
score4hp_fc <- list()
for (fc in cell4DEG_fc){
  score4hp_fcc <- log2(score4hp[,paste0(fc,"_Dox")]/score4hp[,paste0(fc,"_Ctl")])
  score4hp_fc[[fc]] <- score4hp_fcc
}
score4hp_fc_cmp <- purrr::reduce(score4hp_fc,cbind)
colnames(score4hp_fc_cmp) <- names(score4hp_fc)
score4hp_p_fc <- score4hp_p[seq(1,length(cell4DEG_fc)*2,2)]
score4hp_fc_cmp.cor <- score4hp_fc_cmp
score4hp_p_fc.cor <- score4hp_p_fc

colnames(score4hp_fc_cmp)
score4hp_fc_cmp <- score4hp_fc_cmp[,c("Endothelial_cells","Epithelial_cells","Immune_cells",
                   "Perivascular_cells","Smooth_muscle_cells",
                   "Stroma_1","Stroma_2","Unknown",
                   "Overall")]
colnames(score4hp_p_fc) <- gsub("_Dox.*","",colnames(score4hp_p_fc))
score4hp_p_fc <- score4hp_p_fc[,colnames(score4hp_fc_cmp)]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file="SnC-score_ov_cor_cmb_fc.pdf", width=7, height = 5)
ht <- Heatmap(score4hp_fc_cmp,
              col = col_fun,
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

score4hp_fc_cmp.cor <- score4hp_fc_cmp
score4hp_p_fc.cor <- score4hp_p_fc

# for medula
sc_score.ls<-list()
sc_score_p.ls<-list()
for (cell in cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]){
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
sc_score <- sc_score.ls %>% reduce(left_join, by = "SC_gene_set")
sc_score_p <- sc_score_p.ls %>% reduce(left_join, by = "SC_gene_set")

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
pdf(file="SnC-score_ov_med_cmb.pdf", width=10, height = 5)
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
cell4DEG_fc <- c(as.character(cell4DEG)[!cell4DEG %in% c("Epithelial_cells")],"Overall")
score4hp_fc <- list()
for (fc in cell4DEG_fc){
  score4hp_fcc <- log2(score4hp[,paste0(fc,"_Dox")]/score4hp[,paste0(fc,"_Ctl")])
  score4hp_fc[[fc]] <- score4hp_fcc
}
score4hp_fc_cmp <- reduce(score4hp_fc,cbind)
colnames(score4hp_fc_cmp) <- names(score4hp_fc)

score4hp_p_fc <- score4hp_p[seq(1,length(cell4DEG_fc)*2,2)]


colnames(score4hp_fc_cmp)
score4hp_fc_cmp <- score4hp_fc_cmp[,c("Endothelial_cells","Immune_cells",
                                      "Perivascular_cells","Smooth_muscle_cells",
                                      "Stroma_1","Stroma_2","Unknown",
                                      "Overall")]
colnames(score4hp_p_fc) <- gsub("_Dox.*","",colnames(score4hp_p_fc))
score4hp_p_fc <- score4hp_p_fc[,colnames(score4hp_p_fc)]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file="SnC-score_ov_med_cmb_fc.pdf", width=7, height = 5)
ht <- Heatmap(score4hp_fc_cmp,
              col = col_fun,
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

score4hp_fc_cmp.med <- score4hp_fc_cmp
score4hp_p_fc.med <- score4hp_p_fc

## IPA pathway comparison
# load IPA pathway comparison z-score
library(readxl)
PW_Z <- read_excel("Ov_6days.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
colnames(PW_Z) <- gsub("_DoxvsCt.*","",colnames(PW_Z)) # truncate the column names
# load IPA pathway comparison adj.p
PW_P <- read_excel("Ov_6days_P.xls",skip=1)
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
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
column_split <- c(rep("Cortex", 9), rep("Medulla", 8))
tissue.anno <- HeatmapAnnotation(empty = anno_empty(border = FALSE),
                                 foo = anno_block(gp = gpar(fill = 2:3),labels = c("Cortex","Medulla")),
                                 annotation_name_side = "left")

pdf(file="IPA_heatmap_p.pdf", width=10, height = 13)
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
                                          at=c(-6,-4,-2,0,2,4,6),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(PW_P[i, j] < 1.3) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht, column_title = "Pathway Comparison",column_title_gp = gpar(fontsize = 16))
dev.off()


# 6 days vs 10 days
# cortex
PW_Z <- read_excel("Ov_cor_6vs10days.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
PW_Z <- PW_Z[colnames(PW_Z)!="SC_only_Overall_DEGs_cor_dox_vs_ct - 2023-09-01 10:02 AM"]
colnames(PW_Z) <- gsub("- 2023-11-12.*","6Days",colnames(PW_Z)) # truncate the column names
colnames(PW_Z) <- gsub("_cor_DoxvsCt 6Days","_6Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("- 2023-09-27.*","10Days",colnames(PW_Z)) # truncate the column names
colnames(PW_Z) <- gsub("- 2023-08-24.*","10Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("_cor_DoxvsCt 10Days","_10Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("Overall_DEGs_cor_dox_vs_ct 10Days","Overall_10Days",colnames(PW_Z))
PW_Z <- PW_Z[c("Overall_6Days","Overall_10Days","Endothelial_cells_6Days","Endothelial_cells_10Days",
               "Epithelial_cells_6Days","Epithelial_cells_10Days","Immune_cells_6Days","Immune_cells_10Days",
               "Perivascular_cells_6Days","Perivascular_cells_10Days","Smooth_muscle_cells_6Days","Smooth_muscle_cells_10Days",
               "Stromal_cells_6Days","Stromal_cells_10Days","Theca_Stroma_6Days","Theca_Stroma_10Days")] # reorder colnames
# load IPA pathway comparison adj.p
PW_P <- read_excel("Ov_cor_6vs10days_P.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")
PW_P <- PW_P[colnames(PW_P)!="SC_only_Overall_DEGs_cor_dox_vs_ct - 2023-09-01 10:02 AM"]
colnames(PW_P) <- gsub("- 2023-11-12.*","6Days",colnames(PW_P)) # truncate the column names
colnames(PW_P) <- gsub("_cor_DoxvsCt 6Days","_6Days",colnames(PW_P))
colnames(PW_P) <- gsub("- 2023-09-27.*","10Days",colnames(PW_P)) # truncate the column names
colnames(PW_P) <- gsub("- 2023-08-24.*","10Days",colnames(PW_P))
colnames(PW_P) <- gsub("_cor_DoxvsCt 10Days","_10Days",colnames(PW_P))
colnames(PW_P) <- gsub("Overall_DEGs_cor_dox_vs_ct 10Days","Overall_10Days",colnames(PW_P))

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
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
# tissue.anno <- HeatmapAnnotation(empty = anno_empty(border = FALSE),
#                                  foo = anno_block(gp = gpar(fill = 2:9),
#                                                   labels = c("Overall","Endothelial","Epithelial","Immune",
#                                                                                    "Perivascular","Smooth_muscle","Stromal","Theca_Stroma")),
#                                  annotation_name_side = "left")

pdf(file="IPA_heatmap_Cor_6Dvs10D.pdf", width=10, height = 13)
ht <- Heatmap(PW_Z,
              col = col_fun,
              column_split = column_split,
              #top_annotation = tissue.anno,
              column_title = NULL,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-6,-4,-2,0,2,4,6),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(PW_P[i, j] < 1.3) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht, column_title = "Pathway Comparison",column_title_gp = gpar(fontsize = 16))
dev.off()

# 6 days vs 10 days
# medulla
PW_Z <- read_excel("Ov_med_6vs10days.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
colnames(PW_Z) <- gsub("- 2023-11-12.*","6Days",colnames(PW_Z)) # truncate the column names
colnames(PW_Z) <- gsub("_med_DoxvsCt 6Days","_6Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("- 2023-09-27.*","10Days",colnames(PW_Z)) # truncate the column names
colnames(PW_Z) <- gsub("- 2023-08-24.*","10Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("_med_DoxvsCt 10Days","_10Days",colnames(PW_Z))
colnames(PW_Z) <- gsub("Overall_DEGs_med_dox_vs_ct 10Days","Overall_10Days",colnames(PW_Z))
PW_Z <- PW_Z[c("Overall_6Days","Overall_10Days","Endothelial_cells_6Days","Endothelial_cells_10Days",
               "Immune_cells_6Days","Immune_cells_10Days",
               "Perivascular_cells_6Days","Perivascular_cells_10Days","Smooth_muscle_cells_6Days","Smooth_muscle_cells_10Days",
               "Stromal_cells_6Days","Stromal_cells_10Days","Theca_Stroma_6Days","Theca_Stroma_10Days")] # reorder colnames
# load IPA pathway comparison adj.p
PW_P <- read_excel("Ov_med_6vs10days_P.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")
colnames(PW_P) <- gsub("- 2023-11-12.*","6Days",colnames(PW_P)) # truncate the column names
colnames(PW_P) <- gsub("_med_DoxvsCt 6Days","_6Days",colnames(PW_P))
colnames(PW_P) <- gsub("- 2023-09-27.*","10Days",colnames(PW_P)) # truncate the column names
colnames(PW_P) <- gsub("- 2023-08-24.*","10Days",colnames(PW_P))
colnames(PW_P) <- gsub("_med_DoxvsCt 10Days","_10Days",colnames(PW_P))
colnames(PW_P) <- gsub("Overall_DEGs_med_dox_vs_ct 10Days","Overall_10Days",colnames(PW_P))

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
col_fun = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
column_split <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7))

pdf(file="IPA_heatmap_med_6Dvs10D.pdf", width=10, height = 13)
ht <- Heatmap(PW_Z,
              col = col_fun,
              column_split = column_split,
              #top_annotation = tissue.anno,
              column_title = NULL,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-6,-4,-2,0,2,4,6),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5),
              cell_fun = function(j, i, x, y, w, h, f) {
                if(PW_P[i, j] < 1.3) {
                  grid.text("·", x, y)
                }})
ht<-draw(ht, column_title = "Pathway Comparison",column_title_gp = gpar(fontsize = 16))
dev.off()

# compare senescent score changes 6days vs 10days
# cortex
score4hp_fc_cmp <- read.csv("score4hp_fc_cmp_cor_10D.csv")
score4hp_fc_cmp <- column_to_rownames(score4hp_fc_cmp,"X")
colnames(score4hp_fc_cmp) <- paste0(colnames(score4hp_fc_cmp),"_10D")
score4hp_fc_cmp.p <- read.csv("score4hp_p_fc_cor_10D.csv")
score4hp_fc_cmp.p <- column_to_rownames(score4hp_fc_cmp.p,"X")
colnames(score4hp_fc_cmp.p) <- paste0(colnames(score4hp_fc_cmp.p),"_10D")

colnames(score4hp_fc_cmp.cor) <- paste0(colnames(score4hp_fc_cmp.cor),"_6D")
colnames(score4hp_p_fc.cor) <- paste0(colnames(score4hp_p_fc.cor),"_6D")

comp_6vs10 <- cbind(score4hp_fc_cmp.cor,score4hp_fc_cmp)
comp_6vs10.p <- cbind(score4hp_p_fc.cor,score4hp_fc_cmp.p)
colnames(comp_6vs10.p) <- gsub("Dox_p_value_|Dox_","",colnames(comp_6vs10.p))

colnames(comp_6vs10) <- gsub("Stromal_cells","Stroma_1",colnames(comp_6vs10))
colnames(comp_6vs10) <- gsub("Theca_Stroma","Stroma_2",colnames(comp_6vs10))

colnames(comp_6vs10.p) <- gsub("Stromal_cells","Stroma_1",colnames(comp_6vs10.p))
colnames(comp_6vs10.p) <- gsub("Theca_Stroma","Stroma_2",colnames(comp_6vs10.p))

comp_6vs10 <- comp_6vs10[c("Overall_6D","Overall_10D","Endothelial_cells_6D","Endothelial_cells_10D",
               "Epithelial_cells_6D","Epithelial_cells_10D","Immune_cells_6D","Immune_cells_10D",
               "Perivascular_cells_6D","Perivascular_cells_10D","Smooth_muscle_cells_6D","Smooth_muscle_cells_10D",
               "Stroma_1_6D","Stroma_1_10D","Stroma_2_6D","Stroma_2_10D")] # reorder colnames
# reorder the p value matrix by the z-score matrix
comp_6vs10.p <- comp_6vs10.p[rownames(comp_6vs10),colnames(comp_6vs10)]

celltypes.6D <- gsub("_6D","",grep("_6D",colnames(comp_6vs10),value=T)) # cell types in 6D
celltypes.10D <- gsub("_10D","",grep("_10D",colnames(comp_6vs10),value=T)) # cell types in 10D
celltypes <- intersect(celltypes.6D,celltypes.10D)

library(ggsci)
for (c in celltypes){
  sc0d <- data.frame(SnC=rownames(comp_6vs10),log2FC=0,p_value=0,days=0)
  sc6d <- data.frame(SnC=rownames(comp_6vs10),log2FC=comp_6vs10[,paste0(c,"_6D")],p_value=comp_6vs10.p[,paste0(c,"_6D")],days=6)
  sc10d <- data.frame(SnC=rownames(comp_6vs10),log2FC=comp_6vs10[,paste0(c,"_10D")],p_value=comp_6vs10.p[,paste0(c,"_10D")],days=10)
  scdays <- rbind(sc0d,sc6d,sc10d)
  scdays$logP <- -log10(scdays$p_value)
  scdays$logP[is.infinite(scdays$logP)] <- 100
  scdays$logP[is.na(scdays$logP)] <- 100
  ggplot(scdays, aes(x = days, y = log2FC, group = SnC, color = SnC)) +
    geom_line() +
    geom_point(data=scdays, aes(size=ifelse(logP>=100, 100, logP))) +
    scale_x_continuous(breaks = c(0, 6, 10)) +
    theme_minimal() + 
    scale_color_igv() +
    scale_size_continuous(name = "-log10(p.value)", 
                          breaks = c(1.3, 25, 50, 100), 
                          labels = c("1.3", "25", "50", "100+")) +
    labs(title = paste0("Senescence of cortex ",c," over Days"),
         color = "SnC gene sets",
         x = "Days",
         y = "log2FC")
  ggsave(paste0("CompareDays_Cortex_",c,".pdf"),width = 7,height=6)
}

# calculate paired t-test of overall senescence score changes
c="Overall"
sc6d <- data.frame(SnC=rownames(comp_6vs10),log2FC_6D=comp_6vs10[,paste0(c,"_6D")],p_value_6D=comp_6vs10.p[,paste0(c,"_6D")])
sc10d <- data.frame(SnC=rownames(comp_6vs10),log2FC_10D=comp_6vs10[,paste0(c,"_10D")],p_value_10D=comp_6vs10.p[,paste0(c,"_10D")])
scdays <- merge(sc6d,sc10d,by="SnC")
Overall_SC_6vs10D.cor <- t.test(scdays$log2FC_6D,scdays$log2FC_10D, paired=T)

scdays.cor <- scdays



# medulla
score4hp_fc_cmp <- read.csv("score4hp_fc_cmp_med_10D.csv")
score4hp_fc_cmp <- column_to_rownames(score4hp_fc_cmp,"X")
colnames(score4hp_fc_cmp) <- paste0(colnames(score4hp_fc_cmp),"_10D")
score4hp_fc_cmp.p <- read.csv("score4hp_p_fc_med_10D.csv")
score4hp_fc_cmp.p <- column_to_rownames(score4hp_fc_cmp.p,"X")
colnames(score4hp_fc_cmp.p) <- paste0(colnames(score4hp_fc_cmp.p),"_10D")

colnames(score4hp_fc_cmp.med) <- paste0(colnames(score4hp_fc_cmp.med),"_6D")
colnames(score4hp_p_fc.med) <- paste0(colnames(score4hp_p_fc.med),"_6D")

comp_6vs10 <- cbind(score4hp_fc_cmp.med,score4hp_fc_cmp)
comp_6vs10.p <- cbind(score4hp_p_fc.med,score4hp_fc_cmp.p)
colnames(comp_6vs10.p) <- gsub("Dox_p_value_|Dox_","",colnames(comp_6vs10.p))

colnames(comp_6vs10) <- gsub("Stromal_cells","Stroma_1",colnames(comp_6vs10))
colnames(comp_6vs10) <- gsub("Theca_Stroma","Stroma_2",colnames(comp_6vs10))

colnames(comp_6vs10.p) <- gsub("Stromal_cells","Stroma_1",colnames(comp_6vs10.p))
colnames(comp_6vs10.p) <- gsub("Theca_Stroma","Stroma_2",colnames(comp_6vs10.p))

comp_6vs10 <- comp_6vs10[c("Overall_6D","Overall_10D","Endothelial_cells_6D","Endothelial_cells_10D",
                           "Immune_cells_6D","Immune_cells_10D",
                           "Perivascular_cells_6D","Perivascular_cells_10D","Smooth_muscle_cells_6D","Smooth_muscle_cells_10D",
                           "Stroma_1_6D","Stroma_1_10D","Stroma_2_6D","Stroma_2_10D")] # reorder colnames
# reorder the p value matrix by the z-score matrix
comp_6vs10.p <- comp_6vs10.p[rownames(comp_6vs10),colnames(comp_6vs10)]

celltypes.6D <- gsub("_6D","",grep("_6D",colnames(comp_6vs10),value=T)) # cell types in 6D
celltypes.10D <- gsub("_10D","",grep("_10D",colnames(comp_6vs10),value=T)) # cell types in 10D
celltypes <- intersect(celltypes.6D,celltypes.10D)

library(ggsci)
for (c in celltypes){
  sc0d <- data.frame(SnC=rownames(comp_6vs10),log2FC=0,p_value=0,days=0)
  sc6d <- data.frame(SnC=rownames(comp_6vs10),log2FC=comp_6vs10[,paste0(c,"_6D")],p_value=comp_6vs10.p[,paste0(c,"_6D")],days=6)
  sc10d <- data.frame(SnC=rownames(comp_6vs10),log2FC=comp_6vs10[,paste0(c,"_10D")],p_value=comp_6vs10.p[,paste0(c,"_10D")],days=10)
  scdays <- rbind(sc0d,sc6d,sc10d)
  scdays$logP <- -log10(scdays$p_value)
  scdays$logP[is.infinite(scdays$logP)] <- 100
  scdays$logP[is.na(scdays$logP)] <- 100
  ggplot(scdays, aes(x = days, y = log2FC, group = SnC, color = SnC)) +
    geom_line() +
    geom_point(data=scdays, aes(size=ifelse(logP>=100, 100, logP))) +
    scale_x_continuous(breaks = c(0, 6, 10)) +
    theme_minimal() + 
    scale_color_igv() +
    scale_size_continuous(name = "-log10(p.value)", 
                          breaks = c(1.3, 25, 50, 100), 
                          labels = c("1.3", "25", "50", "100+")) +
    labs(title = paste0("Senescence of medulla ",c," over Days"),
         color = "SnC gene sets",
         x = "Days",
         y = "log2FC")
  ggsave(paste0("CompareDays_Medulla_",c,".pdf"),width = 7,height=6)
}

# calculate paired t-test of overall senescence score changes
c="Overall"
sc6d <- data.frame(SnC=rownames(comp_6vs10),log2FC_6D=comp_6vs10[,paste0(c,"_6D")],p_value_6D=comp_6vs10.p[,paste0(c,"_6D")])
sc10d <- data.frame(SnC=rownames(comp_6vs10),log2FC_10D=comp_6vs10[,paste0(c,"_10D")],p_value_10D=comp_6vs10.p[,paste0(c,"_10D")])
scdays <- merge(sc6d,sc10d,by="SnC")
Overall_SC_6vs10D.med <- t.test(scdays$log2FC_6D,scdays$log2FC_10D, paired=T)

scdays.med <- scdays

scdays.cor$SnC == scdays.med$SnC

Overall_SC_cor_vs_med.6D <- t.test(scdays.cor$log2FC_6D,scdays.med$log2FC_6D, paired=T)
Overall_SC_cor_vs_med.10D <- t.test(scdays.cor$log2FC_10D,scdays.med$log2FC_10D, paired=T)



# score4hp_fc_cmp.cor.6d <- as.data.frame(score4hp_fc_cmp.cor)
# score4hp_p_fc.cor.6d <- as.data.frame(score4hp_p_fc.cor)
# score4hp_fc_cmp.cor.10d <- score4hp_fc_cmp
# score4hp_p_fc.cor.10d <- score4hp_fc_cmp.p
#colnames(score4hp_p_fc.cor.6d) <- gsub("_Dox.*","",colnames(score4hp_p_fc.cor.6d))
colnames(score4hp_fc_cmp.cor.6d) <- gsub("_6D","",colnames(score4hp_fc_cmp.cor.6d))
colnames(score4hp_p_fc.cor.6d) <- gsub("_6D","",colnames(score4hp_p_fc.cor.6d))
#colnames(score4hp_fc_cmp.p) <- gsub("_Dox.*","",colnames(score4hp_fc_cmp.p))
colnames(score4hp_fc_cmp.p) <- gsub("_Dox.*","",colnames(score4hp_fc_cmp.p))
colnames(score4hp_fc_cmp.cor.10d) <- gsub("_10D.*","",colnames(score4hp_fc_cmp.cor.10d))

# plot a line graph to show the senescent score changes
library(ggsci)
celltypes <- intersect(colnames(score4hp_fc_cmp.cor.6d),colnames(score4hp_fc_cmp.cor.10d))
line <- list()
for (c in celltypes){
  sc6d <- score4hp_fc_cmp.cor.6d[c]
  sc6d$D <- "6 days"
  sc6d <- tibble::rownames_to_column(sc6d)
  sc6d.p <- score4hp_p_fc.cor.6d[c]
  sc6d.p$D <- "6 days"
  sc6d.p <- tibble::rownames_to_column(sc6d.p)
  sc10d <- score4hp_fc_cmp.cor.10d[c]
  sc10d$D <- "10 days"
  sc10d <- tibble::rownames_to_column(sc10d)
  sc10d.p <- score4hp_fc_cmp.p[c]
  sc10d.p$D <- "10 days"
  sc10d.p <- tibble::rownames_to_column(sc10d.p)
  sc0d <- data.frame(rowname=sc6d$rowname,cell=0,D="0 days")
  colnames(sc0d)[2] <- c
  scdays <- rbind(sc0d,sc6d,sc10d)
  # Convert the 'D' column to a numeric format
  scdays$D <- as.numeric(gsub(" days", "", scdays$D))
  # Plotting
  line[[c]] <- ggplot(scdays, aes(x = D, y = get(c), group = rowname, color = rowname)) +
    geom_line() +
    geom_point(data=scdays, size=2) +
    scale_x_continuous(breaks = c(0, 6, 10)) +
    theme_minimal() + scale_color_igv() +
    labs(title = paste0("Senescence of ",c," over Days"),
         color = "SnC gene sets",
         x = "Days",
         y = c)
  ggsave(plot=line[[c]],paste0(c,"_days_Cortex.pdf"),width = 7,height=6)
}

library(patchwork)
combined_plot <- wrap_plots(line, ncol = floor(length(line)/2))
print(combined_plot)
ggsave(plot=combined_plot,"LinePlot_Cor_Dox_Days.pdf",width=28,height=12)

# Medulla
score4hp_fc_cmp <- read.csv("score4hp_fc_cmp_med_10D.csv")
score4hp_fc_cmp <- column_to_rownames(score4hp_fc_cmp,"X")
colnames(score4hp_fc_cmp) <- paste0(colnames(score4hp_fc_cmp),"_10D")
score4hp_fc_cmp.p <- read.csv("score4hp_p_fc_med_10D.csv")
score4hp_fc_cmp.p <- column_to_rownames(score4hp_fc_cmp.p,"X")
colnames(score4hp_fc_cmp.p) <- paste0(colnames(score4hp_fc_cmp.p),"_10D")

colnames(score4hp_fc_cmp.med) <- paste0(colnames(score4hp_fc_cmp.med),"_6D")
colnames(score4hp_p_fc.med) <- paste0(colnames(score4hp_p_fc.med),"_6D")

comp_6vs10 <- cbind(score4hp_fc_cmp.med,score4hp_fc_cmp)
comp_6vs10.p <- cbind(score4hp_p_fc.med,score4hp_fc_cmp.p)
colnames(comp_6vs10.p) <- gsub("Dox_p_value_","",colnames(comp_6vs10.p))
colnames(comp_6vs10.p) <- gsub("Dox_","",colnames(comp_6vs10.p))



comp_6vs10 <- comp_6vs10[c("Overall_6D","Overall_10D","Endothelial_cells_6D","Endothelial_cells_10D",
                           "Immune_cells_6D","Immune_cells_10D",
                           "Perivascular_cells_6D","Perivascular_cells_10D","Smooth_muscle_cells_6D","Smooth_muscle_cells_10D",
                           "Stromal_cells_6D","Stromal_cells_10D","Theca_Stroma_6D","Theca_Stroma_10D")] # reorder colnames
# reorder the p value matrix by the z-score matrix
comp_6vs10.p <- comp_6vs10.p[rownames(comp_6vs10),colnames(comp_6vs10)]

score4hp_fc_cmp.med.6d <- as.data.frame(score4hp_fc_cmp.med)
score4hp_p_fc.med.6d <- as.data.frame(score4hp_p_fc.med)
score4hp_fc_cmp.med.10d <- score4hp_fc_cmp
#colnames(score4hp_p_fc.cor.6d) <- gsub("_Dox.*","",colnames(score4hp_p_fc.cor.6d))
colnames(score4hp_fc_cmp.med.6d) <- gsub("_6D","",colnames(score4hp_fc_cmp.med.6d))
colnames(score4hp_p_fc.med.6d) <- gsub("_6D|_Dox|_p_value","",colnames(score4hp_p_fc.med.6d))
#colnames(score4hp_fc_cmp.p) <- gsub("_Dox.*","",colnames(score4hp_fc_cmp.p))
colnames(score4hp_fc_cmp.p) <- gsub("_10D.*","",colnames(score4hp_fc_cmp.p))
colnames(score4hp_fc_cmp.p) <- gsub("_Dox|_p_value","",colnames(score4hp_fc_cmp.p))
colnames(score4hp_fc_cmp.med.10d) <- gsub("_10D.*","",colnames(score4hp_fc_cmp.med.10d))

# plot a line graph to show the senescent score changes
celltypes <- intersect(colnames(score4hp_fc_cmp.med.6d),colnames(score4hp_fc_cmp.med.10d ))
line <- list()
for (c in celltypes){
  sc6d <- score4hp_fc_cmp.med.6d[c]
  sc6d$D <- "6 days"
  sc6d <- tibble::rownames_to_column(sc6d)
  sc6d.p <- score4hp_p_fc.med.6d[c]
  sc6d.p$D <- "6 days"
  sc6d.p <- tibble::rownames_to_column(sc6d.p)
  sc10d <- score4hp_fc_cmp.med.10d[c]
  sc10d$D <- "10 days"
  sc10d <- tibble::rownames_to_column(sc10d)
  sc10d.p <- score4hp_fc_cmp.p[c]
  sc10d.p$D <- "10 days"
  sc10d.p <- tibble::rownames_to_column(sc10d.p)
  sc0d <- data.frame(rowname=sc6d$rowname,cell=0,D="0 days")
  colnames(sc0d)[2] <- c
  scdays <- rbind(sc0d,sc6d,sc10d)
  # Convert the 'D' column to a numeric format
  scdays$D <- as.numeric(gsub(" days", "", scdays$D))
  # Plotting
  line[[c]] <- ggplot(scdays, aes(x = D, y = get(c), group = rowname, color = rowname)) +
    geom_line() +
    geom_point(data=scdays, size=2) +
    scale_x_continuous(breaks = c(0, 6, 10)) +
    theme_minimal() + scale_color_igv() +
    labs(title = paste0("Senescence of ",c," over Days"),
         color = "SnC gene sets",
         x = "Days",
         y = c)
  ggsave(plot=line[[c]],paste0(c,"_days_medulla.pdf"),width = 7,height=6)
}

library(patchwork)
combined_plot <- wrap_plots(line, ncol = floor(length(line)/2))
print(combined_plot)
ggsave(plot=combined_plot,"LinePlot_Med_Dox_Days.pdf",width=28,height=12)


# top15 pathways of each cell type
# load IPA pathway comparison z-score
library(readxl)
PW_Z <- read_excel("Ov_6days.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
colnames(PW_Z) <- gsub("_DoxvsCt.*","",colnames(PW_Z)) # truncate the column names
# load IPA pathway comparison adj.p
PW_P <- read_excel("Ov_6days_P.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")
colnames(PW_P) <- gsub("_DoxvsCt.*","",colnames(PW_P)) # truncate the column names

PW_Z_ls <- list()
for (c in names(PW_Z)){
  PW_P_c <- PW_P %>% dplyr::select(c) %>% filter(get(c)>1.3)
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

# PW_c_ud$Pathway <- factor(row.names(PW_c_ud),levels=rev(row.names(PW_c_ud)))
# library(ggplot2)
# library(viridis)
# # Calculate the minimum and maximum values of Adj_P
# min_val <- min(PW_c_ud$Adj_P, na.rm = TRUE)
# max_val <- max(PW_c_ud$Adj_P, na.rm = TRUE)
# # Create an evenly spaced sequence of break points
# breaks <- seq(from = min_val, to = max_val, length.out = 5)
# 
# library(scales)
# ggplot(PW_c_ud, aes(x=Z_Score,y=Pathway,fill=Adj_P))+
#   geom_bar(stat="identity") +
#   ggtitle(paste0("Top 20 in ",c)) +
#   labs(x="Activation z-scores",y="Pathways",fill="-log10(adj.P)") +
#   theme_bw() +
#   scale_fill_viridis_c(breaks = breaks,
#                        labels=label_number(accuracy=0.1))
# ggsave(paste0("top20_",c,".pdf"))

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
                              subtitle = bquote(italic("6 days")))
  vol.plot.ls[[d]] <- vol.plot
}
library(patchwork)
combined_plot <- wrap_plots(vol.plot.ls, ncol = floor(length(vol.plot.ls)/2))
print(combined_plot)
ggsave("Volcano_cor_DoxvsCt.pdf",width=40,height=30)

vol.plot.ls <- list()
cell4vol <- c("Overall",as.character(cell4DEG[!cell4DEG %in% c("Epithelial_cells","Tissue_stem_cells")]))
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
                              subtitle = bquote(italic("6 days")))
  vol.plot.ls[[d]] <- vol.plot
}
library(patchwork)
combined_plot <- wrap_plots(vol.plot.ls, ncol = floor(length(vol.plot.ls)/2))
print(combined_plot)
ggsave("Volcano_med_DoxvsCt.pdf",width=40,height=20)


