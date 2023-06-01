# following the ov_snSeq_01092023.R

setwd("~/ovary_SC")
# import 10X results
# 2nd batch
ov_sample_list <- c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox","RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-ct-1",
                    "RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2","RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1")
Sample.list.cmb <- list()
for (s in ov_sample_list){
  if (s %in% c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox")){
    ov_sample <- Read10X(data.dir = paste0("/bigrock/FurmanLab/Fei/ovary_SC/",s,"/outs/filtered_feature_bc_matrix"))
  } else {
    ov_sample <- Read10X(data.dir = paste0("/data/array2/fwu/ovary_SC_2nd/",s,"/outs/filtered_feature_bc_matrix"))
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

# transfer annotation
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

# DefaultAssay(ov_sub2)<-"SCT"
# Sample.list.cmb[["fist_batch"]] <- ov_sub2
# integrated_variable_features.test <- Sample.list.cmb[["fist_batch"]]@assays[["integrated"]]@var.features
# VariableFeatures(Sample.list.cmb[["fist_batch"]]) <- Sample.list.cmb[["fist_batch"]]@assays[["integrated"]]@var.features

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
                                .default = "batch2")

Sample.cmb$donor <- case_when(Sample.cmb$orig.ident %in%
                                  c("ov-cor-ct","ov-cor-dox","ov-med-ct","ov-med-dox") ~ "donor1",
                              Sample.cmb$orig.ident %in%
                                c("RTL_856_ov-cor-ct-1","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-ct-1","RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2") ~ "donor2",
                              Sample.cmb$orig.ident %in%
                                c("RTL_857_ov-cor-ct-1","RTL_857_ov-cor-dox-1") ~ "donor3")

Sample.cmb$treatment <- case_when(Sample.cmb$orig.ident %in%
                                 c("ov-cor-dox","ov-med-dox","RTL_856_ov-cor-dox-1","RTL_856_ov-cor-dox-2","RTL_856_ov-med-dox-1","RTL_856_ov-med-dox-2","RTL_857_ov-cor-dox-1") ~ "dox",
                               .default = "ctl")

# find markers for every cluster compared to all remaining cells
markers.all <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(markers.all,"markers_all_clusters.csv")

Idents(Sample.cmb) <- "predicted.id"
markers.all.cellType <- FindAllMarkers(Sample.cmb, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="donor")
write.csv(markers.all.cellType,"markers_all_cellTypes.csv")

# Count the frequency of each term in the "cluster" column
cluster_frequency <- count(Sample.cmb@meta.data,predicted.id, sort=T)
# reorder the cell type column by frequency
Sample.cmb$predicted.id <- factor(Sample.cmb$predicted.id, levels = cluster_frequency$predicted.id)
# 
# markers.all.cellType$cluster <- factor(markers.all.cellType$cluster, levels = cluster_frequency$predicted.id)
# markers.all.cellType <- markers.all.cellType %>% arrange(cluster)

# plot top markers in heatmap
markers.all.cellType %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(Sample.cmb, features=top10$gene, group.by="predicted.id",size =4)
ggsave("ov_cell_markers_top10.pdf",width=8,height=16)

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

# plot UMAP based on cell types for tissues
Idents(Sample.cmb) <- "predicted.id"
DimPlot(Sample.cmb, reduction = "umap",split.by = "tissue",label =T,repel =T)
ggsave("ov_cellTypes_combined.pdf",width=12,height=6)


# ------ find DEG for all cell types after Dox treatment ------ #
# with default filter min.pct=0.1 and logfc=0.25
Sample.cmb$celltype_treatment <- paste(Sample.cmb$predicted.id, Sample.cmb$tissue_treat,
                                                   sep = "_")
Idents(Sample.cmb)<-"celltype_treatment"
plan("multisession", workers = 4)
cell4DEG<-unique(Sample.cmb$predicted.id)
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

# ------ Cell type ratio change stack barplot ------ #
sub.prop.all<-data.frame()
for (l in unique(Sample.cmb$tissue_treat)){
  sub.treat<-Sample.cmb@meta.data[Sample.cmb$tissue_treat==l,]
  sub.prop<-data.frame(table(sub.treat$predicted.id)/sum(table(sub.treat$predicted.id)))
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
  group_by(tissue_treat,predicted.id) %>%
  summarise(cell_num = length(predicted.id)) %>%
  group_by(tissue_treat) %>%
  summarise(predicted.id = predicted.id,
            cell_num = cell_num,
            sum = sum(cell_num),
            Freq_mean = cell_num/sum(cell_num))

Cell_Results_donor <- Sample.cmb@meta.data %>%
  group_by(tissue_treat,predicted.id,donor) %>%
  summarise(cell_num_donor = length(donor)) %>%
  group_by(tissue_treat,donor) %>%
  summarise(predicted.id,
            cell_num_donor,
            sum_donor = sum(cell_num_donor),
            Freq_donor = cell_num_donor/sum_donor) %>%
  group_by(predicted.id, tissue_treat) %>%
  summarise(donor,cell_num_donor,sum_donor,Freq_donor,
            se = sd(Freq_donor) / sqrt(n()))
  

Cell_Results4plot_summary <- left_join(Cell_Results4plot_summary,Cell_Results_donor,
                                       by=c("tissue_treat","predicted.id"))

# Create the bar plot with T-shaped error bars
bar_plot <- ggplot(Cell_Results4plot_summary, aes(x = predicted.id, y = Freq_mean, fill = tissue_treat)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Freq_mean - se, ymax = Freq_mean + se), position = position_dodge(0.9)) +
  ggtitle("Ovary Cell Type Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Add the individual "Freq" values as dots
dot_plot <- ggplot(Cell_Results4plot_summary, aes(x = predicted.id, y = Freq_donor, fill = tissue_treat)) +
  geom_jitter(size = 1, width = 0.2 / length(unique(Cell_Results4plot_summary$tissue_treat))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Combine the bar plot with error bars and the dot plot
combined_plot <- bar_plot + 
  geom_point(data = Cell_Results4plot_summary, aes(x = predicted.id, y = Freq_donor, fill = tissue_treat),
             position = position_dodge(0.9), size = 0.5, alpha=0.6)
print(combined_plot)
ggsave("Cell_type_frequency_bar.pdf",width=8,height=6)

# cell type ratio change: log2FC comparision unstim; remove average celltype ratio <1%
# Cortex
sub.prop.cortex_ctl<-subset(Cell_Results4plot_summary, tissue_treat=="cortex_ctl" & Freq_mean>0.01)
sub.prop.cortex_dox<-subset(Cell_Results4plot_summary, tissue_treat=="cortex_dox" & Freq_mean>0.01)
sub.prop.cortex.comm<-intersect(sub.prop.cortex_ctl$predicted.id,sub.prop.cortex_dox$predicted.id)
sub.prop.cortex_ctl <- subset(sub.prop.cortex_ctl, predicted.id %in% sub.prop.cortex.comm)
sub.prop.cortex_dox <- subset(sub.prop.cortex_dox, predicted.id %in% sub.prop.cortex.comm)
sub.prop.cor_DoxvsCtl <- data.frame(Cell_type=sub.prop.cortex_ctl$predicted.id,
                              group="Dox vs Ctl Cortex",
                              Fold_change=sub.prop.cortex_dox$Freq_mean/sub.prop.cortex_ctl$Freq_mean)
sub.prop.cor_DoxvsCtl$log2FC<-log2(sub.prop.cor_DoxvsCtl$Fold_change)

sub.prop.cor_DoxvsCtl<-sub.prop.cor_DoxvsCtl%>% distinct()

# add t-test p-value
sub.prop.cortex_ctl.full <- subset(sub.prop.cortex_ctl, tissue_treat=="cortex_ctl")
sub.prop.cortex_dox.full <- subset(sub.prop.cortex_dox, tissue_treat=="cortex_dox")

t.test.num <- numeric()
for (c in sub.prop.cortex.comm){
  sub.prop.cortex_ctl.temp <- sub.prop.cortex_ctl.full %>% filter(predicted.id == c) %>% select(Freq_donor)
  sub.prop.cortex_dox.temp <- sub.prop.cortex_dox.full %>% filter(predicted.id == c) %>% select(Freq_donor)
  t.test.temp <- t.test(sub.prop.cortex_ctl.temp$Freq_donor, sub.prop.cortex_dox.temp$Freq_donor)
  t.test.num <- c(t.test.num,t.test.temp$p.value)
}
sub.prop.cor_DoxvsCtl$p_value <- t.test.num
# bar plot
ggplot(sub.prop.cor_DoxvsCtl, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Dox vs Ctl Ovary Cortex") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  geom_text(aes(label = ifelse(Cell_type == "CD14.Mono", "*", "")),
            position = position_dodge(0.9), vjust = 1.5,size=5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_Cortex_Dox_vs_Ctl.pdf",width=6,height=4)

# Medulla
sub.prop.medulla_ctl<-subset(Cell_Results4plot_summary, tissue_treat=="medulla_ctl" & Freq_mean>0.01)
sub.prop.medulla_dox<-subset(Cell_Results4plot_summary, tissue_treat=="medulla_dox" & Freq_mean>0.01)
sub.prop.medulla.comm<-intersect(sub.prop.medulla_ctl$predicted.id,sub.prop.medulla_dox$predicted.id)
sub.prop.medulla_ctl <- subset(sub.prop.medulla_ctl, predicted.id %in% sub.prop.medulla.comm)
sub.prop.medulla_dox <- subset(sub.prop.medulla_dox, predicted.id %in% sub.prop.medulla.comm)
sub.prop.med_DoxvsCtl <- data.frame(Cell_type=sub.prop.medulla_ctl$predicted.id,
                                    group="Dox vs Ctl medulla",
                                    Fold_change=sub.prop.medulla_dox$Freq_mean/sub.prop.medulla_ctl$Freq_mean)
sub.prop.med_DoxvsCtl$log2FC<-log2(sub.prop.med_DoxvsCtl$Fold_change)

sub.prop.med_DoxvsCtl<-sub.prop.med_DoxvsCtl%>% distinct()

# add t-test p-value
sub.prop.medulla_ctl.full <- subset(sub.prop.medulla_ctl, tissue_treat=="medulla_ctl")
sub.prop.medulla_dox.full <- subset(sub.prop.medulla_dox, tissue_treat=="medulla_dox")

t.test.num <- numeric()
for (c in sub.prop.medulla.comm){
  sub.prop.medulla_ctl.temp <- sub.prop.medulla_ctl.full %>% filter(predicted.id == c) %>% select(Freq_donor)
  sub.prop.medulla_dox.temp <- sub.prop.medulla_dox.full %>% filter(predicted.id == c) %>% select(Freq_donor)
  t.test.temp <- t.test(sub.prop.medulla_ctl.temp$Freq_donor, sub.prop.medulla_dox.temp$Freq_donor)
  t.test.num <- c(t.test.num,t.test.temp$p.value)
}
sub.prop.med_DoxvsCtl$p_value <- t.test.num
# bar plot
ggplot(sub.prop.med_DoxvsCtl, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Dox vs Ctl Ovary Medulla") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  geom_text(aes(label = ifelse(Cell_type == "CD14.Mono", "*", "")),
            position = position_dodge(0.9), vjust = 1.5,size=5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_medulla_Dox_vs_Ctl.pdf",width=6,height=4)

# normalize and scale RNA assays
DefaultAssay(Sample.cmb) <- "RNA"
Sample.cmb <- NormalizeData(Sample.cmb)
all.genes <- rownames(Sample.cmb)
Sample.cmb <- ScaleData(Sample.cmb, features = all.genes)

# load the senescence gene sets
Senescence_genesets <- read.csv("Senescence_genesets.csv")

sc_module <- Sample.cmb
# scoring the senescence
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
      select(paste0(names(Senescence_genesets[sc]),"1"))
    dox <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_cortex_dox")) %>% 
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
    filter(tissue_treat=="cortex_ctl") %>% 
    select(paste0(names(Senescence_genesets[sc]),"1"))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="cortex_dox") %>% 
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
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.4), c("blue", "white", "red"))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_cor_cmb.pdf", width=10, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Cortex Dox and Control",
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
                  grid.text(round(score4hp[i, j],3), x, y,gp=gpar(fontsize=6.5))
                }})
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
      filter(celltype_treatment==paste0(cell,"_medulla_ctl")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
    dox <- sc_module@meta.data %>% 
      filter(celltype_treatment==paste0(cell,"_medulla_dox")) %>% 
      select(paste0(names(Senescence_genesets[sc]),"1"))
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
    select(paste0(names(Senescence_genesets[sc]),"1"))
  dox <- sc_module@meta.data %>% 
    filter(tissue_treat=="medulla_dox") %>% 
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
score4hp <- column_to_rownames(sc_score,"SC_gene_set")
score4hp_p <- column_to_rownames(sc_score_p,"SC_gene_set")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.4), c("blue", "white", "red"))

score4hp <- as.matrix(score4hp)
pdf(file="SnC-score_ov_med_cmb.pdf", width=9, height = 5)
ht <- Heatmap(score4hp,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              column_title = "Medulla Dox and Control",
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


# compare the 10x single nuclei and Nanostring GeoMx
#GeoMx_Stroma_Cor_P16PosvsP16Neg <-read.csv("~/ovary_SC/GeoMx/DEStroma_Cortex_P16PositivevsP16NegDEALL.csv")

GeoMx_Stroma_Cor_P16PosvsP16Neg <-read.delim("~/ovary_SC/GeoMx/DEStroma_Cortex_P16PositivevsP16NegDEALL.csv",sep = ";")


# p16pos_stromal_cor.sig <- p16pos_stromal_cor %>% filter(p_val_adj<0.05)
# row.names(p16pos_stromal_cor.sig)
# GeoMx_Stroma_Cor_P16PosvsP16Neg$SYMBOL
# 
# intersect(row.names(p16pos_stromal_cor.sig),GeoMx_Stroma_Cor_P16PosvsP16Neg$SYMBOL)
# 
GeoMx_Stroma_Cor_P16PosvsP16Neg<-GeoMx_Stroma_Cor_P16PosvsP16Neg %>% select(c("SYMBOL","logFC","P.Value"))
GeoMx_Stroma_Cor_P16PosvsP16Neg <- column_to_rownames(GeoMx_Stroma_Cor_P16PosvsP16Neg,"SYMBOL")
colnames(GeoMx_Stroma_Cor_P16PosvsP16Neg)<-c("avg_log2FC","p_val")

# overlaps <- inner_join(GeoMx_Stroma_Cor_P16PosvsP16Neg,rownames_to_column(p16pos_stromal_cor), by=join_by(SYMBOL == rowname))




library(RRHO2)
# Create "gene" lists:
p16pos_stromal_cor.rank <- p16pos_stromal_cor %>% 
  select(avg_log2FC,p_val) %>%
  mutate(p_val_adj=ifelse(p_val==0,1e-305,p_val)) %>%
  rownames_to_column("gene_symbol")
p16pos_stromal_cor.rank$weight <- p16pos_stromal_cor.rank$avg_log2FC*-log10((p16pos_stromal_cor.rank$p_val_adj))
p16pos_stromal_cor.rank$rank <- order(p16pos_stromal_cor.rank$weight,decreasing=T)
p16pos_stromal_cor.rank <- p16pos_stromal_cor.rank %>% select(gene_symbol,weight)

GeoMx_Stroma_Cor_P16PosvsP16Neg.rank <- GeoMx_Stroma_Cor_P16PosvsP16Neg %>% 
  select(avg_log2FC,p_val) %>%
  mutate(p_val_adj=ifelse(p_val==0,1e-305,p_val)) %>%
  rownames_to_column("gene_symbol")
GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$weight <- GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$avg_log2FC*-log10((GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$p_val_adj))
GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$rank <- order(GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$weight,decreasing=T)
GeoMx_Stroma_Cor_P16PosvsP16Neg.rank <- GeoMx_Stroma_Cor_P16PosvsP16Neg.rank %>% select(gene_symbol,weight)

overlap_genes<-intersect(p16pos_stromal_cor.rank$gene_symbol,GeoMx_Stroma_Cor_P16PosvsP16Neg.rank$gene_symbol)

p16pos_stromal_cor.rank.o <- filter(p16pos_stromal_cor.rank,gene_symbol%in%overlap_genes)
GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o <- filter(GeoMx_Stroma_Cor_P16PosvsP16Neg.rank,gene_symbol%in%overlap_genes)


RRHO_obj <-  RRHO2_initialize(p16pos_stromal_cor.rank.o,
                              GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,
                              labels = c("sn", "GeoMx"), log10.ind=TRUE,
                              boundary = 0.05)

pdf(file="RRHO2_sn_vs_GeoMx.pdf", width=7, height = 6)
RRHO2_heatmap(RRHO_obj)
dev.off()

pdf(file="RRHO2_sn_vs_GeoMx_DD.pdf", width=6, height = 5)
RRHO2_vennDiagram(RRHO_obj,type="dd")
dev.off()

pdf(file="RRHO2_sn_vs_GeoMx_UU.pdf", width=6, height = 5)
RRHO2_vennDiagram(RRHO_obj,type="uu")
dev.off()

# extract the uu gene expressions
uu<-RRHO_obj[["genelist_uu"]][["gene_list_overlap_uu"]]
# extract the dd gene expressions
dd<-RRHO_obj[["genelist_dd"]][["gene_list_overlap_dd"]]

# uu in sn
uu.sn<-select(p16pos_stromal_cor,c(avg_log2FC,p_val))[uu,]
colnames(uu.sn)<-paste0("sn_",colnames(uu.sn))
# uu in GeoMx
uu.GeoMx <- GeoMx_Stroma_Cor_P16PosvsP16Neg[uu,]
colnames(uu.GeoMx)<-paste0("GeoMx_",colnames(uu.GeoMx))

uu.sn.GeoMx <- cbind(uu.sn,uu.GeoMx)
write.csv(uu.sn.GeoMx,"uu.sn.GeoMx.csv")

# dd in sn
dd.sn<-select(p16pos_stromal_cor,c(avg_log2FC,p_val))[dd,]
colnames(dd.sn)<-paste0("sn_",colnames(dd.sn))
# dd in GeoMx
dd.GeoMx <- GeoMx_Stroma_Cor_P16PosvsP16Neg[dd,]
colnames(dd.GeoMx)<-paste0("GeoMx_",colnames(dd.GeoMx))

dd.sn.GeoMx <- cbind(dd.sn,dd.GeoMx)
write.csv(dd.sn.GeoMx,"dd.sn.GeoMx.csv")

# create dot plot
uu.sn4dot<-select(p16pos_stromal_cor,c(avg_log2FC,p_val))[uu,]
uu.sn4dot$tech <- "sn"
uu.sn4dot$gene <- row.names(uu.sn4dot)
uu.GeoMx4dot <- GeoMx_Stroma_Cor_P16PosvsP16Neg[uu,]
uu.GeoMx4dot$tech <- "GeoMx"
uu.GeoMx4dot$gene <- row.names(uu.GeoMx)
uu4dot <- rbind(uu.sn4dot,uu.GeoMx4dot)

library(dplyr)
library(ggplot2)
library(stringr)
uu4dot %>%
  ggplot(aes(x=tech, y = gene, color = avg_log2FC, size = -log10(p_val))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Gene') +
  xlab('Technology') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("uu_genes_sn_GeoMx.pdf", width = 4, height = 5,limitsize = FALSE)

# # Get the enrichment score
# max_neglog10_pval <- max(-log10(bounded_hypermat))
# 
# cor(p16pos_stromal_cor.rank.o$weight,
#                  GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o$weight)

# 
# inner_join(filter(p16pos_stromal_cor.rank.o,weight>0.4),
#            filter(GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,weight>0.4),by="gene_symbol")
# inner_join(filter(p16pos_stromal_cor.rank.o,weight<(-0.5)),
#            filter(GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,weight<(-0.5)),by="gene_symbol")
# 

library(RRHO)
pth1 <- "/opt/home/buckcenter.org/fwu/ovary_SC/RRHO"
RRHO(p16pos_stromal_cor.rank.o, GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,
     labels=c("sn","GeoMx"),
     alternative = 'enrichment',BY=T,log10.ind = T,plots=T, outputdir = pth1)

RRHO1<-RRHO(p16pos_stromal_cor.rank.o, GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,
     labels=c("sn","GeoMx"),
     alternative = 'enrichment',BY=T,log10.ind = T,plots=F, outputdir = NULL)

correlation.test <- merge(p16pos_stromal_cor.rank.o,GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o,by="gene_symbol")
cor.test(correlation.test$weight.x, correlation.test$weight.y, method = "spearman")

# # Plot density
# density <- density(GeoMx_Stroma_Cor_P16PosvsP16Neg.rank.o$weight)
# plot(density, main = "Density Plot",xlim=c(-0.5,0.5))

# the correlation between reads counts
# counts of sn RNAseq
counts.sn.cor.stroma.ctl <- data.frame(gene_name=row.names(cortex_ctl_stromal@assays[["RNA"]]@counts),sn_means=rowSums(cortex_ctl_stromal@assays[["RNA"]]@counts))

counts.GeoMx.cor.stroma.ctl <- read.csv("/opt/home/buckcenter.org/fwu/ovary_SC/GeoMx/countsStroma.csv")
counts.GeoMx.cor.stroma.ctl <- data.frame(gene_name=counts.GeoMx.cor.stroma.ctl["X"],GeoMx_means=rowSums(counts.GeoMx.cor.stroma.ctl[-1]))
colnames(counts.GeoMx.cor.stroma.ctl)[1]<-"gene_name"

sn.GeoMx <- inner_join(counts.sn.cor.stroma.ctl,counts.GeoMx.cor.stroma.ctl,by="gene_name")

sn.GeoMx$log_sn_means <- log10(sn.GeoMx$sn_means)
sn.GeoMx$log_GeoMx_means <- log10(sn.GeoMx$GeoMx_means)

library(ggpubr)
library(RColorBrewer)
ggscatter(sn.GeoMx, x = "log_sn_means", y = "log_GeoMx_means", 
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "log10 of sn counts", ylab = "log10 of GeoMx counts", title = "sn and GeoMx Cortex Stroma",
          size = 0.5,color = brewer.pal(9, 'Paired')[1]) +
  geom_smooth(method="lm",color = brewer.pal(9, 'Paired')[2]) +
  scale_y_continuous(limits = c(0, 5)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
ggsave("snXGeoMx_cor_stroma_spearman_log.pdf",width = 5, height = 5)


