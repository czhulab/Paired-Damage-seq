library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(sctransform)
library(ggthemes)
library(viridis)
library(Matrix)

set.seed(011024)
theme_set(theme_tufte())
plot_dir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/RNA/Seurat/"
data_outdir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/RNA/"

###########################################################
##################### load the data #######################
###########################################################
data <- Read10X(data.dir="/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/05.20230901_scAnalysis/01.merged_bam_files/02.hg38_RNA/03.gene_matrix/RNA_filtered_matrix/")
seu <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
seu

###########################################################
################# remove unwanted genes ###################
###########################################################
features<-grep("^AC[0-9]", rownames(seu@assays$RNA@counts))
features<-c(features, grep("^MIR[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^FP[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^LINC[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^CR[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^RP[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^AL[0-9]", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^sno", rownames(seu@assays$RNA@counts)))
features<-c(features, grep("^PCDH", rownames(seu@assays$RNA@counts)))
length(features)

All_genes <- rownames(seu)
wanted_genes <- All_genes[-features]
seu <- subset(seu, features = wanted_genes)
seu

###########################################################
#################### QC by each library ###################
###########################################################


pdf(paste0(plot_dir, "Vln_UMI_beforeQC.pdf"))
VlnPlot(seu, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0)
dev.off()

seu_list <- list()
lib <- as.character(unique(seu$orig.ident))

for (l in lib){
    print(l)
    obj <- subset(seu, orig.ident == l)
    low_UMI_cutoff <- quantile(obj$nCount_RNA, 0.01)
    high_UMI_cutoff <- quantile(obj$nCount_RNA, 0.99)
    obj <- subset(obj, nCount_RNA < high_UMI_cutoff & nCount_RNA > low_UMI_cutoff)
    seu_list[[length(seu_list) + 1]] <- obj
}

new_seu <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
new_seu

pdf(paste0(plot_dir, "Vln_UMI_PostQC.pdf"))
VlnPlot(new_seu, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0)
dev.off()

pdf(paste0(plot_dir, "Vln_nGene_PostQC.pdf"))
VlnPlot(new_seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0)
dev.off()

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
pdf(paste0(plot_dir, "Vln_mtCount_PostQC.pdf"))
VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0)
dev.off()


seu <- new_seu

###########################################################
######## add treatment according to barcode ###############
###########################################################

### add treatment according to barcode
bc<-matrix(ncol=3, byrow = T, unlist(strsplit(rownames(seu@meta.data), split=":")))#;head(bc)
Treatment<-rep("Ctrl", dim(bc)[1])
Treatment[bc[,3]=="03" | bc[,3]=="04"]<-"0hr"
Treatment[bc[,3]=="05" | bc[,3]=="06"]<-"2hr"
Treatment[bc[,3]=="07" | bc[,3]=="08"]<-"6hr"
Treatment[bc[,3]=="09" | bc[,3]=="10"]<-"24hr"
Treatment[bc[,3]=="11" | bc[,3]=="12"]<-"48hr"

seu$Treatment<-factor(Treatment, levels=c("Ctrl", "0hr", "2hr", "6hr", "24hr", "48hr"))

###########################################################
##################### pre-processing ######################
###########################################################
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- SCTransform(seu, vst.flavor = "v2", verbose = T)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3, algorithm = 4, method = "igraph")
seu <- FindClusters(seu, resolution = 0.4, algorithm = 4, method = "igraph")
seu <- RunUMAP(seu, dims = 1:30, metric = "cosine", min.dist = 0.1)

#assignment cell cycles
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)




###########################################################
######################## plotting #########################
###########################################################
pdf(paste0(plot_dir, "UMAP_CellPhase.pdf"), width = 7, height = 6)
DimPlot(seu, group.by = "Phase") + theme(
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(), 
  plot.title = element_blank()
)
dev.off()

pdf(paste0(plot_dir, "UMAP_leiden04.pdf"), width = 7, height = 6)
DimPlot(seu, group.by = "SCT_snn_res.0.4", label = TRUE)
dev.off()

pdf(paste0(plot_dir, "UMAP_leiden03.pdf"), width = 7, height = 6)
DimPlot(seu, group.by = "SCT_snn_res.0.3", label = TRUE)
dev.off()

pdf(paste0(plot_dir, "FeatureUMAP_counts.pdf"), width = 16, height = 6)
FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA"), label = TRUE)
dev.off()

colors <- viridis(n = 6, option = "plasma")
pdf(paste0(plot_dir, "UMAP_Treatment.pdf"), width = 7, height = 6)
DimPlot(seu, reduction = "umap", group.by = "Treatment", cols = colors) + theme(
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(), 
  plot.title = element_blank()
)
dev.off()

pdf(paste0(plot_dir, "Vln_leiden03_counts.pdf"), width = 16, height = 6)
VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"), 
        group.by = "SCT_snn_res.0.3", pt.size = 0)
dev.off()
###########################################################
########################### DEG ###########################
###########################################################
Idents(seu) <- seu$SCT_snn_res.0.4
leiden04_markers <- FindAllMarkers(seu, only.pos = F, min.pct = 0.20, logfc.threshold = 0.5)
write.table(leiden04_markers, paste0(plot_dir, "leiden04_marker.csv"), row.names = F, quote = F, sep = ",")


Idents(seu) <- seu$SCT_snn_res.0.3
leiden03_markers <- FindAllMarkers(seu, only.pos = F, min.pct = 0.20, logfc.threshold = 0.5)
write.table(leiden03_markers, paste0(plot_dir, "leiden03_marker.csv"), row.names = F, quote = F, sep = ",")



top_markers <- leiden04_markers %>% 
               group_by(cluster) %>%
               arrange(desc(avg_log2FC))


###########################################################
####################### DEG Plotting#######################
GOI <- c("NEK7", "FGF2", "EPHA7", "PIK3R3", "AKR1C2", "ZC3HAV1")

pdf(paste0(plot_dir, "Dot_markers.pdf"), width = 7, height = 6)
DotPlot(seu, features = GOI, group.by = "RNA_anno") + 
        theme(axis.text.x = element_text(angle = 30, hjust=1))
dev.off()



###########################################################
###################### Annotation #########################
###########################################################
RNA_anno<-rep("Unknown", dim(seu@meta.data)[1])
sc<-seu$seurat_clusters
RNA_anno[sc==1] <- "Early_NEK7"
RNA_anno[sc==2] <- "Early_FGF2"
RNA_anno[sc==3] <- "Late_PIK3R3"
RNA_anno[sc==4] <- "Early_EPHA7"
RNA_anno[sc==5] <- "Mid"
RNA_anno[sc==6] <- "Late_EPHA7"
RNA_anno[sc==7] <- "Late_AKR1C2"
RNA_anno[sc==8] <- "ZC3HAV1"

table(RNA_anno)

seu$RNA_anno<-factor(RNA_anno, levels=c("Early_NEK7","Early_FGF2","Early_EPHA7","Mid",
                                        "Late_PIK3R3","Late_EPHA7","Late_AKR1C2","ZC3HAV1"))


colors <- brewer.pal(8, "Paired")
colors[3] <- "#CAB2D6"
colors[4] <- "#FDBF6F"
colors[7] <- "#FF7F00"
colors[8] <- "#33A02C"
colors 
# "#A6CEE3" "#1F78B4" "#CAB2D6" "#FDBF6F" "#FB9A99" "#E31A1C" "#FF7F00" "#33A02C"

pdf(paste0(plot_dir, "UMAP_leiden03_annotated.pdf"), width = 7, height = 6)
DimPlot(seu, group.by = "RNA_anno", label = F, cols = colors) + theme(
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(), 
  plot.title = element_blank()
)
dev.off()

###########################################################
######################## export ###########################
###########################################################
meta <- seu@meta.data
meta$barcodes <- rownames(meta)
write.table(meta, paste0(data_outdir, "metadata_011024_seu.csv"), 
            sep = ",", row.names = F, quote = F)

umap_cor <- as.data.frame(Embeddings(seu, "umap"))
umap_cor$barcodes <- rownames(umap_cor)
write.table(umap_cor, paste0(data_outdir, "UMAP_cor_011024_seu.csv"), 
            sep = ",", row.names = F, quote = F)

genes_used <- rownames(seu)
write.table(genes_used, paste0(data_outdir, "gene_used_011024_seu.csv"), 
            sep = ",", row.names = F, quote = F)

saveRDS(seu, paste0(data_outdir, "011024_seu.rds"))
seu <- readRDS(paste0(data_outdir, "011024_seu.rds"))



ExprMatrix <- seu@assays$SCT@counts
ExprMatrix <- t(ExprMatrix)
writeMM(ExprMatrix, paste0(data_outdir, "SCTMatrix_030424_seu.mtx"))
SCT_genes <- colnames(ExprMatrix)
SCT_cells <- rownames(ExprMatrix)
write.table(SCT_genes, paste0(data_outdir, "SCTGenes_011024_seu.csv"), 
            sep = ",", row.names = F, quote = F)
write.table(SCT_cells, paste0(data_outdir, "SCTCells_011024_seu.csv"), 
            sep = ",", row.names = F, quote = F)

ExprMatrix <- seu@assays$SCT@data
ExprMatrix <- t(ExprMatrix)
writeMM(ExprMatrix, paste0(data_outdir, "SCTMatrixNorm_030424_seu.mtx"))

ExprMatrix <- seu@assays$SCT@scale.data
ExprMatrix <- t(ExprMatrix)
writeMM(ExprMatrix, paste0(data_outdir, "SCTMatrixScale_030424_seu.mtx"))


