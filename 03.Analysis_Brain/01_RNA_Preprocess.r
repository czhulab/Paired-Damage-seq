library(Seurat)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(reticulate)
use_condaenv(condaenv = "~/.conda/envs/SR/bin/python", conda = "auto", required = NULL)
set.seed(042424)



#deal with the sample list 
sample_list <- read.table("/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/11.Brain_dataset/sample_list.xls", 
                          sep = "\t", header = T, fill = T)

lib_path <- c()
lib_dirs <- list.dirs("/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/11.Brain_dataset/Jupyter_analysis/filtered_matrices/")
for (RNA_lib_name in sample_list$RNA){
    dir <- grep(RNA_lib_name, lib_dirs, value = T)
    lib_path <- c(lib_path, dir)
}

sample_list$RNA_lib_path <- lib_path

# read each individual library
seu_list <- list()
UMI_cutoff <- c()
for (i in 1:nrow(sample_list)){
    lib <- sample_list[i,2]
    print(lib)
    data <- Read10X(sample_list[i,4])
    seu <- CreateSeuratObject(counts = data)
    seu$lib <- lib
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

    top_quantile <- quantile(seu$nCount_RNA, 0.985)
    print(paste("UMI_cutoff:", top_quantile))
    UMI_cutoff <- c(UMI_cutoff, top_quantile)
    seu <- subset(seu, nCount_RNA < top_quantile)


    pdf(paste0("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/RNA/Seurat/individual/", 
                lib, "_3QC.pdf"))
    p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
    print(p)
    dev.off()

    print(paste("finish:", lib))
    seu_list[length(seu_list) + 1] <- seu
}

sample_list$UMI_985 <- UMI_cutoff



#f2 <- rownames(seu)[grepl("^Gm\\d{3}.*", rownames(seu))]
#writeLines(f2, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/RNA/Seurat/merged/GmGenes.txt")

#merge objects and perform joint analysis
plot_dir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/RNA/Seurat/merged/"

seu <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
seu <- subset(seu, nFeature_RNA > 200)

features <- rownames(seu)[(!grepl("^Pcdh", rownames(seu))) & (!grepl("^Gm\\d{3}.*", rownames(seu))) & (!grepl("^Ugt1a", rownames(seu)))]

rownames(seu)[(grepl("^Gm\\d{3}.*", rownames(seu)))]


seu <- subset(seu, features = features)

pdf(paste0(plot_dir, "3QC.pdf"), width = 16, height = 6)
for (f in c("nFeature_RNA", "nCount_RNA", "percent.mt")){
    print(
        VlnPlot(seu, group.by = "lib", features = f, pt.size = 0)
    )
}
dev.off()

rm(seu_list)

# add the age groups 
seu$lib_groups <- gsub(".*:([0-9]+)_.*", "\\1", colnames(seu))
seu$lib_groups <- as.numeric(seu$lib_groups)
seu <- subset(seu, lib_groups <= 6)
seu

# run the preprocessing
seu <- SCTransform(seu, verbose=TRUE, variable.features.n = 5000, conserve.memory=TRUE)
seu <- RunPCA(seu, verbose=TRUE)
seu <- RunUMAP(seu, dims = 1:30, min.dist = 0.1)

pdf(paste0(plot_dir, "UMAP_libs.pdf"))
DimPlot(seu, reduction = "umap", group.by = "lib_groups", raster=FALSE)
dev.off()

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3, algorithm = 4, method = "igraph")


pdf(paste0(plot_dir, "UMAP_Leiden03.pdf"))
DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.3", raster=FALSE, label = TRUE)
dev.off()

DefaultAssay(seu) <- "SCT"
pdf(paste0(plot_dir, "UMAP_InN.pdf"), width = 8, height = 5)
FeaturePlot(seu, features = c("Sst", "Pvalb", "Adarb2"))
dev.off()





genes <- c("Snap25","Slc17a7","Cux2","Rorb","Parm1","Fezf2","Foxp2","Gad1","Gad2","Vip", "Adarb2", "Sst","Pvalb","Pdgfra","Mbp","Mobp","Slc1a2","Flt1","Csf1r")
pdf(paste0(plot_dir, "DotPlot_markers.pdf"))
DotPlot(seu, features=genes) + RotatedAxis()
dev.off()

genes <- c("Snap25","Slc17a7","Cux2","Rorb","Parm1","Fezf2","Foxp2","Gad1","Gad2","Vip", "Adarb2", "Sst","Pvalb","Pdgfra","Mbp","Mobp","Slc1a2","Flt1","Csf1r")
pdf(paste0(plot_dir, "DotPlot_markers_L1.pdf"))
DotPlot(seu, features=genes, group.by = "L1") + RotatedAxis()
dev.off()

genes <- c("Snap25","Slc17a7","Cux2","Rorb","Parm1","Fezf2","Foxp2","Gad1","Gad2","Vip", "Adarb2", "Sst","Pvalb","Pdgfra","Mbp","Mobp","Slc1a2","Flt1","Csf1r")
pdf(paste0(plot_dir, "DotPlot_markers_Anno.pdf"))
DotPlot(seu, features=genes, group.by = "Annotation") + RotatedAxis()
dev.off()

DefaultAssay(seu) <- "RNA"
seu

pdf(paste0(plot_dir, "Vln_leiden03.pdf"))
VlnPlot(seu, group.by = "SCT_snn_res.0.3", features = "nCount_RNA", pt.size=0, log = T)
dev.off()

saveRDS(seu, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524.rds")
seu <- readRDS("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524.rds")


seu <- subset(seu, SCT_snn_res.0.3 != "15")
seu


# run the preprocessing
#seu <- SCTransform(seu, verbose=TRUE, variable.features.n = 5000, conserve.memory=TRUE)
#seu <- RunPCA(seu, verbose=TRUE)
#seu <- RunUMAP(seu, dims = 1:30, min.dist = 0.1)
#seu <- FindNeighbors(seu, dims = 1:30)
#seu <- FindClusters(seu, resolution = 0.3, algorithm = 4, method = "igraph")


#pdf(paste0(plot_dir, "UMAP_Leiden03_2.pdf"))
#DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.3", raster=FALSE, label = TRUE)
#dev.off()


#pdf(paste0(plot_dir, "Vln_leiden03_2.pdf"))
#VlnPlot(seu, group.by = "SCT_snn_res.0.3", features = "nCount_RNA", pt.size=0, log = T)
#dev.off()


#pdf(paste0(plot_dir, "DotPlot_markers_2.pdf"))
#DotPlot(seu, features=genes) + RotatedAxis()
#dev.off()

saveRDS(seu, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524_2.rds")
seu <- readRDS("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524_2.rds")
DefaultAssay(seu) <- "SCT"

markers <- FindAllMarkers(seu)
write.table(markers, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/RNA/Seurat/merged/markers_leiden03_SCT.csv", sep = ",", row.names = F)


pdf(paste0(plot_dir, "DotPlot_markers_3.pdf"))
DotPlot(seu, features=genes) + RotatedAxis()
dev.off()


##Annotation
table(seu$SCT_snn_res.0.3)
Annotation <- rep("Unknown",dim(seu@meta.data)[1])
sc <- seu$SCT_snn_res.0.3

Annotation[sc==1] <-"ExN_L23"
Annotation[sc==2] <-"AST"
Annotation[sc==3] <-"ODC"
Annotation[sc==4] <-"InN_CGE1"
Annotation[sc==5] <-"InN_CGE2"
Annotation[sc==6] <-"ExN_L6"
Annotation[sc==7] <-"OPC"
Annotation[sc==8] <-"MiG_Inpp5d"
Annotation[sc==9] <-"VLMC"
Annotation[sc==10] <-"Endo"
Annotation[sc==11] <-"ExN_L5"
Annotation[sc==12] <-"ExN_L23_2"
Annotation[sc==13] <-"MiG_Stab1"
Annotation[sc==14] <-"MiG_Adap2"
Annotation[sc==16] <-"MiG_Abi3"

seu$Annotation <- Annotation

pdf(paste0(plot_dir, "UMAP_Annotation.pdf"))
DimPlot(seu, reduction = "umap", group.by = "Annotation", raster=FALSE, label = TRUE)
dev.off()


Annotation <- rep("Unknown",dim(seu@meta.data)[1])
sc <- seu$SCT_snn_res.0.3

Annotation[sc==1] <-"ExN"
Annotation[sc==2] <-"AST"
Annotation[sc==3] <-"ODC"
Annotation[sc==4] <-"InN"
Annotation[sc==5] <-"InN"
Annotation[sc==6] <-"ExN"
Annotation[sc==7] <-"OPC"
Annotation[sc==8] <-"MiG"
Annotation[sc==9] <-"VLMC"
Annotation[sc==10] <-"Endo"
Annotation[sc==11] <-"ExN"
Annotation[sc==12] <-"ExN"
Annotation[sc==13] <-"MiG"
Annotation[sc==14] <-"MiG"
Annotation[sc==16] <-"MiG"

seu$L1 <- Annotation

pdf(paste0(plot_dir, "UMAP_L1.pdf"))
DimPlot(seu, reduction = "umap", group.by = "L1", raster=FALSE, label = TRUE)
dev.off()

saveRDS(seu, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524_Anno.rds")
#seu <- readRDS("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/RNA/seu_042524_Anno.rds")

#now, make the barcodes annotation stuff
sample_list$good_prefix  <- sprintf("%02d", sample_list$Prefix)
meta <- seu@meta.data
meta$barcodes = rownames(meta)
meta_merged <- merge(meta, sample_list, by.x = "lib", by.y = "RNA")

barcodes_raw <- do.call(rbind, strsplit(meta$barcodes, "_"))[,1]

meta_merged$good_barcodes <- paste0(barcodes_raw, "_", meta_merged$good_prefix)

assignment <- meta_merged[c("good_barcodes", "L1")]

write.table(assignment, "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/RNA/Seurat/merged/split_bam_dir.txt", sep = "\t", col.names = F, quote = F, row.names = F)



group_colors <- c(
  ExN = "#6CA6CD",
  InN = "#912CEE",
  Endo = "#CD8500",
  VLMC = "#EEDD82",
  MiG = "#EE9A00",
  AST = "#CD0000",
  ODC = "#8B4726",
  OPC = "#8B8989"
)

# Create the violin plot with the specified colors
pdf(paste0(plot_dir, "Vln_UMI_SCT_L1.pdf"), width = 6, height = 4)
VlnPlot(seu, group.by = "L1", features = "nCount_RNA", pt.size = 0, log = TRUE, cols = group_colors) +  xlab("")
dev.off()