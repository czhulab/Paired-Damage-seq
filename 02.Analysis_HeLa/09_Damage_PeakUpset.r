library(dplyr)
library(ComplexHeatmap)
library(cowplot)
library(ggthemes)
library(ggplot2)
library(circlize)
library(GenomicRanges)


# read the files and merge them all together 
files <- list.files("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MACS3_Celltype_Specific", full.names = TRUE)
files  <- grep("narrowPeak", files, value = TRUE)

files_list <- list()
for(f in files){
    temp <- read.table(f, sep = "\t")
    files_list[[length(files_list) + 1]] <- temp
}

############################################################################
####################now, do the peak overlap analysis ###################
celltype_list = c("Early_EPHA7", 
                "Early_FGF2", 
                "Early_NEK7", 
                "Mid", 
                "Late_AKR1C2",
                "Late_EPHA7",
                "Late_PIK3R3", 
                "ZC3HAV1")

peak_list <- lapply(celltype_list, function(c){
    path <- grep(c, files, value = TRUE)
    file <- read.table(path, sep = "\t")
    return(file)
})

peak_list <- lapply(peak_list, function(df) GRanges(seqnames = df[, 1], 
    ranges = IRanges(df[, 2], df[, 3])))

names(peak_list) <- celltype_list
union_matrix <- make_comb_mat(peak_list, mode = "union")
distinct_matrix <- make_comb_mat(peak_list, mode = "distinct")
intersect_matrix <- make_comb_mat(peak_list, mode = "intersect")


pdf("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ComplexHeatmap/upset_union.pdf", 
    height = 5, width = 20)
UpSet(union_matrix)
dev.off()

pdf("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ComplexHeatmap/upset_distinct.pdf", 
    height = 5, width = 30)
UpSet(distinct_matrix)
dev.off()

pdf("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ComplexHeatmap/upset_intersect.pdf", 
    height = 5, width = 30)
UpSet(intersect_matrix)
dev.off()

