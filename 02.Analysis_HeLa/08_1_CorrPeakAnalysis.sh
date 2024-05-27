
chrom_sizes="/gpfs/commons/groups/zhu_lab/shared_tools/reference_genome/hg38/hg38.chrom.sizes"
peak_path="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MACS3_Treatment/filtered_blacklist/multiinner_merged.bed"
outdir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/AllPeak"
infile="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MACS3_Treatment/filtered_blacklist/multiinner_merged_sorted.bed"

# Sort the peaks file
bedtools sort -i "$peak_path" -g "$chrom_sizes" > "$infile"

# Array of input BED files
# Array of input BED files
bed_files=(
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF495VDP_H3K27ac_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF232FMD_H3K27me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/Enhancer/Hela_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF712ATO_H3K9me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF360CQR_H3K4me1_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/EnhancerWithSTR_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE_ENCFF248WXB_H3K36me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE_ENCFF310XFO_H3K4me3_sorted.bed"
)


# Array of output file names
output_names=(
    "UnionPeak_H3K27ac.txt"
    "UnionPeak_H3K27me3.txt"
    "UnionPeak_enhancer.txt"
    "UnionPeak_H3K9me3.txt"
    "UnionPeak_H3K4me1.txt"
    "UnionPeak_STREnhancer.txt"
    "UnionPeak_H3K36me3.txt"
    "UnionPeak_H3K4me3.txt"
)

# Loop through the arrays
for i in "${!bed_files[@]}"; do
    bed_file="${bed_files[$i]}"
    output_name="${output_names[$i]}"
    output="${outdir}/${output_name}"
    bedtools fisher -a "$infile" -b "$bed_file" -g "$chrom_sizes" > "$output"
done

# Sort the peaks file
bedtools sort -i "$peak_path" -g "$chrom_sizes" > "$infile"
mkdir -p $outdir

# Array of output file names
output_names=(
    "UnionPeak_H3K27ac.txt"
    "UnionPeak_H3K27me3.txt"
    "UnionPeak_enhancer.txt"
    "UnionPeak_H3K9me3.txt"
    "UnionPeak_H3K4me1.txt"
    "UnionPeak_STREnhancer.txt"
    "UnionPeak_H3K36me3.txt"
    "UnionPeak_H3K4me3.txt"
)

# Loop through the arrays
for i in "${!bed_files[@]}"; do
    bed_file="${bed_files[$i]}"
    output_name="${output_names[$i]}"
    output="${outdir}/${output_name}"
    bedtools fisher -a "$infile" -b "$bed_file" -g "$chrom_sizes" > "$output"
done



############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################

############################################################################################################################
############################################################################################################################
# Array of input BED files
bed_files=(
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF495VDP_H3K27ac_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF232FMD_H3K27me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/Enhancer/Hela_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF712ATO_H3K9me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE4_ENCFF360CQR_H3K4me1_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/EnhancerWithSTR_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE_ENCFF248WXB_H3K36me3_sorted.bed"
    "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/DNA_only/Peaks_analysis/ChromatinMarks/ENCODE_ENCFF310XFO_H3K4me3_sorted.bed"
)


#sort the peaks after downsampling
chrom_sizes="/gpfs/commons/groups/zhu_lab/shared_tools/reference_genome/hg38/hg38.chrom.sizes"
peak_path="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/DS/TopLinkFixPeaksTopFDR_042324.bed"
outdir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/FDR/Positive/"
infile="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/DS/TopLinkFixPeaksTopFDR_042324_sorted.bed"

bedtools sort -i "$peak_path" -g "$chrom_sizes" > "$infile"
mkdir -p $outdir
# Array of output file names
output_names=(
    "H3K27ac.txt"
    "H3K27me3.txt"
    "enhancer.txt"
    "H3K9me3.txt"
    "H3K4me1.txt"
    "STREnhancer.txt"
    "H3K36me3.txt"
    "H3K4me3.txt"
)

# Loop through the arrays
for i in "${!bed_files[@]}"; do
    bed_file="${bed_files[$i]}"
    output_name="${output_names[$i]}"
    output="${outdir}/${output_name}"
    bedtools fisher -a "$infile" -b "$bed_file" -g "$chrom_sizes" > "$output"
done

#do the negative 
chrom_sizes="/gpfs/commons/groups/zhu_lab/shared_tools/reference_genome/hg38/hg38.chrom.sizes"
peak_path="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/DS/TopLinkFixPeaksBottomFDR_042324.bed"
outdir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/FDR/Negative/"
infile="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/DS/TopLinkFixPeaksBottomFDR_042324_sorted.bed"

bedtools sort -i "$peak_path" -g "$chrom_sizes" > "$infile"
mkdir -p $outdir

# Loop through the arrays
for i in "${!bed_files[@]}"; do
    bed_file="${bed_files[$i]}"
    output_name="${output_names[$i]}"
    output="${outdir}/${output_name}"
    bedtools fisher -a "$infile" -b "$bed_file" -g "$chrom_sizes" > "$output"
done

