# Directory containing the BED files
#input_dir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/bed_celltype/"
#input_dir="/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/11.Brain_dataset/08.split_bams_20240312_Zhenkun/01.DNA_split_Zhenkun/"
input_dir="/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/11.Brain_dataset/11.split_bams_20240429_Zhenkun_4/01.DNA_split_Zhenkun/"

# Output directory for the peaks
#output_dir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/peaks_celltype/"
output_dir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/MACS_L1_BAM_2mo/"

# Loop through each BED file in the directory
for BAM_file in "$input_dir"/Brain_DNA_*_sorted.bam; do
    # Extract the base name of the file (cell type name)
    echo $BAM_file
    Celltype=$(basename "$BAM_file" .bam)
    # Call macs3 callpeak for each file
    macs3 callpeak \
        -t "$BAM_file" \
        --outdir "$output_dir" \
        -B \
        -q 0.01 \
        -f BAM \
        -n "$Celltype" \
        --shift -200 \
        --nolambda \
        --nomodel \
        --extsize 500 \
        --keep-dup all

done


# do the blacklist 
blacklist="/gpfs/commons/groups/zhu_lab/shared_tools/reference_genome/mm10/mm10-blacklist.v2.bed"
#baseDir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/MACS_L1_BAM"
baseDir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/MACS_L1_BAM_2mo/"

outputDir="${baseDir}/blacklist/"
mkdir $outputDir
cd "$baseDir"

for peakFile in *.narrowPeak; do
    basename=$(basename "$peakFile" _peaks.narrowPeak)
    echo "Processing $peakFile..."
    bedtools intersect -a "$peakFile" \
                       -b /gpfs/commons/groups/zhu_lab/shared_tools/reference_genome/mm10/mm10-blacklist.v2.bed \
                       -v \
                       > "${outputDir}/${basename}_blacklist.narrowPeak"
done

# fix the value
#cd /gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/MACS_L1_BAM/blacklist/
cd $outputDir
cell_types=("OPC" "ODC" "MiG" "InN" "ExN" "Endo" "AST" "VLMC")

for c in "${cell_types[@]}"
do
    echo $c
    awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' Brain_DNA_${c}_merged_sorted_blacklist.narrowPeak > ./${c}_valueFixed_blacklist.narrowPeak
done