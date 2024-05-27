
#downsample peak calling expt 
# Directory containing the BAM files
input_dir="/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/06.20230928_DNA_reproce/03.mapping_hg38/"

# Output directory for the downsampled BAM files
downsample_dir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MACS3_downsample/"
# Output directory for the peaks
output_dir="/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MACS3_downsample/MACS3/"
mkdir -p "$downsample_dir"
mkdir -p $output_dir
# Target number of reads to downsample to
target_reads=50000000  # Example: 10 million reads

# Loop through each BAM file in the directory
for BAM_file in "$input_dir"/Damage_DNA_*_sorted_srmdup.bam; do
    echo "Processing $BAM_file"
    
    # Extract the base name of the file (cell type name)
    cell_type=$(basename "$BAM_file" .bam)
    
    # Downsample BAM file
    downsampled_bam="$downsample_dir/${cell_type}_downsampled.bam"
    echo "Downsampling to $target_reads reads."
    samtools view -bs $(echo "$target_reads $(samtools view -c "$BAM_file")" | awk '{print $1/$2}') "$BAM_file" > "$downsampled_bam"
    samtools sort -o "$downsampled_bam"
    samtools index "$downsampled_bam"

    # Call macs3 callpeak for the downsampled file
    echo "Calling peaks for $cell_type"
    macs3 callpeak \
        -t "$downsampled_bam" \
        --outdir "$output_dir" \
        -B \
        -q 0.01 \
        -f BAM \
        -nomodel \
        -n "$cell_type" \
        --shift -200 \
        --extsize 400 \
        --keep-dup all
done

for narrowPeak_file in "$output_dir"/*_peaks.narrowPeak; do
    # Extract the base name of the file (cell type name)
    cell_type=$(basename "$narrowPeak_file" _peaks.narrowPeak)

    # Count the number of lines (peaks) in the narrowPeak file
    echo -n "Number of peaks for $cell_type: "
    wc -l < "$narrowPeak_file"
done
