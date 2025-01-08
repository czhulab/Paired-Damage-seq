d=$1
r=$2

path_to_data="path-to-your-fastq-files-folder/"
path_to_ref="path-to-the-prepared-cell-id-reference-folder/"
path_to_preproc="path-to-the-compiled-preproc-folder/"
path_to_reachtools="path-to-the-compiled-reachtools-folder/"
path_to_gene_annotation_file="path-to-the-gene-annotation-file-from-resources/****.annotation"

### create directories

if [ ! -d "02.trimmed" ]; then
	mkdir -p "02.trimmed"
fi

if [ ! -d "03.mapping_mm10" ]; then
	mkdir -p "03.mapping_mm10"
fi

if [ ! -d "04.matrices" ]; then
	mkdir -p "04.matrices"
fi


### preproc DNA fastq file

cd 01.rawdata
${path_to_preproc}/preproc combine_384plex ${d}
zcat ${d}_combined_TN5.fq.gz | bowtie ${path_to_ref}/Paired_Tag3_384_ID_ref - --norc -m 1 -v 1 -p 4 -S ${d}_BC.sam
${path_to_preproc}/preproc convert ${d}_BC.sam
rm ${d}_BC.sam

trim_galore ${d}_BC_cov.fq.gz
mv ${d}_BC_cov_trimmed.fq.gz ../02.trimmed/

### preproc RNA fastq file


${path_to_preproc}/preproc combine_384plex ${r}
zcat ${r}_combined_RNA.fq.gz | bowtie ${path_to_ref}/Paired_Tag3_384_ID_ref - --norc -m 1 -v 1 -p 4 -S ${r}_BC.sam
${path_to_preproc}/preproc convert ${r}_BC.sam
rm ${r}_BC.sam

trim_galore ${r}_BC_cov.fq.gz
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNNNNNNN ${r}_BC_cov_trimmed.fq.gz
mv ${r}_BC_cov_trimmed_trimmed.fq.gz ../02.trimmed/
#
#### mapping
cd ../02.trimmed
#
bowtie2 -x /gpfs/commons/groups/zhu_lab/czhu/genome_references/mm10/mm10 -U ${d}_BC_cov_trimmed.fq.gz --no-unal -p 4 -S ${d}_mm10_bwt2.sam
samtools sort ${d}_mm10_bwt2.sam -o ${d}_mm10_bwt2_sorted.bam -@ 4
rm ${d}_mm10_bwt2.sam

${path_to_reachtools}/reachtools rmdup2 ${d}_mm10_bwt2_sorted.bam
mv ${d}_mm10_bwt2_sorted.bam ../03.mapping_mm10/
mv ${d}_mm10_bwt2_sorted_rmdup.bam ../03.mapping_mm10/
#
STAR  --runThreadN 4 --genomeDir /gpfs/commons/groups/zhu_lab/czhu/genome_references/refdata-cellranger-mm10-3.0.0/star --readFilesIn ${r}_BC_cov_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${r}_mm10_ --outSAMtype BAM Unsorted
#
samtools view -h -F 256 ${r}_mm10_Aligned.out.bam -b > ${r}\_clean.bam
samtools sort ${r}\_clean.bam -o ${r}_mm10_sorted.bam
${path_to_reachtools}/reachtools rmdup2 ${r}_mm10_sorted.bam

#mv ${r}_mm10_Aligned.out.bam ../03.mapping_mm10/
mv ${r}\_clean.bam ../03.mapping_mm10/
mv ${r}_mm10_sorted.bam ../03.mapping_mm10/
mv ${r}_mm10_sorted_rmdup.bam ../03.mapping_mm10/

cd ../03.mapping_mm10/

${path_to_reachtools}/reachtools bam2Mtx2 ${r}_mm10_sorted_rmdup.bam ${path_to_gene_annotation_file}

mv ${r}*mtx2 ../04.matrices/




#
