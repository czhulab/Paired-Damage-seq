d=$1
r=$2

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
../preproc combine_384plex ${d}
zcat ${d}_combined_TN5.fq.gz | bowtie /gpfs/commons/groups/zhu_lab/czhu/genome_references/cell_id_384/Paired_Tag3_384_ID_ref - --norc -m 1 -v 1 -p 4 -S ${d}_BC.sam
../preproc convert ${d}_BC.sam
rm ${d}_BC.sam

trim_galore ${d}_BC_cov.fq.gz
mv ${d}_BC_cov_trimmed.fq.gz ../02.trimmed/

### preproc RNA fastq file


../preproc combine_384plex ${r}
zcat ${r}_combined_RNA.fq.gz | bowtie /gpfs/commons/groups/zhu_lab/czhu/genome_references/cell_id_384/Paired_Tag3_384_ID_ref - --norc -m 1 -v 1 -p 4 -S ${r}_BC.sam
../preproc convert ${r}_BC.sam
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

../reachtools rmdup2 ${d}_mm10_bwt2_sorted.bam
mv ${d}_mm10_bwt2_sorted.bam ../03.mapping_mm10/
mv ${d}_mm10_bwt2_sorted_rmdup.bam ../03.mapping_mm10/
#
STAR  --runThreadN 4 --genomeDir /gpfs/commons/groups/zhu_lab/czhu/genome_references/refdata-cellranger-mm10-3.0.0/star --readFilesIn ${r}_BC_cov_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${r}_mm10_ --outSAMtype BAM Unsorted
#
samtools view -h -F 256 ${r}_mm10_Aligned.out.bam -b > ${r}\_clean.bam
samtools sort ${r}\_clean.bam -o ${r}_mm10_sorted.bam
../reachtools rmdup2 ${r}_mm10_sorted.bam

#mv ${r}_mm10_Aligned.out.bam ../03.mapping_mm10/
mv ${r}\_clean.bam ../03.mapping_mm10/
mv ${r}_mm10_sorted.bam ../03.mapping_mm10/
mv ${r}_mm10_sorted_rmdup.bam ../03.mapping_mm10/

cd ../03.mapping_mm10/

../reachtools bam2Mtx2 ${r}_mm10_sorted_rmdup.bam /gpfs/commons/groups/zhu_lab/czhu/genome_references/annotations/mm10.annotation

mv ${r}*mtx2 ../04.matrices/




#