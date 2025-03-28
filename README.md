Analysis of Paired-Damage-seq datasets
=====

Paired-Damage-seq is a single-cell multiomics sequencing technique for joint analysis of transcriptome with oxidative DNA damage and single-stranded DNA breaks.

Based on _in situ_ labeling with the "base excision repair" proteins, damaged DNA sites, including 8-oxoG, AP sites, nicks and gaps, can be labelled with biotinylated dUTP and enabled the downstream targeted tagmentation with anti-biotin antibodies and protein A-Tn5 fusion protein. As a new member of the "Paired series" multiomics techniques, Paired-Damage-seq is also based on ligation-based combinatorial barcoding (first introduced by [SPLiT-seq](https://www.science.org/doi/10.1126/science.aam8999)), offering ultra-high-throughput, cost-effective, single-cell indexing without the requirements for specific instruments.

<img width="1000" alt="image" src="https://github.com/user-attachments/assets/1a73d64f-63fd-4e24-8059-a5051a2056a3" />

Using paired transcriptome, we can perform computational "sorting" of cells and analyze the regions most frequently damaged across different cell types and states. Such selective genome vulnerability displays associations with loss of epigenetic memory over time, and could contribute to disease risks.

We used customized barcodes designs and the codes here are specifically for pre-processing of Paired-Damage-seq datasets. If you are using different sets of barcodes sequences, you may need to prepare your own barcode whitelist files.

Preparation
-----
The codes and barcodes whitelist files used for analysis of Paired-Damage-seq were organized into four parts.
- The [pre-processing pipeline](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing) to extract cell barcodes, map to reference genome, and generate RNA cell-to-gene count matrix.
- The codes for [analysis of HeLa Paired-Damage-seq datasets](https://github.com/czhulab/Paired-Damage-seq/tree/main/02.Analysis_HeLa).
- The codes for [analysis of mouse brain Paired-Damage-seq datasets](https://github.com/czhulab/Paired-Damage-seq/tree/main/03.Analysis_Brain).
- The [cell barcode reference and gene annotation reference](https://github.com/czhulab/Paired-Damage-seq/tree/main/resources) (for generation of RNA count matrix) are also uploaded to the "resource" folder.

For additional resources, please find: [Additional resources.](#additional-resources)


> [!NOTE]
> Please have the following softwares installed before running the preprocessing pipeline.

  <details>

  _**<summary> Package Requirements </summary>**_
   If you have previously set up the environment for analysis of [Paired-Tag](https://github.com/cxzhu/Paired-Tag/tree/master) or [SIMPLE-seq](https://github.com/cxzhu/SIMPLE-seq), you may not need to re-install all of them.
  
  Name | Link
  --- | ---
  bowtie 1.x  | http://bowtie-bio.sourceforge.net/index.shtml
  bowtie 2.x  | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  samtools (>1.3) | https://www.htslib.org
  STAR | https://github.com/alexdobin/STAR
  Trim_galore | https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
  FastQC (Optional) | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  > - Generally we do not limit to specific versions of the above softwares, as long as the parameters format in bash scripts here match with the installed versions.
  > - An updated GCC complier is required. 
  </details>

A customized code was prepared for extracting & converting the cell barcodes from the "Paired" series datasets. To compile the code, plesae follow the following steps:
```bash
# Download and uncompressed the scripts
wget https://github.com/czhulab/Paired-Damage-seq/archive/refs/heads/main.zip
unzip main.zip

# Compile the "preprocessing"
cd Paired-Damage-seq-main/01.Preprocessing/preproc
sh make.sh

# Compile the "reachtools":
cd ../reachtools
sh make.sh

# Build the Cellular Barcodes reference for bowtie 1.x:
cd ../../resources
unzip Paired_Tag3_384_ID_ref.fa.zip
bowtie-build ./Paired_Tag3_384_ID_ref.fa ./Paired_Tag3_384_ID_ref

```



> [!IMPORTANT]
> - We have updated the length and numbers of barcode combinations for the split-and-pool ligation, which is different from [SPLiT-seq](https://www.science.org/doi/10.1126/science.aam8999) and Paired-Tag (2021). A [preproc](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing/preproc) tool is provided for correctly extracting the cellular barcodes.
> - The [reachtools](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing/reachtools) is only used for generating cell-to-genes/bins count matrices. 

Pre-processing of Paired-Damage-seq datasets
---
**Step 1.**  Initial QC (_Optional_).

Similar to the other combinatorial barcoding-based "Paired series" datasets, a quick QC can be done with FastQC software. FastQC will give the summary for several key quality metrics of fastq files generated from Illumina bcl2fastq program. The key metrics are similar to the previous Paired-Tag datasets.

```bash
# If installed fastQC software.
fastqc Sample_1_R1.fq.gz
fastqc Sample_1_R2.fq..gz
```

  <details>

  _**<summary> Representative FastQC report from Paired-Tag dataset </summary>**_
  
  The image below shows the "Per base sequence content" and "Adapter Concent" sections of the FastQC output file from a representative Paired-Tag library.
  <img width="877" alt="image" src="https://github.com/user-attachments/assets/fcc61ff1-3f06-4452-bb2e-91734a7da8e4" />
  
   - As shown in Read1 report, there is a high fraction of G base in the 2nd base, that is expected from the library construction (no trimming are needed for this part as Bowtie2 will handle it properly).
   - For read2, there are 3 base balanced regions (UMI and barcode) between 3 base-inbalenced linkers (as indicated in the image).
     - If the linker regions does not show high fluctuation as in the representative image, that may indicate a low ligation efficiency.
   - For "Adapter Content" section, typically we will expect a low fraction of Nextera adaptor sequence (at 100th bp, 5%-20%; expect higher percentage if sequenced to 150 bp or longer) in Read2 library and negligible adaptor content from Read1 library.
   - Higher fraction of adaptors:
     - In RNA library indicates: amount of N5-Tn5 is too high in 2nd adaptor tagging step of RNA library.
     - In DNA library: tagmentation efficiency (antibody efficiency) are low, may expect low library complexity.
  </details>

**Step 2.**  Barcodes extraction, mapping, and matrix generation.

During this step, the shell script will perform barcodes extraction, reads cleaning & mapping, PCR duplicates removal, and generating the cell-counts matrix with the environment prepared in [Preparation](#preparation) section.
> [!IMPORTANT]
> - Don't forget to update the your paths to ```fastq files```, ```barcode references```, ```genome_reference``` ,```preprocessing```, ```reachtools``` folders, and ```gene_annotation``` file in the "01.Preprocessing/per_run.sh" bash script.
> - Cell barcodes whitelist ```barcode references```, and ```gene_annotation``` for mm10 and hg38 are available in "resources" folder.


```bash
# Run this for individual sub-library

# DNA_ID and RNA_ID corresponding to the prefix of fastq files, for example:

# Sample_01_DNA_R1.fq.gz, Sample_01_DNA_R2.fq.gz, Sample_01_RNA_R1.fq.gz, Sample_01_RNA_R2.fq.gz
# DNA_ID = "Sample_01_DNA"
# RNA_ID = "Sample_01_RNA"


sh per_run.sh ${DNA_ID} ${RNA_ID}
```

For batch job submission, we prepared a simple perl script (```01.Preprocessing/01.submit_run.pl```) and an example sample table (```sample_list.txt```). Please modify [this line](https://github.com/czhulab/Paired-Damage-seq/blob/9fec3da3f4a80261b2281da11c3e20b433a914ff/01.Preprocessing/01.submit_run.pl#L12) with your batch submission script.

  <details>
    
  _**<summary>The key output files after this step includes: </summary>**_

   > - ```01.rawdata/*combined_DNA.fq.gz```: Extracted barcode reads with DNA restriction cutting sites, which are derived from PAT tagmentation.
   > - ```01.rawdata/*combined_RNA.fq.gz```: Extracted barcode reads with RNA restriction cutting sites, which are derived from reverse transcription.
   > - ```01.rawdata/*combined_UND.fq.gz```: Extracted barcode reads that cannot be assigned to DNA and RNA modalities, possible due to PCR/sequencing errors, and dimer fragments.
   > - ```01.rawdata/*BC_cov.fq.gz```: Converted fastq files with barcode IDs attached to ReadName lines.
   > - ```02.trimmed/*BC_cov_trimmed.fq.gz```: Cleaned reads that will be used for mapping. An additional optional QC can be performed on them.
   > - ```03.mapping_mm10/*_mm10_sorted.bam```: Mapped DNA and RNA bam files, before PCR duplicates removal.
   > - ```03.mapping_mm10/*_mm10_sorted_rmdup.bam```: Mapped DNA and RNA bam files, after PCR duplicates removal.
   > - ```04.*mtx2/```: Cell-counts matrix for individual sub-libraries, in 10X format.
    
  </details>

**Step 3.**  Pre-filtering, and merging matrices for sub-libraries.

In Paired-Damage-seq, we will aliquote the barcoded nuclei into sub-libraries containing 2-10k cells for library preparation and sequencing. The best approach is to QC & filtering sub-library pairs individually, and then merge them for downstream analyses.
- The matrices files are in standard 10X format, you can use your own scripts to perform this task.
- You can also use our previous [Paired-Tag](https://github.com/cxzhu/Paired-Tag/tree/master?tab=readme-ov-file#3-merge-sub-libraries-for-downstream-analysis) scripts for this step.


<details>
    
_**<summary> Filtering & merging matrices </summary>**_
  
We recommend to filter barcodes with low reads numbers before maerging sub-libraries. The same scripts in [Paired-Tag](https://github.com/cxzhu/Paired-Tag/tree/master?tab=readme-ov-file#3-merge-sub-libraries-for-downstream-analysis) are compatible with Paired-Damage-seq here.
  
- Count & plot reads counts using R: [plot_reads_numbers.R]([rscripts/plot_reads_numbers.R](https://github.com/cxzhu/Paired-Tag/blob/b2f367391aba77b17c833cb5671058ee397a19af/rscripts/plot_reads_numbers.R))

  <img width="250" alt="image" src="https://github.com/user-attachments/assets/36ff403f-9b9d-4bb5-83cd-d8afba04c9e4" />
- Filter matrices pairs using perl: perl [filt_mtx.pl](https://github.com/cxzhu/Paired-Tag/blob/b2f367391aba77b17c833cb5671058ee397a19af/perlscripts/filt_mtx.pl)
  - Do not forget to modify the variables to specific files/prefix in the perl code.
  - Metadata file is generated from the R code above.
 
- Merge matrices using perl: perl [merge_mtx.pl](https://github.com/cxzhu/Paired-Tag/blob/b2f367391aba77b17c833cb5671058ee397a19af/perlscripts/merge_mtx.pl) merge_list.txt
  - Example ```merge_list.txt``` format is annotated in the script.

    <img width="350" alt="image" src="https://github.com/user-attachments/assets/4f21afd9-cb91-48ed-a838-d31108d92790" />
    
</details>
  
> [!CAUTION]
> - When merging sub-libraries, always use **unique prefix** for each DNA-RNA library pairs.
>   - The cells in different sub-libraries may have the same barcodes combinations (BC#3:#2:#1). Merging them without adding sub-library pair-specific prefix (BC#4) may results in barcodes conflicts.
>   - The PCR index will be used as the 4th barcode combination (adding BC#4 -> BC#4:#3:#2:#1) to give sufficient #s of barcodes.


**Step 4.**  Downstream analyses.

We recommend to perform cell clustering on transcriptome profiles of Paired-Damage-seq datasets and then pesudobulk the DNA damage signals.

The computational tools for single-cell genomics are rapidly evoloving and it is impossible to recommend the best ones. Here are some of the softwares we used for these downstream analyses:
- Seurat: https://satijalab.org/seurat/
- Scanpy: https://scanpy.readthedocs.io/en/stable/
- SnapATAC2: https://github.com/kaizhang/SnapATAC2
- SEACell: https://github.com/dpeerlab/SEACells

The code we used to produce the presented figures are deposited to [HeLa cell data](https://github.com/czhulab/Paired-Damage-seq/tree/main/02.Analysis_HeLa) and [Mouse brain data](https://github.com/czhulab/Paired-Damage-seq/tree/main/03.Analysis_Brain). We do not have recommendations for specific versions of the packages and these codes and notebooks are for reference purpose only.
> [!NOTE]
> - Please check with the official documentations for the packages/softwares used in the above analysis.
> - The paths to our original files were kept for record purpose. 


Additional resources
-----
  
  If you are interested in applying our methods or datasets, here are the links for some useful resources:
  - Read out publication: [Paper](https://www.nature.com/articles/s41592-025-02632-3)
  - Download the dataset: [GSE268567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268567)
  - We regularly maintain the step-by-step protocol: [here](Protocol/Protocol_Paired-Damage-seq_Jan2025.pdf).

<details>
  
  _**<summary>Our other techniques: </summary>**_
  
  - Paired-seq: [Paper](https://www.nature.com/articles/s41594-019-0323-x), [Protocol](https://link.springer.com/protocol/10.1007/978-1-0716-2899-7_10), [Codes](https://github.com/cxzhu/Paired-seq), [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4874906)
  - Paired-Tag: [Paper](https://www.nature.com/articles/s41592-021-01060-3), [Protocol](https://github.com/cxzhu/Paired-Tag/tree/master/protocol), [Codes](https://github.com/cxzhu/Paired-Tag), [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152020)
  - Droplet Paired-Tag: [Paper](https://www.nature.com/articles/s41594-023-01060-1), [Protocol](https://protocolexchange.researchsquare.com/article/pex-2310/v1), [Codes](https://github.com/Xieeeee/Droplet-Paired-Tag), [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156683)
  - SIMPLE-seq: [Paper](https://www.nature.com/articles/s41587-024-02148-9), [Protocol](https://github.com/cxzhu/SIMPLE-seq/tree/main/Protocol), [Codes](https://github.com/cxzhu/SIMPLE-seq), [Data](https://www.ncbi.xyz/geo/browse/?view=samples&series=197740)
    
</details>

  
  Please feel free to [contact us](https://czhulab.github.io/contact-us.html) if you have any questions or need anything else.
  

[<img width="85" alt="image" src="https://github.com/user-attachments/assets/66e848d8-b72f-423f-8b9c-1f9eda8038dd"/>](https://czhulab.github.io/index.html) @ New York Genome Center & Weill Cornell Medicine




