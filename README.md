Analysis of Paired-Damage-seq datasets
=====

Single-cell joint analysis of transcriptome with oxidative DNA damage and single-stranded DNA breaks.

![R01_Fig5_Prelim_Damage-seq](https://github.com/user-attachments/assets/666bb112-70e0-49fa-85d8-736a6f13de9a)



Preparation
-----
The codes and whitelist files used for analysis of Paired-Damage-seq were organized into four parts.
- The [pre-processing pipeline](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing) to extract cell barcodes, map to reference genome, and generate RNA cell-to-gene count matrix.
- The codes for [analysis of HeLa Paired-Damage-seq datasets](https://github.com/czhulab/Paired-Damage-seq/tree/main/02.Analysis_HeLa).
- The codes for [analysis of mouse brain Paired-Damage-seq datasets](https://github.com/czhulab/Paired-Damage-seq/tree/main/03.Analysis_Brain).
- The [cell barcode reference and gene annotation reference](https://github.com/czhulab/Paired-Damage-seq/tree/main/resources) (for generation of RNA count matrix) are also uploaded to the "resource" folder.

For additional resources, please refer to: [Additional resources.](#additional-resources)


> [!TIP]
> Please have the following softwares installed before running the preprocessing pipelines.

  <details>

  <summary> Package Requirements </summary>
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
  > - An updated GCC complier is also required. 
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
> - We have updated the length and numbers of combinations of barcodes for the split-and-pool ligation, which is different from Paired-Tag (2021). A [preproc](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing/preproc) tool is provided for correctly extracting the cellular barcodes.
> - The [reachtools](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing/reachtools) is only used for generating cell-to-genes/bins count matrices. 

Pre-processing of Paired-Damage-seq datasets
---
**Step 1.**  Initial QC (_Optional_).

Similar to the other combinatorial barcoding-based "Paired" series datasets, a quick initial QC can be done will FastQC software. FastQC will give a QC summary for several key quality metrics of fastq files generated from Illumina bcl2fastq program. The key metrics are similar to the previous Paired-Tag datasets.

```bash
# If installed fastQC software.
fastqc Sample_1_R1.fq.gz
fastqc Sample_1_R2.fq..gz
```

  <details>

  <summary> Representative QC report from Paired-Tag dataset </summary>
  
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

**Step 2.**  Barcodes extration, mapping, and matrix generation.

During this step, the shell scripts will perform barcodes extraction, reads cleaning & mapping, PCR duplicates removal, and generating the cell-counts matrix with the environment prepared in [Preparation](#preparation) section.
> [!IMPORTANT]
> - Don't forget to update the your paths to ```fastq files```, ```barcode references```, ```preprocessing```, and ```reachtools``` in the "01.Preprocessing/per_run.sh" bash script.

```bash
# Run this for individual sub-library

# DNA_ID and RNA_ID corresponding to the prefix of fastq files, for example:

# Sample_01_DNA_R1.fq.gz, Sample_01_DNA_R2.fq.gz, Sample_01_RNA_R1.fq.gz, Sample_01_RNA_R2.fq.gz
# DNA_ID = "Sample_01_DNA"
# RNA_ID = "Sample_01_RNA"


sh per_run.sh ${DNA_ID} ${RNA_ID}
```





Additional resources
-----
  
  If you are interested in applying Paired-Damage-seq technique or interested in Paired-Damage-seq dataset, here are the links for some useful resources:
  - Read out publication: [Coming soon!]()
  - Downloada the dataset: [GSE268567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268567)
  - We regularly maintain the step-by-step protocol: [here](Protocol/Protocol_Paired-Damage-seq_Jan2025.pdf).
  
  Please feel free to [contact us](https://czhulab.github.io/contact-us.html) if you have any questions or need anything else.
  


