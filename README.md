Analysis of Paired-Damage-seq datasets
=====

Single-cell joint analysis of transcriptome with oxidative DNA damage and single-stranded DNA breaks.

![R01_Fig5_Prelim_Damage-seq](https://github.com/user-attachments/assets/666bb112-70e0-49fa-85d8-736a6f13de9a)



Preparation
-----
The codes and whitelist files used for analysis of Paired-Damage-seq are summarized into three parts.
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

**Step 1:** Compile reachtools for barcode extraction:


> [!IMPORTANT]
> We have updated the length and numbers of combinations of barcodes for the split-and-pool ligation, which is different from Paired-Tag (2021). Thus, the updated [reachtools](https://github.com/czhulab/Paired-Damage-seq/tree/main/01.Preprocessing/reachtools) is required for correctly extracting the cellular barcodes.



Additional resources
-----
  
  If you are interested in applying Paired-Damage-seq technique or interested in Paired-Damage-seq dataset, here are the links for some useful resources:
  - Read out publication: [Coming soon!]()
  - Downloada the dataset: [GSE268567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268567)
  - We regularly maintain the step-by-step protocol: [here](Protocol/Protocol_Paired-Damage-seq_Jan2025.pdf).
  
  Please feel free to [contact us](https://czhulab.github.io/contact-us.html) if you have any questions or need anything else.
  


