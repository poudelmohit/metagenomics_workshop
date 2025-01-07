* I apologize, I had to move the contents of this project to [this repo](https://github.com/poudelmohit/env-metagenomics-analysis)*
* You may also visit the [project website](https://poudelmohit.github.io/env-metagenomics-analysis/) for details.*


# Metagenomic Analysis of Environmental Samples
  This repository contains the pipeline for analyzing metagenomic data from environmental samples using various bioinformatics tools on a Linux system. The project involves quality control, read trimming, assembly, binning, taxonomic profiling, and visualization using open-source tools like FastQC, Trimmomatic, SPAdes, MaxBin, MetaPhlAn, and Krona.


## Notes
Due to the large file size (>100 MB), some files (such as raw reads and assembly results) are not available in this repository. You can generate these files by following the steps outlined in the workflow.

## Workflow Overview
1. Quality Control: FastQC was used to check the quality of the raw FASTQ files.
2. Read Trimming: Trimmomatic was applied to remove low-quality bases and adapter sequences.
3. Assembly: SPAdes was used to assemble the reads into contigs and scaffolds.
4. Binning: MaxBin was used for metagenomic binning of the assembled sequences.
5. Taxonomic Assignment: Kraken2 was used for taxonomic profiling of the metagenomic samples.
6. Visualization: Krona was used to visualize the taxonomic profiles.
7. Diversity calculation: Ecological and statistical anlayses were done with R
   
*All analyses were performed using a Linux HPC cluster environment, ensuring efficient processing of the data*

## Contact
For any questions or further details about the analysis, please contact [**Mohit Poudel**](https://poudelmohit.github.io)

## Special Thanks:
[**Alvarado-Serrano Lab, Ohio University**](https://alvarado-s.weebly.com)
![Logo](https://github.com/poudelmohit/portfolio/blob/main/assets/lablogo-small.png)

## Additional Resources:
paper to follow: https://elifesciences.org/articles/49816

data available at: https://zenodo.org/records/7010950







