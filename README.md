# Schizosaccharomyces-pombe-Chip-seq-analysis

After getting the bed file from the Peak calling, we can annotate and create all downstream analysis from two sets of data

## Overview
This repository contains code and files related to ChIP-seq analysis. It includes scripts, peak files, and result images.

## Files
- `.gitignore`: Specifies files and directories to be ignored by Git.
- `C1_summits.bed`: ChIP-seq peak summit file for condition C1.
- `C2_summits.bed`: ChIP-seq peak summit file for condition C2.
- `Distribution of TF binding loci relative to TSS.png`: Image showing the distribution of transcription factor (TF) binding loci relative to the transcription start site (TSS).
- `Distribution of transcription factor-binding loci.png`: Image showing the distribution of transcription factor-binding loci.
- `Feature Distribution.png`: Image showing the distribution of features.
- `Feature Distribution1.png`: Additional image showing the distribution of features.
- `Feature Distribution2.png`: Additional image showing the distribution of features.
- `LICENSE`: License file specifying the terms and conditions for using the code and files in this repository.
- `Nanog-idr.png`: Image showing the Nanog-idr peaks.
- `Overlap of peaks and annotated genes.png`: Image showing the overlap of peaks and annotated genes.
- `PCF1.png`: Image showing PCF1.
- `PCF2.png`: Image showing PCF2.
- `PCF3.png`: Image showing PCF3.
- `README.md`: This file, providing an overview and instructions for the repository.
- `chip_seq_annotation.R`: R script for ChIP-seq annotation and analysis.
- `upset and vennypie.png`: Image showing the upset plot and venn diagram using pie charts.
- `upset.png`: Image showing the upset plot.
- `vennypie.png`: Image showing the venn diagram using pie charts.

## Usage
1. Clone or download the repository.
2. Install the necessary R packages specified in the code.
3. Place the ChIP-seq summit files (`C1_summits.bed` and `C2_summits.bed`) in the appropriate location.
4. Execute the `chip_seq_annotation.R` script for ChIP-seq annotation and analysis.
5. Refer to the generated images and results for further analysis and interpretation.

## License
The code and files in this repository are licensed under [LICENSE](LICENSE).

## Note
This repository is provided as a reference for ChIP-seq analysis and may require customization or additional files/data for specific use cases.
