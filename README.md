# wanglab-pipelines
Pipelines for the Wang Lab

## Key:
:white_check_mark: - Done and running  
:running: - In Progress  
:finnadie: - Trying, but to no avail  
:shipit: - Low priority  


## Current Pipelines Available:
### NGS alignment + specifics

| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| atac_alignment.sh | Aligns ATAC-Seq data to hg38 genome and calls peaks using the Genrich peakcaller | :white_check_mark: |
| atac_qc.sh | Quality control metrics for ATAC data such as TSS, PT, and NFR score | :white_check_mark: |
| cr_alignment | Aligns CUT&RUN data to hg38 genome and calls peaks | :white_check_mark: |
| cr_epic2 | Calls peaks using the epic2 peakcaller normalized to IgG | :running: |
| rna_alignment.sh | Aligns RNA-Seq data to hg38 genome using the STAR aligner, and calls peaks using epic2 | :white_check_mark: |
| qc.sh | Fastqc for fastq and bam files | :white_check_mark: |
  
### Downstream analysis
| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| ChIPseeker.sh | Peak annotation and pathway analysis | :white_check_mark: |
| ChromHMM.sh | Annotates chromatin states | :shipit: |
| diffBind_createObject | Creates DiffBind object and saves it for faster downstream processing | :white_check_mark: |
| homerMotif.sh | Generates homer motifs | :white_check_mark: |
| heatmaps.sh | Generates heatmaps but only from certain annotated peaks (e.g. only distal intergenic regions) | :running: |
| ROSE.sh | Annotates super enhancers | :shipit: |

### Miscellaneous
| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| move.sh | Moves fastq files from their own subfolder into main folder | :white_check_mark: |
