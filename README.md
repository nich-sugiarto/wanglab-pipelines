# wanglab-pipelines
Pipelines for the Wang Lab

## Key:
:white_check_mark: - Done and running  
:running: - In Progress  
:hammer: - Working, but not needs to be fixed  
:finnadie: - Trying, but to no avail  
:shipit: - Low priority  


## Current Pipelines Available:
### NGS alignment + specifics

| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| atac_alignment.sh | Aligns ATAC-Seq data to hg38 genome and calls peaks using the Genrich peakcaller | :white_check_mark: |
| atac_qc.sh | Quality control metrics for ATAC data such as TSS, PT, and NFR score | :white_check_mark: |
| cr_alignment | Aligns CUT&RUN data to hg38 genome and calls peaks | :white_check_mark: |
| cr_epic2 | Calls peaks using the epic2 peakcaller normalized to IgG | :white_check_mark: |
| rna_alignment.sh | Aligns RNA-Seq data to hg38 genome using the STAR aligner, and calls peaks using epic2 | :white_check_mark: |
| qc.sh | Fastqc for fastq and bam files | :white_check_mark: |
  
### Downstream analysis
| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| ChIPseeker.sh | Peak annotation and pathway analysis | :hammer: |
| ChromHMM.sh | Annotates chromatin states | :shipit: |
| diffBind_analyze | Finds and reports differential peak regions using both DESEQ2 and EdgeR | :white_check_mark: |
| diffBind_ChIPseeker | Annotates regions with ChIPseeker, and performs GO analysis | :hammer: |
| diffBind_createObject | Creates DiffBind object and saves it for faster downstream processing | :white_check_mark: |
| diffBind_HomerMotif | Reports motifs associated with differentially bound regions | :white_check_mark: |
| EChO | Fragment analysis for CUT&RUN data | :white_check_mark: |
| homerMotif.sh | Generates homer motifs | :white_check_mark: |
| heatmaps.sh | Generates heatmaps but only from certain annotated peaks (e.g. only distal intergenic regions) | :shipit: |
| ROSE.sh | Annotates super enhancers | :shipit: |
| rna_multi.sh | DESeq on all RNA samples | :white_check_mark: |
| rna_twoThree.sh | DESeq on two RNA samples, each with three replicates | :white_check_mark: |

### Miscellaneous
| Pipeline name | Description | Status |
| ------------- | ----------- | ------ |
| move.sh | Moves fastq files from their own subfolder into main folder | :white_check_mark: |
