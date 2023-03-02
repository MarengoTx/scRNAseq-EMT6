## Table of Contents

1. [Main Description](#main-description)
2. [File Descriptions](#file-descriptions)
3. [Linked Files](#linked-files)
4. [Installation and Instructions](#installation-and-instructions)

&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;



### **Main Description**
------------------------

This is the GitHub repository for the manuscript titled "A TCR β chain-directed antibody-fusion molecule that activates and expands subsets of T cells and promotes antitumor activity.". The code included in the file titled `marengo_code_for_paper_jan_2023.R` was used to generate the figures from the single-cell RNA sequencing data. 

The following libraries are required for script execution:

```
Seurat
scReportoire
ggplot2
stringr
dplyr
ggridges
ggrepel
ComplexHeatmap

```
&nbsp;
&nbsp;
&nbsp;




### **File Descriptions**
---------------------------

The **code** can be downloaded and opened in **RStudios**. </br>
The `marengo_code_for_paper_jan_2023.R` contains all the code needed to reproduce the figues in the paper.</br>
The `sc_data.rds` file is available at the following address: https://zenodo.org/deposit/7566113 </br>
The `all_res_deg_for_heat.txt` file contains the unfiltered results from DGE anlaysis, also used to create the heatmap with DGE and volcano plots.</br>
The `genes_for_heatmap_fig5F.xlsx` contains the genes included in the heatmap in figure 5F.</br>



&nbsp;
&ensp;
&ensp;

### **Linked Files**
---------------------

This repository contains code for the analysis of single cell RNA-seq dataset. The dataset contains raw FASTQ files, as well as, the aligned files that were deposited in GEO. The **Rdata** or `.Rds` file was deposited in Zenodo. Provided below are descriptions of the linked datasets:

1. Gene Expression Omnibus (GEO) ID: [GSE223311](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223311)
   - **Title**: Gene expression profile at single cell level of CD4+ and CD8+ tumor infiltrating lymphocytes (TIL) originiating from the EMT6 tumor model from mSTAR1302 treatment.
   - **Description**: This submission contains the `matrix.mtx`, `barcodes.tsv`, and `genes.tsv` files for each replicate and condition, corresponding to the aligned files for single cell sequencing data. 
   - **Submission type**: Private. In order to gain access to the repository, you must use a [reviewer token](https://www.ncbi.nlm.nih.gov/geo/info/reviewer.html).

&ensp;

2. Sequence read archive (SRA) repository ID: SRX19088718 and SRX19088719
    - **Title**: Gene expression profile at single cell level of CD4+ and CD8+ tumor infiltrating lymphocytes (TIL) originiating from the EMT6 tumor model from mSTAR1302 treatment 
   - **Description**:  This submission contains the **raw sequencing** or `.fastq.gz` files, which are tab delimited text files. 
   - **Submission type**: Private. In order to gain access to the repository, you must use a [reviewer token](https://www.ncbi.nlm.nih.gov/geo/info/reviewer.html).
   
&ensp;


3. Zenodo DOI: [10.5281/zenodo.7566113](https://zenodo.org/deposit/7566113)
   - **Title**: A TCR β chain-directed antibody-fusion molecule that activates and expands subsets of T cells and promotes antitumor activity. 
   - **Description**:  This submission contains the **Rdata** or `.Rds` file, which is an R object file. This is a necessary file to use the code. 
   - **Submission type**: Restricted Acess. In order to gain access to the repository, you must contact the author.

&nbsp;
&ensp;


&nbsp;
&ensp;




### **Installation and Instructions**
--------------------------------------
The code included in this submission requires several essential packagee, as listed above. Please follow these instructions to install these packages:

&nbsp;

> Ensure you have R version 4.1.2 or higher for compatibility. 

> Although it is not essential, you can use R-Studios (Version 2022.12.0+353 (2022.12.0+353)) for accessing and executing the code. 

&nbsp;

1. Download the **RData** or `.Rds` file from [Zenodo](https://zenodo.org/deposit/7566113) (Zenodo DOI: 10.5281/zenodo.7566113).
2. Open [R-Studios](https://www.rstudio.com/tags/rstudio-ide/) or a similar integrated development environment (IDE) for R. 
3. Set your working directory to where the following files are located:
   - `marengo_code_for_paper_jan_2023.R`
   - `Install_Packages.R`
   - `sc_data.rds`
   - `genes_for_heatmap_fig5F.xlsx`
   - `all_res_deg_for_heat.txt`

&nbsp;

You can use the following code to set working directory in R:

> setwd(directory)

&nbsp;

4. Open the file titled `Install_Packages.R` and execute it in R IDE. This script will attempt to install all the necessary pacakges, and its dependencies in order to set up an environment where the code in `marengo_code_for_paper_jan_2023.R` can be exectued. 
5. Once the `Install_Packages.R` script has been successfully executed, re-start R-Studios or your IDE of choice. 
6. Open the file `marengo_code_for_paper_jan_2023.R` file in R-studios or your IDE of choice. 
7. Execute commands in the file titled `marengo_code_for_paper_jan_2023.R` in R-Studios or your IDE of choice to generate the plots. 


&nbsp;
&nbsp;
&nbsp; 
&nbsp;
&nbsp;
&ensp;
&ensp;



