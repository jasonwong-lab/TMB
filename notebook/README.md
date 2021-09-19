<<<<<<< HEAD
## Panel tumor mutational burden (TMB) correction
<<<<<<< HEAD

__TMBpredict__ is an R package can optimize TMB that is derived from target region sequencing data. Linear regression is employed to modeling the relationship between whole coding regions derived TMB and panel derived mutations in this tool, and it is based on pan-cancer data of 10,179 samples across 33 cancer types from The Cancer Genome Atlas (TCGA). 

## Installation 

library(devtools) # Make sure that the devtools library is loaded
install_github("jasonwong-lab/TMB")

## Usage guide

tmb_pre(ttype, mut, panel.bed, ftype)
ttype: <strings> Tumor type to be evaluated, should be one of these 33 cancer types.
mut: <strings> Input mutation data in VCF format. Can be multiple VCF compressed in tar.gz.
panel.bed: <strings> Input panel region data, BED (Browser Extensible Data) format is accepted and should at least include three columns: chromosome, start position and end position.
ftype: <strings> The type of input mutation file. "s" for single VCF file and "m" for multiple VCF files compressed in tar.gz.

## Examples

tmb_pre("COAD","sample.vcf","panel.bed","s")
tmb_pre("COAD","vcf.tar.gz","panel.bed","m")
Example data can be found at https://github.com/jasonwong-lab/TMB/tree/master/test
=======
__TMBpredict__ is an R package can optimize TMB that is derived from target region sequencing data. Linear regression is employed to modeling the relationship between whole coding regions derived TMB and panel derived mutations in this tool, and it is based on pan-cancer data of 10,179 samples across 33 cancer types from The Cancer Genome Atlas (TCGA).

## Installation
```R
require(devtools) # Make sure that the devtools library is loaded  
install_github("jasonwong-lab/TMB")`  
```
## Usage guide
```R
library(TMBpredict)
library(GenomicRanges)
library(VariantAnnotation)
TMBpredict(ttype, mut, panel.bed, ftype)
```
* ttype: \<strings\> Tumor type to be evaluated, should be one of these [33 cancer types](https://github.com/jasonwong-lab/TMB/blob/main/Cancer_type.txt).  
* mut: \<strings\> Input mutation data in VCF format. Can be multiple VCF compressed in tar.gz.  
* panel.bed: \<strings\> Input panel region data, BED (Browser Extensible Data) format is accepted and should at least include three columns: chromosome, start position and end position.  
* ftype: \<strings\> The type of input mutation file. "s" for single VCF file and "m" for multiple VCF files compressed in tar.gz.  

## Examples
```R
TMBpredict("COAD","COAD_test_sap.vcf","msk_coding.bed","s")  
TMBpredict("COAD","vcf.tar.gz","msk_coding.bed","m")
```
Example data can be downloaded here:<br>
<a id="raw-url" href="https://raw.githubusercontent.com/jasonwong-lab/TMB/master/test/single_file/COAD_test_sap.vcf">COAD_test_sap.vcf</a><br>
<a id="raw-url" href="https://raw.githubusercontent.com/jasonwong-lab/TMB/master/test/multiple_file/vcf.tar.gz">vcf.tar.gz</a><br>
<a id="raw-url" href="https://raw.githubusercontent.com/jasonwong-lab/TMB/master/test/single_file/msk_coding.bed">msk_coding.bed</a>

## Note
Web-based shiny App can be found at: https://cancergenomics-explore.shinyapps.io/shiny_tmb/
>>>>>>> 1e00abcc465be1436f85e656dadce6db38174510
=======

## TMB prediction App

## Introduction

This TMB (tumor mutational burden) prediction tool can harmonize TMB that is derived from target region and whole exome sequencing data. Linear regression is employed to modeling the relationship between whole coding regions derived TMB and panel derived mutations in this tool, and it is based on pan-cancer data of 10,179 samples across 33 cancer types from The Cancer Genome Atlas (TCGA). From the TCGA data module, you can get the landscape of TMB across different cancer types. By uploading mutation file and panel bed file, you can obtain the normalised TMB prediction. 
 
<p align="center"><img src="ui.png"/></p>

## TCGA data module

Mutation burden of each cancer type from TCGA pan-cancer data can be checked in this module. You can get a landscape of TMB across various panels. The region of WES (whole exome sequencing) is defined as all coding regions that are obtained from UCSC table browser. The region of MSK (MSK-IMPACT) and F1CDX (FoundationOne CDx) are defined as coding regions of the genes that are officially released. Please note that the default y-axis limit is 100 mut/Mb, and it can be tuned as required. You also can input a TMB value of interest to see where it is among all the samples. 

## My data module

VCF (variant call format) format is accepted for the input mutation data. For the input panel region data, BED (Browser Extensible Data) format is accepted and should at least include three columns: chromosome, start position and end position. The coordinate is 0-based. After uploading required files, the correlation of panel mutation and whole exome TMB based on TCGA data is provided and the adjusted TMB is reported.

## Availability

The app is available at: https://cancergenomics-explore.shinyapps.io/shiny_tmb/

## Note

The stand-alone R package can be downloaded at: https://github.com/jasonwong-lab/TMB/tree/master

>>>>>>> 48c7d6356ec6c3e89c3284f69f6f60c7fa4741ce
