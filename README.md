## Panel tumor mutational burden (TMB) correction
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
Example data can be downloaded here:

[COAD_test_sap.vcf](https://raw.githubusercontent.com/jasonwong-lab/TMB/master/test/single_file/COAD_test_sap.vcf)

[vcf.tar.gz](https://github.com/jasonwong-lab/TMB/blob/master/test/multiple_file/vcf.tar.gz)

[msk_coding.bed](https://github.com/jasonwong-lab/TMB/blob/master/test/multiple_file/msk_coding.bed)

