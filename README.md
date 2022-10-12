# TMB conversion
Use two layer Poisson model to convert TMB between genomic regions.
## Required input files
Three input files are required to run the model:
* input_mutations: Detected mutations in VCF format.
* sequencing_panel: Input panel regions in BED (Browser Extensible Data) format.
* expected_panel: Expected genomic regions in BED (Browser Extensible Data) format.

## Step 1: Prepare data to train the model
```
intersectBed -a GDC-PANCAN.mutect2_snv_hg19_all.bed -b sequencing_panel | perl get_mut.pl - > panel_mut.txt
intersectBed -a GDC-PANCAN.mutect2_snv_hg19_nosilent.bed -b expected_panel | perl get_mut.pl - > expected_panel_mut.txt
nucBed -fi hg19.fasta -bed sequencing_panel|sed 1d |perl -lane 'BEGIN{$a=0;$t=0;$g=0;$c=0;}{$a+=$F[5];$c+=$F[6];$g+=$F[7];$t+=$F[8]}END{print "$a\t$c\t$g\t$t"}' > panel_nuc.txt
nucBed -fi hg19.fasta -bed expected_panel|sed 1d |perl -lane 'BEGIN{$a=0;$t=0;$g=0;$c=0;}{$a+=$F[5];$c+=$F[6];$g+=$F[7];$t+=$F[8]}END{print "$a\t$c\t$g\t$t"}' > expected_panel_nuc.txt
```

## Step 2: Run the main script
```R
Rscript two_layer_Poisson_model.R input_mutations panel_mut.txt panel_nuc.txt expected_panel_mut.txt expected_panel_nuc.txt
```

## Note
* Several expected genomic regions have been deposited in the data folder.
* Scripts get_mut.pl and two_layer_Poisson_model.R are in the scripts folder.
* bedtools need to be installed beforehand. 



