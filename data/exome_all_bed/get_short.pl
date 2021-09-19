#!/usr/bin/perl -w 
use strict;
open IN,"/home/fanghu/Projects/tmb/shiny/package/TMBpredict/data/ttype_33";
while(<IN>){
	chomp;
	my $type=$_;
	open BED,"/home/fanghu/Projects/tmb/shiny/package/TMBpredict/data/exome_all_bed_v1/$type.bed";
	open OUT,">$type.bed";
	while(<BED>){
		chomp;
		my @field=split;
		print OUT "$field[0]\t$field[1]\t$field[2]\t$field[5]\n";
	}
	close BED;
	close OUT;
}
close IN;
