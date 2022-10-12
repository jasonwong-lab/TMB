#!/usr/bin/perl -w 
use strict;
my ($in,$out)=@ARGV;
open IN,$in;
#open IN,"test.txt";
open OUT,">$out";
my %hash;
while(<IN>){
	chomp;
	my @field=split;
	if(($field[3] eq "C" && $field[4] eq "A")||($field[3] eq "G" && $field[4] eq "T")){
		$hash{$field[5]}{"CA"}++;
	}elsif(($field[3] eq "C" && $field[4] eq "G")||($field[3] eq "G" && $field[4] eq "C")){
		$hash{$field[5]}{"CG"}++;
	}elsif(($field[3] eq "C" && $field[4] eq "T")||($field[3] eq "G" && $field[4] eq "A")){
		$hash{$field[5]}{"CT"}++;
	}elsif(($field[3] eq "T" && $field[4] eq "A")||($field[3] eq "A" && $field[4] eq "T")){
		$hash{$field[5]}{"TA"}++;
	}elsif(($field[3] eq "T" && $field[4] eq "G")||($field[3] eq "A" && $field[4] eq "C")){
		$hash{$field[5]}{"TG"}++;
	}elsif(($field[3] eq "T" && $field[4] eq "C")||($field[3] eq "A" && $field[4] eq "G")){
		$hash{$field[5]}{"TC"}++;
	}else{
		$hash{$field[5]}{"indel"}++;
	}
}
close IN;
foreach my $sap (keys %hash){
	$hash{$sap}{"CA"}=0 if(!exists $hash{$sap}{"CA"});
	$hash{$sap}{"CG"}=0 if(!exists $hash{$sap}{"CG"});
	$hash{$sap}{"CT"}=0 if(!exists $hash{$sap}{"CT"});
	$hash{$sap}{"TA"}=0 if(!exists $hash{$sap}{"TA"});
	$hash{$sap}{"TC"}=0 if(!exists $hash{$sap}{"TC"});
	$hash{$sap}{"TG"}=0 if(!exists $hash{$sap}{"TG"});
	$hash{$sap}{"indel"}=0 if(!exists $hash{$sap}{"indel"});
	print OUT "$sap\t$hash{$sap}{\"CA\"}\t$hash{$sap}{\"CG\"}\t$hash{$sap}{\"CT\"}\t$hash{$sap}{\"TA\"}\t$hash{$sap}{\"TC\"}\t$hash{$sap}{\"TG\"}\t$hash{$sap}{\"indel\"}\n";
}
close OUT;
