#!/usr/bin/perl

# Parse genotype in gwas.gen format and adapt to sQTLseekeR input requirements

use strict;
use warnings;
use Switch;

# Iterate over each line and parse it
while(defined (my $line = <>)){

  chomp $line;
  $line =~ m/^(.*) [ACTGDI] [ACTGDI] (.*)$/ or die "There is at least one line that has not the standard structure\n";
  
  my $head = $1;
  my $tail = parser($2);
  
  my $output = "$head $tail\n"; 
  
  $output =~ s/ /\t/g;
  
  print $output;

}


sub parser{

  my($genotypes) = @_;
  $genotypes =~ s/\s//g;
  my $l= length $genotypes;
  $l == 114 or die "The number of genotypes differs from the number of samples\n";
  my $gen;
  my $newgen;
  my $newgenotypes;

  for (my $i=0; $i<=$l-3; $i+=3){
    $gen = substr($genotypes,$i,3);
    
    switch ($gen) {
      case "001" {$newgen="2"}
      case "010" {$newgen="1"}
      case "100" {$newgen="0"}
      else {$newgen= "9"}
    }
    $newgenotypes.="$newgen ";
  }
  $newgenotypes =~ s/\s+$//;
  return $newgenotypes; 
}
