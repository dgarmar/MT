#!/usr/bin/perl

## Retrieve the position of exon and UTR boundaries for each transcript ID from an annotation file.
## -a : the annotation file to use
## -o : the output file

use warnings;
use strict;
use Getopt::Long;

my $annoF;
my $outfile;

GetOptions('Annotation|a=s'=>\$annoF,'Outfile|o=s'=>\$outfile);

if($annoF =~ /\.gz/){
    open(ANNO,"gunzip -c $annoF |") or die "Can't annotation file: $!";
} else {
    open(ANNO,"$annoF") or die "Can't annotation file: $!";
}
open(OUTF,">$outfile")  or die "Can't output file: $!";
## Print output headers
print OUTF "txId\tstrand\ttxStart\ttxEnd\tcdsCount\tcdsStarts\tcdsEnds\tutrCount\tutrStarts\tutrEnds\n";

my $tx_id_cour = "";
my @cdsStarts = ();
my @cdsEnds = ();
my $cdsCount = 0;
my $txStart = -1;
my $txEnd = -1;
my @utrStarts = ();
my @utrEnds = ();
my $utrCount = 0;
my $strand = "NA";
while(defined(my $line = <ANNO>)){
    chomp $line;
    my @temp = split(/\t/,$line);
    if(defined($temp[2]) && $temp[2] eq "CDS"){
	$line =~ /transcript_id "([^\"]+)";/;
	my $tx_id = $1;
	if($tx_id_cour eq $tx_id || $tx_id_cour eq ""){
	    $tx_id_cour = $tx_id;
	    ## Update values
	    $strand = $temp[6];
	    $cdsCount++;
	    if($txStart > $temp[3] || $txStart == -1){
		$txStart = $temp[3];
	    }
	    if($txEnd < $temp[4] || $txEnd == -1){
		$txEnd = $temp[4];
	    }
	    push(@cdsStarts,$temp[3]);
	    push(@cdsEnds,$temp[4]);
	} else {
	    ## Print previous transcript
	    $" = ",";
	    print OUTF "$tx_id_cour\t$strand\t$txStart\t$txEnd\t$cdsCount\t@cdsStarts\t@cdsEnds\t$utrCount\t@utrStarts\t@utrEnds\n";
	    ## Initialisation
	    $cdsCount = 0;
	    @cdsStarts = ();
	    @cdsEnds = ();
	    $utrCount = 0;
	    @utrStarts = ();
	    @utrEnds = ();
	    $txStart = -1;
	    $txEnd = -1;
	    $tx_id_cour = $tx_id;	
	    ## Update values
	    $strand = $temp[6];
	    $cdsCount++;
	    $txStart = $temp[3];
	    $txEnd = $temp[4];
	    push(@cdsStarts,$temp[3]);
	    push(@cdsEnds,$temp[4]);
	}
    }
    if(defined($temp[2]) && $temp[2] eq "UTR"){
	$line =~ /transcript_id "(\w*\.?\w*)";/;
	my $tx_id = $1;
	if($tx_id_cour eq $tx_id || $tx_id_cour eq ""){
	    $tx_id_cour = $tx_id;
	    ## Update values
	    $utrCount++;
	    if($txStart > $temp[3] || $txStart == -1){
		$txStart = $temp[3];
	    }
	    if($txEnd < $temp[4] || $txEnd == -1){
		$txEnd = $temp[4];
	    }
	    push(@utrStarts,$temp[3]);
	    push(@utrEnds,$temp[4]);
	} else {
	    ## Print previous transcript
	    $" = ",";
	    print OUTF "$tx_id_cour\t$strand\t$txStart\t$txEnd\t$cdsCount\t@cdsStarts\t@cdsEnds\t$utrCount\t@utrStarts\t@utrEnds\n";
	    ## Initialisation
	    $cdsCount = 0;
	    @cdsStarts = ();
	    @cdsEnds = ();
	    $utrCount = 0;
	    @utrStarts = ();
	    @utrEnds = ();
	    $txStart = -1;
	    $txEnd = -1;
	    $tx_id_cour = $tx_id;	
	    ## Update values
	    $strand = $temp[6];
	    $utrCount++;
	    $txStart = $temp[3];
	    $txEnd = $temp[4];
	    push(@utrStarts,$temp[3]);
	    push(@utrEnds,$temp[4]);
	}
    }
}
close(OUTF);
close(ANNO);
