#!/usr/bin/perl -w

use strict;
use warnings;

my $truepos = 0;
my $falsepos = 0;

my ($gold, $fn) = @ARGV;

open I, "$fn";

while (<I>)
{
    if (/^H/)
    {
	chomp;
	my @col = split /\t/;
	
	$col[8] =~ /^[^;]*;([^;]*);[^;]*/;
	my $x = $1;

	$col[9] =~ /^[^;]*;([^;]*);[^;]*/;
	my $y = $1;
	
	if ($x eq $y)
	{
	    $truepos++;
	}
	else
	{
	    $falsepos++;
	}
    }
}

close I;

my $pos = $truepos + $falsepos;
my $falseneg = $gold - $ truepos;
my $prec = 100.0 * $truepos / $pos;
my $fdr  = 100.0 * $falsepos / $pos;
my $recall = 100.0 * $truepos / $gold;
my $fnr = 100.0 * $falseneg / $gold;
my $fscore = 2.0 * ($prec * $recall) / ($prec + $recall);

printf "Total queries:        %7d\n", $gold;
printf "True positives:       %7d\n", $truepos;
printf "False positives:      %7d\n", $falsepos;
printf "False negatives:      %7d\n", $falseneg;
printf "Recall:               %10.2f%%\n", $recall;
printf "Precision:            %10.2f%%\n", $prec;
printf "False negative rate:  %10.2f%%\n", $fnr;
printf "False discovery rate: %10.2f%%\n", $fdr;
printf "F-score:              %10.2f%%\n", $fscore;
