#!/usr/bin/perl -w

use strict;
use warnings;

my $truepos = 0;
my $falsepos = 0;
my $matched = 0;

my ($gold, $fn) = @ARGV;

open I, "$fn";

while (<I>)
{
    chomp;
    my @col = split /\t/;
    
    $col[0] =~ /(^|;)(RF\d+);/;
    my $x = $2;
    
    $col[1] =~ /(^|;)(RF\d+);/;
    my $y = $2;
    
    if (($col[2] >= 75.0) && ($col[3] >= 90))
    {
	$matched++;
    }
    
    if ($x eq $y)
    {
	$truepos++;
    }
    else
    {
	$falsepos++;
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
my $matchedrate = 100.0 * $matched / $gold;

printf "Total queries:             %7d\n", $gold;
printf "True positives:            %7d\n", $truepos;
printf "False positives:           %7d\n", $falsepos;
printf "False negatives:           %7d\n", $falseneg;
printf "\n";
printf "Recall:                    %10.2f%%\n", $recall;
printf "Precision:                 %10.2f%%\n", $prec;
printf "False negative rate:       %10.2f%%\n", $fnr;
printf "False discovery rate:      %10.2f%%\n", $fdr;
printf "F-score:                   %10.2f%%\n", $fscore;
printf "\n";
printf "Matched (id>0.7;cov>=0.9): %7d\n", $matched;
printf "Matched percentage:        %10.2f%%\n", $matchedrate;
