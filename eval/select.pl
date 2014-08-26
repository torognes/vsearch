#!/usr/bin/perl -w

# Print the first representative of each Rfam familiy member
# with more than one member to one file.
# Print the others to another file.

use strict;
use warnings;

my $rfam;
my %count = ();

die "Not enough arguments" if (scalar @ARGV < 3);

my ($input, $first, $second) = @ARGV;

# first pass, just count

open I, "$input";

while(<I>)
{
    if (/[>;](RF\d+);/)
    {
	$rfam = $1;
	if (defined $count{$rfam})
	{
	    $count{$rfam}++;
	}
	else
	{
	    $count{$rfam} = 1;
	}
    }
}

close I;

# second pass, print

open I, "$input";
open F, ">$first";
open S, ">$second";

my $select = 0;

while(<I>)
{
    if (/[>;](RF\d+);/)
    {
	$rfam = $1;
	if ($count{$rfam} > 1)
	{
	    $select = 1;
	    $count{$rfam} = 0;
	}
	else
	{
	    $select = 0;
	}
    }

    if ($select)
    {
	print F;
    }
    else
    {
	print S;
    }
}

close S;
close F;
close I;
