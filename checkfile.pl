#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;

#simple script to check if one query text file has any sentences not present in a reference text file
#case and punctuation agnostic.

# made to help Kalpesh Upadhye on May 12, 2020

#Usage : checkfile.pl queryfile reffile
# output is a file called reffile-queryfile which has sentences from queryfile not found in reffile

my $query=read_file($ARGV[0]);
my $ref=read_file($ARGV[1]);

my %ref;

foreach my $sen(split(/\./,$ref)){
	$sen=lc($sen);
	$sen=~ s/[[:punct:]]//g;
	($sen)=($sen=~/^\s*(\S.*\S)\s*$/);
	$ref{$sen}=0;
	
}

open(OUT,">$ARGV[1]-$ARGV[0]");

foreach my $sen(split(/\./,$query)){
	$sen=lc($sen);
	$sen=~ s/[[:punct:]]//g;
	($sen)=($sen=~/^\s*(\S.*\S)\s*$/);
	if(!exists $ref{$sen}){print OUT "$sen.\n";}
}

close OUT;

