#!/usr/bin/perl
use warnings;
use strict;

#finds best reciprocal hits, assuming blast tabular output with proper E_value threshold is available. 

if(@ARGV!=2){

print "Usage: bbh.pl blastout1 blastout2\n\n";
exit;
}

my (%fwd,%fscore);
my (%rev,%rscore);

open(IN1, $ARGV[0]);

while(<IN1>){
	chomp;
	my @f=split(/\t/);
	if(exists $fscore{$f[0]}){
		if($fscore{$f[0]} < $f[11]){
			$fwd{$f[0]}=$f[1]; $fscore{$f[0]}=$f[11];
		}elsif($fscore{$f[0]} == $f[11]){$fwd{$f[0]}.=":$f[1]";}
	}else{$fwd{$f[0]}=$f[1]; $fscore{$f[0]}=$f[11];}
}
close IN1;

open(IN2, $ARGV[1]);

while(<IN2>){
	chomp;
	my @f=split(/\t/);
	if(exists $rscore{$f[0]}){
		if($rscore{$f[0]} < $f[11]){
			$rev{$f[0]}=$f[1]; $rscore{$f[0]}=$f[11];
		}elsif($rscore{$f[0]} == $f[11]){$rev{$f[0]}.=":$f[1]";}
	}else{$rev{$f[0]}=$f[1]; $rscore{$f[0]}=$f[11];}
}
close IN2;

foreach(keys %fwd){
	my @bbh;
	foreach my $r(split(/:/,$fwd{$_})){
		my $flag=0;
		if(exists $rev{$r}){
			foreach my $f(split(/:/,$rev{$r})){
				if($f eq $_){ $flag=1;}
			}
		}
		if($flag==1){push(@bbh,$r);}
	}
	if((scalar @bbh) >= 1){print "$_\t".join(";",@bbh)."\n";}
}

exit;

