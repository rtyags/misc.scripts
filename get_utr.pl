#!/gsc/bin/perl

use warnings;
use strict;

#open(FAS,$ARGV[1]);
#open(FAS,"temp" );
#my %seqhash;
#my($head,$seq,$len)=("","",-1);
#my $flag=0;
#while(<FAS>){
#	chomp;
#	if(/^>/){
#		if($flag!=0){
#			$seqhash{$head}=$seq;
#			$seq="";
#		}
#		$flag++; 
#		#($head)=(/^>(.*)/);
#		($head)=(/^>(\S+)/);
#	}else{$seq.=$_;}
#}
#$seqhash{$head}=$seq;
#
#close FAS;
#
open(IN,$ARGV[0]);
my @gff=<IN>;
close IN;

chomp(@gff);

my($p1,$pdir,$p5,$p9)=("","+",0,"");
my $flag=0;
for(my $i=0;$i<@gff;$i++){
	
	my ($end,$start)=(-1,-1);
	my @fld=split(/\t/,$gff[$i]);
	if($fld[3]>$p5 || $p1 ne $fld[0]){
		if($pdir eq "-" && $flag==0){
			if($p1 eq $fld[0] && $fld[3]<($p5+2000)){
				$end=$fld[3]-1;
			}else{
				$end=$p5+2000;
			}
			print "$p9\t$p1\t".($p5+1)."\t$end\t$pdir\n";
		}
		if($fld[6] eq "+"){
			if($flag==0 || $p1 ne $fld[0]){
				if($p1 eq $fld[0] && $fld[3]<($p5+2000)){
					$start=$p5+1;
				}elsif($fld[3]>2000){
					$start=$fld[3]-2000;
				}else{$start=1;}
				print "$fld[8]\t$fld[0]\t$start\t".($fld[3]-1)."\t$fld[6]\n";
			}elsif($i>1){
				my @fld1=split(/\t/,$gff[$i-2]);
				if($fld[3]>$fld1[4]){
					if($fld[3]<($fld1[4]+2000)){
						$start=$fld1[4]+1;
					}elsif($fld[3]>2000){
						$start=$fld[3]-2000;
					}else{$start=1;}
					print "$fld[8]\t$fld[0]\t$start\t".($fld[3]-1)."\t$fld[6]\n";
				}elsif($fld[4]>$fld1[4]){
					print "$fld[8]\t$fld[0]\tNA\tNA\t$fld[6]\n";
				}else{print "ERROR! 2 intronic genes under one gene! current=$fld[8]:previous=$p9:end_for_i-2:$fld1[4]\n";}
			}else{print "ERROR! index less than 0!\n";}
		}
		$flag=0;
	}else{
		$flag=1;
		if($pdir eq "-"){
			if($i<(@gff-1)){
				my @fld1=split(/\t/,$gff[$i+1]);
				if($p5 eq $fld1[0]){
					if($fld1[3]<$p5){
						if($fld1[4]>$p5){
							print "$p9\t$p1\tNA\tNA\t$pdir\n";
						}else{ print "ERROR! 2 intronic genes under one gene! current=$fld[8]:next=$fld1[8]:previous_end=$p5\n";}
					}else{
						if($p1 eq $fld1[0] && $fld1[3]<($p5+2000)){
							$end=$fld1[3]-1;
						}else{
							$end=$p5+2000;
						}
						print "$p9\t$p1\t".($p5+1)."\t$end\t$pdir\n";
					}
				}else{
					print "$p9\t$p1\t".($p5+1)."\t".($p5+2000)."\t$pdir\n";
				}
			}else{
				print "$p9\t$p1\t".($p5+1)."\t".($p5+2000)."\t$pdir\n";
			}
		}
		my $res= `grep -P \"\\tCDS\\t.*$p9;\" $ARGV[1]|sort -k 4,4g|cut -f 4,5`;
		my @cds=split(/\n/,$res);
		if($fld[6] eq "+"){
			my $ovlp=0;
			foreach my $cds(@cds){
				my ($st,$en)=split(/\t/,$cds);
				if($en<$fld[3]-1){
					$start=int($en)+1;
				}elsif($st<$fld[3]){$ovlp=1;}
			}
			if($ovlp==0){
				print "$fld[8]\t$fld[0]\t".(($start<($fld[3]-2000))?($fld[3]-2000):$start)."\t".($fld[3]-1)."\t$fld[6]\n";
			}else{
				print "$fld[8]\t$fld[0]\tNA\tNA\t$fld[6]\n";
			}
		}else{
			my $ovlp=0;
			foreach my $cds(@cds){
				my ($st,$en)=split(/\t/,$cds);
				if($st>$fld[4]+1){
					$end=$st-1;
					last;
				}elsif($en>$fld[3]){$ovlp=1;}
			}
			if($ovlp==0){
				print "$fld[8]\t$fld[0]\t".($fld[4]+1)."\t".(($end>($fld[4]+2000))?($fld[4]+2000):$end)."\t$fld[6]\n";
			}else{
				print "$fld[8]\t$fld[0]\tNA\tNA\t$fld[6]\n";
			}
		}
	}


	$p5=$fld[4];$p1=$fld[0];$p9=$fld[8];$pdir=$fld[6];
}
if($pdir eq "-" && $flag==0){
	print "$p9\t$p1\t".($p5+1)."\t".($p5+2000)."\t$pdir\n";
}

#	my $dir=$[4];
#	my $seq="";
#	if(exists $seqhash{$b[1]}){
#		$seq=substr($seqhash{$b[1]},$b[2]-1,$b[3]-$b[2]+1);
#	}else{print "ERROR!! No fasta record for $b[1]\n";}
#
#	if($dir eq "-"){
#		$seq=~tr/actgnACTG/TGACNTGAC/;
#		$seq=reverse $seq;
#	}
#	print ">$b[0]\n$seq\n";
exit;
