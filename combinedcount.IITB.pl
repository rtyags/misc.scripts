#!/usr/bin/perl
use warnings;
use strict;


if(@ARGV != 4){
print "\nUsage: combinedcount.IITB.pl transcript_file old_post_stats_file old_wordcount_stats_file directory_file\n\n";
exit;
}


open(IN,$ARGV[0]);

my %post;
my %word;

my $name="";
my $post="";

my $last="";
while(<IN>){
	next if ((/:\d\d [AP]M - / && !/:.*:/) || !/\S/ );
	chomp;
	my @text_words;
	if(!/:\d\d [AP]M - /){
		@text_words = split(/\s+/);
	}else{
		($name,$post)=(/:\d\d [AP]M - ([^:]+):\s*(.*)$/);
		next if $name=~/(joined using this group)|added|left|changed|created|You|group|(disappearing messages)/;
		if($name ne $last){$post{$name}++;}
		@text_words = split(/\s+/, $post);
		$last=$name;
	}
	$word{$name}+=scalar(@text_words);
}

my (%oldp,%oldpr);
my (%oldw,%oldwr);

my %dir;
open(DIR,$ARGV[3]);
while(<DIR>){
	my($key,$val)=(/^(.*)\t(.*)\s*$/);
	$dir{$key}=$val;
}
close DIR;

open(OLD1,$ARGV[1]);
my $i=0;
while(<OLD1>){
	$i++;
	($name,$post)=(/^\s*([^|]+)\|\s*(\d+)\s*\|\s*\S+\s*\|\s*\S+\s*$/);
	$name=~s/\s+$//;
	my $key=$name;
	if(exists $dir{$name}){$key=$dir{$name};}
	$oldp{$key}=$post;
	$oldpr{$key}=$i;
}
close OLD1;

open(OLD2,$ARGV[2]);
$i=0;
while(<OLD2>){
	$i++;
	($name,$post)=(/^\s*([^|]+)\|\s*(\d+)\s*\|\s*\S+\s*\|\s*\S+\s*$/);
	$name=~s/\s+$//;
	my $key=$name;
	if(exists $dir{$name}){$key=$dir{$name};}
	$oldw{$key}=$post;
	$oldwr{$key}=$i;
}
close OLD2;


open(OUT1,">combinedpost.sort.tsv");
open(OUT2,">words.sort.tsv");

$i=0;
foreach my $key ( sort { $post{$b} <=> $post{$a} } ( keys(%post))) {
	$i++;
	$name=$key;
	if(exists $dir{$key}){$name=$dir{$key};}
	printf OUT1 ("%s | %4d | ",$name, $post{$key});
	print OUT1 (exists($oldpr{$name})?(($oldpr{$name}>$i)?("+".($oldpr{$name}-$i)):(($oldpr{$name}==$i)?"=":($oldpr{$name}-$i))):"NA")." | ".(exists($oldp{$name})?(($post{$key}>$oldp{$name})?("+".($post{$key}-$oldp{$name})):($post{$key}-$oldp{$name})):"NA")."\n";
}


$i=0;
foreach my $key ( sort { $word{$b} <=> $word{$a} } ( keys(%word))) {
	$i++;
	$name=$key;
	if(exists $dir{$key}){$name=$dir{$key};}
	printf OUT2 ("%s | %4d | ",$name, $word{$key});
	print OUT2 (exists($oldwr{$name})?(($oldwr{$name}>$i)?("+".($oldwr{$name}-$i)):(($oldwr{$name}==$i)?"=":($oldwr{$name}-$i))):"NA")." | ".(exists($oldw{$name})?(($word{$key}>$oldw{$name})?("+".($word{$key}-$oldw{$name})):($word{$key}-$oldw{$name})):"NA")."\n";
}
close IN;
close OUT1;
close OUT2;

exit;
