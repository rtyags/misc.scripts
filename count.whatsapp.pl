#!/usr/bin/perl
use warnings;
use strict;

open(IN,$ARGV[0]);

my %post;
my %word;

my $name="";
my $post="";

while(<IN>){
	next if ((/:\d\d [AP]M - / && !/:.*:/) || !/\S/ );
	chomp;
	my @text_words;
	if(!/:\d\d [AP]M - /){
		@text_words = split(/\s+/);
	}else{
		($name,$post)=(/:\d\d [AP]M - ([^:]+):\s*(.*)$/);
		next if $name=~/(joined using this group)|added|left|changed|created|You|group/;
		$post{$name}++;
		@text_words = split(/\s+/, $post);
	}
	$word{$name}+=scalar(@text_words);
}

open(OUT1,">post.sort.tsv");
open(OUT2,">words.sort.tsv");

foreach my $key ( sort { $post{$b} <=> $post{$a} } ( keys(%post))) {
      printf OUT1 ("%4d |\t%s\n", $post{$key}, $key);
}

foreach my $key ( sort { $word{$b} <=> $word{$a} } ( keys(%word))) {
      printf OUT2 ("%4d |\t%s\n", $word{$key}, $key);
}
close IN;
close OUT1;
close OUT2;

exit;
