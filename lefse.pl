use strict;
use warnings;

my $mapFile=$ARGV[0];
my $tableFile=$ARGV[1];

my %mapHash=();
open(MAP,"$mapFile")or die $!;
while(my $line=<MAP>)
{
	chomp($line);
	my @arr=split(/\s+/,$line);
	$mapHash{$arr[0]}=$arr[3];
}
close(MAP);

my %hash=();
my @samples=();
my %sumHash=();
open(RF,"$tableFile") or die $!;
while(my $line=<RF>)
{
	next if($.==1);
	my @arr=split(/\t/,$line);
	if($.==2)
	{
		@samples=@arr;
		next;
	}
	my $taxonomy=$arr[$#arr];
	$taxonomy=~s/;/|/g;
	$taxonomy=~s/(k|p|c|o|f|g|s)__//g;
	$taxonomy=~s/\s+|\[|\]//g;
	$taxonomy=~s/\|*$//g;
	my @taxoArr=split(/\|/,$taxonomy);
	my $taxSetp="";
	for(my $tax=0;$tax<=$#taxoArr;$tax++)
	{
		if($tax==0)
		{
			$taxSetp=$taxoArr[$tax];
		}
		else
		{
			$taxSetp=$taxSetp . "|" . $taxoArr[$tax];
		}
		for(my $i=1;$i<$#samples;$i++)
		{
			${$hash{$taxSetp}}{$samples[$i]}+=$arr[$i];
		}
	}
	for(my $i=1;$i<$#samples;$i++)
	{
		$sumHash{$samples[$i]}+=$arr[$i];
	}
}
close(RF);

open(WF,">lefse_input.txt")or die $!;
#print WF "Taxonomy";
#for(my $i=1;$i<$#samples;$i++)
#{
#	print WF "\t" .$samples[$i];
#}
#print WF"\n";
print WF "sampleType";
for(my $i=1;$i<$#samples;$i++)
{
	print WF "\t" . $mapHash{$samples[$i]};
}
print WF "\n";

foreach my $key(sort (keys %hash))
{
	next if($key=~/Unassigned/);
	print WF $key;
	for(my $i=1;$i<$#samples;$i++)
	{
		print WF "\t" . ${$hash{$key}}{$samples[$i]}/$sumHash{$samples[$i]}*100;
	}
	print WF "\n";
}
close(WF);
