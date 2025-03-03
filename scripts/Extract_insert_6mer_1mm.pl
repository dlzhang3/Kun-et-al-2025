#!/usr/bin/perl

die "usage: $0 <raw data> <3linker> > output\n" if @ARGV<1;

my $linkerpreflength=6; #length of prefix of linker that's checked
my $linkerstr=$ARGV[1]; 
my $linker=substr($linkerstr,0,$linkerpreflength);
my $mlinker=substr($linkerstr,0,$linkerpreflength+1);
my @linker;
for my $i (0..length($mlinker)-2)
{
	my $s=$mlinker;
	substr($s,$i,1,".");
	push @linker,$s;
}
print "Insert\tInsertLength\tLinkerscore\tLinker\n";
open SEQ,$ARGV[0];
while(<SEQ>)
{
	chomp;
	my ($seq)=(split/\t/)[0];
	#perfect match
	if(my ($rna)=$seq=~/(.*)$linker/)
	{
		print $rna,"\t",length($rna),"\t1\t",substr($seq,length($rna)),"\n" if $rna ne "";
	}
	else
	{
		foreach my $l (@linker)
		{
			if(my ($rna)=$seq=~/(.*)$l/)
			{
				print $rna,"\t",length($rna),"\t0\t",substr($seq,length($rna)),"\n" if $rna ne "";
				last;
			}	
		}
	}
}		

