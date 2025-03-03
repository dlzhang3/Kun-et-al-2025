#!/usr/bin/perl
#use warnings;
die "perl $0 <uniq.reads> <bowtie out>\n" if @ARGV<2;
my %nreads;
open unr,$ARGV[0] || die "file could not be opened";
while(<unr>)
{
	chomp;
	my $uniq = (split)[0];
	my $num = (split)[1];
	$nreads{$uniq}=$num;
}
close unr;

my %nloc;
open in,$ARGV[1];
while(<in>)
{
	chomp;
	my ($seq)=(split)[0];
	$nloc{$seq}++;
}
close in;

open in,$ARGV[1];
while(<in>)
{
    chomp;
    my ($seq,$strand,$chr,$pos,$len)=(split)[0..4];
	print $chr,"\t",$pos,"\t",($pos+length($len)),"\t",$seq,"\t",$nreads{$seq}/$nloc{$seq},"\t",$strand,"\n";
}
close in;

