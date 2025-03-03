#!/usr/bin/perl

my %fa;
while(<>)
{
	chomp;
	my $a=split/\t/;
	my @fa=(split/\t/)[0..($a-1)/2-1];
	my @fb=(split/\t/)[($a-1)/2..$a-2];
	my $bp=(split/\t/)[$a-1];
	
	push @{$fa{join("\t",@fa)}{$bp}},[@fb];
}

foreach my $loc (keys %fa)
{
	my $best=shift @{[reverse sort {$a<=>$b} keys %{$fa{$loc}}]};
	if(@{$fa{$loc}{$best}}>1)
	{
		my @fa_info=split/\t/,$loc;
		$fa_info[4]=$fa_info[4]/@{$fa{$loc}{$best}};
		foreach my $l (@{$fa{$loc}{$best}})
		{
			print join("\t",@fa_info),"\t",join("\t",@{$l}),"\n";
		}	
	}
	else
	{
		print $loc,"\t",join("\t",@{$fa{$loc}{$best}[0]}),"\n";
	}
}
