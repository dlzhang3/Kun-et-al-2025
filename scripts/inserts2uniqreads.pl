#!/usr/bin/perl
die "perl $0 <inserts> <min length> <max lenght>\n" if @ARGV<3;
open IN, $ARGV[0];
while(<IN>){
 chomp;
 my ($seq)=(split/\t/)[0];
 next if($seq!~/\S/);
 next if($seq=~/\./);
 if(length($seq)>=$ARGV[1] && length($seq)<=$ARGV[2]){
  $hash{$seq}+=1;
 }
}

for my $s (keys %hash){
 print "$s\t$hash{$s}\n";
}
