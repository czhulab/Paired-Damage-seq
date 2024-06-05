#!/usr/bin/perl
use strict;
use warnings;

open IN, "sample_list.txt" or die $!;
<IN>;
while(<IN>){
  chomp;
  my @tmp = split/\s+/, $_;
  my $dna = $tmp[0];
  my $rna = $tmp[1];
  system("sh per_run.sh $dna $rna"); ### or change to your job submission script

}
close IN;
