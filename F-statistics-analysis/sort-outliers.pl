#!/usr/bin/perl -w

use strict;

my %f_statistics;
while (<>) {
  next unless /TP/;
  chomp;
  my ($locus, $fst, $kld) = split;
  $f_statistics{$locus}=$kld;
}
foreach my $locus (sort numerically keys %f_statistics) {
  printf("%10s: %5.3f\n", $locus, $f_statistics{$locus});
}

sub numerically {
  $f_statistics{$b} <=> $f_statistics{$a}
}
