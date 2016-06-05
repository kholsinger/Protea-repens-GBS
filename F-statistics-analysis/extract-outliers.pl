#!/usr/bin/perl -w

use strict;

my %f_statistics;
while (<>) {
  next unless /TP/;
  chomp;
  my ($locus, $fst, $kld) = split;
  $locus =~ s/^(TP.*)\[.*\]:.*$/$1/;
  print $locus, "\n";
}
