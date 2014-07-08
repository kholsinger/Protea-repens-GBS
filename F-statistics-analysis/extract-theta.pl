#!/usr/bin/perl -w

use strict;

while (<>) {
  next unless m/theta/;
  if (m/locus/) {
    my @fields = split;
    print $fields[1], "\n";
  }
}
