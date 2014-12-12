#!/usr/bin/perl -w

use strict;

my $base_command = "treemix -i loci_20.treemix.gz -bootstrap -o ";

my $n_boot = shift
  or die "Specify number of bootstrap replicated on command line: $!";

for (my $i = 0; $i < $n_boot; $i++) {
  my $command = sprintf("%s repens-%04d", $base_command, $i+1);
  print $command, "\n";
  system($command) == 0
    or die "system $command failed: $?";
}
