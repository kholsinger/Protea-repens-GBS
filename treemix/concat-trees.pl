#!/usr/bin/perl -w

use strict;

use File::Find;
use File::Temp qw/tempfile/;
use IO::Uncompress::Gunzip qw/gunzip $GunzipError/;

my $concat_file = shift
  or die "Enter output file name on commandline";

find(\&wanted, ".");

sub wanted {
  if (m/treeout/) {
    print "$_\n";
    my ($fh, $filename) = tempfile();
    gunzip $_ => $fh
      or die "gunzip $_ failed: $GunzipError\n";
    close $fh;
    open(INFILE, "< $filename")
      or die "Could not open $filename for input: $!";
    open(OUTFILE, ">> $concat_file")
      or die "Could not open $concat_file for output: $!";
    while (<INFILE>) {
      print OUTFILE;
    }
    close(OUTFILE);
    close(INFILE);
  }
}
