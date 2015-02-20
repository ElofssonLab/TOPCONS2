#!/usr/bin/perl -w

if($#ARGV != 1) {
  printf "Usage: splitseqfile.pl <seqnamefile> <outdir>\n";
  exit;
}

my $fastafile = $ARGV[0];
my $outdir = $ARGV[1]."/";

my $protname = "";
my $seq = "";
my $first = 'YES';
my $outfile = "";
open(FASTAFILE, "$fastafile")
  or die "could not open $fastafile\n";
while(<FASTAFILE>) {
  chomp $_;
  my $protname = $_;
  $outfile = "$outdir"."$protname".".seq";
  open(OUTFILE, ">"."$outfile")
      or die "could not open $outfile\n";
  print OUTFILE "$protname\n";
  close OUTFILE;
}


