#!/usr/bin/perl -w


if($#ARGV != 6) {
  printf "Usage: generate_matlab_score_m_file.pl <outfile> <netdir> <inputfiledir><inputfile> <netfile> <predoutfile> <netlabdir>\n";
  exit;
}

my $outfile = $ARGV[0];
my $netdir = $ARGV[1]."/";
my $inputfiledir = $ARGV[2]."/";
my $inputfile = $ARGV[3];
my $netfile = $ARGV[4];
my $predoutfile = $ARGV[5];
my $netlabdir = $ARGV[6]."/";

my $HIDDENLAYER = "8";
my $OUTPUTLAYER = "1";




open(OUTFILE, ">"."$outfile")
    or die "could not open $outfile";

print OUTFILE "cd $netdir\n";
print OUTFILE "load -mat $netfile\n";
print OUTFILE "cd $inputfiledir\n";
print OUTFILE "load $inputfile\n";
print OUTFILE "cd $netlabdir\n";
$inputfile =~ s/\.dat//g;
print OUTFILE "pred = mlpfwd(net,$inputfile);\n";
print OUTFILE "save $predoutfile pred -ascii;\n";


close OUTFILE;
