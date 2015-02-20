#!/usr/bin/perl -w


if($#ARGV != 7) {
  printf "Usage: generate_matlab_score_m_file.pl <outfile> <basedir> <inputfile> <netfile> <inputlayer> <predoutfile> <letter> <scorefiledir>\n";
  exit;
}

my $outfile = $ARGV[0];
my $basedir = $ARGV[1]."/";
my $inputfile = $ARGV[2];
my $netfile = $ARGV[3];
my $inputlayer = $ARGV[4];
my $predoutfile = $ARGV[5];
my $letter = $ARGV[6];
my $scorefiledir = $ARGV[7]."/";

my $HIDDENLAYER = "8";
my $OUTPUTLAYER = "2";




open(OUTFILE, ">"."$outfile")
    or die "could not open $outfile";

print OUTFILE "cd $basedir/OPM_101_TRAINED_NETS/$letter\n";
print OUTFILE "load -mat $netfile\n";
print OUTFILE "cd $scorefiledir\n";
print OUTFILE "load $inputfile\n";
print OUTFILE "cd $basedir/netlab/\n";
$inputfile =~ s/\.dat//g;
print OUTFILE "pred = mlpfwd(net,$inputfile);\n";
print OUTFILE "save $predoutfile pred -ascii;\n";


close OUTFILE;
