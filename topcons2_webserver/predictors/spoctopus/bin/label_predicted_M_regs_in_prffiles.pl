#!/usr/bin/perl -w

if($#ARGV != 4) {
  printf "Usage: label_predicted_M_regs_in_prffiles.pl <protnamefile> <prfdir> <topodir> <outdir> <prf alphabet size>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $topodir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";
my $a_size = $ARGV[4];



open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;


  my $outfile = "$outdir"."$seqname".".prf";

  open(OUTFILE, ">"."$outfile")
      or die "could not open $outfile";

  my $topofile = "$topodir"."$seqname".".top";
  $topo = "";
  open(TOPOFILE,"$topofile") or die "Could not open $topofile";
  while(<TOPOFILE>) {
    chomp;
    if($_ =~ '>') {
    }
    else {
      $topo .= $_;
    }
  }
  $topo =~ s/[^M]/./g;
  @topo = split //, $topo;
  
  
  my $topopos = 0;
  open(PRFFILE, "$prfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      my $row = $_;
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      if($topopos > $#topo) {
	  die "$seqname";
      }
      if($rowlist[$a_size+4] ne $topo[$topopos]) {
	$rowlist[$a_size+4] = $topo[$topopos];
	#print "@rowlist\n";
	$topopos++;
	printf OUTFILE "COL %4d", $topopos;
	print OUTFILE ":  ";
	for $b (2 .. $#rowlist ) {
	  if($b >= ($#rowlist - 1)) {
	    print OUTFILE "$rowlist[$b]       ";
	  }
	  else {
	    printf OUTFILE "%5.03f   ", $rowlist[$b];
	  }
	}
	print OUTFILE "\n";
      }
      else {
	$topopos++;
	print OUTFILE "$row\n";
      }
      
    }
    else {
      print OUTFILE "$_\n";
    }
  }
  close(PRFFILE);
  close OUTFILE;

}
