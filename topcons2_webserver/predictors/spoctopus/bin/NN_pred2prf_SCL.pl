#!/usr/bin/perl -w

if($#ARGV != 3) {
  printf "Usage: NN_pred2prf.pl <protnamefile> <S_prfdir> <C_prfdir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $Sprfdir = $ARGV[1]."/";
my $Cprfdir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";



my @preds  =();


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";





while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my @Scols = ();
  open(PRFFILE, "$Sprfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @Scols, [@rowlist];

    }
    else {
      
    }
  }
  close(PRFFILE);


  my @Ccols = ();
  my $intro = "";
  my $extro = "";
  my $at_intro = 1;
  open(PRFFILE, "$Cprfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @Ccols, [@rowlist];

    }
    elsif($at_intro > 0)  {
      if($_ =~ 'ALPHABET') {
	$at_intro = 0;
      }
      else {
	$intro .= "$_"."\n";
      }
    }
    else {
      $extro .= "$_"."\n";
    }
  }
  close(PRFFILE);
  
  

  
  open(OUTFILE, ">"."$outdir"."$seqname".".prf")
      or die "could not open outfile for $seqname";
  print OUTFILE "$intro";
  print OUTFILE "ALPHABET:       S       C      L       -     <SPACE> <LABEL> <QUERY> \n";

  for(my $i = 0; $i <= $#Scols; $i++) {
    my $col_nr = $i + 1;
    my @Scol = @{$Scols[$i]};
    my @Ccol = @{$Ccols[$i]};
    if(!$Scol[2]) {
	die "@Scol";
    }
    if($Scol[2] == 0) {
      $Scol[2] = 0.005;
    }
    if($Ccol[2] == 0) {
      $Ccol[2] = 0.005;
    }

    if($i >= 70) {
	$Scol[2] = 0.005;
    }

    printf OUTFILE "COL %4d".":\t", $col_nr;
    printf OUTFILE "%.3f\t%.3f\t%.3f\t", $Scol[2], $Ccol[2], 1.0 - $Scol[2];
    printf OUTFILE "0.00   0.00   .       $Scol[6]       \n";
  }




  print OUTFILE "$extro";
  #print "done running $seqname\n";

  close OUTFILE;
  
}
