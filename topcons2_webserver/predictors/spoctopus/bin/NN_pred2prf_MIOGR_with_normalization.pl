#!/usr/bin/perl -w

if($#ARGV != 5) {
  printf "Usage: NN_pred2prf.pl <protnamefile> <M_prfdir> <IO_prfdir> <G_prfdir> <R_prfdir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $Mprfdir = $ARGV[1]."/";
my $IOprfdir = $ARGV[2]."/";
my $Gprfdir = $ARGV[3]."/";
my $Rprfdir = $ARGV[4]."/";
my $outdir = $ARGV[5]."/";



my @preds  =();


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";





while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my @Mcols = ();
  open(PRFFILE, "$Mprfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @Mcols, [@rowlist];

    }
    else {
      
    }
  }
  close(PRFFILE);

  
   
  
  my @IOcols = ();
  open(PRFFILE, "$IOprfdir"."$seqname".".prf")
    or die "could not open prffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @IOcols, [@rowlist];

    }
    else {
      
    }
  }
  close(PRFFILE);

 
  
  my @Gcols = ();
  open(PRFFILE, "$Gprfdir"."$seqname".".prf")
    or die "could not open prffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @Gcols, [@rowlist];

    }
    else {
      
    }
  }
  close(PRFFILE);


  my @Rcols = ();
  my $intro = "";
  my $extro = "";
  my $at_intro = 1;
  open(PRFFILE, "$Rprfdir"."$seqname".".prf")
    or die "could not open prffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @Rcols, [@rowlist];

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
    or die "could not open outfile for $seqname\n";
  print OUTFILE "$intro";
  print OUTFILE "ALPHABET:       M       IO      G       R       -     <SPACE> <LABEL> <QUERY> \n";

  for(my $i = 0; $i <= $#Mcols; $i++) {
    my $col_nr = $i + 1;
    my @Mcol = @{$Mcols[$i]};
    my @IOcol = @{$IOcols[$i]};
    my @Gcol = @{$Gcols[$i]};
    my @Rcol = @{$Rcols[$i]};
    if(!$Mcol[2]) {
      die "@Mcol\n";
    }
    

    $Mcol[2] += 0.1;
    $Mcol[2] /= 1.1;
    $IOcol[2] += 0.1;
    $IOcol[2] /= 1.1;
    $Gcol[2] += 0.1;
    $Gcol[2] /= 1.1;
    $Rcol[2] += 0.1;
    $Rcol[2] /= 1.1;

    #$Mcol[2] = exp($Mcol[2]) / exp(1);
    #$IOcol[2] = exp($IOcol[2]) / exp(1);
    #$Gcol[2] = exp($Gcol[2]) / exp(1);
    #$Rcol[2] = exp($Rcol[2]) / exp(1);

    printf OUTFILE "COL %4d".":\t", $col_nr;
    printf OUTFILE "%.3f\t%.3f\t%.3f\t%.3f\t", $Mcol[2], $IOcol[2], $Gcol[2], $Rcol[2];
    printf OUTFILE "0.00   0.00   $IOcol[5]       $Mcol[6]       \n";
  }




  print OUTFILE "$extro";
  #print "done running $seqname\n";

  close OUTFILE;
  
}
