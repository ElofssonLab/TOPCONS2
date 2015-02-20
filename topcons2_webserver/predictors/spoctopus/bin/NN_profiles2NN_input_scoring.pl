#!/usr/bin/perl -w

if($#ARGV != 6) {
  printf "Usage: NN_profiles2NN_input.pl.pl <protnamefile> <M_prfdir> <G_prfdir> <IO_prfdir> <R_prfdir> <inoutdir> <size>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $M_prfdir = $ARGV[1]."/";
my $G_prfdir = $ARGV[2]."/";
my $IO_prfdir = $ARGV[3]."/";
my $R_prfdir = $ARGV[4]."/";
my $outdir = $ARGV[5]."/";
my $size = $ARGV[6];

my $WING = int($size/2);

my $COLSIZE = 4;

open(PROTNAMEFILE, "$protnamefile")
  or die "could not open $protnamefile\n";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my $inoutfile = "$outdir"."$seqname".".input";
  
  
  open(INOUTFILE, ">"."$inoutfile")
    or die "could not open $inoutfile\n";
 
  my @Mcols = ();
  my @labels = ();
  open(PRFFILE, "$M_prfdir"."$seqname".".prf")
    or die "could not open Mprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @Mcols, $rowlist[2];
    }
    else {
	
    }
  }
  close(PRFFILE);
  
  my @Gcols = ();
  open(PRFFILE, "$G_prfdir"."$seqname".".prf")
    or die "could not open Gprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @Gcols, $rowlist[2];
      push @labels, $rowlist[5];
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  my @IOcols = ();
  open(PRFFILE, "$IO_prfdir"."$seqname".".prf")
    or die "could not open IOprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @IOcols, $rowlist[2];
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  my @Rcols = ();
  open(PRFFILE, "$R_prfdir"."$seqname".".prf")
    or die "could not open Rprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @Rcols, $rowlist[2];
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  if($#Mcols != $#Gcols || $#Gcols != $#IOcols || $#Rcols != $#IOcols || $#Rcols != $#Mcols) {
      die "Unequal no. cols: $#Mcols $#Gcols $#IOcols $#Rcols \n";
  }


  my @wholeinput = ();
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput, (0.0, 0.0, 0.0, 0.0);
  }
  for(my $i = 0; $i <= $#Mcols; $i++) {
      push @wholeinput, ($Mcols[$i], $Gcols[$i], $IOcols[$i], $Rcols[$i]);  
  }
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput,(0.0, 0.0, 0.0, 0.0);
  }
  
  my @output = splice @wholeinput, 0, $size * $COLSIZE;
  my $inoutrow = "";
  my $outoutrow = "";
  for(my $i = 0; $i <= $#output; $i++) {
      $inoutrow .= "$output[$i]\t";
  }
  ##print " << $#wholeinput   $#labels\n";
  
  for(my $i = 0;;$i++) {
      if($i > $#labels) {
	  print "$i  $#labels\n";
	  die "$seqname\n";
      }
      my $label = $labels[$i];
      
      for(my $j = 0; $j <= $#output; $j++) {
	  printf INOUTFILE "%.4f\t", $output[$j];
      }
      print INOUTFILE "\n";
      
      if($#wholeinput < 0) {
	  last;
      }
      #print "@output\n";
      #print "@wholeoutput\n";
      splice @output, 0, $COLSIZE;
      push @output, (splice @wholeinput, 0, $COLSIZE);
      #print "@output\n";
  }
  
  close INOUTFILE;
  
  #print "$seqname\n";
}
