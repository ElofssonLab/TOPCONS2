#!/usr/bin/perl -w

if($#ARGV != 9) {
  printf "Usage: NN_profiles2NN_input_I_vs_O.pl.pl <protnamefile> <topodir> <M_prfdir> <G_prfdir> <IO_prfdir> <IO_prfdir> <inoutdir> <outoutdir> <scoring/training> <size>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $topodir = $ARGV[1]."/";
my $M_prfdir = $ARGV[2]."/";
my $G_prfdir = $ARGV[3]."/";
my $I_prfdir = $ARGV[4]."/";
my $O_prfdir = $ARGV[5]."/";
my $inoutdir = $ARGV[6]."/";
my $outoutdir = $ARGV[7]."/";
my $mode = $ARGV[8];
my $size = $ARGV[9];

my $WING = int($size/2);

my $COLSIZE = 2;

open(PROTNAMEFILE, "$protnamefile")
  or die "could not open $protnamefile\n";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my $outoutfile = "$outoutdir"."$seqname".".desired_output";
  my $inoutfile = "$inoutdir"."$seqname".".input";
  
  open(OUTOUTFILE, ">"."$outoutfile")
    or die "could not open $outoutfile\n";
  open(INOUTFILE, ">"."$inoutfile")
    or die "could not open $inoutfile\n";
  
  my $topo = "";
  open(TOPOFILE, "$topodir"."$seqname".".top")
    or die "could not open topofile for $seqname\n";
  while(<TOPOFILE>) {
    chomp;
    if($_ =~ '>') {
      
    }
    else {
	$topo .= $_;
    }
  }
  close(TOPOFILE);
  my @topo = split //, $topo;
 
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
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  my @Icols = ();
  open(PRFFILE, "$I_prfdir"."$seqname".".prf")
    or die "could not open Iprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @Icols, $rowlist[2];
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  my @Ocols = ();
  open(PRFFILE, "$O_prfdir"."$seqname".".prf")
    or die "could not open Oprffile for $seqname\n";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      push @Ocols, $rowlist[3];
      push @labels, $rowlist[6];
    }
    else {
      
    }
  }
  close(PRFFILE);
  
  if($#Mcols != $#Gcols || $#Gcols != $#Icols || $#Ocols != $#Icols || $#Ocols != $#Mcols) {
      die "$#Mcols $#Gcols $#Ocols $#Icols \n";
  }


  my @wholeinput = ();
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput, (0.0, 0.0);
  }
  for(my $i = 0; $i <= $#Icols; $i++) {
      push @wholeinput, ($Icols[$i], $Ocols[$i]);
  }
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput,(0.0, 0.0);
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
      my $outoutrow = "";
      my $use_col = "YES";
      my $label = $labels[$i];
      
      $outoutrow = "0 0\n";
      if($label eq '.' || $label eq 'T') {
	  if($mode eq 'training') {
	      $use_col = "NO";
	  }
	  elsif($mode eq 'scoring') {
	      $use_col = "YES";
	  }
      }
      elsif($label eq 'i') {
	  $outoutrow = "1 0\n";
      }
      elsif($label eq 'o') {
	  $outoutrow = "0 1\n";
      }
      else {
	  if($mode eq 'training') {
	      $use_col = "NO";
	  }
	  elsif($mode eq 'scoring') {
	      $use_col = "YES";
	  }
      }
      #print "$use_col\n";
      if($use_col eq 'YES') {
	  for(my $j = 0; $j <= $#output; $j++) {
	      printf INOUTFILE "%.4f\t", $output[$j];
	  }
	  print INOUTFILE "\n";
	  print OUTOUTFILE "$outoutrow";
      }
      
      if($#wholeinput < 0) {
	  last;
      }
      #print "@output\n";
      #print "@wholeoutput\n";
      splice @output, 0, $COLSIZE;
      push @output, (splice @wholeinput, 0, $COLSIZE);
      #print "@output\n";
  }
  
  close OUTOUTFILE;
  close INOUTFILE;
  
  print "$seqname\n";
}
