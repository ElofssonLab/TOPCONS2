#!/usr/bin/perl -w

if($#ARGV != 3) {
  printf "Usage: prf2NN_input_fast_scoring.pl <protnamefile> <prfdir> <outdir> <winsize> \n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $outdir = $ARGV[2]."/";;
my $winsize = $ARGV[3];


my $COLSIZE = 20;
my $WING = int($winsize/2); # = window size *2+1

my @background_freqs = (0.075520, 0.016973, 0.053029, 0.063204, 0.040762, 0.068448, 0.022406, 0.057284, 0.059398, 0.093399, 0.023569, 0.045293, 0.049262, 0.040231, 0.051573, 0.072214, 0.057454, 0.065252, 0.012513, 0.031985);


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my $outfile = "$outdir"."$seqname".".input";
  
  open(INOUTFILE, ">"."$outfile")
      or die "could not open $outfile";
  

 
  my @cols = ();
  open(PRFFILE, "$prfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
      chomp;
      if($_ =~ 'COL ') {
	  $_ =~ s/\s+/ /g;
	  my $rowsum = 0;
	  my @rowlist = split / /, $_;
	  #print "@rowlist\n";
	  splice @rowlist, 0, 2;
	  splice @rowlist, 20;
	  for(my $i = 0; $i <= $#rowlist; $i++) {
	      $rowsum += $rowlist[$i];
	  }
	  if($rowsum > 0) {
	      for(my $i = 0; $i <= $#rowlist; $i++) {
		  $rowlist[$i] /= $rowsum;
	      }
	  }
	  push @cols, @rowlist;
      }
      else {
      
      }
  }
  close(PRFFILE);

  
  my @wholeinput = ();
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput, @background_freqs;
  }
  push @wholeinput, @cols;
  for(my $i = 0; $i < $WING; $i++) {
      push @wholeinput, @background_freqs;
  }
  
  my @output = splice @wholeinput, 0, $winsize * $COLSIZE;
  my $inoutrow = "";
  for(my $i = 0; $i <= $#output; $i++) {
      $inoutrow .= "$output[$i]\t";
  }
  for(my $i = 0;;$i++) {
      #print "$i\n";
      for(my $j = 0; $j <= $#output; $j++) {
	  printf INOUTFILE "%.4f\t", $output[$j];
      }
      print INOUTFILE "\n";
      
      if($#wholeinput <= 0) {
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
