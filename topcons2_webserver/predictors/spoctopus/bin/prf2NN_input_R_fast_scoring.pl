#!/usr/bin/perl -w

if($#ARGV != 8) {
  printf "Usage: z_input.pl <protnamefile> <zdir> <prfdir> <inoutdir> <outoutdir> <scoring/training> <winsize> <inlimit> <outlimit>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $zdir = $ARGV[1]."/";
my $prfdir = $ARGV[2]."/";
my $inoutdir = $ARGV[3]."/";
my $outoutdir = $ARGV[4]."/";
my $mode = $ARGV[5];
my $winsize = $ARGV[6];
my $inlimit = $ARGV[7];
my $outlimit = $ARGV[8];


my $COLSIZE = 20;
my $WING = int($winsize/2); # = window size *2+1
my $NET_CHARGE_WING = int($winsize/2); # = window size *2+1
my $CHARGE_FRAC_WING = int($winsize/2); # = window size *2+1

my @background_freqs = (0.075520, 0.016973, 0.053029, 0.063204, 0.040762, 0.068448, 0.022406, 0.057284, 0.059398, 0.093399, 0.023569, 0.045293, 0.049262, 0.040231, 0.051573, 0.072214, 0.057454, 0.065252, 0.012513, 0.031985);


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my $outoutfile = "$outoutdir"."$seqname".".desired_output";
  my $inoutfile = "$inoutdir"."$seqname".".input";
  
  open(OUTOUTFILE, ">"."$outoutfile")
      or die "could not open $outoutfile";
  open(INOUTFILE, ">"."$inoutfile")
      or die "could not open $inoutfile";
  

 
  my @cols = ();
  my @labels = ();
  open(PRFFILE, "$prfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
    chomp;
    if($_ =~ 'COL ') {
      $_ =~ s/\s+/ /g;
      my @rowlist = split / /, $_;
      #print "@rowlist\n";
      push @labels, $rowlist[24];
      splice @rowlist, 0, 2;
      splice @rowlist, 20;
      push @cols, @rowlist;
      
      

    }
    else {
      
    }
  }
  close(PRFFILE);

  my @z_coords = ();
  if($mode ne 'scoring') {
    open(ZFILE, "$zdir"."$seqname".".zco")
	or die "could not open zfile for $seqname";
    while(<ZFILE>) {
      chomp;
      if($_ =~ '>') {
	
      }
      else {
	my $row = "$_";
	$row =~ s/\s+/ /g;
	@z_coords = split / /, $row;
      }
    }
    close(ZFILE);
  }
  else {
    for(my $i = 0; $i <= $#cols; $i++) {
      $z_coords[$i] = 0.0;
    }
  }

  
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
  my $outoutrow = "";
  for(my $i = 0; $i <= $#output; $i++) {
      $inoutrow .= "$output[$i]\t";
  }
  for(my $i = 0;;$i++) {
      #print "$i\n";
      my $outoutrow = "";
      my $use_col = "YES";
      my $label = $labels[$i];
      
      $outoutrow = "0\n";
      if($z_coords[$i] eq 'X') {
	  if($mode eq 'training') {
	      $use_col = "NO";
	  }
	  elsif($mode eq 'scoring') {
	      $use_col = "YES";
	  }
      }
      elsif($z_coords[$i] >= $inlimit && $z_coords[$i] <= $outlimit && ($label eq 'R' || $label eq 'r')) {
	  $outoutrow = "1\n";
      }
      elsif($z_coords[$i] >= $inlimit && $z_coords[$i] <= $outlimit && $label eq 'M') {
	  $outoutrow = "0\n";
      }
      elsif($z_coords[$i] < $inlimit || $z_coords[$i] > $outlimit) {
	  $outoutrow = "0\n";
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
      
      if($#wholeinput <= 0) {
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
