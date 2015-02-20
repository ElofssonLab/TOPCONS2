#!/usr/bin/perl -w

if($#ARGV != 5) {
  printf "Usage: prf2KR_WY_NN_input_scoring.pl <protnamefile> <prfdir> <inoutdir> <outoutdir> <first_winsize> <second_winsize>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $outdir = $ARGV[2]."/";
my $outoutdir = $ARGV[3]."/";
my $winsize = $ARGV[4];
my $s_winsize = $ARGV[5];


my $COLSIZE = 20;
my $WING = int($winsize/2); # = window size *2+1
my $S_WING = int($s_winsize/2); # = window size *2+1

my @background_freqs = (0.075520, 0.016973, 0.053029, 0.063204, 0.040762, 0.068448, 0.022406, 0.057284, 0.059398, 0.093399, 0.023569, 0.045293, 0.049262, 0.040231, 0.051573, 0.072214, 0.057454, 0.065252, 0.012513, 0.031985);
my @nullfreqs = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

my @use_aas = (1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1);
open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
  chomp;
  $seqname = $_;

  my $outfile = "$outdir"."$seqname".".input";
  my $outoutfile = "$outdir"."$seqname".".desired_output";
  
  open(INOUTFILE, ">"."$outfile")
      or die "could not open $outfile";
  open(OUTOUTFILE, ">"."$outoutfile")
      or die "could not open $outoutfile";
  

 
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
      my $tot = 0;
      for(my $i = 0; $i <= $#rowlist; $i++) {
	  $tot += $rowlist[$i];
      }
      for(my $i = 0; $i <= $#rowlist;$i++) {
	  $rowlist[$i] /= $tot;
      }
      push @cols, @rowlist;
     
      

    }
    else {
      
    }
  }
  close(PRFFILE);


 
  my @seq_background_freqs = (0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001);

  my $nr_cols = 0;
  for(my $i = 0; $i <= $#cols; $i++) {
      $seq_background_freqs[$i % $COLSIZE] += $cols[$i];
      if($i % $COLSIZE == 0) {
	  $nr_cols++;
      }
  }
  for(my $i = 0; $i <= $#seq_background_freqs; $i++) {
      $seq_background_freqs[$i] /= $nr_cols;
  }
  $seq_background_freqs[8] += $seq_background_freqs[14]; # add R to K
  $seq_background_freqs[19] += $seq_background_freqs[18]; # add W to Y
  #print "@seq_background_freqs\n";

  my @wholeinput = ();
  for(my $i = 0; $i < $WING; $i++) {
      #push @wholeinput, @background_freqs;
      push @wholeinput, @nullfreqs;
  }
  push @wholeinput, @cols;
  for(my $i = 0; $i < $WING; $i++) {
      #push @wholeinput, @background_freqs;
      push @wholeinput, @nullfreqs;
  }
  


  my @aa_avg_cols = ();
  my $max_aa_avg_KR = 0;
  my $max_aa_avg_WY = 0;
  my $max_aa_avg_D = 0;
  my $max_aa_avg_E = 0;
  my $max_aa_avg_G = 0;
  my $max_aa_avg_N = 0;
  for(my $i = 0; $i <= $#wholeinput - ($winsize * $COLSIZE) + 1; $i = $i + $COLSIZE) {
      my @avg_col = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
      for(my $j = 0; $j < $winsize; $j++) {
	  for(my $k = 0; $k < $COLSIZE; $k++) {
	      my $index = $i + $k + ($COLSIZE * $j);
	      #print "$i $j $k $index\n";
	      $avg_col[$k] += $wholeinput[$i + ($COLSIZE * $j) + $k];
	  }
      }
      $avg_col[8] += $avg_col[14]; # add R to K
      $avg_col[19] += $avg_col[18]; # add W to Y
      $avg_col[0] =  $seq_background_freqs[8];
      for(my $j = 1; $j <= $#avg_col; $j++) {
	  $avg_col[$j] /= $winsize;
	  #$avg_col[$j] /= $seq_background_freqs[$j];
	  if($use_aas[$j] == 0) {
	      $avg_col[$j] = 0;
	  }
	  #printf "%.3f ", $avg_col[$j];
      }
      #print "\n";
      push @aa_avg_cols, @avg_col;
      if($avg_col[8] > $max_aa_avg_KR) {
	  $max_aa_avg_KR = $avg_col[8];
      }
      if($avg_col[19] > $max_aa_avg_WY) {
	  $max_aa_avg_WY = $avg_col[19];
      }
      if($avg_col[2] > $max_aa_avg_D) {
	  $max_aa_avg_D = $avg_col[2];
      }
      if($avg_col[3] > $max_aa_avg_E) {
	  $max_aa_avg_E = $avg_col[3];
      }
      if($avg_col[5] > $max_aa_avg_G) {
	  $max_aa_avg_G = $avg_col[5];
      }
      if($avg_col[11] > $max_aa_avg_N) {
	  $max_aa_avg_N = $avg_col[11];
      }
  }
  
  if($max_aa_avg_KR > 0) {
      for(my $i = 8; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_KR;
      }
  }
  if($max_aa_avg_WY > 0) {
      for(my $i = 19; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_WY;
	  $aa_avg_cols[$i] = 0;
      }
  }
  if($max_aa_avg_D > 0) {
      for(my $i = 2; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_D;
      }
  }
  if($max_aa_avg_E > 0) {
      for(my $i = 3; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_E;
      }
  }
  if($max_aa_avg_G > 0) {
      for(my $i = 5; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_G;
      }
  }
  if($max_aa_avg_N > 0) {
      for(my $i = 11; $i <= $#aa_avg_cols; $i = $i + 20) {
	  $aa_avg_cols[$i] /= $max_aa_avg_N;
      }
  }
  @wholeinput = @aa_avg_cols;


  my @second_cols = ();
  for(my $i = 0; $i < $S_WING; $i++) {
      push @second_cols, [0,0];
  }
  my @output = splice @wholeinput, 0, $COLSIZE;
  for(my $i = 0;;$i++) {
      #print "$i\n";
      push @second_cols, [$output[8],$output[19]];
      
      if($#wholeinput <= 0) {
	  last;
      }
      splice @output, 0, $COLSIZE;
      push @output, (splice @wholeinput, 0, $COLSIZE);
  }
  for(my $i = 0; $i < $S_WING; $i++) {
      push @second_cols, [0,0];
  }
  
  for(my $i = 0; $i <= $#labels; $i++) {
      my $use_col = "YES";
      $outoutrow = "0 0\n";
      if($labels[$i] eq 'i') {
	  $outoutrow = "1 0\n";
      }
      elsif($labels[$i] eq 'o') {
	  $outoutrow = "0 1\n"; 
      }
      else {
	  $use_col = "NO";
      }
      if($use_col eq 'YES') {
	  for(my $j = $i; $j < $i + $s_winsize; $j++) {
	      printf INOUTFILE "%.4f\t", $second_cols[$j][0];
	      printf INOUTFILE "%.4f\t", $second_cols[$j][1];
	  }
	  print INOUTFILE "\n";
	  print OUTOUTFILE "$outoutrow";
      }
  }

  close INOUTFILE;
  close OUTOUTFILE;
  
  #print "$seqname\n";
  #exit;
}
