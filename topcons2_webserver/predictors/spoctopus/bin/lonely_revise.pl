#!/usr/bin/perl -w

if($#ARGV != 3) {
  printf "Usage: lonely_revise.pl <protnamefile> <prfdir> <topodir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $topodir = $ARGV[2]."/";
my $prfoutdir = $ARGV[3]."/";

my $GLOBCUTOFF = 60;


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";
while(<PROTNAMEFILE>) {  
  chomp;
  $seqname = $_;
  
 

  my @cols = ();
  my @labels = ();
  my $start_prf = "";
  my $end_prf = "";
  my $instart = "YES";
  my $inend = "NO";
  open(PRFFILE, "$prfdir"."$seqname".".prf")
      or die "could not open prffile for $seqname";
  while(<PRFFILE>) {
      chomp;
      if($_ =~ 'START 1') {
	  $start_prf .= "$_\n";
	  $instart = "NO";
      }
      if($_ =~ 'END 1') {
	  $inend = "YES";
      }

      if($instart eq 'YES') {
	  $start_prf .= "$_\n";
      }
      if($inend eq 'YES') {
	  $end_prf .= "$_\n";
      }
      
      if($_ =~ 'COL ' && $instart eq "NO" && $inend eq "NO") {
	  $_ =~ s/\s+/ /g;
	  my @rowlist = split / /, $_;   
	  push @cols, [@rowlist];
      }
      else {
	  
      }
  }
  close(PRFFILE);
  
  
  my $topofile = "$topodir"."$seqname".".top";
  open(TOPOFILE, "$topofile")
      or die "could not open $topofile";
  my $topo = "";
  while(<TOPOFILE>) {
      chomp;
      if($_ =~ '>') {

      }
      else {
	  $topo .= $_;
      }
  }
  my @topo = split //, $topo;
  close TOPOFILE;

 

  my @TM_starts = ();
  for(my $i = 1; $i <= $#topo; $i++) {
      if($topo[$i] eq 'M' && $topo[$i-1] ne 'M') {
	  push @TM_starts, $i;
      }
  }



 


  for(my $i = 0; $i <= $#TM_starts; $i++) {
      my $TM_dist = 10000;
      if(($i-1 < 0 || $TM_starts[$i] - $TM_starts[$i-1] > $GLOBCUTOFF) && ($i+1 > $#TM_starts || $TM_starts[$i+1] - $TM_starts[$i] > $GLOBCUTOFF)) {
	  if($#TM_starts <= 0) {
	      next;
	  }
	  if($i > 0) {
	      $TM_dist =  $TM_starts[$i] - $TM_starts[$i-1];
	  }
	  if($i < $#TM_starts && $TM_starts[$i+1] - $TM_starts[$i] < $TM_dist) {
	      $TM_dist =  $TM_starts[$i+1] - $TM_starts[$i];
	  }
	  if(($i == 0 || $i == $#TM_starts) && $TM_dist < $GLOBCUTOFF * 2) {
	      next;
	  }
      }
      else {
	  next;
      }

      my $TM_end = 0;
      my $TM_start = $TM_starts[$i];
      for(my $j = $TM_start; $j <= $#topo; $j++) {
	  if($topo[$j-1] eq 'M' && $topo[$j] ne 'M') {
	      $TM_end = $j-1;
	      last;
	  }
      }
      
      my $TM_length = $TM_end - $TM_start + 1;
      my $avg_I = 0;
      my $avg_M = 0;
      my $avg_L = 0;
      my $avg_G = 0;
      for(my $k = $TM_start; $k <= $TM_end; $k++) {
	  my $I_val = $cols[$k][5];
	  my $M_val = $cols[$k][2];
	  my $L_val = $cols[$k][3];
	  my $G_val = $cols[$k][4];
	  $I_val += 0.000001;
	  $M_val += 0.000001;
	  $L_val += 0.000001;
	  $G_val += 0.000001;
	  $I_val = log $I_val;
	  $M_val = log $M_val;
	  $L_val = log $L_val;
	  $G_val = log $G_val;
	  $avg_I += $I_val;
	  $avg_M += $M_val;
	  $avg_L += $L_val;
	  $avg_G += $G_val;
      }
      #$avg_I /= $TM_length - 6;
      #$avg_M /= $TM_length - 6;
      

      $avg_I = $avg_I / ($TM_length);
      $avg_I = exp $avg_I;
      $avg_L = $avg_L / ($TM_length);
      $avg_L = exp $avg_L;
      $avg_M = $avg_M / ($TM_length);
      $avg_M = exp $avg_M;
      $avg_G = $avg_G / ($TM_length);
      $avg_G = exp $avg_G;
      
      my $relabel = "NO";
      if($avg_M <= 1.0) {
	  if(((1.33 * $avg_M) - 0.61) < $avg_L) {
	      printf "$seqname ($TM_start - $TM_end): %.3f\t%.3f\n", $avg_M, $avg_L;
	      $relabel = "YES";
	  }
      }
      else {
	  
      }
      
      
      if($relabel eq 'YES') {
	  #print "$seqname ($TM_start - $TM_end)\n";
	  for(my $j = $TM_start; $j <= $TM_end; $j++) {
	      $cols[$j][2] = 0.005;
	      #$cols[$j][5] = 0.500;
	  }
      }
      
      
      #printf "$seqname\t\t$TM_start\t$TM_end\t%.3f\t%.3f\t$corr_pred\t$wrote\n", $avg_M, $avg_G;
      
      
      
  }
  




  open(PRFOUTFILE, ">"."$prfoutdir"."$seqname".".prf")
      or die "could not open prfoutfile for $seqname";
  
  print PRFOUTFILE "$start_prf";
   printf PRFOUTFILE "ALPHA:   M       L       G       I       -      <SPACE> <LABEL> <QUERY>\n";      
  for(my $i = 0; $i <= $#cols; $i++) {
      my @col = @{$cols[$i]};
      printf PRFOUTFILE "$col[0]  $col[1]  %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   $col[8]       $col[9]\n", $col[2],$col[3],$col[4],$col[5],$col[6],$col[7];
  }
  print PRFOUTFILE "$end_prf";
  close PRFOUTFILE;
}
