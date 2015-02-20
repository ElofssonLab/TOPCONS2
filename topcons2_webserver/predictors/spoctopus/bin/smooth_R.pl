#!/usr/bin/perl -w

if($#ARGV != 2) {
  printf "Usage: get_M_vs_R_stats.pl <protnamefile> <prfdir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $prfoutdir = $ARGV[2]."/";


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";
while(<PROTNAMEFILE>) {  
  chomp;
  $seqname = $_;
  
 

  my @cols = ();
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
      
      if($_ =~ 'COL ' && $instart eq "NO" && $inend eq "NO" ) {
	  $_ =~ s/\s+/ /g;
	  my @rowlist = split / /, $_;   
	  push @cols, [@rowlist];
      }
      else {
	  
      }
  }
  close(PRFFILE);
  
  
 


  my @new_R_cols = ();


  
  for(my $i = 0; $i <= $#cols; $i++) {
      my $high_R = 0.0;
      for(my $j = $i - 1; $j <= $i + 1; $j++) {
	  if($j >= 0 && $j <= $#cols) {
	      if($cols[$j][5] > $high_R) {
		  $high_R = $cols[$j][5];
	      }
	  }
      }
      push @new_R_cols, $high_R;
  }
  

  open(PRFOUTFILE, ">"."$prfoutdir"."$seqname".".prf")
      or die "could not open prfoutfile for $seqname";
  
  print PRFOUTFILE "$start_prf";
  for(my $i = 0; $i <= $#cols; $i++) {
      my @col = @{$cols[$i]};
      printf PRFOUTFILE "$col[0]  $col[1]  %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   $col[8]       $col[9]\n", $col[2],$col[3],$col[4],$new_R_cols[$i],$col[6],$col[7];
  }
  print PRFOUTFILE "$end_prf";
  close PRFOUTFILE;
  

}
