#!/usr/bin/perl -w

if($#ARGV != 2) {
  printf "Usage: make_pred_output.pl <protnamefile> <prfdir> <prfoutdir>\n";
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
  my $io_cols = ();
  my @labels = ();
  my $start_prf = "";
  my $end_prf = "";
  my $instart = "YES";
  my $inend = "NO";
  my $start_prf_2 = "";
  my $end_prf_2 = "";
  my $instart_2 = "YES";
  my $inend_2 = "NO";
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
      if($_ =~ 'START 2') {
	  $start_prf_2 .= "$_\n";
	  $instart_2 = "NO";
      }
      if($_ =~ 'END 2') {
	  $inend_2 = "YES";
      }
      if($instart_2 eq 'YES') {
	  $start_prf_2 .= "$_\n";
      }
      if($inend_2 eq 'YES') {
	  $end_prf_2 .= "$_\n";
      }
      if($_ =~ 'COL ' && $instart_2 eq "NO" && $inend_2 eq "NO") {
	  $_ =~ s/\s+/ /g;
	  my @rowlist = split / /, $_;   
	  push @io_cols, [@rowlist];
      }

      else {
	  
      }
  }
  close(PRFFILE);

  open(PRFOUTFILE, ">"."$prfoutdir"."$seqname".".nnprf")
      or die "could not open prfoutfile for $seqname";
  
  #print PRFOUTFILE "$start_prf";
  printf PRFOUTFILE "Pos\tM       L       G       R       S       i       o      AA\n";
  for(my $i = 0; $i <= $#cols; $i++) {
      my @col = @{$cols[$i]};
      $col[1] =~ s/://g;
      my @io_col = @{$io_cols[$i]};
      printf PRFOUTFILE "%d\t%.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   $col[10]\n", $col[1], $col[2],$col[3],$col[4],$col[5], $col[6], $io_col[2], $io_col[3];
  }
  #print PRFOUTFILE "$end_prf";
  close PRFOUTFILE;
}
