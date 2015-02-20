#!/usr/bin/perl -w

if($#ARGV != 4) {
  printf "Usage: glob_filter.pl <protnamefile> <M_prfdir> <IO_prfdir> <G_prfdir> <R_prfdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $Mprfdir = $ARGV[1]."/";
my $IOprfdir = $ARGV[2]."/";
my $Gprfdir = $ARGV[3]."/";
my $Rprfdir = $ARGV[4]."/";



my @preds  =();

my $M_CUTOFF = 0.9;

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
      or die "could not open prffile for $seqname";
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
      or die "could not open prffile for $seqname";
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
      or die "could not open prffile for $seqname";
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
  
  
  my $is_glob = "YES";
  my $tot_M_over_L = 0;
  my $L = 14;
  my $max_avg_over_L = 0;
  for(my $i = 0; $i <= $#Mcols; $i++) {
      my @col = @{$Mcols[$i]};
      if($i < $L) {
	  $tot_M_over_L += $col[2];
	  my $avg = $tot_M_over_L / $L;
	  if($avg > $max_avg_over_L) {
	      $max_avg_over_L = $avg;
	  }
      }
      elsif($i > $#Mcols - $L) {
	  my @old_col =  @{$Mcols[$i-$L]};
	  $tot_M_over_L -= $old_col[2];
	  my $avg = $tot_M_over_L / $L;
	  if($avg > $max_avg_over_L) {
	      $max_avg_over_L = $avg;
	  }
      }
      else {
	  my @old_col =  @{$Mcols[$i-$L]};
	  $tot_M_over_L += $col[2];
	  $tot_M_over_L -= $old_col[2];
	  my $avg = $tot_M_over_L / $L;
	  if($avg > $max_avg_over_L) {
	      $max_avg_over_L = $avg;
	  }
      }

     
  }

  if($max_avg_over_L > $M_CUTOFF) {
      $is_glob = "NO";
  }

  print "$seqname: $is_glob\n";

}
