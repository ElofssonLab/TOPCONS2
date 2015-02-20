#!/usr/bin/perl -w

if($#ARGV != 4) {
  printf "Usage: get_NN_pred_stats_from_prf.pl <protnamefile> <Sdir> <length> <cutoff> <outdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $LENGTH = $ARGV[2];
my $PRED_CUTOFF = $ARGV[3];
my $outdir = $ARGV[4]."/";




open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
    chomp;
    my $prffile = "$prfdir"."$_".".prf";
    my $outfile = "$outdir"."$_".".sigpeppa";

    
    open(PRFFILE, "$prffile")
	or die "could not open $prffile";
    @cols = ();
    while(<PRFFILE>) {
	chomp;
	if($_ =~ 'COL ') {
	    $_ =~ s/\s+/ /g;
	    my @rowlist = split / /, $_;
	    push @cols, [@rowlist];
	    
	}
	else {
	    
	}
    }
    close(PRFFILE);
    
 
    





    my $lsum = 0;
    for(my $i = 0; $i < $LENGTH; $i++) {
	my @col = @{$cols[$i]};
	if($col[2] == 0) {
	    $col[2] = 0.0001;
	}
	my $logval = log $col[2];
	$lsum += $logval;
    }
    
    my $max_lsum = $lsum;
    my $max_lpos = 0;
    for(my $i = $LENGTH; $i <= $#cols; $i++) {
	if($i > 35 + $LENGTH) {
	    last;
	}
	my @col = @{$cols[$i]};
	my @lastcol = @{$cols[$i-$LENGTH]};
	if($col[2] == 0) {
	    $col[2] = 0.0001;
	}
	if($lastcol[2] == 0) {
	    $lastcol[2] = 0.0001;
	}
	$lsum += log $col[2];
	$lsum -= log $lastcol[2];
	if($lsum > $max_lsum) {
	    $max_lsum = $lsum;
	    $max_lpos = $i - $LENGTH;
	}
    }
    
    $max_lsum /= $LENGTH;
    my $max_gavg = exp $max_lsum;
    
    
    if($max_gavg > $PRED_CUTOFF) {
	open OUTFILE, ">"."$outfile"
	    or die "Could not open $outfile";
	print OUTFILE "YES\n";
	close OUTFILE;
    }
}
