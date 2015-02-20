#!/usr/bin/perl -w

if($#ARGV != 3) {
  printf "Usage: NN_pred2prf_2.pl <protnamefile> <predfile> <prfdir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $predfile = $ARGV[1];
my $prfdir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";




my @preds  =();
open(PREDFILE, "$predfile")
  || die "could not open $predfile";

while(<PREDFILE>) {
  chomp;
  my $row = $_;
  $row =~ s/\s+/ /g;
  $row =~ s/^\s+//g;
  $row =~ s/\s+$//g;
  my @row_list = split / /, $row;
  push @preds, [@row_list];

}
close PREDFILE;

open(PROTNAMEFILE, "$protnamefile")
    || die "could not open $protnamefile";




my $predpos = 0;
while(<PROTNAMEFILE>) {
    chomp;
    $seqname = $_;
    my $seq_length = 0;
    
    open(PRFOUTFILE, ">"."$outdir"."$seqname".".prf")
	or die "could not open prfoutfile for $seqname";
    
    open(PRFFILE, "$prfdir"."$seqname".".prf")
	or die "could not open prffile for $seqname";
    while(<PRFFILE>) {
	chomp;
	my $row = $_;
	
	if($row =~ 'Sequence') {
	    print PRFOUTFILE "$row\n";
	}
	elsif($row =~ 'NR of aligned') {
	    print PRFOUTFILE "$row\n";
	}
	elsif($row =~ 'Length of') {
	    print PRFOUTFILE "$row\n";
	    (my $slask, $seq_length) = split /:/, $row;
	    $seq_length =~ s/^s+//g;
	    #print "$seqname $seq_length\n";
	}
	elsif($row =~ 'ALPHABET') {
	    print PRFOUTFILE "ALPHABET:   i       o     <SPACE> <LABEL> <QUERY> \n";
	}
	elsif($row =~ 'COL') {
	    $row =~ s/\s+/ /g;
	    my @col = split / /, $row;
	    my $col_nr = $col[1];
	    $col_nr =~ s/://g;
	    
	    if(!$preds[$predpos][0]) {
		die "no $seqname";
	    }
	    
	    printf PRFOUTFILE "COL %4d".":   ", $col_nr;
	    if( $preds[$predpos][0] < 0) {
		$preds[$predpos][0] = 0;
	    }
	    printf PRFOUTFILE "%.2f    %.2f    ", $preds[$predpos][0], $preds[$predpos][1];
	    printf PRFOUTFILE "0.00   0.00   $col[24]       $col[25]       \n";
	    $predpos++;
	}
	else {
	    print PRFOUTFILE "$row\n";
	}
    }
    close PRFFILE;
    close PRFOUTFILE;
  
  
    
}
