#!/usr/bin/perl -w

if($#ARGV != 2) {
  printf "Usage: normalize_io_prffiles.pl <protnamefile> <prfdir> <prfoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prfdir = $ARGV[1]."/";
my $outdir = $ARGV[2]."/";


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";


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
	    
	    
	    printf PRFOUTFILE "COL %4d".":   ", $col_nr;
	    
	    my $icol = $col[2];
	    my $ocol = $col[3];
	    
	    $icol = exp($icol) / exp(1);
	    $ocol = exp($ocol) / exp(1);

	    printf PRFOUTFILE "%.2f    %.2f    ", $icol, $ocol;
	    printf PRFOUTFILE "0.00   0.00   $col[4]       $col[5]       \n";
	}
	else {
	    print PRFOUTFILE "$row\n";
	}
    }
    close PRFFILE;
    close PRFOUTFILE;
  
  
    
}
