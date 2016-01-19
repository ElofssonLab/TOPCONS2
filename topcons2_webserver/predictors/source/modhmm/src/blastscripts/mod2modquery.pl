if($#ARGV != 2) {
  printf "Usage: mod2modquery.pl <protnamefile> <modfiledir> <modfileoutdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $modfiledir = $ARGV[1]."/";
my $modfileoutdir = $ARGV[2]."/";

open(PROTNAMEFILE,"$protnamefile") or die "Could not open $protnamefile\n";
while(<PROTNAMEFILE>) {
  chomp;
  my $protname = $_;
  my $modfile = "$modfiledir"."$protname".".mod";
  my $modoutfile = "$modfileoutdir"."$protname".".mod";

  my @first_row = ();
  my $first = "YES";

  open(MODFILE,"$modfile") or die "Could not open $modfile\n";
  while(<MODFILE>) {
    chomp;
    my $row = $_;
    $row =~ s/<//g;
    $row =~ s/>//g;
    my @row_list = split /;/, $row;

    for(my $i = 0; $i <= $#row_list; $i++) {
      if($first eq 'YES') {
	@first_row = @row_list;
	last;
      }
    }
    $first = "NO";
  }

  close MODFILE;
  open(MODFILE,"$modfile") or die "Could not open $modfile\n";
  open(OUTFILE,">"."$modoutfile") or die "Could not open $modoutfile\n";
  while(<MODFILE>) {
    chomp;
    my $row = $_;
    $row =~ s/<//g;
    $row =~ s/>//g;
    my @row_list = split /;/, $row;
    
    print OUTFILE "<";
    for(my $i = 0; $i <= $#row_list; $i++) {
      if($first_row[$i] eq ' ' || $first_row[$i] eq '-' || $first_row[$i] eq '.' || $first_row[$i] eq '_') {
      }
      else {
	print OUTFILE "$row_list[$i]".";";
      }
    }
    print OUTFILE ">\n";
  }

  print "finished running $protname\n";
}



