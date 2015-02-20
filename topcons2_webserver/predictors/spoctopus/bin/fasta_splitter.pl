#!/usr/bin/perl -w

if($#ARGV != 1) {
  printf "Usage: fasta_splitter.pl <fastafile> <outdir>\n";
  exit;
}

my $fastafile = $ARGV[0];
my $outdir = $ARGV[1]."/";

my $protname = "";
my $seq = "";
my $first = 'YES';
my $outfile = "";
open(FASTAFILE, "$fastafile")
    || die "could not open $fastafile";
while(<FASTAFILE>) {
  chomp $_;
  my $row = $_;
  #$row =~ s/\s//g;
  #$row =~ s/\t//g;
  #print "row = $row\n";
  if($row =~ ">" && $first eq 'YES') {
    
    my @row_list = split / /, $row;
    $protname = substr $row_list[0], 1;
    $protname =~ s/\s//g;
    #print "$protname\n";
    $first = 'NO';
  }
  elsif($row =~ ">") {
    $outfile = "$outdir"."$protname".".fa";
    open(OUTFILE, ">"."$outfile")
	|| die "could not open $outfile";
    print OUTFILE ">"."$protname\n";
    my @seq_list = split //, $seq;
    for(my $j = 0; $j <= $#seq_list; $j++) {
      print OUTFILE "$seq_list[$j]";
      if(($j+1) != 0 && ($j + 1) % 60 == 0) {
	print OUTFILE "\n";
      }
    }
    print OUTFILE "\n";
    close OUTFILE;
    $seq = "";

    #print "$row\n";
    my @row_list = split / /, $row;
    $protname = substr $row_list[0], 1;
    $protname =~ s/\s//g;
    #print "$protname\n";
  }
  elsif($row ne "" && $row ne "\n") {
    $seq .= $row;
  }
}

$outfile = "$outdir"."$protname".".fa";

#print "outfile = $outfile\n";
#print "seq = $seq\n";
open(OUTFILE, ">"."$outfile")
    || die "could not open $outfile";
print OUTFILE ">"."$protname\n";
my @seq_list = split //, $seq;
for(my $j = 0; $j <= $#seq_list; $j++) {
  print OUTFILE "$seq_list[$j]";
  if(($j+1) != 0 && ($j + 1) % 60 == 0) {
    print OUTFILE "\n";
  }
}
print OUTFILE "\n";
close OUTFILE;
