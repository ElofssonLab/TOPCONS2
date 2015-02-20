#!/usr/bin/perl -w

if($#ARGV != 2) {
  printf "Usage: ioM29_state.pl <protnamefile> <intopodir> <outtopodir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $intopodir = $ARGV[1]."/";
my $outtopodir = $ARGV[2]."/";


open(PROTNAMEFILE, "$protnamefile")
    or die "could not open $protnamefile";


while(<PROTNAMEFILE>) {
  chomp;
  my $prot = $_;
  
  my $topofile = "$intopodir"."$prot".".top";
  open(TOPOFILE, "$topofile")
      or die "could not open $topofile";
  
  my $outfile = "$outtopodir"."$prot".".top";
  open(OUTFILE, ">"."$outfile")
      or die "could not open $outfile";

  my $topo = "";
  my $seqname = "";
  while(<TOPOFILE>) {
    chomp;
    if($_ =~ ">") {
      $seqname = $_;
    }
    else{
      $topo .= $_;
    }
  }
  
  ### if no TM reg, relabel with all g ###
  if($topo !~ 'M' && $topo !~ 'S') {
      $topo =~ s/./g/g;
  }

  
  my @topo_list = split //, $topo;
  print OUTFILE "$seqname\n";
  for (my $i = 0; $i <= $#topo_list; $i++) {
    if($i % 60 == 0 && $i != 0) {
      print OUTFILE "\n";
    }
    print OUTFILE "$topo_list[$i]";
  }
  print OUTFILE "\n";

  close TOPOFILE;
  close OUTFILE;
}
