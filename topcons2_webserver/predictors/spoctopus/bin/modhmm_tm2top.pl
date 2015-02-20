#!/usr/bin/perl -w

# reads a resultfile from modhmms_tm and outputs a number of
# .shorttop and .top files

my $infile = $ARGV[0];
my $shorttopdir;
my $topdir;
if($#ARGV > 0) {
  $shorttopdir = $ARGV[1]."/";
}
else {
  $shorttopdir = "./";
}

if($#ARGV > 1) {
  $topdir = $ARGV[2]."/";
}
else {
  $topdir = "./";
}

if($#ARGV != 2 && $#ARGV != 1 && $#ARGV != 0) {
  printf "Usage: modhmm_tm2top.pl <modhmm_tm resfile>";
  printf " [.shorttop outputdir] [.top outputdir]\n";
  exit;
}

open(INFILE,$infile) or die "Could not open $infile";

my $nr_regions = 0;;
my $N_term = "N";
my $seqfile;
my $seqname;
my $seq_length;
my $N_term_short;
my $relscore;

my $i;


my $inside = -1;
my $seq = "";

while(<INFILE>) {
  if($_ =~ '#') {
    next;
  }
  elsif($_ =~ '.fa' || $_ =~ '.mod' || $_ =~ '.mseq' || $_ =~ '.prf') {
    close(OUT);
    $inside = 0;
    $seqfile = $_;
    my @seqfile = split /\./, $seqfile;
    $seqname = "";
    for(my $i = 0; $i < $#seqfile; $i++) {
	$seqname .= $seqfile[$i];
	if($i < $#seqfile - 1) {
	    $seqname .= ".";
	}
    }
    
    open(OUT,">"."$topdir"."$seqname".".top");
    open(OUTSHORT,">"."$shorttopdir"."$seqname".".shorttop");
  }
  elsif($_ =~ 'Seq length') {
    (my $stuff, my $length) = split /:/, $_;
    $seq_length = $length;
  }
  elsif($_ =~ 'Labeling:') {
    $inside = 1;
  }
  elsif($inside == 1) {
    chomp;
    $seq .= $_;
    if(length $seq == $seq_length) {
      $N_term_short = substr($seq,0,1);
      my @seq_list = split //, $seq;
      for(my $j = 1; $j <= $#seq_list; $j++) {
	if(($seq_list[$j] eq 'M' || $seq_list[$j] eq 'W') && ($seq_list[$j-1] ne 'M' && $seq_list[$j-1] ne 'W')) {
	  $nr_regions++;
	}
      }
      
      #print short file
      printf OUTSHORT ">"."$seqname\n";
      printf OUTSHORT "$nr_regions";
      printf OUTSHORT "|$N_term_short\n";
      close(OUTSHORT);
      
      printf OUT ">"."$seqname\n";
      for(my $j = 0; $j <= $#seq_list; $j++) {
	if($j % 60 == 0 && $j != 0) {
	  printf OUT "\n";
	}
	printf OUT "$seq_list[$j]";
      }
      printf OUT "\n";
      close(OUT);
      
      $seq = "";
      $nr_regions = 0;
      $N_term_short = "N";
      $inside = -1;
    }
    elsif(length $seq > $seq_length) {
      print "sequence is longer than it should be\n";
      exit(0);
    }
  }
}
