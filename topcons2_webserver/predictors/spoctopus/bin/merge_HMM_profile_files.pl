#!/usr/bin/perl -w


if($#ARGV != 3) {
  printf "Usage: merge_HMM_profile_files.pl <protnamefile> <MLGRdir> <iodir> <outdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $prffiledir = $ARGV[1]."/";
my $ssprffiledir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";




open(PROTNAMEFILE, "$protnamefile")
    || die "could not open $protnamefile";

while(<PROTNAMEFILE>) {
  chomp;
  my $seqname = $_;
  
  


  if(-e "$prffiledir"."$seqname".".prf") {
    open(PRFFILE, "$prffiledir"."$seqname".".prf")
	|| die "could not open $prffiledir"."$seqname".".prf";
  }
  else {
      die "could not open $prffiledir"."$seqname".".prf";
  }

  if(-e "$ssprffiledir"."$seqname".".prf") {
    open(SSPRFFILE, "$ssprffiledir"."$seqname".".prf")
	|| die "could not open $ssprffiledir"."$seqname".".prf";
  }
  else {
      die "could not open $ssprffiledir"."$seqname".".prf";
  }
 
 
 
 
  open(OUTFILE, ">"."$outdir"."$seqname".".prf")
      or die "could not open "."$outdir"."$seqname".".prf"."";
  
 
  while(<PRFFILE>) {
    print OUTFILE $_;
  }

  while(<SSPRFFILE>) {
    chomp;
    if($_ =~ 'Sequence' || $_ =~ 'NR' || $_ =~ 'Length' || length $_ == 0) {

    }
    elsif( $_ =~ 'START') {
      print OUTFILE "\nSTART 2\n";
    }
    elsif( $_ =~ 'END') {
      print OUTFILE "END 2\n";
    }
    else {
      print OUTFILE "$_\n";
    }
  }

  print "done merging $seqname\n";
  close OUTFILE;
  close SSPRFFILE;
  close PRFFILE;
}
 
