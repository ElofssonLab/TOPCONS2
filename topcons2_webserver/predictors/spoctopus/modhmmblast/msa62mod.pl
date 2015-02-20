#!/usr/bin/perl -w

#use strict;

if($#ARGV != 4) {
  printf "Usage: msa2mod.pl <protnamefile> <round> <msafiledir>";
  printf " <modfiledir> <fastafiledir>\n";
  exit;
}

my $infile = $ARGV[0];
my $round = $ARGV[1];
my $indir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";
my $fastadir = $ARGV[4]."/";
my $default = ' ;' x 60;

open(INFILE,$infile) or die "Could not open $infile\n";
while(<INFILE>) {
  chomp;
  my @seqs = ();
  my $query_start = 0;
  my $query_end = 0;
  my $query = "";
  my $msaname = $_;
  my $msafile = "$indir"."$msaname".".msa";
  my $fasta_length;
  my $msa_length;
  my $cur_block = 0;

  if(-e "$fastadir"."$msaname".".seq") {
    open(FASTAFILE, "$fastadir"."$msaname".".seq")
      || die "could not open $fastadir"."$msaname".".seq\n";
  }
  elsif(-e "$fastadir"."$msaname".".fa") {
    open(FASTAFILE, "$fastadir"."$msaname".".fa")
      || die "could not open $fastadir"."$msaname".".fa\n";
  }
  else {
    die "could not open $fastadir"."$msaname"."\n";
  }

  <FASTAFILE>;
  while(<FASTAFILE>) {
    chomp;
    $query .= $_;
  }
  $fasta_length = length $query;
  close(FASTAFILE);
  
  if (!open(MSAFILE,$msafile)) {
    print "Could not open $msafile\n";
    next;
  }
  
  #check if file is ok
  my $OK = -1;
  while(<MSAFILE>) {
    if($_ =~ 'Query=') {
      $OK = 1;
    }
    if($_ =~ ' letters\)') {
      @_ = split;
      $msa_length = substr $_[0],1;
      $msa_length =~ s/,//g;
      last;
    }
  }
  if($OK < 0) {
    print "Error: msa_file $msafile is corrupt\n";
    next;
  }
  
  if($msa_length ne $fasta_length) {
    print "Error: non equal lengths of msa and fasta sequence\n";
    print "fasta_length = $fasta_length, msa_length = $msa_length \n";
    next;
  }
  
  #find beginning of msa info
  #print "$round\n";
  while(<MSAFILE>) {
    chomp;
    if($round > 0) {
      if($_ =~ "Results from round $round" || $_ =~ 'CONVERGED!') {
	#print "$_\n";
	last;
      }
      if($_ =~ 'No hits found') {
	print "No hits found for $msafile\n";
	$OK = -1;
	open(OUT,">".$outdir.$msaname.".mod");
	$query =~ s/(.)/$+;/g;
	printf OUT "<${query}>\n";
	last;
      }
    }
    else {
      if($_ =~ 'No hits found') {
	print "No hits found for $msafile\n";
	$OK = -1;
	open(OUT,">".$outdir.$msaname.".mod");
	$query =~ s/(.)/$+;/g;
	printf OUT "<${query}>\n";
	last;
      }
      elsif($_ =~ 'Sequences producing significant alignments:') {
	last;
      }
    }
  }

  if($OK < 0) {
    print "Error: strange format in msafile, could not find start of alignment info\n";
    next;
  }

  


  #find first sequence
  while(<MSAFILE>) {
    chomp;
    #print "$_";
    if($_ =~ 'QUERY' || $_ =~ '^1\_0' ) {
      my $id;
      my $begin;
      my $end;
      my $seq;
      
      ($id, $begin, $seq, $end) = split /\s+/, $_;
      $query_start = $begin;
      $query_end = $end;
      
      @seqs = (@seqs, $seq);
      last;
    }
  }

  my $seq_index = 1; 
  while(<MSAFILE>) {
    if($_ eq "\n") {
      $cur_block++;
    }
    elsif($_ =~ 'Database') {
      last;
    }
    else {
      chomp;

      my $id;
      my $begin;
      my $end;
      my $seq;
      my @list = split /\s+/, $_;
      if($#list == 3) {
	($id, $begin, $seq, $end) = @list;
      }
      elsif($#list == 1 && (substr($_,0,6) ne 'QUERY ' || substr($_,0,4) ne '1_0 ')) {
	($id, $seq) = @list;
      }
      elsif($#list == 2 && ( substr($_,0,6) ne 'QUERY ' ||substr($_,0,4) ne '1_0 ')) {
	
	($id, $begin, $seq) = @list;
	if(length $seq < 50) {
	  ($id, $seq, $end) = @list;
	}
	if(length $seq < 50) {
	  print "strange row in msafile: $_\n";
	  exit;
	}
      }
      else {
	print "strange row in msafile: $_\n";
	exit;
	  
      }
	    
      if(substr($_,0,6) eq 'QUERY '  || substr($_,0,4) eq '1_0 ' ) {
	$query_end = $end;
	$seqs[0] .= $seq;
	$seq_index = 1;
      }
      else {
	if($cur_block == 0) {
	  @seqs = (@seqs, $seq);
	}
	else {
	  if($seq_index > $#seqs) {
	    print "Error: msafile does not have the same number of aligned sequences in all blocks\n";
	    print "No output produced for $msaname\n";
	    exit;
	  }
	  #print "whole $_\n";
	  #print "id $id\n";
	  #print "begin $begin\n";
	  #print "end $end\n";
	  #print "seq $seq\n";
	  $seqs[$seq_index] .= $seq;
	  $seq_index++;
	}
      }
    }
  }

  my $outfile = "$outdir"."$msaname".".mod";
  open(OUT,">"."$outfile")
    or die "could not open $outfile\n";
  
  my $pos = $#seqs;
  my $pream = "";
  if($query_start > 1) {
    printf OUT "<";
    my $pre = substr($query,0,$query_start-1);
    #$pre =~ s/(.)/$+;/g;
    my @query_pre_list = split //, $pre;
    for(my $k = 0; $k <= $#query_pre_list; $k++) {
      printf OUT "$query_pre_list[$k]".";";
    }
    my @query_seq_list = split //, $seqs[0];
    for(my $k = 0; $k <= $#query_seq_list; $k++) {
      printf OUT "$query_seq_list[$k]".";";
    }
    $pream = " ;" x ($query_start-1);
  }
  else {
    my @query_seq_list = split //, $seqs[0];
    printf OUT "<";
    for(my $k = 0; $k <= $#query_seq_list; $k++) {
      printf OUT "$query_seq_list[$k]".";";
    }
  }
  if($query_end < length($query)) {
    my $rest = substr($query,$query_end);
    #$rest =~ s/(.)/$+;/g;
    my @query_end_list = split //, $rest;
    for(my $k = 0; $k <= $#query_end_list; $k++) {
      printf OUT "$query_end_list[$k]".";";
    }
  }
  printf OUT ">\n";
  for(my $j = 1; $j <= $pos; $j++) {
    my @seq_list = split //, $seqs[$j];
    for(my $k = 0; $k <= $#seq_list; $k++) {
      if($seq_list[$k] ne '-') {
	last;
      }
      else {
	$seq_list[$k] = ' ';
      }
    }
    for(my $k = $#seq_list; $k >= 0; $k--) {
      if($seq_list[$k] ne '-') {
	last;
      }
      else {
	$seq_list[$k] = ' ';
      }
    }
    printf OUT "<${pream}";
    for(my $k = 0; $k <= $#seq_list; $k++) {
      printf OUT "$seq_list[$k]".";";
    }
    printf OUT ">\n";
  }

  close(OUT);
  print "Done modding $msaname\n";
}
