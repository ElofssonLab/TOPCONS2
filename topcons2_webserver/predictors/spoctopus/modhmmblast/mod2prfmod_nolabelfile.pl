#!/usr/bin/perl -w


if($#ARGV != 3) {
  printf "Usage: mod2prfmod.pl <protnamefile> <fastafiledir> <modfiledir>";
  printf " <outdir>\n";
  exit;
}

my $protnamefile = $ARGV[0];
my $fastafiledir = $ARGV[1]."/";
my $modfiledir = $ARGV[2]."/";
my $outdir = $ARGV[3]."/";

my @alphabet = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
		'B', 'U', 'X', 'Z');
my @init_row = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);


my $init_info = "";
for(my $l = 0; $l <= 19; $l++) {
  $init_info .= "$alphabet[$l]"." " x 7;
}
$init_info .= "-"." " x 5;
$init_info .= "<SPACE> ";
$init_info .= "<LABEL> ";
$init_info .= "<QUERY> \n";



open(PROTNAMEFILE, "$protnamefile")
  || die "could not open $protnamefile\n";


while(<PROTNAMEFILE>) {
  chomp;
  my $seqname = $_;
  my $seq;
  my $top;
  my @seq_list;
  my @top_list;
  my $row;
  my @row_list;

  my $nr_alignseqs = 0;
  my $seq_length = 0;

  my @profile_info;

  if(-e "$fastafiledir"."$seqname".".seq") {
    open(FASTAFILE, "$fastafiledir"."$seqname".".seq")
      || die "could not open $fastafiledir"."$seqname".".seq\n";
  }
  elsif(-e "$fastafiledir"."$seqname".".fa") {
    open(FASTAFILE, "$fastafiledir"."$seqname".".fa")
      || die "could not open $fastafiledir"."$seqname".".fa\n";
  }
  else {
    die "could not open $fastafiledir"."$seqname".".fa\n";
  }

  
  <FASTAFILE>;
  while(<FASTAFILE>) {
    chomp;
    if($_ =~ '>') {
    }
    else {
      $seq .= $_;
    }
  }
  @seq_list = split //, $seq;
  $seq_length = $#seq_list + 1;
  close FASTAFILE;

  

  $top = "." x $seq_length;
  @top_list = split //, $top;

  open(MODFILE, "$modfiledir"."$seqname".".mod")
  || die "could not open $modfiledir"."$seqname".".mod\n";

  $row = readline MODFILE;
  chomp $row;
  @row_list = split /<|>|;/ , $row;

  my $seq_index = 0;
  for(my $i = 1; $i <= $#row_list; $i++) {
    #print "$row_list[$i]";
    push @profile_info, [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    if($row_list[$i] eq ' ') {
      $profile_info[$i-1][21] += 1;
      $profile_info[$i-1][23] = '-';
      $profile_info[$i-1][22] = '-';
    }
    elsif($row_list[$i] eq '-') {
      $profile_info[$i-1][20] += 1;
      $profile_info[$i-1][22] = '-';
      $profile_info[$i-1][23] = '-';
    }
    else {
      if($#top_list >= $seq_index) {
	$profile_info[$i-1][22] = $top_list[$seq_index];
      }
      else {
	$profile_info[$i-1][22] = 'd';
      }
      $profile_info[$i-1][23] = $seq_list[$seq_index];
      $seq_index++;
      my $a_index;
      for($a_index = 0; $a_index <= $#alphabet; $a_index++) {
	if($row_list[$i] eq $alphabet[$a_index]) {
	  last;
	}
      }
      if($a_index >= 24) {
	die "error: letter $row_list[$i] not in sequence\n";
      }
      elsif($a_index == 23) {
	#Z
	$profile_info[$i-1][3] += 0.5;
	$profile_info[$i-1][13] += 0.5;
      }
      elsif($a_index == 22) {
	#X	
	for(my $j = 0; $j < 20; $j++) {
	  $profile_info[$i-1][$j] += 0.05;
	}
      }
      elsif($a_index == 21) {
	#U
	$profile_info[$i-1][10] += 1;
      }
      elsif($a_index == 20) {
	#B
	$profile_info[$i-1][2] += 0.5;
	$profile_info[$i-1][11] += 0.5;
      }
      else {
	$profile_info[$i-1][$a_index] += 1;
      }
    }
  }

  while(<MODFILE>) {
    $nr_alignseqs++;
    $row = $_;
    chomp $row;
    if($row =~ '/') {
      $row = substr $row, 0, 1;
      $row = substr $row, (length($row) - 1);
      @top_list = split //, $row;
    }
    else {
      @row_list = split /<|>|;/ , $row;
      for(my $i = 1; $i <= $#row_list; $i++) {
	if($row_list[$i] eq ' ') {
	  $profile_info[$i-1][21] += 1;	
	}
	elsif($row_list[$i] eq '-') {
	  $profile_info[$i-1][20] += 1;
	}
	else {
	  my $a_index;
	  for($a_index = 0; $a_index <= $#alphabet; $a_index++) {
	    if($row_list[$i] eq $alphabet[$a_index]) {
	      last;
	    }
	  }
	  if($a_index >= 24) {
	    #die "error: letter $row_list[$i] not in sequence\n";
	  }
	  elsif($a_index == 23) {
	    #Z
	    $profile_info[$i-1][3] += 0.5;
	    $profile_info[$i-1][13] += 0.5;
	  }
	  elsif($a_index == 22) {
	    #X	
	    for(my $j = 0; $j < 20; $j++) {
	      $profile_info[$i-1][$j] += 0.05;
	    }
	  }
	  elsif($a_index == 21) {
	    #U
	    $profile_info[$i-1][10] += 1;
	  }
	  elsif($a_index == 20) {
	    #B
	    $profile_info[$i-1][2] += 0.5;
	    $profile_info[$i-1][11] += 0.5;
	  }
	  else {
	    $profile_info[$i-1][$a_index] += 1;
	  }
	}
      }
    }
  }
  close MODFILE;

  my $a;
  my $b;

  open(OUTFILE, ">"."$outdir"."$seqname".".prf")
    or die "could not open "."$outdir"."$seqname".".prf"."\n";

  print OUTFILE "Sequence: $seqname \n";
  print OUTFILE "NR of aligned sequences: $nr_alignseqs \n";
  print OUTFILE "Length of query sequence: $seq_length \n\n";
  print OUTFILE "START 1\n";

  print OUTFILE "ALPHABET:   $init_info";
  for $a (0 .. $#profile_info) {
    printf OUTFILE "COL %4d", $a+1;
    print OUTFILE ":  ";
    for $b (0 .. $#{$profile_info[$a]} ) {
      #print "b = $b\n";
      #print "profile_info = $#{$profile_info[$a]}\n";
      if($b >= ($#{$profile_info[$a]} - 1)) {
	print OUTFILE "$profile_info[$a][$b]       ";
	#print("a = $a och b = $b\n"); 
      }
      else {
	printf OUTFILE "%5.02f   ", $profile_info[$a][$b];
      }
  }
    print OUTFILE "\n";
  }
  print OUTFILE "END 1\n";

  close OUTFILE;
  print "done profiling $seqname \n";
}

