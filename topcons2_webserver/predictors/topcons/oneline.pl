#!/usr/bin/perl -w

# converts a multi line sequence in a fasta file into one line sequence
if(@ARGV!=1) { die "correct usage:\n\t perl oneline.pl <fasta_file> "};

open(FASTA,$ARGV[0]) || die "can't open fasta file";
$line="";
while($line !~ /^>/){
	$line=<FASTA>;
}
$prev="";
while (1){
	$line=<FASTA>;
	$line=~ s/[\r]//g; # remove carriage return
	while ($line !~ /^>/){
		chomp $line;
		$prev=$prev.$line;
		$prev =~ s/[\r]//g; # remove carriage return
		$line=<FASTA>;
		if(!(defined $line)) { print $prev."\n"; close(FASTA); exit(1);}
	}
	print $prev."\n";
	$line  =~ s/[\r]//g; # remove carriage return
	print $line;
	$prev="";
}


