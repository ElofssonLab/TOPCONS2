#!/usr/bin/perl -w

($file=$ARGV[0]) || die("Give blast output file name.");
open(IN,$file) || die("Couldn't open file $file");

$blast_block_length=60;

$id = $file; $id =~ s/^.*\///; $id =~ s/\..*//;
$pathID = $file; $pathID =~ s/\..*//;
#-----------------------------------------------------------------------------------------------------------

# Look for a sequence file
if ((-e ($seq_file="$id.seq")) || (-e ($seq_file="$id.fasta")) || (-e ($seq_file="$pathID.fa"))) {
    $ref_seq=`tail -n 1 $seq_file`; chomp($ref_seq);
}

# Find out length of query sequence
while(<IN>) { if (/^\s+\((\d+)\sletters\)$/) { $query_length=$1; last; } }

# Find beginning of MSA
while(<IN>) { last if /^Searching\.+done$/; }
while(<IN>) { last unless /^((CONVERGED!)|(\s*))$/; }
seek IN, (-1*length), 1;

# Parse MSA
<IN>; seek IN, (-1*length), 1;
if (/^\s*\*+\sNo\shits\sfound\s\*+$/) {
    if ($ref_seq) {                                                  # If there is a sequence file...
	$sequences{'QUERY'}{'seq'}=$ref_seq;
	$sequences{'QUERY'}{'nr'}=1;
	$sequences{'2'}{'seq'}=$ref_seq;
	$sequences{'2'}{'nr'}=2;
    }
} else {
    if ((/^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)$/) && ($2 != 1)) {
	$align_begin=$2;
    }
    $block=1; $nr=0;
    %block_hits=();               # Keep track of which id:s are already represented in each block
    while(<IN>) {
	last if /^\s+Database:/;
	unless (/^$/) {
	    /^(\S+)\s+(0\s+)?(-+)$/ || /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)$/ || die("Faulty line '$_'");
	    $id=$1; $seq=$3; #$ul=$4; $ll=$2;
	    if ($block_hits{$id}{$block}) {
		$j=0;
		while($block_hits{$id."###".++$j}{$block}) {}
		$id = $id."###".$j;
	    }
	    if (exists($sequences{$id})) {
		$sequences{$id}{'seq'} .= $seq;
	    } else {
		$sequences{$id}{'seq'} .= "-" x ($blast_block_length*($block-1)) . $seq;
		$sequences{$id}{'nr'} = ++$nr;
	    }
	    $block_hits{$id}{$block}=1;
	} else {
	    $block++;
	}
    }
    close(IN);
}

if ($align_begin && $ref_seq) {             # If alignment doesn't begin at position 1 and there is a sequence file
    for $id (sort {$sequences{$a}{'nr'} <=> $sequences{$b}{'nr'}} keys %sequences) {
	if ($sequences{$id}{'nr'} == 1) {
	    $sequences{$id}{'seq'} = substr($ref_seq,0,$align_begin-1) . $sequences{$id}{'seq'};
	} else {
	    $sequences{$id}{'seq'} = "-" x ($align_begin-1) . $sequences{$id}{'seq'};
	}
    }
}    
$query_sequence=$sequences{'QUERY'}{'seq'};
$query_sequence =~ s/-//g;
if ((length($query_sequence) != $query_length) && $ref_seq) { # If !(aln. ends at last pos.) && exists(seq. file)
    for $id (sort {$sequences{$a}{'nr'} <=> $sequences{$b}{'nr'}} keys %sequences) {
	if ($sequences{$id}{'nr'} == 1) {
	    $sequences{$id}{'seq'} = $sequences{$id}{'seq'} . substr($ref_seq,length($query_sequence),($query_length-length($query_sequence)));
	} else {
	    $sequences{$id}{'seq'} = $sequences{$id}{'seq'} . "-" x ($query_length-length($query_sequence));
	}
    }
}    

for $id (sort {$sequences{$a}{'nr'} <=> $sequences{$b}{'nr'}} keys %sequences) {
    $sequences{$id}{'seq'} =~ s/-//g;
    print ">$sequences{$id}{'nr'}\n$sequences{$id}{'seq'}\n";
}

#-----------------------------------------------------------------------------------------------------------

sub spaces {
    ($nr)=@_;
    $str="";
    for ($k=0;$k<8-length($nr);$k++) {
	$str .= " ";
    }
    return $str;
}
