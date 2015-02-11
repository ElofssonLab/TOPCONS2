#!/usr/bin/perl -w

($file=$ARGV[0]) || die("Give blast output file name.");
open(IN,$file) || die("Couldn't open file $file");

$blast_block_length=60;
@aminoacids = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');

$id = $file; $id =~ s/^.*\///; $id =~ s/\..*//;
$query_id=$id;
#-----------------------------------------------------------------------------------------------------------

# Look for a sequence file
#if ((-e ($seq_file="$id.seq")) || (-e ($seq_file="$id.fasta"))) {
if ($seq_file=$ARGV[1]) {
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
	    /^(\S+)\s+(0\s+)?(-+)$/ || /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)$/ || die("Faulty line: '$_'");
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

$query_no_gaps=$sequences{'QUERY'}{'seq'};
$query_no_gaps =~ s/-//g;
$total_length=length($query_no_gaps);

print "Sequence: $query_id\n";
print "NR of aligned sequences: 100\n";
print "Length of query sequence: $total_length\n\n\n";
print "START 1\n";
print "ALPHABET:       ";
for (@aminoacids) {
    print "$_       ";
}
print "-     <SPACE> <LABEL> <QUERY>\n";
$i=0; $col=1;
$ali_length=length($sequences{'QUERY'}{'seq'});
for ($i=1;$i<=$ali_length;$i++) {
    $query_chr=substr($sequences{'QUERY'}{'seq'},$i-1,1);
    ($query_chr eq "-") && next;
    %aa_count=();
    $tot_count=0;
    for $id (sort {$sequences{$a}{'nr'} <=> $sequences{$b}{'nr'}} keys %sequences) {
	$chr=substr($sequences{$id}{'seq'},$i-1,1);
	unless($chr eq "-") {
	    $aa_count{$chr}++;
	    $tot_count++;
	}
    }
    print "COL" . spaces($col,5) . "$col:       ";
    for (@aminoacids) {
	if (exists($aa_count{$_})) {
	    $print=sprintf("%.2f",100*$aa_count{$_}/$tot_count);
	    print $print . spaces($print,8);
	} else {
	    $print="0.00 ";
	    print $print . spaces($print,8);
	}
    }
    $print="0.00 ";
    print $print . spaces($print,8);
    print $print . spaces($print,8);
    print "." . spaces(".",8);
    print "A" . spaces("X",8);
    
    print "\n";
    $col++;
}
print "END 1\n";
	

#-----------------------------------------------------------------------------------------------------------

sub spaces {
    ($nr,$total)=@_;
    $str="";
    for ($k=0;$k<$total-length($nr);$k++) {
	$str .= " ";
    }
    return $str;
}
