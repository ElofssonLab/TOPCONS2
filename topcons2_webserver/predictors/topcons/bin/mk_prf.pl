#!/usr/bin/perl -w

use strict;

#-----------------------------------------------------------------------------------------------------------

my @alphabet=('i','M','o','S','p');

#-----------------------------------------------------------------------------------------------------------

@ARGV || die("Syntax: $0 seqname topologydir1 topologydir2 ...\n");
my %profile=();
my $length=0;
my $name="";
my $alphabet_str=join("",@alphabet);

# Read topologies
for (my $i = 1; $i <= $#ARGV; $i++) {
    my $file = $ARGV[$i]."/"."$ARGV[0]".".top";
    #open(IN,$file) || die("Could not open file $file\n");
    #while(my $id=<IN>) {
#	my $top="";
#	while(my $row = <IN>) {
	    #if ($row =~ /^>/) {
#		seek IN, -1*length($row), 1;
#		last;
	    #}
	#    $top .= $row;
   	#}
   my $top;
    open(my $fh, '<', $file) or die "cannot open file $file";
    {
        $top = <$fh>;
    }
    close($fh);
 
	$top =~ s/R/o/g;
	$top =~ s/r/i/g;
	$top =~ s/O/o/g;
	$top =~ s/I/i/g;
	#$id =~ s/^>//;
	#chomp($id,$top);
	$top =~ s/\s//g;
	((length($top) != $length) && ($length != 0)) && die("Topologies have different lengths, file $file\n");
	$length=length($top);
	#($id ne $name) && ($name ne "") && die("Topologies have different names, file $file\n");
	$name="query";#$id;
	($top =~ /^[$alphabet_str]+/) || die("Topology contains unsupported characters, file $file\n");
	my $pos=1;
	my $chr;
	$profile{$pos++}{$chr}++ while ($chr = substr($top,0,1,""));
    #}
    #close(IN);
}

# Normalize
for my $pos (sort {$a <=> $b} keys %profile) {
    my $sum=0;
    for (@alphabet) {
	exists($profile{$pos}{$_}) || ($profile{$pos}{$_}=0);
	$sum += $profile{$pos}{$_};
    }
    $profile{$pos}{$_} = sprintf("%.2f",100*$profile{$pos}{$_}/$sum) for (@alphabet);
}

# Print out profile
print "Sequence: $name\n";
print "NR of aligned sequences: 100\n";
print "Length of query sequence: $length\n\n";
print "START 1\n";
print "ALPHABET:   ";
for (@alphabet) {
    print "$_       ";
}
print "-     <SPACE> <LABEL> <QUERY>\n";
	
for my $pos (sort {$a <=> $b} keys %profile) {
    print "COL" . spaces($pos,5) . "$pos:   ";
    for (@alphabet) {
	print $profile{$pos}{$_} . spaces($profile{$pos}{$_},8);
    }
    print "0.00    0.00    0.00    0.00    .\n";
}

print "END 1\n";


#-----------------------------------------------------------------------------------------------------------

sub spaces {
    my ($nr,$total)=@_;
    my $str="";
    for (my $k=0;$k<$total-length($nr);$k++) {
	$str .= " ";
    }
    return $str;
}

