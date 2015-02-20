#!/usr/bin/perl -w
use Math::Trig;

if($#ARGV != 2) {
  printf "Usage: score_net.pl <netfile> <datafile> <outfile>\n";
  exit;
}

my $netfile = $ARGV[0];
my $datafile = $ARGV[1];
my $outfile = $ARGV[2];


#read net
my $nin = 0;
my $nhidden = 0;
my $nout = 0;
my $nwts = 0;
my @w1 =();
my @w2 = ();
my @b1 =();
my @b2 = ();


#make sure outfile is empty
open(OUTFILE, ">"."$outfile")
    or die "could not open $outfile";

open(NETFILE, "$netfile")
    or die "could not open $netfile";
my $reading = "";
print "Reading net ...\n";
while(<NETFILE>) {
    chomp;
    $_ =~ s/^\s+//g;
    $_ =~ s/\s+/ /g;
    if($_ =~ 'nin') {
	$reading = "nin";
    }
    elsif($_ =~ 'nhidden') {
	$reading = "nhidden";	
    }
    elsif($_ =~ 'nout') {
	$reading = "nout";
    }
    elsif($_ =~ 'nwts') {
	$reading = "nwts";
    }
    elsif($_ =~ 'w1') {
	$reading = "w1";
    }
    elsif($_ =~ 'w2') {
	$reading = "w2";
    }
    elsif($_ =~ 'b1') {
	$reading = "b1";
    }
    elsif($_ =~ 'b2') {
	$reading = "b2";
    }
    elsif($reading eq 'nin') {
	$nin = $_;
	print "nin = $nin\n";;
    }
    elsif($reading eq 'nhidden') {
	$nhidden = $_;
	print "nhidden = $nhidden\n";;
    }
    elsif($reading eq 'nout') {
	$nout = $_;
	print "nout = $nout\n";;
    }
    elsif($reading eq 'nwts') {
	$nwts = $_;
	print "nwts = $nwts\n";
    }
    elsif($reading eq 'w1') {
	my @weights_row = split /\s/, $_;
	push @w1, [@weights_row];
    }
    elsif($reading eq 'w2') {
	my @weights_row = split /\s/, $_;
	push @w2, [@weights_row];
    }
    elsif($reading eq 'b1') {
	my @bias_row = split /\s/, $_;
	push @b1, @bias_row;
    }
    elsif($reading eq 'b2') {
	my @bias_row = split /\s/, $_;
	push @b2, @bias_row;
    }
}
#print "$#w1  $#w2  $#b1  $#b2\n";
#print "@b1\n";

print "done\n";

#read data
print "scoring net ...\n";
open(DATAFILE, "$datafile")
    or die "could not open $datafile";
while(<DATAFILE>) {
    chomp;
    $_ =~ s/\s+/ /g;
    my @dataline = split /\s/, $_;
    
    #score net
    my @hidden_scores = ();
    for(my $i = 0; $i < $nhidden; $i++) {
	$hidden_score = $b1[$i];
	for(my $j = 0; $j <= $#w1; $j++) {
	    $hidden_score += $w1[$j][$i] * $dataline[$j];
	}
	$hidden_score = tanh_hidden($hidden_score);
	$hidden_scores[$i] = $hidden_score;
	#print "$hidden_score\n";
	
    }
    
    my @outscores = ();
    for(my $i = 0; $i < $nout; $i++) {
	my $outscore =  $b2[$i];
	for(my $j = 0; $j <= $#hidden_scores; $j++) {
	    $outscore += $w2[$j][$i] * $hidden_scores[$j]; 
	}
	$outscore = logistic_out($outscore);
	#$outscore = tanh_hidden($outscore);
	$outscores[$i] = $outscore;
    }
    #print "@outscores\n";
    
    #print data
    #open(OUTFILE, ">>"."$outfile")
    #or die "could not open $outfile";
    for(my $i = 0; $i <= $#outscores; $i++) {
	print OUTFILE "$outscores[$i] ";
    }
    print OUTFILE "\n";
    
    #close OUTFILE;
}

close OUTFILE;
print "done scoring\n";



sub logistic_hidden {
    my $in = shift;
    my $beta = 5.0;
    my $e = 2.7182818;
    my $result = 1 / (1 + $e ** (0.0 - ($beta * $in)));

    return $result;
}
sub tanh_hidden {
    my $in = shift;
    my $result = tanh($in);

    return $result;
}


sub logistic_out {
    my $in = shift;
    my $beta = 1;
    my $e = 2.7182818;
    my $result = 1 / (1 + $e ** (0.0 - ($beta * $in)));
    
}
