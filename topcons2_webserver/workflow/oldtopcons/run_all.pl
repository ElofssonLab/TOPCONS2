#!/usr/bin/perl -w

use strict;
use POSIX qw(:signal_h :errno_h :sys_wait_h);

####################################################
# Parameters:
my $pwd="/bin/pwd";
my $maxtime=36000; # 10 min.
####################################################

my @alphabet=('i','M','o');
$ENV{'PATH'}="/bin:/usr/bin";

$ARGV[0] =~ /^([\w\.\-\/]+)$/ || die("Faulty format for tmpdir: '$ARGV[0]'");

my $tmpdir=$1;

system("/bin/echo Started run_all.pl > $tmpdir/log");

# Run jobs serially! If another topcons process is running, this process will wait until the file "occupied" is released by that process. This works as a que system for topcons processes!

# open(OCCUPIED, ">/server/var/www/topcons.cbr.su.se/lock/occupied") || die("Could not lock file");
# flock(OCCUPIED, 2); # Lock file for reading by anyone else

system("/bin/echo BLASTing > $tmpdir/log");

print "BLASTing...\n";
system("/bin/sh fa2prfs.sh $tmpdir/query");
print "BLAST finished!\n";

system("/bin/echo Finished BLASTing > $tmpdir/log");

system("/bin/echo $tmpdir/query.raw.prf > $tmpdir/query.raw.prf.snf");
system("/bin/echo $tmpdir/query.fa > $tmpdir/query.fa.snf");
system("/bin/echo query > $tmpdir/query.pnf");

my $fix_topo="";
if(open(IN,"$tmpdir/query.fix")) {
    $fix_topo = <IN>;
    close(IN);
    $fix_topo =~ /^[iMo\.]+$/ || die("Failed in fix topo");
    relabel_prf_file("$tmpdir/query.prf",$fix_topo);
    relabel_prf_file("$tmpdir/query.raw.prf",$fix_topo);
    relabel_fa_file("$tmpdir/query.fa",$fix_topo);
}
    

# Run in parallel on 4 processors! The process forks into four children, keeping the parent process
#print "0x" . sprintf("%08x",($$*4)) . "\n";

my $pid1 = open(CHILD1,"-|");
if ($pid1 == 0) { # We are in CHILD1 process
    open(OUT1, ">$tmpdir/child1.out");
    print OUT1 "SCAMPI_multi: ";
    print OUT1 SCAMPI_multi_run($tmpdir,"query.raw.prf",$fix_topo);
    close(OUT1);
    res2topofile("$tmpdir/DGHMM_KR_21_multi.hmg.res","$tmpdir/query.top");
    system("./DGpred_multi_scan $tmpdir/query.raw.prf 21 > $tmpdir/query.DG");
    ZPRED_run($tmpdir, "$tmpdir/query.psi","$tmpdir/query.top","$tmpdir/query.zpred");
    exit(149); # CHILD1 process is finished, exit here.
} 
my $pid2 = open(CHILD2,"-|");
if ($pid2 == 0) { # We are in CHILD2 process 
    open(OUT2, ">$tmpdir/child2.out");
    print OUT2 "PRODIV: ";
    print OUT2 PRODIV_run($tmpdir,"query.raw.prf");
    close(OUT2);
    exit(149); # CHILD2 process is finished, exit here.
} 
my $pid3 = open(CHILD3,"-|");
if ($pid3 == 0) { # We are in CHILD3 process
    open(OUT3, ">$tmpdir/child3.out");
    print OUT3 "PRO: ";
    print OUT3 PRO_run($tmpdir,"query.raw.prf");
    print OUT3 "SCAMPI_single: ";
    print OUT3 SCAMPI_single_run($tmpdir,"query.fa",$fix_topo);
    close(OUT3);
    exit(149); # CHILD3 process is finished, exit here.
} 
my $pid4 = open(CHILD4,"-|");
if ($pid4 == 0) { # We are in CHILD4 process
    open(OUT4, ">$tmpdir/child4.out");
    print OUT4 "OCTOPUS: ";
    print OUT4 OCTOPUS_run($tmpdir,$pwd);
    close(OUT4);
    exit(149); # CHILD4 process is finished, exit here.
}

# From here on, we are in PARENT process.
my $nr_children_finished=0;
my @all_topo=();
my $finished_CHILD1 = 0;
my $finished_CHILD2 = 0;
my $finished_CHILD3 = 0;
my $finished_CHILD4 = 0;
print "\nRESULTS...\n";
while ( $nr_children_finished < 4 && (times)[0] < $maxtime ) {

    if ( $finished_CHILD1 == 0 && my_wait($pid1) ) {
	$finished_CHILD1 = 1;
	open(RES1, "$tmpdir/child1.out");
	while ( my $line = <RES1> ) {
	    $line =~ /^\S+\s(\S+)$/;
	    push(@all_topo,$1);
	}
	close(RES1);
	$nr_children_finished++;
    }

    if ( $finished_CHILD2 == 0 && my_wait($pid2) ) {
	$finished_CHILD2 = 1;
	open(RES2, "$tmpdir/child2.out");
	while ( my $line = <RES2> ) {
	    $line =~ /^\S+\s(\S+)$/;
	    push(@all_topo,$1);
	}
	close(RES2);
	$nr_children_finished++;
    }

    if ( $finished_CHILD3 == 0 && my_wait($pid3) ) {
	$finished_CHILD3 = 1;
	open(RES3, "$tmpdir/child3.out");
	while ( my $line = <RES3> ) {
	    $line =~ /^\S+\s(\S+)$/;
	    push(@all_topo,$1);
	}
	close(RES3);
	$nr_children_finished++;
    }

    if ( $finished_CHILD4 == 0 && my_wait($pid4) ) {
	$finished_CHILD4 = 1;
	open(RES4, "$tmpdir/child4.out");
	while ( my $line = <RES4> ) {
	    $line =~ /^\S+\s(\S+)$/;
	    push(@all_topo,$1);
	}
	close(RES4);
	$nr_children_finished++;
    }

    # Sleep for 0.5 seconds
    select undef, undef, undef, 0.5; 
}

# print STDERR "Topos:\n";
# foreach my $top ( @all_topo ) {
#     print STDERR $top."\n";
# }

close(CHILD1);
close(CHILD2);
close(CHILD3);
close(CHILD4);

open(OUT5, ">$tmpdir/master.out") || die("Couldn't open $tmpdir/master.out");
print OUT5 TOPCONS_run(@all_topo); # Run the TOPCONS HMM.
close(OUT5);

# close(OCCUPIED);
exit;


#------------------------------------------------------------------------------------

sub my_wait {
    my ($waitpid) = @_;
    
    my $pid = waitpid($waitpid, &WNOHANG);
    my $ret = $?;
    
    if ($pid == -1) {
        # no child waiting
    } elsif (WIFEXITED($ret)) {
	#   print "Process $pid exited.\n";
	
	my $exit_value = $ret >> 8;
	
	if ($exit_value == 149) {  # Is CHILD finished?
	    return 1;
	}
	else { return 0; }
    }
    return 0;
}

sub spaces {
    my ($nr,$total)=@_;
    my $str="";
    for (my $k=0;$k<$total-length($nr);$k++) {
	$str .= " ";
    }
    return $str;
}

sub PRO_run {
    my ($tmpdir,$infile) = @_;
    system("/server/modhmm/bin/modhmms -m PRODIV/PRO_TMHMM_0.91.txt -s $tmpdir/$infile.snf -f prf -o $tmpdir -M GM -r PRODIV/replacement_letter_multi.rpl -L --nopostout --viterbi > $tmpdir/PRO_TMHMM_0.91.hmg.res.xml");
    system("/bin/cat  $tmpdir/PRO_TMHMM_0.91.hmg.res.xml | /server/modhmm/bin/modhmmxml2res > $tmpdir/PRO_TMHMM_0.91.hmg.res");
    return res2topo("$tmpdir/PRO_TMHMM_0.91.hmg.res");
}

sub PRODIV_run {
    my ($tmpdir,$infile) = @_;
    system("/server/modhmm/bin/modhmms -m PRODIV/PRODIV_TMHMM_0.91.txt -s $tmpdir/$infile.snf -f prf -o $tmpdir -M GM -r PRODIV/replacement_letter_multi.rpl -L --nopostout --viterbi --max_d > $tmpdir/PRODIV_TMHMM_0.91.hmg.res.xml");
    system("/bin/cat  $tmpdir/PRODIV_TMHMM_0.91.hmg.res.xml | /server/modhmm/bin/modhmmxml2res > $tmpdir/PRODIV_TMHMM_0.91.hmg.res");
    return res2topo("$tmpdir/PRODIV_TMHMM_0.91.hmg.res");
}

sub SCAMPI_multi_run {
    my ($tmpdir,$infile,$fix_topo) = @_;
    if ($fix_topo =~ /[iMo]/) {
	system("/server/modhmm/bin/modhmms_scampi -f prf -s $tmpdir/$infile.snf -m SCAMPI/DGHMM_KR_21_multi.txt -o $tmpdir -r SCAMPI/replacement_letter_multi.rpl --nopostout --viterbi -u -L  > $tmpdir/DGHMM_KR_21_multi.hmg.res.xml");
    } else {
	system("/server/modhmm/bin/modhmms_scampi -f prf -s $tmpdir/$infile.snf -m SCAMPI/DGHMM_KR_21_multi.txt -o $tmpdir -r SCAMPI/replacement_letter_multi.rpl --nopostout --viterbi -u -L -g >  $tmpdir/DGHMM_KR_21_multi.hmg.res.xml");
    }
    system("/bin/cat  $tmpdir/DGHMM_KR_21_multi.hmg.res.xml | /server/modhmm/bin/modhmmxml2res >  $tmpdir/DGHMM_KR_21_multi.hmg.res");
    return res2topo("$tmpdir/DGHMM_KR_21_multi.hmg.res");
}

sub SCAMPI_single_run {
    my ($tmpdir,$infile,$fix_topo) = @_;
    if ($fix_topo =~ /[iMo]/) {
	system("/server/modhmm/bin/modhmms_scampi -f fa -s $tmpdir/$infile.snf -m SCAMPI/DGHMM_KR_21_single.txt -o $tmpdir -r SCAMPI/replacement_letter_multi.rpl --nopostout --viterbi -u -L > $tmpdir/DGHMM_KR_21_single.hmg.res.xml");
    } else {
	system("/server/modhmm/bin/modhmms_scampi -f fa -s $tmpdir/$infile.snf -m SCAMPI/DGHMM_KR_21_single.txt -o $tmpdir -r SCAMPI/replacement_letter_multi.rpl --nopostout --nolabels --viterbi -u -L -g > $tmpdir/DGHMM_KR_21_single.hmg.res.xml");
    }
    system("/bin/cat $tmpdir/DGHMM_KR_21_single.hmg.res.xml | /server/modhmm/bin/modhmmxml2res > $tmpdir/DGHMM_KR_21_single.hmg.res ");
    return res2topo("$tmpdir/DGHMM_KR_21_single.hmg.res");
}

sub ZPRED_run {
    my ($tmpdir, $prof_file,$top_file,$out_file) = @_;
    system("zpred2/scripts/zpred.pl -profile $prof_file -topology $top_file -out $out_file > $tmpdir/zpred.error_log 2>&1");
}

sub OCTOPUS_run {
    my ($tmpdir,$pwd) = @_;
    my $topo = "";
    `$pwd` =~ /^([\w\-\/\.]*)$/ || die;
    my $octopusdir="$1/OCTOPUS";
    system("/bin/sh OCTOPUS/bin/run_OCTOPUS_no_matlab.sh $tmpdir/query.pnf $tmpdir $tmpdir $tmpdir $tmpdir $octopusdir > $tmpdir/octopus.log 2>&1");
    open(IN,"$tmpdir/query.top") || die("Couldn't open $tmpdir/query.top");
    while(<IN>) {
	chomp;
	/^>(.*)/ || ($topo .= $_);
    }
    close(IN);
    return "$topo\n";
}

sub TOPCONS_run {
    my @all_topo = (@_);
    my %profile=();
    for my $top (@all_topo) {
	my $pos=1;
	my $chr;
	$profile{$pos++}{$chr}++ while ($chr = substr($top,0,1,""));
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
    my $length = scalar(keys(%profile));
    open(OUT,">$tmpdir/pred_topo.prf") || die;
    print OUT "Sequence: query\n";
    print OUT "NR of aligned sequences: 100\n";
    print OUT "Length of query sequence: $length\n\n";
    print OUT "START 1\n";
    print OUT "ALPHABET:   ";
    for (@alphabet) {
	print OUT "$_       ";
    }
    print OUT "-     <SPACE> <LABEL> <QUERY>\n";
    
    for my $pos (sort {$a <=> $b} keys %profile) {
	print OUT "COL" . spaces($pos,5) . "$pos:   ";
	for (@alphabet) {
	    print OUT $profile{$pos}{$_} . spaces($profile{$pos}{$_},8);
	}
	print OUT "0.00    0.00    0.00    .\n";
    }
    
    print OUT "END 1\n";
    close(OUT);
    
    # Run TOPCONS
    `/bin/echo $tmpdir/pred_topo.prf > $tmpdir/pred_topo.prf.snf`;
    system("/server/modhmm/bin/modhmms -m TOPCONS/topcons.txt -s $tmpdir/pred_topo.prf.snf -f prf -o $tmpdir -L -M DP --nolabels --viterbi --nopostout > $tmpdir/topcons.hmg.res.xml");
    system("/bin/cat  $tmpdir/topcons.hmg.res.xml | /server/modhmm/bin/modhmmxml2res > $tmpdir/topcons.hmg.res");
    my $topo = res2topo("$tmpdir/topcons.hmg.res");
    chomp($topo);

    # Make reliability curve
    my $pos=1;
    my @support=();
    while(my $chr=substr($topo,$pos-1,1)) {
	push(@support,$profile{$pos++}{$chr}/100);
    }
    open(OUT,">$tmpdir/query.reliability") || die;
    for (my $i=10;$i<scalar(@support)-10;$i++) {
	my $sum=0;
	my $count=0;
	for (my $k=$i-10;$k<=$i+10;$k++) {
	    $sum+=$support[$k];
	    $count++;
	}
	print OUT "" . ($i+1) . " " . sprintf("%.3f",($sum/$count)) . "\n";
    }
    close(OUT);
    
    return "TOPCONS: $topo\n";
}


sub res2topo {
    my ($infile) = @_;
    open(IN,$infile) || die("Could not open file $infile\n");
    <IN>;<IN>;
    my $new=1; my $topo=""; my $TM=-1; my $seq_file=""; my $seq_length="";
    while(<IN>) {
	chomp;
	if ($new) {
	    next if (/^$/);
	    $seq_file=$_;
	    $new=0;
	} elsif (/^Seq\slength:\s(\d+)/) {
	    $seq_length=$1;
	} elsif (/^((Is)|(No))\sTM\sprotein/) {
	    if ($1 eq "Is") {
		$TM=1;
	    } elsif ($1 eq "No") {
		$TM=0;
	    }
	} elsif (/^Labeling:/) {
	} elsif (/^$/) {	
	    $new=1;
	    $seq_file="";
	    $seq_length="";
	    $TM="";
	} else {
	    $topo .= $_;
	}
    }
    close(IN);
    $topo =~ tr/[A-Z]/[a-z]/;
    $topo =~ s/e/m/g;
    $topo =~ s/l/m/g;
    $topo =~ s/j/i/g;
    $topo =~ s/m/M/g;
    return "$topo\n";
}


sub res2topofile {
    my ($infile,$outfile)=@_;
    open(IN,$infile) || die("Could not open file $infile\n");
    <IN>;<IN>;
    my $new=1; my $topo=""; my $TM=-1; my $seq_file=""; my $seq_length="";
    while(<IN>) {
	chomp;
	if ($new) {
	    next if (/^$/);
	    $seq_file=$_;
	    $new=0;
	} elsif (/^Seq\slength:\s(\d+)/) {
	    $seq_length=$1;
	} elsif (/^((Is)|(No))\sTM\sprotein/) {
	    if ($1 eq "Is") {
		$TM=1;
	    } elsif ($1 eq "No") {
		$TM=0;
	    }
	} elsif (/^Labeling:/) {
	} elsif (/^$/) {
	    print_topofile($seq_file, $seq_length, $TM, $topo, $outfile);
	    $new=1;
	    $seq_file="";
	    $seq_length="";
	    $TM="";
	    $topo="";
	} else {
	    $topo .= $_;
	}
    }
    close(IN);
}

sub print_topofile {
    my ($seq_file, $seq_length, $TM, $topo, $outfile) = @_;
    open(OUT,">$outfile") || die("Could not write to file");
    ($seq_length == length($topo)) || die("$seq_file $seq_length " . $topo . " ");
    $topo =~ tr/[A-Z]/[a-z]/;
    $topo =~ s/e/m/g;
    $topo =~ s/l/m/g;
    $topo =~ s/j/i/g;
    $topo =~ s/m/M/g;
    print OUT ">query\n$topo\n";
    close(OUT);
}

sub relabel_prf_file {
    my $prffile = shift;
    my $topology = shift;
    if($prffile eq "" || $topology eq "") {
	return;
    }
    else {
	my @topology = split //, $topology;
	
	
	open(PRFFILE, "$prffile")
	    or die "Could not open $prffile\n";
	
	my $intro = "";
	my @cols = ();
	my $extro = "";
	my $inintro = 1;
	my $inextro = 0;
	my $topopos = 0;
	while(<PRFFILE>) {
	    if($_ =~ "START 1") {
		$intro .= $_;
		$inintro = 0;
	    }
	    if($_ =~ "END 1") {
		$inextro = 1;
	    }
	    if($_ =~ 'ALPHABET' && $inextro <= 0) {
		$intro .= $_;
	    }
	    if($_ =~ 'COL ') {
		chomp;
		$_ =~ s/\s+/ /g;
		my @col = split /\s/, $_;
		$col[$#col - 1] = $topology[$topopos];
		push @cols, [@col];
		$topopos++;
	    }
	    if($inintro > 0) {
		$intro .= $_;
	    }
	    if($inextro > 0) {
		$extro .= $_;
	    }
	    
	}
	close PRFFILE;
	open OUTFILE, ">"."$prffile"
	    or die "Could not open $prffile for writing\n";
	print OUTFILE "$intro";
	for(my $i = 0; $i <= $#cols; $i++) {
	    my @col = @{$cols[$i]};
	    printf OUTFILE "COL %4d", $i+1;
	    print OUTFILE ":  ";
	    for $b (2 .. $#col ) {
		if($b >= ($#col - 1)) {
		    print OUTFILE "$col[$b]       ";
		}
		else {
		    printf OUTFILE "%5.02f   ", $col[$b];
		}
	    }
	    print OUTFILE "\n";
	}
	print OUTFILE "$extro";
	close OUTFILE;
    }
    return;
}

sub relabel_fa_file {
    my ($fa_file,$topo) = @_;
    chomp($topo);
    open(OUT,">>$fa_file") || die;
    print OUT "/$topo/\n";
    close(OUT);
}

