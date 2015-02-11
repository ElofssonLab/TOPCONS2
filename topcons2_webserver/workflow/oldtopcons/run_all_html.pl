#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(gettimeofday tv_interval);

##########################################
# Parameters:
my $convert="/usr/bin/convert";
my $mktemp="/bin/mktemp";
my $sort="/usr/bin/sort";
my $awk="/usr/bin/awk";
my $head="/usr/bin/head";
my $gnuplot="/usr/bin/gnuplot";
##########################################

$ENV{'PATH'}="/bin:/usr/bin";
$ARGV[0] =~ /^([\w\-\/\.]+)$/ || die;
my $tmpdir=$1;



$ARGV[1]=~ /^(\w+)$/ || die;
my $resultdir=$1;


$ARGV[2]=~ /^([\w\-\/:\.]+)$/ || die;
#$ARGV[1] =~ /^(.*)$/ || die;
my $url=$1;

$ARGV[3] =~ /^(.+)$/ || die;
my $workdir=$1;

chdir($workdir); #added 2014-05-26 by nanjiang

$ENV{'GDFONTPATH'}="/server/var/www/modhmm-projects-common/fonts/";


my $start=[gettimeofday];

my $seq=""; 
my %topo=();
my $topcons_pred="";
open(IN,"$tmpdir/query.tmp.fa") || die;
open(OUT,">$tmpdir/query.fa") || die;
my $name=<IN>;
chomp($name);
$name =~ s/^>//;
print OUT ">query\n";
while(<IN>) {
    $seq .= $_;
    print OUT;
}
close(OUT);
close(IN);
chomp($seq);



my $fix_str="";
my $fix_topo="";
if (open(IN,"$tmpdir/query.tmp.fix")) {
    $fix_str=<IN>;
    close(IN);
    chomp($fix_str);
    open(OUT,">$tmpdir/query.fix");
    $fix_topo=fix2topo( $fix_str, length($seq) );
    print OUT $fix_topo . "\n";
    close(OUT);
}


my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date=sprintf("%4d-%02d-%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);


`/usr/bin/perl -w -T run_all.pl $tmpdir > $tmpdir/run_all.pl_result.txt 2> $tmpdir/run_all.pl.error_log`;

system("cat $tmpdir/child1.out $tmpdir/child2.out $tmpdir/child3.out $tmpdir/child4.out $tmpdir/master.out > $tmpdir/children.out");

open(OUTPUT, "$tmpdir/children.out") || die("Couldn't open file '$tmpdir/run_all.pl_result.txt'");
while (<OUTPUT>) {
    ($_ =~ /^OCTOPUS:\s(.*)$/) && ($topo{'OCTOPUS'}=$1);
    /^PRO:\s(.*)$/ && ($topo{'PRO'}=$1);
    /^PRODIV:\s(.*)$/ && ($topo{'PRODIV'}=$1);
    /^SCAMPI_multi:\s(.*)$/ && ($topo{'SCAMPI-msa'}=$1);
    /^SCAMPI_single:\s(.*)$/ && ($topo{'SCAMPI-seq'}=$1);
    /^TOPCONS:\s(.*)$/ && ($topcons_pred=$1);
}
close(OUTPUT);

my $i=0;
while(my $fix_chr=substr($fix_topo,0,1,"")) {
    if ($fix_chr ne ".") {
	for my $method (sort keys %topo) {
	    if ($topo{$method} ne "ERR" && (substr($topo{$method},$i,1) ne $fix_chr)) {
		$topo{$method} = "ERR";
	    }
	}
    }
    $i++;
}

for my $method (sort keys %topo) {
    $topo{$method} =~ s/^[io]*$/NOTM/;
    if ($topo{$method} eq "ERR") {
	$topcons_pred = "ERR";
    }
}
$topcons_pred =~ s/^[io]*$/NOTM/;

my $length=length($seq);
$seq =~  s/(.{60})/$1\n/g;
mk_plot($tmpdir,\%topo, $length);
mk_topcons_plot($tmpdir,$topcons_pred,$length);

my $time=sprintf("%.2f",tv_interval($start));
mk_text($tmpdir,\%topo, $topcons_pred, $length, $name, $seq, $time, $url, $date, $fix_str);



system("/bin/cp $tmpdir/query.png /server/var/www/topcons.cbr.su.se/result/$resultdir/result.png");
system("/bin/cp $tmpdir/query.large.png /server/var/www/topcons.cbr.su.se/result/$resultdir/result.large.png");
system("/bin/cp $tmpdir/query.topcons.png /server/var/www/topcons.cbr.su.se/result/$resultdir/result.topcons.png");
system("/bin/cp $tmpdir/query.topcons.large.png /server/var/www/topcons.cbr.su.se/result/$resultdir/result.topcons.large.png");
system("/bin/cp $tmpdir/query.txt /server/var/www/topcons.cbr.su.se/result/$resultdir/topcons.txt");
system("/bin/cp $tmpdir/query.blast /server/var/www/topcons.cbr.su.se/result/$resultdir/blast.txt");

my $output= "";

print "Sequence name: $name<br>Sequence length: $length aa.<br>\n\n";
$output .= "Sequence name: $name<br>Sequence length: $length aa.<br>\n\n";

if ($fix_str) {
    print "Restraints: $fix_str<br>\n\n";
    $output .= "Restraints: $fix_str<br>\n\n";
}
print "Total request time: $time seconds\n";
$output .= "Total request time: $time seconds\n";

print "<br><br> The output data can be found in the <a href=\"result/$resultdir/topcons.txt\">TOPCONS result file (txt)</a><br>";
print "The multiple sequence alignment can be found in the <a href=\"result/$resultdir/blast.txt\">BLAST result file (txt)</a><br>";
print "<br>Predicted topologies, predicted &Delta;G values and predicted distances to the membrane center (Z=0):<br><br><br><img width=600 src=\"result/$resultdir/result.png\" border=0>\n";
print "<br><div style=\"text-align:right\"><a href=\"result/$resultdir/result.large.png\">High-resolution image&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></div>\n";

$output .= "<br><br> The output data can be found in the <a href=\"result/$resultdir/topcons.txt\">TOPCONS result file (txt)</a><br>";
$output .= "The multiple sequence alignment can be found in the <a href=\"result/$resultdir/blast.txt\">BLAST result file (txt)</a><br>";
$output .= "<br>Predicted topologies, predicted &Delta;G values and predicted distances to the membrane center (Z=0):<br><br><br><img width=600 src=\"result/$resultdir/result.png\" border=0>\n";
$output .= "<br><div style=\"text-align:right\"><a href=\"result/$resultdir/result.large.png\">High-resolution image&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></div>\n";

if (($topcons_pred ne "NOTM") && ($topcons_pred ne "ERR")) {
    print "<br><br>Consensus prediction (TOPCONS):<br><img width=600 src=\"result/$resultdir/result.topcons.png\" border=0>\n";
    print "<br><div style=\"text-align:right\"><a href=\"result/$resultdir/result.topcons.large.png\">High-resolution image&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></div>\n";
    $output .= "<br><br>Consensus prediction (TOPCONS):<br><img width=600 src=\"result/$resultdir/result.topcons.png\" border=0>\n";
    $output .= "<br><div style=\"text-align:right\"><a href=\"result/$resultdir/result.topcons.large.png\">High-resolution image&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></div>\n";
}
print "<br><br>\n";
$output .= "<br><br>\n";

my @sort=sort {$b cmp $a} (keys %topo);
$topo{'TOPCONS'}=$topcons_pred;
push(@sort,'TOPCONS');

print "Predicted TM-helix positions:<br><br>" . print_TMcoordinates(\@sort,\%topo) . "<br><br>\n";
$output .= "Predicted TM-helix positions:<br><br>" . print_TMcoordinates(\@sort,\%topo) . "<br><br>\n";
$seq =~ s/\s//g;
print "Sequence and predicted topologies:<br><br>" . print_seq_topo($seq,\@sort,\%topo) . "<br><br>\n";
$output .= "Sequence and predicted topologies:<br><br>" . print_seq_topo($seq,\@sort,\%topo) . "<br><br>\n";

my $text = "";
$text .= $output;

open(RESULTS, ">$tmpdir/results.html") || die("Couldn't open file $tmpdir/results.html for writing");
print RESULTS $text;
close(RESULTS);

sub fix2topo {
    my ($fix_str,$length) = @_;
    my @topo = split(//,mk_str(".",$length));
    for (split(/;/,$fix_str)) {
	my ($ll,$ul,$label)=split /-/;
	for (my $i=$ll;$i<=$ul;$i++) {
	    $topo[$i-1]=$label;
	}
    }
    return join("",@topo);
}

sub mk_str {
    my ($chr,$length) = @_;
    my $str="";
    for (my $i=0;$i<$length;$i++) {
	$str .= $chr;
    }
    return $str;
}

sub print_seq_topo {
    my ($seq,$sort,$topo) = @_;
    my $rowL=50; # Parameter
    my $colL=10; # Parameter
    my $cols_per_block=$rowL/$colL; ($rowL % $colL) && die; # Must be an integer...
    my $length=length($seq);
    my $nr_complete_blocks=int($length/$rowL);
    my $final_row_L=$length-$nr_complete_blocks*$rowL;
    my $ret = "<table cellspacing=0 cellpadding=0 border=0>\n";
    my $i;
#    for my $method (sort {$b cmp $a} keys %topo) {
    for my $method (@{$sort}) {
	my @BOUNDARIES=get_boundaries_from_topology($topo->{$method});
	if (($topo->{$method} =~ /NOTM/i) || ($topo->{$method} =~ /ERR/i) || scalar(@BOUNDARIES) == 0 || scalar(@BOUNDARIES) == 1) {
	    delete($topo->{$method});
	}
    }
    $topo{'Seq.'}=$seq;
    @{$sort}=("Seq.",@{$sort});
    for ($i=1;$i<=$nr_complete_blocks;$i++) {
	$ret.= "<tr><td>&nbsp;</td><td align=\"left\" style=\"font-family: Courier;\">";
	for (my $j=1;$j<=$cols_per_block;$j++) {
	    my $nr=(($i-1)*$rowL+($j-1)*$colL+1);
	    my $spaces=spaces($nr,$colL+1);
	    (!(($nr-1)%$rowL) || (!(($nr+$colL-1)%$rowL))) || ($nr =~ s/./&nbsp;/g);
	    $ret .= $nr . $spaces;
	}
	$ret .= "</td></tr>";
#	for my $method (sort {$b cmp $a} keys %topo) {
	for my $method (@{$sort}) {
	    if ($method eq "TOPCONS") {
		$ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\"><span class=nowrap>$method&nbsp;&nbsp;&nbsp;</span></td>";
	    } else {
		$ret .= "<tr><td style=\"font-family: Courier;\"><span class=nowrap>$method&nbsp;&nbsp;&nbsp;</span></td>";
	    }
	    my $row = substr($topo->{$method},0,$rowL,"");
	    $row =~ s/(.{$colL})/$1\ /g;
	    if ($method eq "TOPCONS") {
		$ret .= "<td style=\"font-family: Courier; Color:#990000;\">$row</td></tr>\n";
	    } else {
		$ret .= "<td style=\"font-family: Courier;\">$row</td></tr>\n";
	    }
	}
	$ret .= "<tr><td colspan=2>&nbsp;</td></tr>\n";
    }
    $ret.= "</table>\n";
    unless ($length % $rowL == 0) {
	$ret.="<table cellspacing=0 cellpadding=0 border=0>\n"; 
	$ret.= "<tr><td>&nbsp;</td><td align=\"left\" style=\"font-family: Courier;\">";
	for (my $j=1;$j<=int(($length-1-$nr_complete_blocks*$rowL)/$colL)+1;$j++) {
	    my $nr=(($i-1)*$rowL+($j-1)*$colL+1);
	    my $spaces=spaces($nr,$colL+1);
	    ((!(($nr-1)%$rowL)) || ($length-$nr<$colL)) || ($nr =~ s/./&nbsp;/g);
	    $ret .= $nr . $spaces;
	}
	$ret .= "</td></tr>";
#	for my $method (sort {$b cmp $a} keys %topo) {
	for my $method (@{$sort}) {
	    if ($method eq "TOPCONS") {
		$ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\">$method&nbsp;&nbsp;&nbsp;</td>";
	    } else {
		$ret .= "<tr><td style=\"font-family: Courier;\">$method&nbsp;&nbsp;&nbsp;</td>";
	    }
	    my $row = substr($topo->{$method},0,$rowL,"");
	    $row =~ s/(.{$colL})/$1\ /g;
	    if ($method eq "TOPCONS") {
		$ret .= "<td style=\"font-family: Courier; Color:#990000;\">$row</td></tr>\n";
	    } else {
		$ret .= "<td style=\"font-family: Courier;\">$row</td></tr>\n";
	    }
	}
	$ret .= "</table>\n";
    }
    return $ret;
}

sub spaces {
    my ($str,$length)=@_;
    my $nr_spaces=$length-length($str);
    my $space_str="";
    for (my $i=0;$i<$nr_spaces;$i++) {
	$space_str .= "&nbsp;";
    }
    return $space_str;
}

sub print_TMcoordinates {
    my ($sort,$topo) = @_;
    my $result="";
    $result .= "<table cellspacing=0 cellpadding=0 border=0 width=\"100%\">\n";
    $result .= "<tr><td><img src=\"images/0.gif\" width=50 height=1></td><td><img src=\"images/0.gif\" width=50 height=1></td><td><img src=\"images/0.gif\" width=500 height=1></td></tr>\n";
#    for my $method (sort {$b cmp $a} keys %topo) {
    for my $method (@{$sort}) {
	$topo->{$method} =~ s/\s//g;
	if ($method eq "TOPCONS") {
	    $result .= "<tr><td colspan=2 valign=\"top\"><div style=\"color:#990000\">$method</div></td><td valign=\"top\"><div style=\"color:#990000\">&nbsp;";
	} else {
	    $result .= "<tr><td colspan=2 valign=\"top\">$method</td><td valign=\"top\">&nbsp;";
	}
	my @BOUNDARIES=get_boundaries_from_topology($topo->{$method});

	if ($topo{$method} =~ /ERR/i) {
	    $result .= "***No\ topology\ could\ be\ produced,\ check\ prediction\ constraints***";
	} elsif (($topo->{$method} =~ /NOTM/i) || scalar(@BOUNDARIES) == 0 || scalar(@BOUNDARIES) == 1) {
	    $result .= "***No TM-regions predicted***";
	} else {
	    my $TMnr=1;
	    for(my $i=0;$i<scalar(@BOUNDARIES);$i++){
		my ($start,$end,$type)=split(/\t/,$BOUNDARIES[$i]);
		if ($type eq "M") {
		    $result .= ($TMnr++).". ".($start+1)."-$end, ";
		}
	    }
	}
	$result =~ s/,\s$//;
	if ($method eq "TOPCONS") {
	    $result .= "</div></td></tr>\n";
	} else {
	    $result .= "</td></tr>\n";
	}
	$result .= "<tr><td colspan=3><img src=\"images/0.gif\" width=1 height=5></td></tr>\n";
    }
    $result .= "</table>\n";
    return $result;
}


sub mk_text {
    my ($tmpdir,$ext_topo,$topcons_pred,$length,$name,$seq,$time,$url,$date,$fix_str) = @_;
    my $i;
    my %topo = %{$ext_topo};
    for my $method (sort keys %topo) {
	$topo{$method} =~ s/(.{60})/$1\n/g;
	$topo{$method} =~ s/NOTM/***No\ TM-regions\ predicted***/;
	$topo{$method} =~ s/ERR/***No\ topology\ could\ be\ produced,\ check\ prediction\ constraints***/;
    }
    $topcons_pred =~ s/(.{60})/$1\n/g;
    $topcons_pred =~ s/NOTM/***No\ TM-regions\ predicted***/;
    $topcons_pred =~ s/ERR/***No\ topology\ could\ be\ produced,\ check\ prediction\ constraints***/;
    open(OUT,">$tmpdir/query.txt") || die;
    print OUT "##############################################################################\n";
    print OUT "TOPCONS result file\n";
    print OUT "Generated from http://$url at $date\n";
    print OUT "Total request time: $time seconds.\n";
    print OUT "##############################################################################\n\n\n";
    print OUT "Sequence name: $name\n";
    print OUT "Sequence length: $length aa.\n";
    if ($fix_str) {
	print OUT "Restraints: $fix_str\n";
    }
    print OUT "Sequence:\n$seq\n\n";
    print OUT "SCAMPI-seq predicted topology:\n" . $topo{'SCAMPI-seq'} . "\n\n";
    print OUT "SCAMPI-msa predicted topology:\n" . $topo{'SCAMPI-msa'} . "\n\n";
    print OUT "PRODIV predicted topology:\n" . $topo{'PRODIV'} . "\n\n";
    print OUT "PRO predicted topology:\n" . $topo{'PRO'} . "\n\n";
    print OUT "OCTOPUS predicted topology:\n" . $topo{'OCTOPUS'} . "\n\n";
    print OUT "TOPCONS predicted topology:\n" . $topcons_pred . "\n\n";
    if (open(IN,"$tmpdir/query.zpred")) {
	print OUT "Predicted Z-coordinates (left column=sequence position; right column=Z-coordinate)\n";
	while(<IN>) {
	    /^(\S+)/ || die;
	    print OUT "".(++$i)." \t ".sprintf("%.2f",$1)."\n";
	}
	close(IN);
    }
    if (open(IN,"$tmpdir/query.DG")) {
	print OUT "\nPredicted Delta-G-values (kcal/mol) (left column=sequence position; right column=Delta-G)\n";
	while(<IN>) {
	    /^(\S+)\s+(\S+)/ || die;
	    print OUT $1."\t".sprintf("%.2f",$2)."\n";
	}
	close(IN);
    }
    unless ($topcons_pred =~ /^\*\*\*/) {
	if (open(IN,"$tmpdir/query.reliability")) {
	    print OUT "\nPredicted TOPCONS reliability (left column=sequence position; right column=reliability)\n";
	    while(<IN>) {
		/^(\S+)\s+(\S+)/ || die;
		print OUT $1."\t".sprintf("%.3f",$2)."\n";
	    }
	    close(IN);
	}
    }
    close(OUT);
}

sub mk_plot {
    my ($tmpdir,$topo, $length)=@_;
    my $infile="$tmpdir/query";
    my $DG_min=`$sort -gk 2 $tmpdir/query.DG|$head -n 1|$awk '{print \$2}'`;
    my $DG_max=`$sort -grk 2 $tmpdir/query.DG|$head -n 1|$awk '{print \$2}'`;
    my $DG_tick_max=int($DG_max+0.5 * ($DG_max <=> 0))+1;
    my $DG_tick_min=int($DG_min+0.5 * ($DG_min <=> 0))-1;
    my $DG_axis_max=1.8*($DG_tick_max-$DG_tick_min)+$DG_tick_min;
#    print "$DG_min<br>\n$DG_max<br>\n$DG_tick_min<br>\n$DG_tick_max<br>\n$DG_axis_max<br>\n";
    chomp($DG_min,$DG_max);
    my $obj_nr;
    open(OUT,">$infile.gnu") || die("Could not open file");
    print OUT "set encoding iso_8859_1\n";
    print OUT "set yrange [0:45]\n";
    print OUT "set xrange [1:$length]\n";
    print OUT "set y2range [$DG_tick_min:$DG_axis_max]\n";
    print OUT "set autoscale xfix\n";
    print OUT "set ter png enh interlace size 2400,1680 font 'Nimbus,40'\n";
    print OUT "set ylabel \'Absolute value of Z-coord. (Å)                                           \' tc lt 2 offset 6\n";
    print OUT "set ylabel \'Absolute value of Z-coord. (Å)                                           \' tc lt 2\n";
    print OUT "set y2label \'{/Symbol D}G (kcal/mol)                                              \' tc lt 3\n";
    print OUT "set ytics scale 1,0.5 nomirror \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25, \"OCTOPUS\" 30.5 0, \"PRO\" 33.5 0, \"PRODIV\" 36.5 0, \"SCAMPI-msa\" 39.5 0, \"SCAMPI-seq\" 42.5 0\)\n";
    if ($DG_tick_max-$DG_tick_min < 10) {
	print OUT "set y2tics nomirror $DG_tick_min,1,$DG_tick_max\n";
    } else {
	print OUT "set y2tics nomirror $DG_tick_min,2,$DG_tick_max\n";
    }
    print OUT "set out \'$infile.large.png\'\n";

    print OUT "set lmargin 11.5\n";
    print OUT "set rmargin 6.5\n";


    print OUT "set tmargin 1.3\n";
    print OUT "set object " . (++$obj_nr) . " rect from screen 0.19,0.986 to screen 0.21,0.992 fc rgb \"red\" fs noborder\n";
    print OUT "set label 'Inside' font 'Nimbus,30' at screen 0.215,0.982\n";
    print OUT "set object " . (++$obj_nr) . " rect from screen 0.28,0.986 to screen 0.30,0.992 fc rgb \"blue\" fs noborder\n";
    print OUT "set label 'Outside' font 'Nimbus,30' at screen 0.305,0.982\n";
    print OUT "set object " . (++$obj_nr) . " rect from screen 0.38,0.978 to screen 0.40,1 fc rgb \"grey\" fs noborder\n";
    print OUT "set label 'TM-helix (IN->OUT)' font 'Nimbus,30' at screen 0.405,0.982\n";
    print OUT "set object " . (++$obj_nr) . " rect from screen 0.57,0.978 to screen 0.59,1 fc rgb \"white\"\n";
    print OUT "set label 'TM-helix (OUT->IN)' font 'Nimbus,30' at screen 0.595,0.982\n";
    print OUT "set object " . (++$obj_nr) . " rect from screen 0.76,0.978 to screen 0.78,1 fc rgb \"black\"\n";
    print OUT "set label 'Reentrant region' font 'Nimbus,30' at screen 0.785,0.982\n";

    my $offset=0;
    for my $method (sort keys %{$topo}) {
	my @BOUNDARIES=get_boundaries_from_topology($topo->{$method});
	if ($topo->{$method} =~ /ERR/) {
	    print OUT "set label \"***No topology could be produced, check prediction constraints***\" at " . (int($length/20+1)+1) . "," . (30+$offset+0.5) . "\n";
	} elsif (($topo->{$method} =~ /NOTM/i) || scalar(@BOUNDARIES) == 0 || scalar(@BOUNDARIES) == 1) {
	    print OUT "set label \"***No TM-regions predicted***\" at " . (int($length/20+1)+1) . "," . (30+$offset+0.5) . "\n";
	} else {
	    for(my $i=0;$i<scalar(@BOUNDARIES);$i++){
		my @fields=split(/\t/,$BOUNDARIES[$i]);
		my $start=$fields[0]+1-0.5;
		my $end=$fields[1]+0.5;
# 		my ($start,$end);
# 		if (($fields[2] eq "M") || ($fields[2] eq "R")) {
# 		    $start=$fields[0]+1;
# 		    $end=$fields[1];
# 		} else {
# 		    $start=$fields[0];
# 		    $end=$fields[1]+1;
# 		}
		if($fields[2] eq "i"){
		    print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+0.2) . " fc rgb \"red\" fs noborder\n";
		}
		if($fields[2] eq "o"){
		    print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset+1-0.2) . " to $end," . (30+$offset+1) . " fc rgb \"blue\" fs noborder\n";
		}
		if($fields[2] eq "M"){
		    if ((substr($topo->{$method},$start-2,1) eq "i") || (substr($topo->{$method},$end,1) eq "o")) {
			print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+1) . " fc rgb \"grey\" fs noborder\n";# fs pattern 5\n";
		    } else {
			print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+1) . " fc rgb \"white\"\n";# fs pattern 4\n";
		    }
		}
		if($fields[2] eq "R"){
		    if ((substr($topo->{$method},$start-2,1) eq "i") || (substr($topo->{$method},$end,1) eq "i")) {
			print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+0.5) . " fc rgb \"black\" \n";
		    } else {
			print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset+0.5) . " to $end," . (30+$offset+1) . " fc rgb \"black\" \n";
		    }
		}
	    }
	}
	print "\n";
	$offset+=3;
    }
    
    print OUT "plot '$infile.zpred' u (\$0+1):1 w l t '' lt 2 lw 4 ,'$infile.DG' axes x1y2 w l t '' lt 3 lw 4\n";
    print OUT "exit\n";
    close(OUT);
    `$gnuplot $infile.gnu`;
    `$convert -scale 600x420 $infile.large.png $infile.png`;
#    `mv $infile.2.png $infile.png`;
}

sub mk_topcons_plot {
    my ($tmpdir,$topo, $length)=@_;
    my $infile="$tmpdir/query";
    my $rel_min=`$sort -gk 2 $tmpdir/query.reliability|$head -n 1|$awk '{print \$2}'`;
    my $rel_max=`$sort -grk 2 $tmpdir/query.reliability|$head -n 1|$awk '{print \$2}'`;
    my $rel_tick_max=int(10*$rel_max+1)/10;
    ($rel_tick_max > 1) && ($rel_tick_max = 1);
    my $rel_tick_min=int(10*$rel_min)/10;
    my $rel_axis_max=1.5*($rel_tick_max-$rel_tick_min)+$rel_tick_min;
    my $factor=$rel_axis_max-$rel_tick_max;
    my $offset=0.4*$factor;

#    print "$rel_min<br>\n$rel_max<br>\n$rel_tick_min<br>\n$rel_tick_max<br>\n$rel_axis_max<br>\n";
    chomp($rel_min,$rel_max);
    my $obj_nr;
    open(OUT,">$infile.topcons.gnu") || die("Could not open file");
    print OUT "set encoding iso_8859_1\n";
    print OUT "set xrange [1:$length]\n";
    print OUT "set yrange [$rel_tick_min:$rel_axis_max]\n";
    print OUT "set autoscale xfix\n";
    print OUT "set ter png enh interlace size 2400,840 font 'Nimbus,40'\n";
    print OUT "set xlabel \'Position\'\n";
    print OUT "set ylabel \'Reliability           \' \n";
#    print OUT "set ytics scale 1,0.5 nomirror \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25, \"OCTOPUS\" 30.5 0, \"PRO\" 33.5 0, \"PRODIV\" 36.5 0, \"SCAMPI-msa\" 39.5 0, \"SCAMPI-seq\" 42.5 0\)\n";
 #   if ($rel_tick_max-$rel_tick_min < 10) {
    print OUT "set ytics nomirror $rel_tick_min,0.1,$rel_tick_max\n";
#    } else {
#	print OUT "set y2tics nomirror $rel_tick_min,2,$rel_tick_max\n";
#    }
    print OUT "set out \'$infile.topcons.large.png\'\n";

    print OUT "set tmargin 1.3\n";


    print OUT "set lmargin 11.5\n";
    print OUT "set rmargin 6.5\n";
    print OUT "set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775\n";


#     print OUT "set object " . (++$obj_nr) . " rect from screen 0.19,0.986 to screen 0.21,0.992 fc rgb \"red\" fs noborder\n";
#     print OUT "set label 'Inside' font 'Nimbus,30' at screen 0.215,0.982\n";
#     print OUT "set object " . (++$obj_nr) . " rect from screen 0.28,0.986 to screen 0.30,0.992 fc rgb \"blue\" fs noborder\n";
#     print OUT "set label 'Outside' font 'Nimbus,30' at screen 0.305,0.982\n";
#     print OUT "set object " . (++$obj_nr) . " rect from screen 0.38,0.978 to screen 0.40,1 fc rgb \"grey\" fs noborder\n";
#     print OUT "set label 'TM-helix (IN->OUT)' font 'Nimbus,30' at screen 0.405,0.982\n";
#     print OUT "set object " . (++$obj_nr) . " rect from screen 0.57,0.978 to screen 0.59,1 fc rgb \"white\"\n";
#     print OUT "set label 'TM-helix (OUT->IN)' font 'Nimbus,30' at screen 0.595,0.982\n";
#     print OUT "set object " . (++$obj_nr) . " rect from screen 0.76,0.978 to screen 0.78,1 fc rgb \"black\"\n";
#     print OUT "set label 'Reentrant region' font 'Nimbus,30' at screen 0.785,0.982\n";
    my @BOUNDARIES=get_boundaries_from_topology($topo);
    if ($topo =~ /ERR/) {
	print OUT "set label \"***No topology could be produced, check prediction constraints***\" at " . (int($length/20+1)+1) . "," . ($rel_tick_max+$offset+0.5) . "\n";
    } elsif (($topo =~ /NOTM/i) || scalar(@BOUNDARIES) == 0 || scalar(@BOUNDARIES) == 1) {
	print OUT "set label \"***No TM-regions predicted***\" at " . (int($length/20+1)+1) . "," . ($rel_tick_max+$offset+0.5) . "\n";
    } else {
	for(my $i=0;$i<scalar(@BOUNDARIES);$i++){
	    my @fields=split(/\t/,$BOUNDARIES[$i]);
	    my $start=$fields[0]+1-0.5;
	    my $end=$fields[1]+0.5;
# 		my ($start,$end);
# 		if (($fields[2] eq "M") || ($fields[2] eq "R")) {
# 		    $start=$fields[0]+1;
# 		    $end=$fields[1];
# 		} else {
# 		    $start=$fields[0];
# 		    $end=$fields[1]+1;
# 		}
	    my $small=0.0325;
	    my $big=0.175;
	    if($fields[2] eq "i"){
		print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$small*$factor) . " fc rgb \"red\" fs noborder\n";
	    }
	    if($fields[2] eq "o"){
		print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset+$big*$factor-$small*$factor) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"blue\" fs noborder\n";
	    }
	    if($fields[2] eq "M"){
		if ((substr($topo,$start-2,1) eq "i") || (substr($topo,$end,1) eq "o")) {
		    print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"grey\" fs noborder\n";# fs pattern 5\n";
		} else {
		    print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"white\"\n";# fs pattern 4\n";
		}
	    }
#	    if($fields[2] eq "R"){
#		if ((substr($topo,$start-2,1) eq "i") || (substr($topo,$end,1) eq "i")) {
#		    print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+0.5) . " fc rgb \"black\" \n";
#		} else {
#		    print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset+0.5) . " to $end," . ($rel_tick_max+$offset+1) . " fc rgb \"black\" \n";
#		}
#	    }
	}
    }
    print "\n";
    
    print OUT "plot '$infile.reliability' w l t '' lc rgb \"black\" lw 4\n";
    print OUT "exit\n";
    close(OUT);
#    $ENV{'GDFONTPATH'}='./fonts/';
    `$gnuplot $infile.topcons.gnu`;
    `$convert -scale 600x210 $infile.topcons.large.png $infile.topcons.png`;
#    `mv $infile.topcons.2.png $infile.topcons.png`;
}

sub get_boundaries_from_topology {
    my @BOUNDARIES;
    my ($topology)=@_;
    my $temp;
    
    while($topology=~m/i+/g){
	$temp=$-[0]."\t".$+[0]."\ti";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/o+/g){
	$temp=$-[0]."\t".$+[0]."\to";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/O+/g){
	$temp=$-[0]."\t".$+[0]."\to";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/M+/g){
	$temp=$-[0]."\t".$+[0]."\tM";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/r+/g){
	$temp=$-[0]."\t".$+[0]."\tR";
	push @BOUNDARIES, $temp;
    }
    while($topology=~m/R+/g){
	$temp=$-[0]."\t".$+[0]."\tR";
	push @BOUNDARIES, $temp;
    }
    return(@BOUNDARIES);
}
