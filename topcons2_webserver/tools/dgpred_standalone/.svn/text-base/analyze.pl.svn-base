#!/usr/bin/perl -w

#($helix=$ARGV[0]) || die;

# ($with_length=$ARGV[1]) || ($with_length = "");
use File::Basename;
use Cwd qw(realpath);

my $fullpath = realpath($0);
my $programpath=dirname($fullpath);


my $with_length="on";
my $allow_sub="off";
my $outformat="html";
my $outfile="/dev/stdout";

my $usage="
Usage:   analyze.pl [Options] <helix sequence>
Note: the output is in HTML format
Options:
    -lc |--with-length-corr on|off    : whether use length correction, default=$with_length
    -sub|--allow-sub        on|off    : whether allow subsequences, default=$allow_sub
    -o|--outfile            <file>    : output the result to outpath, default=$outfile

Modified from the original perl script by Nanjiang Shu, created 2010-09-01, updated 2010-09-01

Examples:
    ./analyze.pl  YIYLGGAILAEVIGTTLMKF -o rst.html
";

my $numArgs = $#ARGV+1;
if($numArgs < 1)
{
    &PrintHelp;
    exit;
}

my $helix="";

if(@ARGV)#{{{
{
    my $i = 0;
    while($ARGV[$i])
    {
        if($ARGV[$i] eq "-h" || $ARGV[$i] eq "--help" )
        {
            &PrintHelp;
            exit;
        }
        elsif($ARGV[$i] eq "-o" || $ARGV[$i] eq "--outfile"  || $ARGV[$i] eq "-outfile" )
        {   
            $outfile = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-of" || $ARGV[$i] eq "--outformat"  || $ARGV[$i] eq "-outformat" )
        {   
            $outformat = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-lc" || $ARGV[$i] eq "--with-length-corr" )
        {   
            $with_length = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-sub" || $ARGV[$i] eq "--allow-sub" )
        {   
            $allow_sub = $ARGV[$i+1];
            $i +=2;
        }
        else 
        {
            $helix = $ARGV[$i];
            $i ++;
        }
    }
}#}}}

my $outfilepath="";
if ($outfile ne "/dev/stdout")
{
    $outfilepath=dirname($outfile);
    system("cp -f $programpath/*.gif  $outfilepath/");
}

$fpout=STDOUT;
if ($outfile ne "/dev/stdout")
{
    open($fpout,">", $outfile) || die "Can not open the output file $outfile for write\n";
}

my $dG_matrix_length=19;
my $c0=0.27045;
my $c1=9.29274167549645;
my $c2=-0.64513139783394;
my $c3=0.00822196628688;

my %profiles = (  'A' => [  0.1267255,  0.0215152  ],
		  'B' => [  1.5556351,  0.0133523  ],
                  'C' => [ -0.0765051,  0.0994228  ],
	          'D' => [  1.7939795,  0.0172643  ],
   	          'E' => [  1.4193720,  0.0089351  ],
	          'F' => [ -0.2766953,  0.0010297  ],
	          'G' => [  0.4813492,  0.0047210  ],
	          'H' => [  1.1998590,  0.0080127  ],
	          'I' => [ -0.4597384,  0.0181495  ],
	          'K' => [  1.8485768,  0.0218446  ],
	          'L' => [ -0.4282992,  0.0023804  ],
	          'M' => [ -0.0774786,  0.0984413  ],
   	          'N' => [  1.3266132,  0.0092375  ],
	          'P' => [  1.0860888,  0.0100568  ],
	          'Q' => [  1.3336109,  0.0111996  ],
	          'R' => [  1.6492534,  0.0512044  ],
	          'S' => [  0.7023921,  0.0077661  ],
	          'T' => [  0.5266550,  0.0311973  ],
	          'U' => [ -0.0774786,  0.0984413  ],
	          'V' => [ -0.2447218,  0.0979201  ],
	          'W' => [  0.2909390,  0.0189282, -0.5479140,  0.0930222,  6.4736619  ],
	          'X' => [  0.6348884,  0.0180273  ],                                           
	          'Y' => [  0.6275249,  0.0103896, -0.5744404,  0.0947821,  6.9164963  ],
	          'Z' => [  1.3761092,  0.0099898  ]  );

unless ($with_length eq "on") {
    $c1=0; $c2=0; $c3=0;
}

my %aa2nr = ( "A" => 1,"C" => 2,"D" => 3,"E" => 4,"F" => 5,"G" => 6,"H" => 7,"I" => 8,"K" => 9,"L" => 10,
	      "M" => 11,"N" => 12,"P" => 13,"Q" => 14,"R" => 15,"S" => 16,"T" => 17,"V" => 18,"W" => 19,"Y" => 20);
my $piover180=atan2(1,1)/45;
my $totwidth=40;

print $fpout  "<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\" width=\"100%\" class=\"TableGrey\">\n";
print $fpout  "<tr><td><table cellpadding=0 cellspacing=0 width='100%' border=0><tr><td class='orangebg'><img src='0.gif' border=0 width=1 height=1></td></tr></table></td></tr>";
print $fpout  "<tr><td><img src='0.gif' width=1 height=10></td></tr>\n";
print $fpout  "<tr><td colspan=6 align=left><table width=\"100%\" cellspacing=0 cellpadding=0 border=0 class='analyze'>\n";
print $fpout  "<tr><td colspan=6>";

print $fpout  "<table cellspacing=0 cellpadding=0 border=0><tr><td>";
print $fpout  "<big><big>Contributions to predicted &Delta;G<small><span style='vertical-align: sub'>app</span></small></big></big>";
print $fpout  "<table cellpadding=0 cellspacing=0 width='100%' border=0><tr><td class='orangebg'><img src='0.gif' border=0 width=1 height=1></td></tr></table>";
print $fpout  "</td></tr></table>\n";

print $fpout  "</td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6><i><b>Side-chain contributions</b></i></td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td><i>Amino acid</i></td><td><i>Position</i></td><td><i>Rel. position</i></td><td>&nbsp;</td><td colspan=2 align='center'><i>&Delta;G<small><span style='vertical-align: super'>aa(i)</span><span style='vertical-align: sub'>app</span></small></i></td></tr>\n";
print $fpout  "<tr><td><img src='0.gif' height=1 width=100></td><td><img src='0.gif' height=1 width=100></td><td><img src='0.gif' height=1 width=100></td><td><img src='0.gif' height=1 width=10></td><td><img src='0.gif' height=1 width=$totwidth></td><td><img src='0.gif' height=1 width=$totwidth></td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";

($helix =~ /^\s+$/) && (print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n") && exit;
$helix =~ s/\s*$//;
$helix =~ s/^\s*//;
$helix =~ tr/a-z/A-Z/;
$Length=length($helix);
unless ($helix =~ /^[ABCDEFGHIKLMNPQRSTUVWXYZ]+$/ && $Length <= 50) {
    print $fpout  "<tr><td colspan=6>Error. </td></tr>\n";
    exit;
}
$Lmin=$Length; $Lmax=$Length;
%putative_TMs=();
($dg_sum,$hphobmom_dg,$len_dg) = scan_for_best_TM_dGraw(\%putative_TMs, $Length, 0);
for (keys(%putative_TMs)) {
    $helix_dG=sprintf("%.3f",$putative_TMs{$_});
}

$dg_sum = number_print($dg_sum);
$hphobmom_dg = number_print($hphobmom_dg);
$len_dg = number_print($len_dg);
$helix_dG = number_print($helix_dG);

print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=4>&nbsp;</td><td colspan=2 align='center'><b><big>&Sigma;</big>&nbsp;&nbsp;&Delta;G<small><span style='vertical-align: super'>aa(i)</span><span style='vertical-align: sub'>app</span></small> :&nbsp;&nbsp;&nbsp;&nbsp;" . $dg_sum . "</b></td><td>&nbsp;</td>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6><i><b>Hydrophobic moment contribution</b></i></td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td>&Delta;G<small><span style='vertical-align: sub'>hyd.mom.</span></small></td><td colspan=5>=&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>" . $hphobmom_dg . "</b></td></tr>\n";

print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6><i><b>Length contribution</b></i></td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td>&Delta;G<small><span style='vertical-align: sub'>length</span></small></td><td colspan=5>=&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>" . $len_dg . "</b></td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6><table cellpadding=0 cellspacing=0 width='100%' border=0><tr><td class='orangebg'><img src='0.gif' border=0 width=1 height=1></td></tr></table></td></tr>";
print $fpout  "<tr><td colspan=6>&nbsp;</td></tr>\n";
print $fpout  "<tr><td colspan=6><big>&Delta;G<small><span style='vertical-align: super'>pred</span></small><small><span style='vertical-align: sub'>app</span></small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;=&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<big>&Sigma;</big>&nbsp;&nbsp;&Delta;G<small><span style='vertical-align: super'>aa(i)</span><span style='vertical-align: sub'>app</span></small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&Delta;G<small><span style='vertical-align: sub'>hyd.mom.</span></small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&Delta;G<small><span style='vertical-align: sub'>length</span></small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;=&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>$helix_dG</b></big></td></tr>\n";
#print $fpout  "<tr><td><small>$helix</small></td><td>&nbsp;</td><td><small>$helix_dG</small></td>";
#print $fpout  "<td><div onClick=\"javascript:analyze('$helix')\"><small><span style='text-decoration: none;color: #800000;'>[analyze]</span></small></div></td></tr>\n";
print $fpout  "</table><br><br></td></tr>";
print $fpout  "<tr><td><table cellpadding=0 cellspacing=0 width='100%' border=0><tr><td class='orangebg'><img src='0.gif' border=0 width=1 height=1></td></tr></table></td></tr>";
print $fpout  "</table><br>\n";

if ($fpout ne STDOUT)
{
    close($fpout);
}
# --------------------------------------------------------------------------------------------

sub number_print {
    ($nr) = @_;
    $nr = sprintf("%.2f",$nr);
    if ($nr < 0) {
	$ret = "<span style='color: green'>$nr</span>";
    } elsif ($nr > 0) {
	$ret = "<span style='color: red'>+$nr</span>";
    } else {
	$ret = "+-$nr";
    }
    return $ret;
}

sub scan_for_best_TM_dGraw {
    my ($putative_TMs, $seq_length, $print_plot) = @_;
    $print_plot && (open(OUT,">dG_curve.m") || die);
    my $lowest_segment=1000000;
    my ($lowest_k,$lowest_L,$half_step);
    my $dg_sum;
    my $switch=0;
    for (my $L=$Lmin;$L<=$Lmax;$L++) {
	$print_plot && (print OUT "dG_curve_${L} = [\n");
	$flank = int(($L-1)/2+0.5);
	((int($L/2) == ($L/2)) && ($half_step=0.5)) || ($half_step=0); # If L is even, half_steps are introduced
	for (my $k=$flank+1-$half_step;$k<=$seq_length-$flank+$half_step;$k++) {
	    my $segment_dG=0; my $dg=0; $dg_sum=0; my $dg_sin_sum=0; my $dg_cos_sum=0;
	    for (my $i=1;$i<=$L;$i++) {
		$aa = substr($helix,$k-2-$flank+$i+$half_step,1);
		#(!exists($dG{$aa})) && ($aa = "X");   #
		$dg = pos_spec_dG($aa,$i,$L);
		$rel_pos = 9 * (2 * ($i-1)/($L-1) - 1);
		$bgcolor = $switch ? 'greybg' : 'whitebg';
		print $fpout  "<tr><td class='$bgcolor'>$aa</td><td class='$bgcolor'>$i</td><td class='$bgcolor'>" . sprintf("%.2f",$rel_pos) . "</td><td class='$bgcolor'>&nbsp;</td>";
		$dg_print=sprintf("%.2f",$dg);
		if ($dg < 0) {
		    $bar_length = -1 * $totwidth * $dg / 1.9;
		    print $fpout  "<td class='$bgcolor' align='right'>";

		    print $fpout  "<table cellspacing=0 cellpadding=0 border=0 class='$bgcolor'><tr>";
		    print $fpout  "<td><span style='color:green; font-size:75%'>$dg_print&nbsp;&nbsp;&nbsp;</span></td>";
		    print $fpout  "<td valign='center'><img src='green.gif' height=5 width='$bar_length'></td>";
		    print $fpout  "</tr></table>";

		    print $fpout  "</td><td class='$bgcolor'>&nbsp;</td></tr>";
		} elsif ($dg > 0) {
		    $bar_length = $totwidth * $dg / 1.9;
		    print $fpout  "<td class='$bgcolor'>&nbsp;</td>";
		    print $fpout  "<td class='$bgcolor' align='left'>";

		    print $fpout  "<table cellspacing=0 cellpadding=0 border=0 class='$bgcolor'><tr>";
		    print $fpout  "<td valign='center'><img src='red.gif' height=5 width='$bar_length'></td>";
		    print $fpout  "<td><span style='color:red; font-size:75%'>&nbsp;&nbsp;&nbsp;+$dg_print</span></td>";
		    print $fpout  "</tr></table>";

		    print $fpout  "</td></tr>";
		} else {
		    print $fpout  "<td class='$bgcolor' colspan=2>+-0.00</td></tr>";
		}
		$switch = 1 - $switch;
		$dg_sum += $dg;
		$dg_sin_sum += ($dg * sin(100*($i-1)*$piover180));
		$dg_cos_sum += ($dg * cos(100*($i-1)*$piover180));
	    }
	    $hphobmom_dg = $c0*sqrt($dg_sin_sum**2+$dg_cos_sum**2);
	    $segment_dG=$dg_sum+$hphobmom_dg;

	    # Correct for length
	    $len_dg = $c1 + $c2*$L + $c3*$L*$L;
	    $segment_dG += $len_dg;

	    if ($segment_dG < $lowest_segment) {
		$lowest_segment=$segment_dG;
		$lowest_k=$k;
		$lowest_L=$L;
		$lowest_ll=$k-$flank+$half_step;
		$lowest_ul=$lowest_ll+$lowest_L-1;
	    }
	    $print_plot && (print OUT "$k $segment_dG\n");
	}
	$print_plot && (print OUT "];\n");
    }
    $putative_TMs->{"${lowest_ll}_$lowest_ul"}=$lowest_segment;
    $print_plot && close(OUT);
    return($dg_sum,$hphobmom_dg,$len_dg);
}

sub pos_spec_dG {
    my ($aa,$i,$L) = @_;
    my $pos = 9 * (2 * ($i-1)/($L-1) - 1);    
    my $ret;
    if (($aa eq "W") || ($aa eq "Y")) {
	$ret = $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2) 
	    + $profiles{$aa}[2] * (   exp(-1*$profiles{$aa}[3]*($pos-$profiles{$aa}[4])**2) 
                                      + exp(-1*$profiles{$aa}[3]*($pos+$profiles{$aa}[4])**2)   );
    } else {
	$ret = $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2);
    }
    return $ret;
}

sub min {
    my (@nrs) = @_;
    my $min=$nrs[0];
    for (@nrs) { ($_ < $min) && ($min = $_); }
    return $min;
}

sub max {
    my (@nrs)=@_;
    my $max=$nrs[0];
    for (@nrs) { ($_ > $max) && ($max = $_); }
    return $max;
}
sub PrintHelp
{
    print $usage;
}
