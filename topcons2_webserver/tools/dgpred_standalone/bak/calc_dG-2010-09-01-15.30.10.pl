#!/usr/bin/perl -w

#use strict;

my $with_length="on";
my $allow_sub="off";
my $outfile="/dev/stdout";

my $usage="
Usage:   calc_dG.pl [Options] <amino acid fragment file>
Note: the output is in HTML format
Options:
    -lc |--with-length-corr on|off    : whether use length correction, default=$with_length
    -sub|--allow-sub        on|off    : whether allow subsequences, default=$allow_sub
    -o|--outfile            <file>    : output the result to outpath, default=$outfile
    -nc|--not-clean                   : do not clean the temporary dir

Modified from the original perl script by Nanjiang Shu, created 2010-08-27, updated 2010-09-01

Format of the input file: a list of amino acid fragments (assumed to be helices), one line per record
==============================
ASFKSDASD
FJSDFIQFQKLFLSDFA
==============================
Examples:
    calc_dG.pl fragfile.txt -o rst.html
";

my $numArgs = $#ARGV+1;
if($numArgs < 1)
{
    &PrintHelp;
    exit;
}

my $inFile="";
my $isClean=1;

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
            $inFile = $ARGV[$i];
            $i ++;
        }
    }
}#}}}

my @seq=();
open(IN,$inFile) || die "Can not open the input file $inFile\n";
while(<IN>) {
    chomp;
    push(@seq,$_);
}
close(IN);

$fpout=STDOUT;
if ($outfile ne "/dev/stdout")
{
    open($fpout,">", $outfile) || die "Can not open the output file $outfile for write\n";
}

#($fptmp, $tmpfile) = tempfile();

my $max_seq_length=10000;
my $dG_matrix_length=19;
my $c0=0.27045;
my $c1=9.29274167549645;
my $c2=-0.64513139783394;
my $c3=0.00822196628688;

my %biological = (  'A' => [  0.1267255,  0.0215152  ],
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

my %opm = ( 'A'  => [   -0.35564862900004,   0.01332130021828    ],
	    'C'  => [   -0.03501895680136,   0.99999999999436	 ],
	    'D'  => [    1.59956132209071,   0.00947489754213	 ],
	    'E'  => [    1.44456437954037,   0.00991793633818	 ],
	    'F'  => [   -0.58372719807522,   0.01010636364230	 ],
	    'G'  => [   -0.21135828304947,   0.00345594635635	 ],
	    'H'  => [    0.98713634884570,   0.21148435891199 	 ],
	    'I'  => [   -0.53139640795947,   0.02460014333285	 ],
	    'K'  => [    1.65710232354055,   0.01319508610467	 ],
	    'L'  => [   -0.40964877807885,   0.01877474391562	 ],
	    'M'  => [   -0.31236254433492,   0.01063405577742	 ],
	    'N'  => [    0.95392420347656,   0.01345897376904	 ],
	    'P'  => [    0.47894080618646,   0.02083850626139	 ],
	    'Q'  => [    1.13718143293910,   0.01049378516167	 ],
	    'R'  => [    1.45575556617996,   0.01812501297927	 ],
	    'S'  => [    0.16601192335102,   0.00036243232009 	 ],
	    'T'  => [   -0.04699968531037,   0.05979821658957	 ],
	    'V'  => [   -0.45345841058687,   0.02575637636981	 ],
	    'W'  => [   -0.97553693183770,   0.00549871324932,    0.39792640886763,   0.05557263615780,    0.00000000350551 ],                  
	    'X'  => [ 	 0.34521774973501,   0.01890642845475    ],
	    'Y'  => [   -0.47290915514074,   0.01054870058203,    0.51698866596097,   0.06638567043267,   -2.21478619615780 ]  );                 
  
my %kdo      = (  'A' => [ -1.8, 0 ],	      
	     	  'B' => [  3.5, 0 ],	      
                  'C' => [ -2.5, 0 ],	      
	          'D' => [  3.5, 0 ],	      
   	          'E' => [  3.5, 0 ],	      
	          'F' => [ -2.8, 0 ],	      
	          'G' => [  0.4, 0 ],	      
	          'H' => [  3.2, 0 ],	      
	          'I' => [ -4.5, 0 ],	      
	          'K' => [  3.9, 0 ],	      
	          'L' => [ -3.8, 0 ],	      
	          'M' => [ -1.9, 0 ],	      
   	          'N' => [  3.5, 0 ],	      
	          'P' => [  1.6, 0 ],	      
	          'Q' => [  3.5, 0 ],	      
	          'R' => [  4.5, 0 ],	      
	          'S' => [  0.8, 0 ],	      
	          'T' => [  0.7, 0 ],	      
	          'U' => [ -1.9, 0 ],	      
	          'V' => [ -4.2, 0 ],	      
	          'W' => [  0.9, 0, 0, 0, 0 ],
	          'X' => [  0.5, 0 ],                                           
	          'Y' => [  1.3, 0, 0, 0, 0 ],
	          'Z' => [  3.5   , 0 ]  );   

my %zhl = (  'A' => [0.38 , 0 ],	      
	     'B' => [2.45 , 0 ],	      
	     'C' => [0.30 , 0 ],	      
	     'D' => [3.27 , 0 ],	      
	     'E' => [2.90 , 0 ],	      
	     'F' => [1.98 , 0 ],	      
	     'G' => [0.19 , 0 ],	      
	     'H' => [1.44 , 0 ],	      
	     'I' => [1.97 , 0 ],	      
	     'K' => [3.46 , 0 ],	      
	     'L' => [1.82 , 0 ],	      
	     'M' => [1.40 , 0 ],	      
	     'N' => [1.62 , 0 ],	      
	     'P' => [1.44 , 0 ],	      
	     'Q' => [1.84 , 0 ],	      
	     'R' => [2.57 , 0 ],	      
	     'S' => [0.53 , 0 ],	      
	     'T' => [0.32 , 0 ],	      
	     'U' => [1.40 , 0 ],	      
	     'V' => [1.46 , 0 ],	      
	     'W' => [1.53 , 0, 0, 0, 0 ],
	     'X' => [1.55 , 0 ],         
	     'Y' => [0.49 , 0, 0, 0, 0 ],
	     'Z' => [2,37 , 0 ]  );   

my %profiles = %biological;


unless ($with_length eq "on") {
    $c1=0; $c2=0; $c3=0;
}

if ($allow_sub eq "on") {
    $allow_sub=1;
} else {
    $allow_sub=0;
}

my %aa2nr = ( "A" => 1,"C" => 2,"D" => 3,"E" => 4,"F" => 5,"G" => 6,"H" => 7,"I" => 8,"K" => 9,"L" => 10,
	      "M" => 11,"N" => 12,"P" => 13,"Q" => 14,"R" => 15,"S" => 16,"T" => 17,"V" => 18,"W" => 19,"Y" => 20);
my $piover180=atan2(1,1)/45;

print $fpout  "<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\" width=\"100%\" class=\"TableGrey\">\n";
#print $fpout  "<tr><td rowspan=100><img src='0.gif' border=0 width=30 height=1></td></tr>\n";
#print $fpout  "<tr><td colspan=4><img src=\"0.gif\" border=0 width=750 height=30></td></tr>\n";
#print $fpout  "<tr><td align=\"left\" colspan=4><big><big>Results</big></big><br><br><br><br><br></td></tr>\n";

unless (scalar(@seq)) {
    print $fpout  "<tr><td>No sequence entered!</td></tr></table><br><br><br><br>\n";
    exit;
} elsif (scalar(@seq) > $max_seq_length) {
    print $fpout  "<tr><td>Sequence too long (> $max_seq_length aa.)</td></tr></table><br><br><br><br>\n";
    exit;
}    

print $fpout  "<tr><td colspan=4 align=left><table width=\"100%\" cellspacing=0 cellpadding=0 border=0>\n";
if ($allow_sub) {
    print $fpout  "<tr><td valign='top'><i>Sequence</i><br><small>(Subsequence with lowest &Delta;G marked in <span style='color:#FF8800'>orange</span>)</small></td><td>&nbsp;</td><td valign='top'><i>Predicted&nbsp;&Delta;G</i></td><td>&nbsp;</td></tr>\n";
    print $fpout  "<tr><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n";
} else {
    print $fpout  "<tr><td><i>Sequence</i></td><td><img src='0.gif' height=1 width=50></td><td><i>Predicted&nbsp;&Delta;G</i></td><td>&nbsp;</td></tr>\n";
    print $fpout  "<tr><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n";
}

$count=0;
while($helix=shift(@seq)) {
    $count++;
    if ($count > 1000) {
	last;
    }
    ($helix =~ /^\s+$/) && (print $fpout  "<tr><td colspan=4>&nbsp;</td></tr>\n") && next;
    $helix =~ s/\s*$//;
    $helix =~ s/^\s*//;
    $helix =~ tr/a-z/A-Z/;
    unless ($helix =~ /^[ABCDEFGHIKLMNPQRSTUVWXYZ]+$/) {
	print $fpout  "<tr><td>$helix</td><td>&nbsp;</td><td colspan=2>Sequence error. </td></tr>\n";
	next;
    }
    $Length=length($helix);
    if ($Length > 40) {
	print $fpout  "<tr><td>" . substr($helix,0,40) . "... </td><td>&nbsp;</td><td colspan=2>Length error (> 40 residues)</td></tr>\n";
	next;
    }
    if ($Length < 9) {
	print $fpout  "<tr><td>$helix</td><td>&nbsp;</td><td colspan=2>Length error (< 9 residues)</td></tr>\n";
	next;
    }

    if ($allow_sub) {
	$Lmin=min(9,$Length); $Lmax=$Length;
    } else {
	$Lmin=$Length; $Lmax=$Length;
    }
    
    %putative_TMs=();
    scan_for_best_TM_dGraw(\%putative_TMs, $Length, 0);

    for (keys(%putative_TMs)) {
	$helix_dG=sprintf("%.3f",$putative_TMs{$_});
	if ($allow_sub) {
	    ($ll,$ul)=split(/_/,$_);
	    $helix = substr($helix,0,$ll-1) . "<span style='color: #FF8800'>" . substr($helix,$ll-1,$ul-$ll+1) . "</span>" . substr($helix,$ul);
	}
    }
    print $fpout  "<tr><td>$helix</td><td>&nbsp;</td><td>$helix_dG</td>";
    $helix =~ s/^.*?\>//; $helix =~ s/<.*$//;
    print $fpout  "<td><a href=\"analyze.php?with_length=$with_length&seq=$helix\" target=\"_blank\" style='text-decoration:none' onClick=\"javascript:analyze('$helix');return false;\"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>\n";
}
print $fpout  "</table><br><br><br><br></td></tr></table>\n";
if ($fpout ne STDOUT)
{
    close($fpout);
}

# $number = `cat count`;
# $number =~ /^\d+$/ || die;
# $number += $count;
# system("echo $number > count");

# --------------------------------------------------------------------------------------------

sub scan_for_best_TM_dGraw {
    my ($putative_TMs, $seq_length, $print_plot) = @_;
    $print_plot && (open(OUT,">dG_curve.m") || die);
    my $lowest_segment=1000000;
    my ($lowest_k,$lowest_L,$half_step);
    for (my $L=$Lmin;$L<=$Lmax;$L++) {
	$print_plot && (print OUT "dG_curve_${L} = [\n");
	$flank = int(($L-1)/2+0.5);
	((int($L/2) == ($L/2)) && ($half_step=0.5)) || ($half_step=0); # If L is even, half_steps are introduced
	for (my $k=$flank+1-$half_step;$k<=$seq_length-$flank+$half_step;$k++) {
	    my $segment_dG=0; my $dg=0; my $dg_sum=0; my $dg_sin_sum=0; my $dg_cos_sum=0;
	    for (my $i=1;$i<=$L;$i++) {
		$aa = substr($helix,$k-2-$flank+$i+$half_step,1);
		#(!exists($dG{$aa})) && ($aa = "X");   #
		$dg = pos_spec_dG($aa,$i,$L);
		$dg_sum += $dg;
		$dg_sin_sum += ($dg * sin(100*($i-1)*$piover180));
		$dg_cos_sum += ($dg * cos(100*($i-1)*$piover180));
	    }
	    $segment_dG=$dg_sum+$c0*sqrt($dg_sin_sum**2+$dg_cos_sum**2);

	    # Correct for length
	    $segment_dG += $c1 + $c2*$L + $c3*$L*$L;

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
