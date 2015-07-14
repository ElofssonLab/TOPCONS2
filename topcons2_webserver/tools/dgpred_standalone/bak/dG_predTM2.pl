#!/usr/bin/perl -w

my $max_seq_length=10000;

$gnuplot="/var/www/dgpred.cbr.su.se/gnuplot/bin/gnuplot";

($tmpdir=$ARGV[0]) || die;
($tmpdir =~ /^\/tmp\/(\w*)$/) || die;
$tmpname=$1;
open(IN,"$tmpdir/options") || die;
while(<IN>) {
    /^with_length: (.*)/ && ($with_length=$1);
    /^graphics: (.*)/ && ($graphics=$1);
    /^Lmin: (.*)/ && ($Lmin=$1);
    /^Lmax: (.*)/ && ($Lmax=$1);
}
close(IN);
#open(IN,"$tmpdir/sequences") || die;
#while(<IN>) {
#    chomp;
#    push(@seq,$_);
#}
#close(IN);
(($graphics eq "on") && ($print_plot=1)) || ($print_plot=0);

(($Lmin =~ /^\d+$/) && ($Lmax =~ /^\d+$/)) || dieHTML("Something is wrong with the minimum and maximum length.\n");
#(ceil($Lmin/2) == ($Lmin/2)) && ($Lmin++);
#(ceil($Lmax/2) == ($Lmax/2)) && ($Lmax--);
($Lmin < 9) && ($Lmin=9);
($Lmax > 40) && ($Lmax=40);
($Lmin <= $Lmax) || dieHTML("Something is wrong with the minimum and maximum length. (max length must be >= min length and max length must be >=9)\n");

#print "L:$in{'withlength'} G:$in{'graphics'} F:$in{'withflanks'}<br>";
$infile="$tmpdir/sequences";
#$suffix="symm.k_0_lt1p0.expfit.noPolyA";
#$suffix="withRscan.fixed.47par";
#$suffix="GES";

my $dG_matrix_length = 19;

#my $dG_file;
#if ($invivo_invitro eq "invivo") {
#    $dG_file = "pos_spec_dG/pos_spec_dG.symm.k_0_lt1p0.direct.withRscan.final.invivo";
#    $phi = 0.49; $xi=0.049;
#} else {
#    $dG_file = "pos_spec_dG/pos_spec_dG.$suffix";
#$charged_DDG_file="flanking_KR_matrix_additive";

#$c1=9.74843852799261; $c2=-0.69096610396195; $c3=0.00936651460524;
#$c1=9.81504713559865; $c2=-0.69588262140385; $c3=0.00944060947990;

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

my %kdo = (  'A' => [ -1.8, 0 ],
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

#    $c1=9.42399201351392; $c2=-0.65041775608542; $c3=0.00830540592606; 
#    $pa=0.00831; $pb=-0.643; $pc=9.22;
#    $phi = 0.35; $xi=0.044;
#}

#($with_phi eq "on") || ($phi=0);
#($with_xi eq "on") || ($xi=0);
unless ($with_length eq "on") {
    $c1=0; $c2=0; $c3=0;
}

#my $dG_file = "pos_spec_dG/pos_spec_dG.$suffix";
#my $phi=0.54; my $xi=0.056;
#my $phi=0; my $xi=0;

my $maxDG=3.0;

# -----------------------------------------------------------------------------------------------------------------


my (%dG, %charged_DDG, %nr2aa, %list, $aa, $name, $seq, $segment_dG, $segment_dG1, $segment_dG2, $flank);
my %aa2nr = ( "A" => 1,"C" => 2,"D" => 3,"E" => 4,"F" => 5,"G" => 6,"H" => 7,"I" => 8,"K" => 9,"L" => 10,
	      "M" => 11,"N" => 12,"P" => 13,"Q" => 14,"R" => 15,"S" => 16,"T" => 17,"V" => 18,"W" => 19,"Y" => 20);
while( my ($key, $value) = each(%aa2nr) ) { $nr2aa{$value} = $key; }
my $piover180=atan2(1,1)/45;
my $c;
my $plot_cmd; 
my %number=();
my %occupied_names=();
my $global_max_dG; my $global_min_dG;
#read_dG_matrix();
#read_charged_DDG_matrix_additive();
#system("rm plots/plot* plots/*.plot");
#print_dG();
open(FA,$infile) || dieHTML("Couldn't open input file.\n");
$nr_sequences=0;
$nr_helices=0;
$prev_number=0;
mkdir("$tmpdir/plots") || dieHTML("Could not create temporary directory2 $tmpdir/plots.\n");
(-e "plots") || mkdir("plots");
mkdir("plots/$tmpname",0777) || dieHTML("Could not create temporary directory2 plots/$tmpname.\n");
my $break=1;
while (($name,$seq) = read_one_entry(\%occupied_names)) {
    $break=0;
#    print "$name<br>\n";
    last if ($nr_sequences > 10);
    @order=(@order,$name);
    my $seq_length=length($seq);
    my %putative_TMs = ();
    my %TMs = ();
    my %TM_list = ();
    $global_max_dG = -1000; $global_min_dG = 1000;
    $seq =~ tr/a-z/A-Z/;
    $seq_length=length($seq);
    if (($seq_length >= $Lmin) && ($seq_length <= $max_seq_length)) {
	if ($seq =~ /^[ABCDEFGHIKLMNPQRSTUVWXYZ]+$/) {
	    scan_for_putative_TMs_dGraw(\%putative_TMs,$seq_length,$print_plot);
	    resolve_TMs(\%putative_TMs,\%TMs,\%TM_list);
	    mk_plot($name,\%number,$seq_length,\%TM_list,$tmpdir);
	    clean_plot_files($tmpdir);
	} else {
	    $TMs{"0-0"}="seq_err";
	}
    } else {
	$TMs{"0-0"}="length_err";
    }	
    for (keys %TMs) {
	($ll,$ul)=split(/-/,$_);
	$list{$name}{$_}{'dG'} = $TMs{$_};
	$list{$name}{$_}{'sort'} = $ll;
	$list{$name}{$_}{'seq'} = substr($seq,$ll-1,$ul-$ll+1);
	$nr_helices++;
    }
    $nr_sequences++;

#    print_TMs($name,\%TMs);
}
rmdir("$tmpdir/plots");
close(FA);
print "<table cellspacing=\"0\" cellpadding=\"0\" border=\"0\" width=\"720\" class=\"TableGrey\">\n";

if ($break) {
    print "<tr><td>No sequence entered!</td></tr></table><br><br><br><br>\n";
    exit;
}


print "<tr><td align=\"left\">";
print_out_list(\%list, \%number, @order);
print "</td></tr></table>\n";

#print "<br><br></body>\n<head>\n<meta http-equiv=\"pragma\" content=\"nocache\">\n<META HTTP-EQUIV=\"Expires\" CONTENT=\"-1\">\n</head></html>\n";

close(OCCUPIED);

$number = `cat count`;
$number =~ /^\d+$/ || die;
$number += $nr_helices;
system("echo $number > count");

exit;

# -----------------------------------------------------------------------------------------------------------------

sub print_out_list {
    my ($list, $number, @order) = @_;
    for $name (@order) {
	print "<table cellspacing=0 cellpadding=5 border=0 class=\"TableWhite\" width='100%'>\n";
	print "<tr><td colspan=3><big>&gt;$name</big></td></tr>\n";
	$exists_helices=0;
	for (keys(%{$list{$name}})) {
	    if ($list{$name}{$_}{'dG'} < $maxDG) {
		$exists_helices=1;
		last;
	    }
	}
	if ($exists_helices == 0) {
	    print "<tr><td colspan=2 valign=\"top\"><img src=\"0.gif\" height=1 width=200><br><div style=\"color:red\">No predicted TM helices</div></td></tr><tr><td colspan=2>\n";
	} elsif (exists($list{$name}{"0-0"})) {
	    if ($list{$name}{"0-0"}{'dG'} eq "seq_err") {
		print "<tr><td>Sequence error.</td><td>&nbsp;&nbsp;&nbsp;</td><td>&nbsp;</td></tr>\n";
	    } elsif ($list{$name}{"0-0"}{'dG'} eq "length_err") {
		print "<tr><td>Length error. (Sequence shorter than \"Helix min length\" or longer than $max_seq_length aa.)</td><td>&nbsp;&nbsp;&nbsp;</td><td>&nbsp;</td></tr>\n";
	    }		
	} else {	
	    print "<tr><td colspan=2 valign=\"top\"><img src=\"0.gif\" height=1 width=200><br><table cellspacing=0 cellpadding=0 border=0><tr><td colspan=9><i><b>Predicted TM helices:</b><br></i></td></tr><tr><td align=center><i>Position</i></td><td>&nbsp;&nbsp;&nbsp;</td><td align=center><i>Length</i></td><td align=center>&nbsp;&nbsp;&nbsp;</td><td><i>Predicted &Delta;G</i></td><td>&nbsp;&nbsp;&nbsp;</td><td align=left><i>Sequence</i></td><td><img src='0.gif' height=1 width=100></td><td align=left>&nbsp;</td></tr>";
	    for (sort {$list{$name}{$a}{'sort'} <=> $list{$name}{$b}{'sort'}} keys %{$list{$name}}) {
		($list{$name}{$_}{'dG'} > $maxDG) && next;
		($ll,$ul)=split(/-/,$_); $length=$ul-$ll+1;
		$list{$name}{$_}{'dG'}=sprintf("%.3f",$list{$name}{$_}{'dG'});
		print "<tr><td align=center>$_</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>$length</td><td>&nbsp;&nbsp;&nbsp;</td><td align=center>$list{$name}{$_}{'dG'}</td><td>&nbsp;</td><td>$list{$name}{$_}{'seq'}</td><td>&nbsp;</td><td><a href=\"analyze.php?with_length=$with_length&seq=$list{$name}{$_}{'seq'}\" target=\"_blank\" style='text-decoration:none'  onClick=\"javascript:analyze('$list{$name}{$_}{'seq'}');return false;\"><span style='text-decoration: none;color: #800000;'>[analyze]</span></a></td></tr>\n";
	    }
	    print "</table><br>";
	}
	if ($print_plot && !exists($list{$name}{"0-0"})) {
	    print "&nbsp;&nbsp;&nbsp;<img src=\"plots/$tmpname/plot_$number{$name}.gif\" border=\"0\"></td>";
	} else {
	    print "</td>\n";
	}
	print "</tr></table><br><br><br>\n";
    }
}

sub mk_plot {
    ($name,$number,$seq_length,$TM_list,$tmpdir)=@_;

    open(OUT,">$tmpdir/plots/helices.list") || die;
    $exists_helices=0;
    for $L_k (keys %{$TM_list}) {
	($TM_list->{$L_k} > $maxDG) && next;
	$exists_helices=1;
	($L,$k)=split(/_/,$L_k);
#	$L=$L/2;
	print OUT "$k 100 $L\n$k -100 $L\n";
    }
    close(OUT);

#    $prev_number=`cat img_number`; chomp($prev_number);
    $number{$name}=$prev_number+1;
    $prev_number=$number{$name};
#    system("echo $number{$name} > img_number");

    $x="position"; $y="Delta-G";
    open(OUT,">$tmpdir/plots/gnuplot.command") || die;
    $plot_cmd =~ s/,\s*$//;
    print OUT "set rmargin 13\n";
    print OUT "set key outside\n";
    print OUT "set yrange [" . ($global_min_dG-1) . ":" . ($global_max_dG+1) . "] reverse\n";
    print OUT "set ylabel \"Predicted {/Symbol D}G (kcal/mol)\"\n";
    print OUT "set xlabel \"Sequence position (of center of helix)\"\n";
#    print OUT "set ter postscript enhanced color\n";
#    print OUT "set out \"$tmpdir/plots/plot.ps\"\n";
    print OUT "set ter png enh interlace size 2400,1680 font 'Nimbus,40'\n";
    print OUT "set out \"$tmpdir/plots/plot.png\"\n";
    if ($exists_helices) {
	print OUT "set style line 1 lc rgb \"#DDDDDD\"\n";
	print OUT "plot [0:$seq_length] \"$tmpdir/plots/helices.list\" w boxes ls 1 fill solid noborder notitle, 0 notitle w l 0, $plot_cmd \n";
    } else {
	print OUT "plot [0:$seq_length] 0 notitle w l 0, $plot_cmd \n";
    }
    close(OUT);
    $ENV{'GDFONTPATH'}='./fonts/';
    system("$gnuplot $tmpdir/plots/gnuplot.command");
#     open(IN,"$tmpdir/plots/plot.ps") || die;
#     open(OUT,">$tmpdir/plots/plot_fixed.ps") || die;
#     while(<IN>) {
# 	if (/^\/LT0\s/) {
# 	    print OUT "/LT0 { PL [] 0.9 0.9 0.9 DL } def\n";
# 	} else {
# 	    print OUT $_;
# 	}
#     }
#     close(IN);
#     close(OUT);

    `convert -scale 600x420 $tmpdir/plots/plot.png $tmpdir/plots/plot.2.png`;
    `mv $tmpdir/plots/plot.2.png $tmpdir/plots/plot.png`;
##    system("./ps2bmp -f ppm -s 0.9 $tmpdir/plots/plot_fixed.ps $tmpdir/plots/plot.bmp >/dev/null");
##    system("pamflip -cw $tmpdir/plots/plot.bmp > $tmpdir/plots/plot.ppm");
##    system("ppmtogif $tmpdir/plots/plot.ppm > $tmpdir/plots/plot.gif");

#    system("$tmpdir/plots/ps2gif >/dev/null");

##    system("mv $tmpdir/plots/plot.gif plots/$tmpname/plot_$number{$name}.gif");
    system("mv $tmpdir/plots/plot.png plots/$tmpname/plot_$number{$name}.gif");
    system("chmod 0777 plots/$tmpname");
    system("chmod 0777 plots/$tmpname/plot_$number{$name}.gif");
    $plot_number++;
    $plot_cmd="";
}

sub clean_plot_files {
    ($tmpdir) = @_;
#    system("rm $tmpdir/plots/plot.ps $tmpdir/plots/plot.bmp $tmpdir/plots/plot.ppm $tmpdir/plots/dG_curve_*.plot $tmpdir/plots/gnuplot.command $tmpdir/plots/helices.list");
}

sub scan_for_putative_TMs_dGraw {
    my ($putative_TMs, $seq_length, $print_plot) = @_;
    for (my $L=$Lmin;$L<=$Lmax;$L++) {
	($L > $seq_length) && next;
	$print_plot && (open(OUT,">$tmpdir/plots/dG_curve_$L.plot") || die);
	$flank = int(($L-1)/2+0.5);
	$nr_points=0;
	((int($L/2) == ($L/2)) && ($half_step=0.5)) || ($half_step=0); # If L is odd, half_steps are introduced
	for (my $k=$flank+1-$half_step;$k<=$seq_length-$flank+$half_step;$k++) {
	    my $segment_dG=0; my $dg=0; my $dg_sum=0; my $dg_sin_sum=0; my $dg_cos_sum=0;
	    for (my $i=1;$i<=$L;$i++) {
		$aa = substr($seq,$k-2-$flank+$i+$half_step,1);
		$dg = pos_spec_dG($aa,$i,$L);
		$dg_sum += $dg;
		$dg_sin_sum += ($dg * sin(100*($i-1)*$piover180));
		$dg_cos_sum += ($dg * cos(100*($i-1)*$piover180));
	    }
	    $segment_dG=$dg_sum+$c0*sqrt($dg_sin_sum**2+$dg_cos_sum**2);
		
	    # Add delta-delta-G for charged flanking residues upstreams (C-terminal/Cytoplasmic side) of segment
	        my $upstr_seq_length=min(20,$k-$flank-1);
	        my $upstreams_sequence=reverse(substr($seq,$k-$flank-$upstr_seq_length-1,$upstr_seq_length));
	        $segment_dG += charged_DDG_additive($upstreams_sequence);
	    
	    # Correct for length
#	    $segment_dG += (-1*$phi*$L + 19*$phi) + $xi * (max(0,$L-19)**2);
	    $segment_dG += $c1 + $c2*$L + $c3*$L*$L;
	    
	    $putative_TMs->{"${L}_$k"}=$segment_dG;
	    $print_plot && (print OUT "$k $segment_dG\n");
	    if ($segment_dG > $global_max_dG) {
		$global_max_dG=$segment_dG;
	    }
	    if ($segment_dG < $global_min_dG) {
		$global_min_dG=$segment_dG;
	    }
	    $nr_points++;
	}
	if ($print_plot) {
	    if ($nr_points > 1) {
		$plot_cmd.="\"$tmpdir/plots/dG_curve_$L.plot\" w l lw 6 t \"L=$L\", ";
	    } else {
		$plot_cmd.="\"$tmpdir/plots/dG_curve_$L.plot\" t \"L=$L\", ";
	    }
	    close(OUT);
	}
    }
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

sub pos_spec_dG_old {
    my ($aa,$i,$L) = @_;
    my $dG_profile_exact_position=($i-1)/($L-1)*($dG_matrix_length-1)+1;
    my $lower_pos=floor($dG_profile_exact_position);
    my $upper_pos=ceil($dG_profile_exact_position);
    my $ret;
    if ($lower_pos == $upper_pos) {
	$ret = $dG{$aa}{$dG_profile_exact_position};
    } else {
	my $lower_fraction=$upper_pos-$dG_profile_exact_position;
	$ret=$lower_fraction*$dG{$aa}{$lower_pos}+(1-$lower_fraction)*$dG{$aa}{$upper_pos};
    }
    return $ret;
}

sub resolve_TMs {
    my ($putative_TMs,$TMs,$TM_list) = @_;
    for (sort {$putative_TMs->{$a} <=> $putative_TMs->{$b}} keys %{$putative_TMs}) {
	my ($L,$k)=split(/_/,$_);
	my $put_ll = $k - ($L-1)/2;
	my $put_ul = $k + ($L-1)/2;
	my $as=0;
	for my $rgn (keys %{$TMs}) {
	    my ($ll,$ul) = split(/\-/,$rgn);
	    unless (($put_ll > $ul) || ($put_ul < $ll)) {
		$as=1;
		last;
	    }
	}
	unless ($as) {
	    $TMs->{"$put_ll-$put_ul"} = $putative_TMs->{"${L}_$k"};
	    $TM_list->{"${L}_$k"}=$putative_TMs->{"${L}_$k"};
	}
    }
}

sub print_TMs {
    my ($name,$TMs) = @_;
    print "$name ";
    for (sort {$TMs->{$a} <=> $TMs->{$b}} keys %{$TMs}) {
        my ($ll,$ul)=split /-/;
	print " ".($ul-$ll+1).":$ll-$ul $TMs->{$_} | ";
    }
    print "\n";

}

sub read_dG_matrix {
    open(DG,$dG_file) || dieHTML("Couldn't find dG matrix file $dG_file.\n");
    my $aa_nr=0;
    my ($aa, @dG_arr);
    $c=<DG>; $c =~ s/\s.*//;
    while(<DG>) {
	chomp;
	@dG_arr=split(/\s+/,$_);
	$aa_nr++;
	$aa=$nr2aa{$aa_nr};
	(scalar(@dG_arr) != $dG_matrix_length) && dieHTML("Error in dG matrix file.\n".scalar(@dG_arr)."\n");
	for (my $i=1;$i<=$dG_matrix_length;$i++) {
	    $dG{$aa}{$i}=shift(@dG_arr);
	}
    }
    close(DG);
}

sub charged_DDG_additive {
    my ($down_stream_seq)=@_;
    $down_stream_seq =~ s/[^KR]/-/g;
    my $nr_charges=($down_stream_seq =~ s/[KR]/+/g);
    ($nr_charges == 0) && (return 0);
    my $distance=0; my $DDG=0;
    while(my $chr=substr($down_stream_seq,0,1,"")) {
	if ($chr eq "+") {
	    $DDG += $charged_DDG{$distance};
	}
	$distance++;
    }
    unless ($with_flanks eq "on") {
	$DDG=0;
    }
    return $DDG;
}

sub read_charged_DDG_matrix_additive {
    #This is a vector rather than a matrix
    my $distance;
    open(DDG,$charged_DDG_file) || die("Couldn't find DDG matrix file.\n");
    my $line=<DDG>;
    close(DDG);
    $distance=0;
    for (split(/\s/,$line)) {
	$charged_DDG{$distance++}=$_;
    }
}

sub print_dG {
    for (sort keys %dG) {
	$aa=$_;
	print "$aa ";
	for (sort {$a <=> $b} keys %{$dG{$aa}}) {
	    print "$dG{$aa}{$_} ";
	}
	print "<br>\n";
    }
}

sub read_one_entry {
    ($occupied_names)=@_;
    my ($line, $seq);
    defined(my $name = <FA>) || return;
    while ($name =~ /^\s+$/) { $name=<FA>; }
    if ($name =~ /^\s*>(.*)$/) {
	$name=$1;
	$name =~ s/^\s*//; $name =~ s/\s*$//;
    } elsif ($name =~ /^\s*([A-Z]+)\s*$/) {
	$seq=$1;
	$name="***Unnamed sequence";
    } else {
	dieHTML("Error in FASTA format: $name");
    }
#    $name =~ /^\s*>(.*)$/ || dieHTML("Error in FASTA format: $name");
#    $name=$1;
    while($line=<FA>) {
	if ($line =~ /^>/) {
	    seek FA, -1*length($line), 1;
	    last;
	}
	if ($line =~ /^\s*$/) {
	    last;
	}
	$seq .= $line;
    }
    ($seq && ($seq =~ s/\s//g));
    if ($occupied_names{$name}) {
	$i=2;
	while($occupied_names{"$name-$i"}) {$i++;}
	$name="$name-$i";
    }
    $occupied_names{$name}=1;
    return ($name,$seq);
}

sub dieHTML {
    ($str)=@_;
    print "$str<br><br><br><br>\n";
    exit(0);
}

#------------------------------------------------------------------------------------------------------------
sub scan_for_putative_TMs_dGraw_both_dir {                   
    my ($putative_TMs, $seq_length, $print_plot) = @_;
    $print_plot && (open(OUT,">dG_curve.m") || die);
    for (my $L=$Lmin;$L<=$Lmax;$L+=2) {
	$print_plot && (print OUT "dG_curve_${L} = [\n");
	$flank = ($L-1)/2;
	for (my $k=$flank+1;$k<=$seq_length-$flank;$k++) {
	    my $segment_dG1=0; my $segment_dG2=0; my $dg1=0; my $dg2=0; my $dg_sum1=0; my $dg_sum2=0; my $dg_sin_sum1=0; my $dg_sin_sum2=0; my $dg_cos_sum1=0; my $dg_cos_sum2=0;
	    for (my $i=1;$i<=$L;$i++) {
		$aa = substr($seq,$k-2-$flank+$i,1);
		$dg1 = pos_spec_dG($aa,$i,$L);
		$dg2 = pos_spec_dG($aa,$L-$i+1,$L);
		$dg_sum1 += $dg1;
		$dg_sum2 += $dg2;
		$dg_sin_sum1 += ($dg1 * sin(100*($i-1)*$piover180));
		$dg_sin_sum2 += ($dg2 * sin(100*($i-1)*$piover180));
		$dg_cos_sum1 += ($dg1 * cos(100*($i-1)*$piover180));
		$dg_cos_sum2 += ($dg2 * cos(100*($i-1)*$piover180));
	    }
	    $segment_dG1=$dg_sum1+$c*sqrt($dg_sin_sum1**2+$dg_cos_sum1**2);
	    $segment_dG2=$dg_sum2+$c*sqrt($dg_sin_sum2**2+$dg_cos_sum2**2);
	    $segment_dG=min($segment_dG1,$segment_dG2);
#	    print "$segment_dG1 $segment_dG2 ==> $segment_dG\n";
	    $putative_TMs->{"${L}_$k"}=$segment_dG;
	    $print_plot && (print OUT "$k $segment_dG\n");
	}
	$print_plot && (print OUT "];\n");
    }
    $print_plot && close(OUT);
}

