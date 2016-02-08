#!/usr/bin/perl -W
use strict;
use warnings;

my $usage = "
USAGE: $0 DATADIR

Examples:
    $0 rst/seq_0
";
my $numArgs = $#ARGV+1;
if($numArgs < 1) {
    print "$usage\n";
    exit 1;
}

my $gnuplot = "gnuplot";
my $convert = "convert";
my $sort = "sort";
my $awk = "awk";
my $head="head";

$ENV{'GDFONTPATH'}='./fonts/';
my $res_folder = $ARGV[0] . "/";  #provide this from the command line, contains all sub-folders for various prediction methods

if (! -d $res_folder){
    print STDERR "DATADIR $res_folder does not exist. EXIT!\n";
    exit 1;
}

my $seqfile = $res_folder."seq.fa";
my $octopus_res_file = $res_folder."OCTOPUS/query.top";
my $philius_res_file = $res_folder."philius/query.top";
my $polyphobius_res_file = $res_folder."PolyPhobius/query.top";
my $scampi_res_file = $res_folder."SCAMPI_MSA/query.top";
my $spoctopus_res_file = $res_folder."SPOCTOPUS/query.top";
my $topcons_res_file = $res_folder."Topcons/topcons.top";
my $homology_res_file = $res_folder."Homology/query.top";
my $DG_res_file = $res_folder."Topcons/topcons.DG";

my $seq = "";
my $seqid = "";
my $seqanno = "";
my $pdbchain_homo = "";
my $pdbcode_homo = "";

if (! -e $seqfile){
    warn "$seqfile does not exist";
}else{
    ($seqid, $seqanno, $seq) = ReadSingleFastaSeq($seqfile);
}


my %topo=();

sub ReadSeq#{{{
{
    my $file=shift;
    my $seq=`grep -v '>' $file`;
    $seq=~s/\s//g;
    $seq=~s/\n//g;
    return $seq;
}#}}}
sub ReadSingleFastaSeq{#{{{
    my $infile=shift;
    open IN, $infile;
    my $seqid = "";
    my $anno = "";
    my $seq = "";
    while(<IN>) {
        if($_ =~ /^>/) {
            $anno = $_;
            $anno =~ s/^>+//;
            if ($anno ne ""){
                my @fields=split(/\s+/,$anno);
                $seqid = $fields[0];
            }
        }else{
            my $res = $_;
            $res=~s/\s//g;
            $res=~s/\n//g;
            $seq .= $res;
        }
    }
    close IN;
    return ($seqid, $anno, $seq);
}#}}}

$topo{'OCTOPUS'} = "";
$topo{'Philius'} = "";
$topo{'CONSENSUS_TOPCONS'} = "";
$topo{'SPOCTOPUS'} = "";
$topo{'SCAMPI'} = "";
$topo{'PolyPhobius'} = "";
$topo{'VHomology'} = "";
if(-e $octopus_res_file) {
    ($seqid, $seqanno, $topo{'OCTOPUS'}) = ReadSingleFastaSeq($octopus_res_file);
}
if(-e $philius_res_file) {
    ($seqid, $seqanno, $topo{'Philius'}) = ReadSingleFastaSeq($philius_res_file);
}
if(-e $polyphobius_res_file) {
    ($seqid, $seqanno, $topo{'PolyPhobius'}) = ReadSingleFastaSeq($polyphobius_res_file);
}
if(-e $scampi_res_file) {
     ($seqid, $seqanno, $topo{'SCAMPI'}) = ReadSingleFastaSeq($scampi_res_file);
}
if(-e $spoctopus_res_file) {
    ($seqid, $seqanno, $topo{'SPOCTOPUS'}) = ReadSingleFastaSeq($spoctopus_res_file);
}
if(-e $homology_res_file) {
    ($seqid, $seqanno, $topo{'VHomology'}) = ReadSingleFastaSeq($homology_res_file);
    if ($topo{'VHomology'} =~ /Nohomologoushitsfound/){
        $seqid = "";
        $seqanno = "";
        $topo{'VHomology'} = "";
    }
    if ($seqid ne ""){
        $pdbchain_homo = $seqid;
        $pdbcode_homo = substr $pdbchain_homo, 0, 4;
        $pdbcode_homo = lc $pdbcode_homo;
    }
}

my $showtext_homo = "";
my $url_homo = "";
if ($pdbchain_homo ne ""){
    $showtext_homo = $pdbchain_homo;
    $url_homo = "http://www.rcsb.org/pdb/explore/explore.do?structureId=$pdbcode_homo";
}else{
    $showtext_homo = "PDB-homology";
}
# if(-e $octopus_res_file)
# {
#     open IN, $octopus_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $octopus_res=$_;
#             chomp $octopus_res;
#             $topo{'OCTOPUS'} = $octopus_res;
#         }
#     }
#     close IN;
# }
# 
# if(-e $philius_res_file) {
#     open IN, $philius_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $philius_res=$_;
#             chomp $philius_res;
#             $topo{'Philius'} = $philius_res;
#         }
#     }
#     close IN;
# }
# 
# 
# if(-e $polyphobius_res_file) {
#     open IN, $polyphobius_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $polyphobius_res=$_;
#             chomp $polyphobius_res;
#             $topo{'PolyPhobius'} = $polyphobius_res;
#         }
#     }
#     close IN;
# }
# 
# if(-e $scampi_res_file) {
#     open IN, $scampi_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $scampi_res=$_;
#             chomp $scampi_res;
#             $topo{'SCAMPI'} = $scampi_res;
#         }
#     }
#     close IN;
# }
# 
# if(-e $spoctopus_res_file) {
#     open IN, $spoctopus_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $spoctopus_res=$_;
#             chomp $spoctopus_res;
#             $topo{'SPOCTOPUS'} = $spoctopus_res;
#         }
#     }
#     close IN;
# }
# 
# if(-e $homology_res_file) {
#     open IN, $homology_res_file;
#     while(<IN>) {
#         if($_=~/[iMoS]/i) {
#             my $res = $_;
#             chomp $res;
#             $topo{'VHomology'} = $res;
#         }
#     }
#     close IN;
# }

my $topcons_res='';
my $length=0;
if(-e $topcons_res_file)
{
    open TOPCONS, $topcons_res_file;
    while(<TOPCONS>)
    {
        if($_=~/[iMoS]/i)
        {
            $topcons_res=$_;
            chomp $topcons_res;
            $length=length($topcons_res);
            $topo{'CONSENSUS_TOPCONS'} = $topcons_res;
       }
    }
}
close TOPCONS;

mk_topcons_plot($res_folder, $topcons_res, $length);
mk_plot($res_folder,\%topo, $length);

my @sort=sort {$a cmp $b} (keys %topo);

my $nice_top = print_seq_topo($seq,\@sort,\%topo);
my $nice_top_file = "$res_folder"."nicetop.html";
open(OUT, ">$nice_top_file");
print OUT $nice_top;
close(OUT);

sub get_boundaries_from_topology {#{{{
    my @BOUNDARIES;
    my ($topology)=@_;
    my $temp;

    while($topology=~m/i+/g)            #cytoplasm
    {
       $temp=$-[0]."\t".$+[0]."\ti";
       push @BOUNDARIES, $temp;
    }

    while($topology=~m/o+/g)            #extracellular
    {
       $temp=$-[0]."\t".$+[0]."\to";
       push @BOUNDARIES, $temp;
    }

    while($topology=~m/M+/g)            #TM helix
    {
       $temp=$-[0]."\t".$+[0]."\tM";
       push @BOUNDARIES, $temp;
    }
    while($topology=~m/M[M\.-]*M/g)     #TM helix with gaps inside
    {
       $temp=$-[0]."\t".$+[0]."\tM";
       push @BOUNDARIES, $temp;
    }

    while($topology=~m/S+/g)            #signal peptide
    {
       $temp=$-[0]."\t".$+[0]."\tS";
       push @BOUNDARIES, $temp;
    }
    while($topology=~m/\.+/g)           # un-aligned region, this is
                                        # specifically for homologous preidction
    {
       $temp=$-[0]."\t".$+[0]."\t.";
       push @BOUNDARIES, $temp;
    }
    while($topology=~m/-+/g)           # un-aligned region, this is
                                        # specifically for homologous preidction
    {
       $temp=$-[0]."\t".$+[0]."\t-";
       push @BOUNDARIES, $temp;
    }
    return(@BOUNDARIES);
}
#}}}
sub mk_plot {
    my ($outfolder, $topo, $protein_length) = @_;
    my $DG_res_file = $outfolder."dg.txt";
    my $DG_res_used = $outfolder."DG1.txt";
    my $whole_img_gnu_file = $outfolder."Topcons/total_image.gnu"; 
    my $whole_img_large_image = $outfolder."Topcons/total_image.large.png";
    my $whole_img_small_image = $outfolder."Topcons/total_image.png";

    open (DG, ">$DG_res_used");

    open (IN, $DG_res_file);
    while(<IN>)
    {
        if ($_=~/^\d+/)
        {
            print DG $_;
        }
    }
    close IN;
    close DG;
#    unlink ($DG_res_file);
#    unlink ($DG_res_used);

#     print "$sort -gk 2 $DG_res_used|$head -n 1|$awk '{print \$2}'"."\n";
#     print "$sort -grk 2 $DG_res_used|$head -n 1|$awk '{print \$2}'"."\n";
    my $DG_min = -2;
    my $DG_max = 6;

    $DG_min=`$sort -gk 2 $DG_res_used | $head -n 1 | $awk '{print \$2}'`;
    $DG_max=`$sort -grk 2 $DG_res_used | $head -n 1 | $awk '{print \$2}'`;

    my $DG_tick_max=int($DG_max+0.5 * ($DG_max <=> 0))+1;
    my $DG_tick_min=int($DG_min+0.5 * ($DG_min <=> 0))-1;
    #my $DG_axis_max=1.8*($DG_tick_max-$DG_tick_min)+$DG_tick_min;
    # the factor 2.0 controls the height of the DG plot, the larger the factor,
    # the smaller the height the DG plot 
    my $DG_axis_max=2.2*($DG_tick_max-$DG_tick_min)+$DG_tick_min;

    chomp($DG_min,$DG_max);
    my $obj_nr;
    open (OUT,">$whole_img_gnu_file") || die("Could not open file");
    print OUT "set style line 11 lt 1 lw 1 lc rgb \"blue\"\n";
    print OUT "set encoding iso_8859_1\n";
    print OUT "set yrange [0:50]\n";
    print OUT "set xrange [1:$length]\n";
    print OUT "set y2range [-2:$DG_axis_max]\n";

    print OUT "set autoscale xfix\n";
    print OUT "set ter png enh interlace size 2400,1680 font 'Nimbus,40'\n";
    print OUT "set y2label \'{/Symbol D}G (kcal/mol)                                             \' tc lt 3\n";
    #print OUT "set ytics scale 1,0.5 nomirror \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25, \"SPOCTOPUS\" 30.5 0, \"SCAMPI\" 33.5 0, \"PolyPhobius\" 36.5 0, \"Philius\" 39.5 0, \"OCTOPUS\" 42.5 0, \"TOPCONS\" 45.5 0\)\n";

    print OUT "set ytics scale 1,0.0 nomirror \(\"$showtext_homo\" 26.9 0, \"SPOCTOPUS\" 32.9 0, \"SCAMPI\" 35.9 0, \"PolyPhobius\" 38.9 0, \"Philius\" 41.9 0, \"OCTOPUS\" 44.9 0, \"TOPCONS\" 47.9 0\)\n";

    if ($DG_tick_max-$DG_tick_min < 10) {
        print OUT "set y2tics nomirror $DG_tick_min,1,$DG_tick_max\n";
    } else {
        print OUT "set y2tics nomirror $DG_tick_min,2,$DG_tick_max\n";
    }
    print OUT "set out \'$whole_img_large_image\'\n";

    print OUT "set lmargin 13.5\n";
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
    print OUT "set label 'Signal peptide' font 'Nimbus,30' at screen 0.785,0.982\n";
    print OUT "set ytic scale 0\n";

    #my $offset=0;
    my $offset=-3.6;
    for my $method (sort {$b cmp $a} keys %{$topo}) 
    {
        if ($topo->{$method} eq ""){
            my $label = "";
            if ($method eq "VHomology"){
                $label = "***No homologous TM proteins detected***";
            }else{
                $label = "***No signal peptide nor TM-regions predicted***";
            }
            print OUT "set label \"$label\" at 10, " . (30+$offset)  . " font \"Nimbus,35\"\n";
        }else{
            my @BOUNDARIES=get_boundaries_from_topology($topo->{$method});

            for(my $i=0;$i<scalar(@BOUNDARIES);$i++)
            {
                my @fields=split(/\t/,$BOUNDARIES[$i]);
                my $start=$fields[0]+1-0.5;
                my $end=$fields[1]+0.5;

                if($fields[2] eq "S")
                {
                    print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+1) . " fc rgb \"black\" fs noborder\n";
                }

#             if($fields[2] eq ".") # for unaligned region, just not shown
#             {
#                 print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+0.2) . " fc rgb \"red\" fs noborder\n";
#             }

                if($fields[2] eq "i")
                {
                    print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+0.2) . " fc rgb \"red\" fs noborder\n";
                }

                if($fields[2] eq "o")
                {
                    print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset+1-0.2) . " to $end," . (30+$offset+1) . " fc rgb \"blue\" fs noborder\n";
                }

                if($fields[2] eq "M")
                {
                    if ((substr($topo->{$method},$start-2,1) eq "i") || (substr($topo->{$method},$end,1) eq "o")) 
                    {
                        print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+1) . " fc rgb \"grey\" fs noborder\n";
                    } 

                    else 
                    {
                        print OUT "set object " . (++$obj_nr) . " rect from $start," . (30+$offset) . " to $end," . (30+$offset+1) . " fc rgb \"white\"\n";
                    }
                }
            }
        }
        if ($method eq "VHomology") { 
            $offset += 6;
        }else{
            $offset += 3;
        }
    }

    print OUT "plot '$DG_res_used' axes x1y2 w l t '' lt 3 lw 4\n";
    print OUT "exit\n";
    close(OUT);
    `$gnuplot $whole_img_gnu_file`;
    `$convert -scale 580 $whole_img_large_image $whole_img_small_image`;
}

sub mk_topcons_plot 
{
    my ($outfolder, $topcons_topo, $protein_length) = @_;

    my $reliability_res_file_start = $outfolder."Topcons/reliability.txt";
    my $reliability_res_file_used = $outfolder."Topcons/reliability.final";

    open (OUT_RELIAB, ">$reliability_res_file_used");
    open (RELIAB_BEG, $reliability_res_file_start);
    while(<RELIAB_BEG>)
    {
        if($_=~/(.*?)\t(.*)/)
        {
            my $pos_reliab=$1;
            my $reliab_score=$2;

            my $reliab_score_used = ($reliab_score/100);
            print OUT_RELIAB $pos_reliab."\t".$reliab_score_used."\n";
        }
    }
    close RELIAB_BEG;
    close OUT_RELIAB;

    my $topcons_gnu_file = $outfolder."Topcons/topcons.gnu"; 
    my $topcons_large_image = $outfolder."Topcons/topcons.large.png";  
    my $topcons_small_image = $outfolder."Topcons/topcons.png";  

#     print "$sort -gk 2 $reliability_res_file_used|$head -n 1|$awk '{print \$2}'"."\n";
#     print "$sort -grk 2 $reliability_res_file_used|$head -n 1|$awk '{print \$2}'"."\n";
    my $rel_min = 0;
    my $rel_max = 100;
    $rel_min=`$sort -gk 2 $reliability_res_file_used | $head -n 1 | $awk '{print \$2}'`;
    $rel_max=`$sort -grk 2 $reliability_res_file_used | $head -n 1 | $awk '{print \$2}'`;
    my $rel_tick_max=int(10*$rel_max+1)/10;
    ($rel_tick_max > 1) && ($rel_tick_max = 1);
    my $rel_tick_min=int(10*$rel_min)/10;
    my $rel_axis_max=1.5*($rel_tick_max-$rel_tick_min)+$rel_tick_min;
    my $factor=$rel_axis_max-$rel_tick_max;
    my $offset=0.4*$factor;

    chomp($rel_min,$rel_max);
    my $obj_nr;
    open(OUT,">$topcons_gnu_file") || die("Could not open file");
    print OUT "set encoding iso_8859_1\n";
    print OUT "set xrange [1:$length]\n";
    print OUT "set yrange [0.83:$rel_axis_max]\n";
    print OUT "set autoscale xfix\n";
    print OUT "set ter png enh interlace size 2400,840 font 'Nimbus,40'\n";
    print OUT "set xlabel \'Position\'\n";
    print OUT "set ylabel \'Reliability           \' \n";
    print OUT "set ytics nomirror $rel_tick_min,0.1,$rel_tick_max\n";

    print OUT "set out \'$topcons_large_image\'\n";
    print OUT "set tmargin 1.3\n";
    print OUT "set lmargin 11.5\n";
    print OUT "set rmargin 6.5\n";
    print OUT "set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775\n";

    my @BOUNDARIES = get_boundaries_from_topology($topcons_topo);

    for(my $i=0;$i<scalar(@BOUNDARIES);$i++)
    {
        my @fields=split(/\t/,$BOUNDARIES[$i]);
        my $start=$fields[0]+1-0.5;
        my $end=$fields[1]+0.5;

        my $small=0.0325;
        my $big=0.175;
        if($fields[2] eq "i")
        {
            print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$small*$factor) . " fc rgb \"red\" fs noborder\n";
        }

        if($fields[2] eq "S")
        {
            print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"black\" fs noborder\n";
        }

        if($fields[2] eq "o")
        {
            print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset+$big*$factor-$small*$factor) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"blue\" fs noborder\n";
        }

        if($fields[2] eq "M")
        {
            if ((substr($topcons_topo,$start-2,1) eq "i") || (substr($topcons_topo,$end,1) eq "o")) 
            {
                print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"grey\" fs noborder\n";# fs pattern 5\n";
            } 
            else 
            {
                print OUT "set object " . (++$obj_nr) . " rect from $start," . ($rel_tick_max+$offset) . " to $end," . ($rel_tick_max+$offset+$big*$factor) . " fc rgb \"white\"\n";# fs pattern 4\n";
            }
        }
    }
    print OUT "plot '$reliability_res_file_used' w l t '' lc rgb \"black\" lw 4\n";
    print OUT "exit\n";
    close(OUT);
    `$gnuplot $topcons_gnu_file`;
    `$convert -scale 580 $topcons_large_image $topcons_small_image`;
}

sub print_seq_topo {#{{{
    my ($seq,$sort,$topo) = @_;
    my $rowL=50; # Parameter
    my $colL=10; # Parameter
    my $cols_per_block=$rowL/$colL; ($rowL % $colL) && die; # Must be an integer...
    my $length=length($seq);
    my $nr_complete_blocks=int($length/$rowL);
    my $final_row_L=$length-$nr_complete_blocks*$rowL;
    my $ret = "<table cellspacing=0 cellpadding=0 border=0>\n";
    my $i;
    for my $method (@{$sort}) {
        my @BOUNDARIES = get_boundaries_from_topology($topo->{$method});
        if (($topo->{$method} =~ /NOTM/i) || ($topo->{$method} =~ /ERR/i) || scalar(@BOUNDARIES) == 0 || scalar(@BOUNDARIES) == 1) {
            #delete($topo->{$method});
            $topo->{$method} = "";
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
        for my $method (@{$sort}) {
            if ($method eq "CONSENSUS_TOPCONS") {
                $ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\"><span class=nowrap>TOPCONS&nbsp;&nbsp;&nbsp;</span></td>";
            } elsif ($method eq "VHomology") {
                if ($url_homo eq ""){
                    $ret .= "<tr><td style=\"font-family: Courier; Color:#0000FF\"><span class=nowrap>$showtext_homo&nbsp;&nbsp;&nbsp;</span></td>";
                }else{
                    $ret .= "<tr><td style=\"font-family: Courier; Color:#0000FF\"><span class=nowrap><a href=\"$url_homo\">$showtext_homo</a>&nbsp;&nbsp;&nbsp;</span></td>";
                }
            } else {
                $ret .= "<tr><td style=\"font-family: Courier;\"><span class=nowrap>$method&nbsp;&nbsp;&nbsp;</span></td>";
            }
            my $row = substr($topo->{$method},0,$rowL,"");
            #print "DEBUG: topo->method=<".$topo->{$method}.">\n" ;
            $row =~ s/(.{$colL})/$1\ /g;
            if ($method eq "CONSENSUS_TOPCONS") {
                $ret .= "<td style=\"font-family: Courier; Color:#990000;\">$row</td></tr>\n";
            } elsif ($method eq "VHomology") {
                $ret .= "<td style=\"font-family: Courier; Color:#0000FF;\">$row</td></tr>\n";
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
        for my $method (@{$sort}) {
            if ($method eq "CONSENSUS_TOPCONS") {
                $ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\">TOPCONS&nbsp;&nbsp;&nbsp;</td>";
            } elsif ($method eq "VHomology") {
                if ($url_homo eq ""){
                    $ret .= "<tr><td style=\"font-family: Courier; Color:#0000FF\">$showtext_homo&nbsp;&nbsp;&nbsp;</td>";
                }else{
                    $ret .= "<tr><td style=\"font-family: Courier; Color:#0000FF\"><span class=nowrap><a href=\"$url_homo\">$showtext_homo</a>&nbsp;&nbsp;&nbsp;</span></td>";
                }
            } else {
                $ret .= "<tr><td style=\"font-family: Courier;\">$method&nbsp;&nbsp;&nbsp;</td>";
            }
            my $row = substr($topo->{$method},0,$rowL,"");
            $row =~ s/(.{$colL})/$1\ /g;
            if ($method eq "CONSENSUS_TOPCONS") {
                $ret .= "<td style=\"font-family: Courier; Color:#990000;\">$row</td></tr>\n";
            } elsif ($method eq "VHomology") {
                $ret .= "<td style=\"font-family: Courier; Color:#0000FF;\">$row</td></tr>\n";
            } else {
                $ret .= "<td style=\"font-family: Courier;\">$row</td></tr>\n";
            }
        }
        $ret .= "</table>\n";
    }
    return $ret;
}#}}}
sub spaces {
    my ($str,$length)=@_;
    my $nr_spaces=$length-length($str);
    my $space_str="";
    for (my $i=0;$i<$nr_spaces;$i++) {
        $space_str .= "&nbsp;";
    }
    return $space_str;
}
