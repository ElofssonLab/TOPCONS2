#!/usr/bin/perl
use strict;
use warnings;

my $gnuplot = "/usr/bin/gnuplot";
my $convert = "/usr/bin/convert";
my $sort = "/usr/bin/sort";
my $awk = "/usr/bin/awk";
my $head="/usr/bin/head";

$ENV{'GDFONTPATH'}='./fonts/';
my $res_folder = $ARGV[0];  #provide this from the command line, contains all sub-folders for various prediction methods

my $seqfile = $res_folder."seq.fa";
my $octopus_res_file = $res_folder."OCTOPUS/query.top";
my $philius_res_file = $res_folder."philius/query.top";
my $polyphobius_res_file = $res_folder."PolyPhobius/query.top";
my $scampi_res_file = $res_folder."SCAMPI_MSA/query.top";
my $spoctopus_res_file = $res_folder."SPOCTOPUS/query.top";
my $topcons_res_file = $res_folder."Topcons/topcons.top";
my $DG_res_file = $res_folder."Topcons/topcons.DG";

my $seq = "";
if (! -e $seqfile){
    warn "$seqfile does not exist";
}else{
    $seq = ReadSeq($seqfile);
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

if(-e $octopus_res_file)
{
    open OCTOPUS, $octopus_res_file;
    while(<OCTOPUS>)
    {
        if($_=~/[iMoS]/i)
        {
            my $octopus_res=$_;
            chomp $octopus_res;
            $topo{'OCTOPUS'} = $octopus_res;
        }
    }
    close OCTOPUS;
}
else
{
    $topo{'OCTOPUS'} = ' ';
}

if(-e $philius_res_file)
{
    open PHILIUS, $philius_res_file;
    while(<PHILIUS>)
    {
        if($_=~/[iMoS]/i)
        {
            my $philius_res=$_;
            chomp $philius_res;
            $topo{'PHILIUS'} = $philius_res;
        }
    }
    close PHILIUS;
}
else
{
    $topo{'PHILIUS'} = ' ';
}


if(-e $polyphobius_res_file)
{
    open POLYPHOBIUS, $polyphobius_res_file;
    while(<POLYPHOBIUS>)
    {
        if($_=~/[iMoS]/i)
        {
            my $polyphobius_res=$_;
            chomp $polyphobius_res;
            $topo{'PolyPhobius'} = $polyphobius_res;
        }
    }
    close POLYPHOBIUS;
}
else
{
    $topo{'POLYPHOBIUS'} = ' ';
}

if(-e $scampi_res_file)
{
    open SCAMPI, $scampi_res_file;
    while(<SCAMPI>)
    {
        if($_=~/[iMoS]/i)
        {
            my $scampi_res=$_;
            chomp $scampi_res;
            $topo{'SCAMPI'} = $scampi_res;
        }
    }
    close SCAMPI;
}
else
{
    $topo{'SCAMPI'} = ' ';
}

if(-e $spoctopus_res_file)
{
    open SPOCTOPUS, $spoctopus_res_file;
    while(<SPOCTOPUS>)
    {
        if($_=~/[iMoS]/i)
        {
            my $spoctopus_res=$_;
            chomp $spoctopus_res;
            $topo{'SPOCTOPUS'} = $spoctopus_res;
        }
    }
    close SPOCTOPUS;
}
else
{
    $topo{'SPOCTOPUS'} = ' ';
}

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

sub get_boundaries_from_topology 
{
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

    while($topology=~m/M+/g)            #membrane
    {
       $temp=$-[0]."\t".$+[0]."\tM";
       push @BOUNDARIES, $temp;
    }

    while($topology=~m/S+/g)            #signal peptide
    {
       $temp=$-[0]."\t".$+[0]."\tS";
       push @BOUNDARIES, $temp;
    }
    return(@BOUNDARIES);
}

sub mk_plot 
{
    my ($outfolder, $topo, $protein_length) = @_;
    my $DG_res_file = $outfolder."dg.txt";
    my $DG_res_used = $outfolder."DG1.txt";
    my $whole_img_gnu_file = $outfolder."Topcons/total_image.gnu"; 
    my $whole_img_large_image = $outfolder."Topcons/total_image.large.png";
    my $whole_img_small_image = $outfolder."Topcons/total_image.png";

    open (DG, ">$DG_res_used");

    open (TMP, $DG_res_file);
    while(<TMP>)
    {
        if ($_=~/^\d+/)
        {
            print DG $_;
        }
    }
    close TMP;
    close DG;
    #unlink ($DG_res_file);

    my $DG_min=`$sort -gk 2 $DG_res_used|$head -n 1|$awk '{print \$2}'`;
    my $DG_max=`$sort -grk 2 $DG_res_used|$head -n 1|$awk '{print \$2}'`;
    
    my $DG_tick_max=int($DG_max+0.5 * ($DG_max <=> 0))+1;
    my $DG_tick_min=int($DG_min+0.5 * ($DG_min <=> 0))-1;
    my $DG_axis_max=1.8*($DG_tick_max-$DG_tick_min)+$DG_tick_min;

    chomp($DG_min,$DG_max);
    my $obj_nr;
    open (OUT,">$whole_img_gnu_file") || die("Could not open file");
    print OUT "set encoding iso_8859_1\n";
    print OUT "set yrange [0:50]\n";
    print OUT "set xrange [1:$length]\n";
    print OUT "set y2range [-2:$DG_axis_max]\n";

    print OUT "set autoscale xfix\n";
    print OUT "set ter png enh interlace size 2400,1680 font 'Nimbus,40'\n";
    print OUT "set y2label \'{/Symbol D}G (kcal/mol)                                             \' tc lt 3\n";
    #print OUT "set ytics scale 1,0.5 nomirror \(\"0\" 0, \"5\" 5, \"10\" 10, \"15\" 15, \"20\" 20, \"25\" 25, \"SPOCTOPUS\" 30.5 0, \"SCAMPI\" 33.5 0, \"PolyPhobius\" 36.5 0, \"Philius\" 39.5 0, \"OCTOPUS\" 42.5 0, \"TOPCONS\" 45.5 0\)\n";
    print OUT "set ytics scale 1,0.5 nomirror \(\"SPOCTOPUS\" 30.5 0, \"SCAMPI\" 33.5 0, \"PolyPhobius\" 36.5 0, \"Philius\" 39.5 0, \"OCTOPUS\" 42.5 0, \"TOPCONS\" 45.5 0\)\n";
     
    if ($DG_tick_max-$DG_tick_min < 10) 
    {
        print OUT "set y2tics nomirror $DG_tick_min,1,$DG_tick_max\n";
    } 

    else 
    {
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

    my $offset=0;
    for my $method (sort {$b cmp $a} keys %{$topo}) 
    {
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
        $offset+=3;
    }
    
    print OUT "plot '$DG_res_used' axes x1y2 w l t '' lt 3 lw 4\n";
    print OUT "exit\n";
    close(OUT);
    `$gnuplot $whole_img_gnu_file`;
    `$convert -scale 600x420 $whole_img_large_image $whole_img_small_image`;
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

    my $rel_min=`$sort -gk 2 $reliability_res_file_used|$head -n 1|$awk '{print \$2}'`;
    my $rel_max=`$sort -grk 2 $reliability_res_file_used|$head -n 1|$awk '{print \$2}'`;
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

    my @BOUNDARIES=get_boundaries_from_topology($topcons_topo);
    
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
    `$convert -scale 600x210 $topcons_large_image $topcons_small_image`;
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
    for my $method (@{$sort}) {
        if ($method eq "CONSENSUS_TOPCONS") {
        $ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\"><span class=nowrap>TOPCONS&nbsp;&nbsp;&nbsp;</span></td>";
        } else {
        $ret .= "<tr><td style=\"font-family: Courier;\"><span class=nowrap>$method&nbsp;&nbsp;&nbsp;</span></td>";
        }
        my $row = substr($topo->{$method},0,$rowL,"");
        $row =~ s/(.{$colL})/$1\ /g;
        if ($method eq "CONSENSUS_TOPCONS") {
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
    for my $method (@{$sort}) {
        if ($method eq "CONSENSUS_TOPCONS") {
        $ret .= "<tr><td style=\"font-family: Courier; Color:#990000;\">TOPCONS&nbsp;&nbsp;&nbsp;</td>";
        } else {
        $ret .= "<tr><td style=\"font-family: Courier;\">$method&nbsp;&nbsp;&nbsp;</td>";
        }
        my $row = substr($topo->{$method},0,$rowL,"");
        $row =~ s/(.{$colL})/$1\ /g;
        if ($method eq "CONSENSUS_TOPCONS") {
        $ret .= "<td style=\"font-family: Courier; Color:#990000;\">$row</td></tr>\n";
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
