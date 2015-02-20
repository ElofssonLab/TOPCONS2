#!/usr/bin/perl -w


if($#ARGV != 0) {
  printf "Usage: make_gnuplot_input_from_nnprf_output.pl <nnprffile>\n";
  exit;
}

my $nnprffile = $ARGV[0];
#my $outfile = $ARGV[1];


open NNPRFFILE, "$nnprffile"
    or die "could not open $nnprffile";


my $gnufile = $nnprffile;
$gnufile =~ s/nnprf/gnu/;
my $pngfile = $nnprffile;
$pngfile =~ s/nnprf/png/; 
my $nr_cols = 0;
my @cols = ();
while(<NNPRFFILE>) {
    chomp;
    if($_ =~ 'Pos') {
	
    }
    else {
	$nr_cols++;
	$_ =~ s/\s+/ /g;
	my @row = split /\s/, $_;
	my $tot_MLGR = 0;
	my $tot_io = 0;
	for(my $i = 1; $i <= 4; $i++) {
	    $tot_MLGR += $row[$i];
	}
	if($tot_MLGR > 0) {
	    for(my $i = 1; $i <= 4; $i++) {
		$row[$i] /= $tot_MLGR;
	    } 
	}
	for(my $i = 5; $i <= 6; $i++) {
	    $tot_io += $row[$i];
	}
	if($tot_io > 0) {
	    for(my $i = 5; $i <= 6; $i++) {
		$row[$i] /= $tot_io;
	    } 
	}
	push @cols, @row;
    }
}

my $obj_nr;

open GNUFILE, ">"."$gnufile"
    or die "could not open $gnufi";
print GNUFILE "set terminal png enh interlace size 600,450 font \"Nimbus\" 8\n";
print GNUFILE "set style line 2\n";
print GNUFILE "set output \"$pngfile\"\n";
print GNUFILE "set multiplot\n";
## plot MLGR ##
print GNUFILE "set size 1,0.4\n";  
print GNUFILE "set origin 0.0,0.4\n";
print GNUFILE "set ylabel \'Network value\'\n";
print GNUFILE "set xlabel \'Position\' textcolor rgb \"white\"\n";
print GNUFILE "plot [1:${nr_cols}][0.0:1.0] 'test.nnprf' using 1:2 notitle smooth unique with lines lt 1,'test.nnprf' using 1:5 notitle smooth unique with lines lt 4, '/scratch/test.nnprf' using 1:3 notitle smooth unique with lines lt 5, '/scratch/test.nnprf' using 1:4 notitle smooth unique with lines lt 3\n";

## plot io ##
print GNUFILE "set origin 0.0,0.0\n";
print GNUFILE "unset ylabel\n";
print GNUFILE "set ylabel \'Network value\'\n";
print GNUFILE "unset xlabel\n";
print GNUFILE "set xlabel \'Position\' textcolor rgb \"black\"\n";
print GNUFILE "plot [1:${nr_cols}][0.0:1.0] '/scratch/test.nnprf' using 1:6 notitle smooth unique with lines lt 2,'/scratch/test.nnprf' using 1:7 notitle smooth unique with lines lt 6\n";
## plot topology ##
print GNUFILE "set size 1,0.2\n";
print GNUFILE "set origin 0.0,0.8\n";
print GNUFILE "unset ylabel\n";
print GNUFILE "set ylabel \'Topology\'\n";
print GNUFILE "set ytics 0.0,.1,0.1 textcolor rgb \"white\"\n";
print GNUFILE "unset xlabel\n";
print GNUFILE "set xlabel \'Position\' textcolor rgb \"white\"\n";
print GNUFILE "set label '\***No TM-regions predicted***\' at graph 0.3,0.5\n";
print GNUFILE "plot [1:$nr_cols][0.0:0.1] 2 w l t \'\' lc rgb \"white\"\n";
print GNUFILE "unset label\n";

print GNUFILE "set size 1,1\n";
print GNUFILE "set origin 0.0,0.0\n";
print GNUFILE "unset ytics\n";
print GNUFILE "unset ylabel\n";
print GNUFILE "unset xtics\n";
print GNUFILE "unset xlabel\n";
print GNUFILE "unset border\n";

# topology labels
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.10,0.987 to screen 0.12,0.993 fc rgb \"green\" fs noborder\n";
print GNUFILE "set label 'Inside' at screen 0.13,0.99\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.22,0.987 to screen 0.24,0.993 fc rgb \"brown\" fs noborder\n";
print GNUFILE "set label 'Outside' at screen 0.25,0.99\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.35,0.980 to screen 0.37,1.0 fc rgb \"red\" fs noborder\n";
print GNUFILE "set label 'TM-helix' at screen 0.38,0.99\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.49,0.980 to screen 0.51,1.0 fc rgb \"magenta\" fs noborder\n";
print GNUFILE "set label 'Reentrant region' at screen 0.52,0.99\n";
# MLGR labels
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.10,0.785 to screen 0.12,0.795 fc rgb \"red\" fs noborder\n";
print GNUFILE "set label 'Membrane (0-13)' at screen 0.13,0.79\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.30,0.785 to screen 0.32,0.795 fc rgb \"magenta\" fs noborder\n";
print GNUFILE "set label 'Reentrant (0-13)' at screen 0.33,0.79\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.51,0.785 to screen 0.53,0.795 fc rgb \"cyan\" fs noborder\n";
print GNUFILE "set label 'Close loop (13-23)' at screen 0.54,0.79\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.75,0.785 to screen 0.77,0.795 fc rgb \"blue\" fs noborder\n";
print GNUFILE "set label 'Globular loop (23-)' at screen 0.78,0.79\n";


#io labels
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.10,0.385 to screen 0.12,0.395 fc rgb \"green\" fs noborder\n";
print GNUFILE "set label 'Inside' at screen 0.13,0.39\n";
print GNUFILE "set object " . (++$obj_nr) . " rect from screen 0.22,0.385 to screen 0.24,0.395 fc rgb \"brown\" fs noborder\n";
print GNUFILE "set label 'Outside' at screen 0.25,0.39\n";

print GNUFILE "plot [0.0:1.0][0.0:1.0] 2 w l t \'\' lc rgb \"white\"\n";
#print GNUFILE "plot [1:${nr_cols}][0.0:1.0] '/scratch/test.nnprf' using 1:6 title 'in' smooth unique with lines lt 2,'/scratch/test.nnprf' using 1:7 title 'out' smooth unique with lines lt 6\n";
close GNUFILE;
