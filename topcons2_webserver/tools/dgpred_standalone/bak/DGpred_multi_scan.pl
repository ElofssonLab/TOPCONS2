#!/usr/bin/perl -w

use strict;

#my $ncbi_dir="/afs/pdc.kth.se/home/a/abernsel/.vol02/align/blast/";
my $ncbi_dir="/afs/pdc.kth.se/misc/pdc/volumes/pdc/prj.pdc.nanjiang.1/usr/share/blast/blast-2.2.17/";
#my $blastdb="/afs/pdc.kth.se/home/a/abernsel/.vol02/align/nr/uniref90.fasta";
#my $blastdb="/afs/pdc.kth.se/home/a/abernsel/.vol10/psi_pmembr/sprot.fas";
my $blastdb="$ENV{'BLASTDB'}/uniprotWG100Filtered";
my $scampi_msa_bin="/scampi-msa";

my ($Lmin,$Lmax,$seq_or_file,$seq,$infile);
(($seq_or_file=$ARGV[0]) &&
 ($Lmin=$ARGV[1]) &&
 ($Lmin =~ /^\d+$/))
    || die("Syntax: DGpred_multi_scan prffile Lmin [Lmax]\n");

unless (($Lmax=$ARGV[2]) && ($Lmax =~ /^\d+$/)) {
    $Lmax=$Lmin;
}

my @aminoacids = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
my %query=();

if ($seq_or_file =~ /^[A-Z]+$/) {     # Seems like a sequence
    $seq=$seq_or_file;
    ($seq =~ /[JO]/) && die("Sequence error: $seq\n");
    read_modhmm_prf($seq,\%query);
} else {                              # Seems like a file name
    $infile=$seq_or_file;
    open(my $IN,$infile) || die("Could not open file $infile\n");
    read_modhmm_prf($IN,\%query);
    close($IN);
}

my $seq_length=scalar(keys(%query));
$seq_length || die("Error in file format: $infile\n");
($Lmax > $seq_length) && ($Lmax = $seq_length);
($Lmin > $Lmax) && ($Lmin = $Lmax);

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


my $piover180=atan2(1,1)/45;
my ($lowest_ll,$lowest_ul,$lowest_segment);
my $first=1;
for (my $L=$Lmin;$L<=$Lmax;$L++) {
    my $flank = int(($L-1)/2+0.5);
    my $half_step;
    ((int($L/2) == ($L/2)) && ($half_step=0.5)) || ($half_step=0); # If L is even, half_steps are introduced
    for (my $k=$flank+1-$half_step;$k<=$seq_length-$flank+$half_step;$k++) {
	my $DG=0; my $DG_sum=0; my $DG_sin_sum=0; my $DG_cos_sum=0;
	for (my $i=1;$i<=$L;$i++) {
	    my $DG=0;
	    for my $aa (@aminoacids) {
		my $aa_freq=$query{$k-1-$flank+$i+$half_step}{$aa};
		$DG += $aa_freq*pos_spec_contrib($aa,$i,$L);
	    }
	    $DG_sum += $DG;
	    $DG_sin_sum += ($DG * sin(100*($i-1)*$piover180));
	    $DG_cos_sum += ($DG * cos(100*($i-1)*$piover180));
	}
	my $helix_DG=$DG_sum+$c0*sqrt($DG_sin_sum**2+$DG_cos_sum**2)+$c1+$c2*$L+$c3*$L*$L;
	print "$k $helix_DG\n";
	if (($first) || ($helix_DG < $lowest_segment)) {
	    $lowest_segment=$helix_DG;
	    $lowest_ll=$k-($L-1)/2;
	    $lowest_ul=$k+($L-1)/2;
	}
	$first && ($first=0);
    }
}

#print "$lowest_ll-$lowest_ul: " . sprintf("%.4f",$lowest_segment) . "\n";


# --------------------------------------------------------------------------------------------

sub pos_spec_contrib {
    my ($aa,$i,$L) = @_;
    my $pos = 9 * (2 * ($i-1)/($L-1) - 1);    
    if (($aa eq "W") || ($aa eq "Y")) {
	return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2) 
	    + $profiles{$aa}[2] * (   exp(-1*$profiles{$aa}[3]*($pos-$profiles{$aa}[4])**2) 
                                      + exp(-1*$profiles{$aa}[3]*($pos+$profiles{$aa}[4])**2)   );
    } else {
	return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2);
    }
}

sub read_modhmm_prf {
    my ($seq_or_fh,$profile)=@_;
    my $fh;
    if ($seq_or_fh =~ /^[A-Z]+$/) {     # Seems like a sequence
	$seq=$seq_or_fh;
	open(OUT,">/tmp/DGpred_tmpseq.fasta.2");
	print OUT ">temp\n$seq\n";
	close(OUT);
	open($fh,"/tmp/DGpred_tmpseq.fasta.2");
    } else {                            # Seems like a file handle
	$fh=$seq_or_fh;
    }
    my $i=0;
    while(<$fh>) {
	if (/^>(.+)$/) {    # Looks like a FASTA file!
	    my $name=$1;
	    my $seq = "";
	    while(<$fh>) {
		/^[A-Z]+$/ || last;
		$seq .= $_;
	    }
	    $seq =~ /[JO]/ && die("Error in FASTA format: $seq\n");
	    chomp($seq);
	    open(OUT,">/tmp/DGpred_tmpseq.fasta") || die("Could not open file: /tmp/DGpred_tmpseq.fasta");
	    print OUT ">$name\n$seq\n";
	    close(OUT);
	    `$ncbi_dir/bin/blastpgp -j 1 -h 1e-5 -e 1e-5 -b 100 -v 0 -d $blastdb -m 6 -i /tmp/DGpred_tmpseq.fasta > /tmp/DGpred_tmpseq.blast`;
	    `$scampi_msa_bin/msa62mod_oneround.pl /tmp/DGpred_tmpseq.blast /tmp/DGpred_tmpseq.fasta > /tmp/DGpred_tmpseq.prf`;
	    close($fh);
	    open($fh,"/tmp/DGpred_tmpseq.prf");
	    return read_modhmm_prf($fh,$profile);     # Recursive
	} else {
	    /^ALPHABET/ && last;
	}
    }
    while(<$fh>) {
	/^END/ && last;
	/^COL\s+(\d+):\s+(((\d|\.)+\s+){20})/ || die("Error in profile format: $_\n");
	my $col=$1; my @aa_count=split(/\s+/,$2);
	($col == ++$i) || die("Error in profile format: $_\n");
	my $total_count=0;
	$total_count += $_ for (@aa_count);
	$total_count || ($total_count=1); # If non-aligned X...
	my $j=0;
	for (@aminoacids) {
	    $profile->{$col}->{$_} = $aa_count[$j++]/$total_count;
	}
    }
}
