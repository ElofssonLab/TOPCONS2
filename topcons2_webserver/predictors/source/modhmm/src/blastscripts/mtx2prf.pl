#!/usr/bin/perl -w

if($#ARGV < 1) {
    print "Usage: mtx2prf.pl <mtxfile> <outdir>\n";
    exit;
}

my $mtx_file=$ARGV[0];
my $outdir = $ARGV[1]."/";
($mtx_file =~ /^((.+\/)?(.+))\.mtx$/) or die;
$seq_id=$3;
my $prf_file = "$outdir"."$seq_id".".prf";
#die "$3\n";

%aa_freq=('A',7.588652e-02,
	  'C',1.661404e-02,
	  'D',5.288085e-02,
	  'E',6.371048e-02,
	  'F',4.091684e-02,
	  'G',6.844532e-02,
	  'H',2.241501e-02,
	  'I',5.813866e-02,
	  'K',5.958784e-02,
	  'L',9.426655e-02,
	  'M',2.372756e-02,
	  'N',4.455974e-02,
	  'P',4.908899e-02,
	  'Q',3.974479e-02,
	  'R',5.169278e-02,
	  'S',7.125817e-02,
	  'T',5.677277e-02,
	  'V',6.585309e-02,
	  'W',1.235821e-02,
	  'Y',3.188801e-02);

print "$mtx_file\n";
open(IN,$mtx_file) or die;
open(OUT,">$prf_file") or die;
$length=<IN>;
$seq=<IN>;
for ($i=0;$i<12;$i++) { <IN>; }
$lambda=0.3176;
print OUT "Sequence: $seq_id\nNR of aligned sequences: 100\nLength of query sequence: $length\n\nSTART 1\nALPHABET:";
for (sort keys %aa_freq) { print OUT "\t$_"; }
print OUT "\t-\t<SPACE>\t<LABEL>\t<QUERY>\n";
$col=1;
while(<IN>) {
    /^((\d|-)+\s+){26}$/ or die($_);
    @log_list=split(/\s+/,$_);
    $query_aa=substr($seq,$col-1,1);
    $i=0;
    for $aa (sort keys %aa_freq) {
	(($i == 0) || ($i == 2) || ($i == 21)) && $i++;  # Some columns are junk
	$log_vector{$aa}=$log_list[$i];
	$i++;
    }
    $sum=0;
    while (($aa,$log_nr) = each(%log_vector)) {
	$prob_vector{$aa} = $aa_freq{$aa}*exp($log_nr/100 * $lambda);
	$sum += $prob_vector{$aa};
    }
    # Normalize and round
    for $aa (keys %prob_vector) {
	$norm_perc=sprintf("%.4f",sprintf("%.3e",100*$prob_vector{$aa}/$sum));
	$decimals=4;
	(substr($norm_perc,0,1) eq "0") && $decimals++;
	$prob_vector{$aa}=substr($norm_perc,0,$decimals+1);
	$prob_vector{$aa} = $prob_vector{$aa} / 100;
    }
    $spaces = " ";
    for ($i=1; $i < 5 - length(sprintf("%s",$col)); $i++) { $spaces .= " "; }
    print OUT "COL$spaces$col:";
    for $aa (sort keys %prob_vector) {
      printf OUT "\t%.3f", $prob_vector{$aa};
      #print OUT "\t$prob_vector{$aa}";
    }
    print OUT "\t0.0000\t0.0000\t.\t$query_aa\n";
    $col++;
    
}
print OUT "END 1\n";
close(OUT);
