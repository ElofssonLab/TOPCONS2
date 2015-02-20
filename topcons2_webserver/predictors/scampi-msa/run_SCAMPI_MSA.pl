#!/usr/bin/perl

use strict;
use warnings;
use File::Temp "tempdir";
use File::Basename;

#ämy $fastafile=$ARGV[0];
my $seqfolder = $ARGV[0];#dirname($fastafile);
my $outDir = $ARGV[1]; 
my $tmpdir = $seqfolder; #$ARGV[2];#tempdir("/tmp/scampi-msa_tmpdir_XXXXXX")
my $outTmpFile = $tmpdir."SCAMPI_MSA.top";

#my $bindir = "./";

$ENV{'PATH'}="/bin:/usr/bin";
system("/bin/echo DGHMM_KR_21_multi.hmg > DGHMM_KR_21_multi.txt");
system("/bin/echo $seqfolder/query.raw.prf > $seqfolder/query.raw.prf.snf");
system("./modhmms_scampi -f prf -s $seqfolder/query.raw.prf.snf -m DGHMM_KR_21_multi.txt -r replacement_letter_multi.rpl --nopostout --viterbi -u -L -g > $tmpdir/scampi_modhmmres.xml");
system("./modhmmxml2top < $tmpdir/scampi_modhmmres.xml > $outTmpFile");

my $final_pred_scampi_msa;
mkdir("$outDir/SCAMPI_MSA");
open FINAL, ">$outDir".'SCAMPI_MSA/query.top';
open OUT, $outTmpFile;
while(<OUT>)
{
	if($_=~/^>/)
	{
		$final_pred_scampi_msa = <OUT>;
		$final_pred_scampi_msa=~s/[rTJI]/i/ig;
		$final_pred_scampi_msa=~s/[nhsRCGO]/o/ig;
		$final_pred_scampi_msa=~s/[le]/M/ig;
	}

    print FINAL $final_pred_scampi_msa;
}
close OUT;
close FINAL;


