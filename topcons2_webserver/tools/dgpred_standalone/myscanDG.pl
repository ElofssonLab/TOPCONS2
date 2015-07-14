#!/usr/bin/perl -w

# 2010-09-16 
# modified from Bernsel's code "DGpred_multi_scan.pl"
# this code should work on both single sequence and profile 
# Nanjiang 
# ChangeLog 2012-02-23 
#   Add '//' to the end of each record 
use strict;

my $ncbi_dir="/afs/pdc.kth.se/misc/pdc/volumes/pdc/prj.pdc.nanjiang.1/usr/share/blast/blast-2.2.17/";
my $blastdb="$ENV{'BLASTDB'}/uniprotWG100Filtered";
my $scampi_msa_bin="/scampi-msa";
my $Lmin = 19;
my $Lmax = 19;
my $informat="";
my $infile = "";
my $listFile = "";
my $outfile="/dev/stdout";
my $isClean=1;
my $method = 0; # the fast_version with the method 1 does not work very well by perl

##other global variables
my @aminoacids = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
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


my $usage="
Usage:   myscanDG.pl [Options] <fasta amino acid file> or <prf file>
Note: the output is in HTML format
Options:
    -l|--list    <file> : set the list file, one line per file name
    -o|--outfile <file> : output the result to outpath, default=$outfile
    -informat    fa|prf : set the format of the input file, if not set, it will be detected automatically
    -multi              : using profile to scan the DG, in that case,
                        : blastgpg will be run if the input is fasta sequecne
    -lmin        <int>  : set the minimum length of the sliding window, default=$Lmin
    -lmax        <int>  : set the maximum length of the sliding inwdow, default=$Lmax
    -ncbidir     <str>  : set the path which contains the program blastpgp, default=$ncbi_dir
    -blastdb     <str>  : set the blast database, default=$blastdb
    -scampipath  <str>  : set the scampi-msa bin path, default=$scampi_msa_bin
    -nc|--not-clean               : do not clean the temporary dir

Modified from the original perl script 'DGpred_multi_scan.pl', 
by Nanjiang Shu, created 2010-09-16, updated 2010-09-17

==============================
Examples:
# 1. calculate DG values for 'test.fa' based on single sequence
    myscanDG.pl test.fa
# 2. calculate DG values for 'test.fa' based on multiple sequence alignment
    myscanDG.pl -multi test.fa
# 3. calculate DG values from the prf file 
    myscanDG.pl test.prf
";


# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
my $numArgs = $#ARGV+1;
if($numArgs < 1)
{
    &PrintHelp;
    exit;
}

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
        elsif($ARGV[$i] =~ /^-/)
        {
            if($ARGV[$i] eq "-o" || $ARGV[$i] eq "--outfile"  || $ARGV[$i] eq "-outfile" )
            {   
                $outfile = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-multi" || $ARGV[$i] eq "--multi" )
            {   
                $informat = "prf";
                $i +=1;
            }
            elsif($ARGV[$i] eq "-lmin" || $ARGV[$i] eq "--lmin" )
            {   
                $Lmin = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-lmax" || $ARGV[$i] eq "--lmax" )
            {   
                $Lmax = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-informat" || $ARGV[$i] eq "--informat" )
            {   
                $informat = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-l" || $ARGV[$i] eq "--list" )
            {   
                $listFile = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-ncbidir" || $ARGV[$i] eq "--ncbidir" )
            {   
                $ncbi_dir = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-blastdb" || $ARGV[$i] eq "--blastdb" )
            {   
                $blastdb = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-scampipath" || $ARGV[$i] eq "--scampipath" )
            {   
                $scampi_msa_bin = $ARGV[$i+1];
                $i +=2;
            }
            elsif($ARGV[$i] eq "-nc" || $ARGV[$i] eq "--not-clean" )
            {   
                $isClean = 0;
                $i +=1;
            }
            else
            {
                die "wrong argument $ARGV[$i]";
            }
        }
        else 
        {
            $infile = $ARGV[$i];
            $i ++;
        }
    }
}#}}}

if (($infile eq "" ) && ($listFile eq ""))
{
    die("Error! Input file not set. type myscanDG.pl -h to see syntax.");
}
my $fpout;
if ($outfile ne "/dev/stdout")
{
    open($fpout,">", $outfile) || die "Can not open the output file $outfile for write\n";
}
else #output to STDOUT
{
    open($fpout,'>-' );
}

my %EXPValue=();
for (my $L = $Lmin; $L <= $Lmax; $L++)
{
    for (my $i=1;$i<=$L;$i++)
    {
        my $pos = 9 * (2 * ($i-1)/($L-1) - 1);    
        for my $aa (@aminoacids) 
        {
            if (($aa eq "W") || ($aa eq "Y")) 
            {
            } 
            else
            {
                $EXPValue{$L}{$i}{$aa} = exp(-1*$profiles{$aa}[1]*$pos**2); 
            }
        }
    }
}

if ($infile ne "")
{
    MyScanDG($infile, $fpout);
}

if ($listFile ne "")
{
    my @filelist;
    open(IN, "<", $listFile) || die "Can not open the input list file $listFile\n";
    while (<IN>)
    {
        chomp;
        push(@filelist,$_);
    }
    close(IN);

    my $file;
    for $file (@filelist)
    {
        MyScanDG($file, $fpout);
    }
}

close($fpout);

# --------------------------------------------------------------------------------------------
# subroutines
# --------------------------------------------------------------------------------------------

sub GetSeqIDFromAnnotation#{{{
{
    my ($annoLine) = @_;
    my $seqID ="";
    my @strs = split(' ', $annoLine);
    $annoLine = $strs[0];
    if ($annoLine && $annoLine =~ '|')
    {
        @strs = split('\|', $annoLine);
        if ($strs[0] eq "sp" || $strs[0] eq "lcl" || $strs[0] eq "tr")
        {
            $seqID= $strs[1];
        }
        else
        {
            $seqID = $strs[0];
        }
    }
    else
    {
        $seqID = $annoLine;
    }
    return $seqID;
}#}}}
sub ReadFasta#{{{
{
    my ($infile) = @_;
    my $line;
    my @seqList=();
    my @seqIDList=();
    my @annoLineList=();
    open(IN, "<", $infile) || die "Cannot open FASTA file $infile.\n";
    my $seqnum = 0;
    my $first = 0;
    while (<IN>)
    {
        $line = $_;
        chomp($line);
        if ($line =~ /^>/)
        {
            $annoLineList[$seqnum] = substr($line,1);
            $seqList[$seqnum] = "";
            $seqIDList[$seqnum] = GetSeqIDFromAnnotation($annoLineList[$seqnum]);
            if ($seqIDList[$seqnum] eq "" )
            {
                $seqIDList[$seqnum]=$infile . '_' . $seqnum;
            }
            $seqnum +=  1;
            if ($first == 0)
            {
                $first = 1;
            }
            next;
        }
        if ($first == 0)
        {
            die "Not a standard FASTA file. Stop.\n";
        }
        $seqList[$seqnum - 1] .= $line;
    }
    close(IN);

    #print "$#seqIDList+1\n";
    #print "$seqIDList[0]\n";
    #print "$annoLineList[0]\n";
    #print "$seqList[0]\n";

    return (\@seqIDList, \@annoLineList, \@seqList);
} #}}}
sub ReadSingleFasta#{{{
{
    my ($infile) = @_;
    open(IN, "<", $infile) || die "Can not open the file $infile\n";
    my $seq = "";
    my $annoLine = "";
    my $seqID="";
    while(<IN>) 
    {
        if (/^>(.+)$/) 
        {   
            chomp;
            $annoLine = substr($_,1);
            $seqID = GetSeqIDFromAnnotation($annoLine);
            if ($seqID eq "" )
            {
                $seqID=$infile;
            }
            while(<IN>) 
            {
                chomp;
                /^[A-Z]+$/ || last;
                $seq .= $_;
            }
        }
    }
    close(IN);
    return ($seqID, $annoLine, $seq);
}#}}}
sub GetFileFormat#{{{
{
    my ($filename) = @_;
    my $fileformat = "";
    my $ext = ($filename =~ m/([^.]+)$/)[0];
    if ($ext eq "fa" || $ext eq "fasta" || $ext eq "aa")
    {
        $fileformat = "fa";
    }
    elsif ($ext eq "prf")
    {
        $fileformat = "prf";
    }
    if ($fileformat eq "")
    {
        open(IN, "<$filename") || die "Can not open the file $filename\n";
        while(<IN>) 
        {
            if (/^>/ )
            {
                $fileformat = "fa" ;
                last;
            }
            if (/^Sequence/ )
            {
                $fileformat = "prf"; 
                last; 
            }
        }
        close(IN);
    }
    # if still the format is not recognized 
    if ($fileformat eq "")
    {
        die "The file $filename has wrong file format\n";
    }
    return $fileformat;
}#}}}
sub CalcDG#{{{
{
    my ($query_ref, $seq, $seq_length, $Lmin, $Lmax, $informat, $fpout) = @_;
    my %query = %$query_ref;
    my ($lowest_ll,$lowest_ul,$lowest_segment);
    my $first=1;
    for (my $L=$Lmin;$L<=$Lmax;$L++) 
    {
        my $flank = int(($L-1)/2+0.5);
        my $half_step;
        ((int($L/2) == ($L/2)) && ($half_step=0.5)) || ($half_step=0); # If L is even, half_steps are introduced
        print $fpout "#L: $L\n";
        printf $fpout ("#Number of sliding windows: %d\n",$seq_length-$flank+$half_step-($flank+1-$half_step)+1);
        for (my $k=$flank+1-$half_step;$k<=$seq_length-$flank+$half_step;$k++) 
        {
            my $DG=0; my $DG_sum=0; my $DG_sin_sum=0; my $DG_cos_sum=0;
            for (my $i=1;$i<=$L;$i++) 
            {
                my $DG=0;
                my $aa;
                my $idx;
                if ($informat eq "fa")
                {                          
                    $idx = $k-1-$flank+$i+$half_step-1;
                    $aa = substr($seq, $idx, 1);
                    #print STDERR "$idx  $aa\n";
                    if ($method == 0)
                    {
                        $DG = pos_spec_contrib($aa,$i,$L);
                    }
                    else
                    {
                        $DG = pos_spec_contrib_fast_version($aa,$i,$L);
                    }
                }
                else #prf format
                {
                    for my $aa (@aminoacids) 
                    {
                        my $aa_freq=$query{$k-1-$flank+$i+$half_step}{$aa};
                        if ($method == 0)
                        {
                            $DG += $aa_freq*pos_spec_contrib($aa,$i,$L);
                        }
                        else
                        {
                            $DG += $aa_freq*pos_spec_contrib_fast_version($aa,$i,$L);
                        }
                    }
                }

                $DG_sum += $DG;
                $DG_sin_sum += ($DG * sin(100*($i-1)*$piover180));
                $DG_cos_sum += ($DG * cos(100*($i-1)*$piover180));
            }
            my $helix_DG=$DG_sum+$c0*sqrt($DG_sin_sum**2+$DG_cos_sum**2)+$c1+$c2*$L+$c3*$L*$L;

            ##print out the data
#            print $fpout "$k $helix_DG\n";
            printf $fpout ("%g %.3f\n",$k, $helix_DG);

            if (($first) || ($helix_DG < $lowest_segment)) 
            {
                $lowest_segment=$helix_DG;
                $lowest_ll=$k-($L-1)/2;
                $lowest_ul=$k+($L-1)/2;
            }
            $first && ($first=0);
        }
    }
}#}}}
sub MyScanDG #{{{
{
    my ($infile, $fpout)  = @_;

    my $tmp_informat = $informat;
    if ($informat eq "" )
    {
        $tmp_informat = GetFileFormat($infile);
    }

    my %query=();
    my $seq = "";
    #my $annoLine = "";
    #my $seqID = "";
    my $seq_length;
    if ($tmp_informat eq "prf" )
    {
        open(my $IN,$infile) || die("Could not open file $infile\n");
        read_modhmm_prf($IN,\%query);
        close($IN);
        $seq_length=scalar(keys(%query));
        $seq_length || die("Error in file format: $infile\n");
        print $fpout "#DG calculated from profile\n";
        print $fpout "#filename: $infile\n";
        print $fpout "#SeqLength: $seq_length\n";
        CalcDG(\%query, $seq, $seq_length, $Lmin, $Lmax, $tmp_informat, $fpout)
    }
    elsif ($tmp_informat eq "fa") #fasta file can include multiple sequences
    {
        #($seqID, $annoLine, $seq)= ReadSingleFasta($infile);
        my ($seqIDList_ref, $annoLineList_ref, $seqList_ref)= ReadFasta($infile);
        my @seqIDList = @$seqIDList_ref;
        my @annoLineList = @$annoLineList_ref;
        my @seqList = @$seqList_ref;

        print $fpout "#DG calculated from single sequence\n";
        my $numseq=@seqIDList;
        for (my $i = 0; $i < $numseq; $i ++)
        {
            $seq_length=length($seqList[$i]);
            print $fpout "#SeqID: $seqIDList[$i]\n";
            print $fpout "#>$annoLineList[$i]\n";
            print $fpout "#Seq: $seqList[$i]\n";
            print $fpout "#SeqLength: $seq_length\n";
            ($Lmax > $seq_length) && ($Lmax = $seq_length);
            ($Lmin > $Lmax) && ($Lmin = $Lmax);
            printf $fpout ("#Number of window sizes: %d\n",$Lmax-$Lmin+1);
            CalcDG(\%query, $seqList[$i], $seq_length, $Lmin, $Lmax, $tmp_informat, $fpout);
            print $fpout "//\n";  
        }
    }
}#}}}
sub pos_spec_contrib  #{{{
{
    my ($aa,$i,$L) = @_;
    my $pos = 9 * (2 * ($i-1)/($L-1) - 1);    
    if (($aa eq "W") || ($aa eq "Y")) {
	return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2) 
	    + $profiles{$aa}[2] * (   exp(-1*$profiles{$aa}[3]*($pos-$profiles{$aa}[4])**2) 
                                      + exp(-1*$profiles{$aa}[3]*($pos+$profiles{$aa}[4])**2)   );
    } else {
	return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2);
    }
}#}}}
sub pos_spec_contrib_fast_version  #{{{
# does not work very well
{
    my ($aa,$i,$L) = @_;
    my $pos = 9 * (2 * ($i-1)/($L-1) - 1);    
    #printf STDERR ("%g\n",$pos); #debug
    if (($aa eq "W") || ($aa eq "Y")) 
    {
        return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2) 
        + $profiles{$aa}[2] * (   exp(-1*$profiles{$aa}[3]*($pos-$profiles{$aa}[4])**2) 
            + exp(-1*$profiles{$aa}[3]*($pos+$profiles{$aa}[4])**2)   );
    } 
    else 
    {
        return $profiles{$aa}[0] * $EXPValue{$L}{$i}{$aa}; 
        #return $profiles{$aa}[0] * exp(-1*$profiles{$aa}[1]*$pos**2);
    }
}#}}}
sub read_modhmm_prf  #{{{
{
    my ($seq_or_fh,$profile)=@_;
    my $fh;
    if ($seq_or_fh =~ /^[A-Z]+$/) {     # Seems like a sequence
	my $seq=$seq_or_fh;
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
}#}}}
sub PrintHelp#{{{
{
    print $usage;
}#}}}
#END --------------------------------------------------------------------------------------------
