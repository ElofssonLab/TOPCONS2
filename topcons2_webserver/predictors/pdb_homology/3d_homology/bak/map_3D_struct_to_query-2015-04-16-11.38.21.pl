#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( max );
use Cwd 'abs_path';
use File::Basename;
my $rundir = dirname(abs_path(__FILE__));

my $usage = "USAGE: $0 seqfile DB_fasta DB_3line outpath";
my $numArgv = $#ARGV+1;
if ($numArgv < 4){
    print "$usage\n";
    exit 1;
}
my ($query_file, $DB_fasta, $DB_3line, $outpath) = @ARGV;

my $basename_queryfile = basename($query_file);

if (! -d $outpath){
    `mkdir -p $outpath`;
}

my ($query_id, $query_seq, $query_len, $struct_id, $struct_seq, $seq_length, $struct_top);
my %hash_seq=();
my %hash_top=();
my %hash_len=();
my @hit_lengths=();

my $alignment_file = "$outpath/$basename_queryfile.aln";
my $alignment_file_one_hit = "$outpath/$basename_queryfile.aln_single";
my $alignment_AFA_file = "$outpath/$basename_queryfile.afa";
my $alignment_AFA_oneline_file = "$outpath/$basename_queryfile.afa.oneline";
my $final_topo_file = "$outpath/$basename_queryfile.homo_top";

open (HOMO, ">$final_topo_file");
open (ALN_SIN, ">$alignment_file_one_hit");

open (IN, $query_file);
while(<IN>)
{
	if($_=~/^>(.*)/)
	{
		$query_id = $1;
		$query_seq = <IN>;
			chomp $query_seq;
		$query_len = length($query_seq);
	}
}
close IN;

open (STRUCTURES, $DB_3line);
while(<STRUCTURES>)
{
	if($_=~/^>(.*)/)
	{
		$struct_id=$1;
		$struct_seq=<STRUCTURES>;
			$seq_length=length($struct_seq);
		$struct_top=<STRUCTURES>;
			chomp $struct_seq;
			chomp $struct_top;

		$hash_seq{$struct_id} = $struct_seq;
		$hash_top{$struct_id} = $struct_top;
		$hash_len{$struct_id} = $struct_seq;
	}
}
close STRUCTURES;

#run JACKHMMER
system("jackhmmer -E 0.001 -A $alignment_file --noali --acc -o /dev/null $query_file $DB_fasta");

#esl-reformat to transform the alignment from jackhmmer into AFA format
system ("$rundir/esl-reformat afa $alignment_file > $alignment_AFA_file");

#AFA-formatted alignment to oneline
multi2oneline ($alignment_AFA_file, $alignment_AFA_oneline_file);

#keep only the best hit (top row) + query
my $best_hit_lines = `head -2 $alignment_AFA_oneline_file`;
my $query_lines = `tail -2 $alignment_AFA_oneline_file`;

print ALN_SIN $best_hit_lines, $query_lines;
close ALN_SIN;


#map the topology to the HIT protein first
my $final_top_3D='';
my $final_query = '';
my $pdb_code = "";

open (FINAL_ALN, $alignment_file_one_hit);
while(<FINAL_ALN>)
{

    if($_=~/^>(3D.*?\|(.*?))\/(\d+)\-(\d+)/)
	{
		my $hit_id_wanted = $1;
		$pdb_code = $2;
		my $hit_id_total = $_;
		my $hit_seq_start = $3-1;
		my $hit_seq_end = $4;

		my $hit_seq_len = $hit_seq_end - $hit_seq_start;

		#keep only substring of initial hit sequence
		my $respective_hit_seq = substr($hash_top{$hit_id_wanted}, $hit_seq_start, $hit_seq_len );
		my $hit_seq=<FINAL_ALN>;
			chomp $hit_seq;

		my $query_id_fin=<FINAL_ALN>;
		my $query_seq=<FINAL_ALN>;
			chomp $query_seq;

		my @split_hit_seq=();
		my @split_topo_hit_seq=();
		my @final_topo=();
		my ($a, $b, $c);

		@split_hit_seq = split(//, $hit_seq);
		@split_topo_hit_seq = split(//,$respective_hit_seq);

		my $count_pos = -1;
		my $count_top = -1;
                                                          
		foreach $b(@split_hit_seq)
		{
			if ($b=~/-/)
			{
				$count_pos++;
				$final_topo[$count_pos]=$b;
			}
			elsif ($b!~/-/)	 
			{			
				$count_pos++;
				$count_top++;
				$final_topo[$count_pos]=$split_topo_hit_seq[$count_top];
			}
		}
		
		my $final_3D = join('', @final_topo); 	#this is the final version of 3D topology (with "-") | subject to change the terminals
		$final_3D=~s/I/i/g;
		$final_3D=~s/O/o/g;
		
		#keep only "significant gaps" and ignore everything which is mapped to a gap in the query seq
		my @split_query = split(//, $query_seq);
		my @split_final_top = split(//, $final_3D);
		my @final_hmg_array_seq = ();
		my @final_hmg_array_top = ();

		for ($c=0; $c<=$#split_final_top; $c++)
		{
			if($split_query[$c] ne '-' )
			{
				push @final_hmg_array_seq, $split_query[$c];
				push @final_hmg_array_top, $split_final_top[$c];
			}
		}

		$final_query = join('', @final_hmg_array_seq);
		$final_top_3D = join('', @final_hmg_array_top);
	}

	else
	{
		$final_top_3D = "No homologous hits found";
	}

}
#print HOMO $final_query, "\n", $final_top_3D."\n";
print HOMO ">$pdb_code\n".$final_top_3D."\n";
close FINAL_ALN;
close HOMO;

#remove unused files
#unlink($alignment_AFA_oneline_file);
unlink($alignment_AFA_file);

sub multi2oneline 
{
	my ($input, $output)=@_;
	open (ONELINE_OUT, ">$output");
	open(IN,"$input");

	my $line = <IN>; 
	print ONELINE_OUT $line;

	while ($line = <IN>)
	{
		chomp $line;
		if ($line=~m/^>/) { print ONELINE_OUT "\n",$line,"\n"; }
		else { print ONELINE_OUT $line; }
	}	

	print ONELINE_OUT "\n";
	close IN;
	close ONELINE_OUT;
}

=pod
In the end we have:
1) [query_file].aln : all proteins that had a significant hit and their alignment with the query
2) [query_file].aln_single : only the top hit aligned with our query seq
3) [query_file].homo_top : only the topology inferred from homology
