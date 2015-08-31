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

if (! -d $outpath)
{
    `mkdir -p $outpath`;
}

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

my ($query_id, $query_seq, $query_len, $struct_id, $struct_seq, $seq_length, $struct_top);
my %hash_seq=();
my %hash_top=();
my %hash_len=();
my @hit_lengths=();

my $query_file_to_use = "$outpath/$basename_queryfile.fasta_to_use";
my $temp_fasta = "$outpath/$basename_queryfile.temp_fasta";
my $alignment_file = "$outpath/$basename_queryfile.aln";
my $hits_table = "$outpath/$basename_queryfile.hits_table";
my $temp_file_aln = "$outpath/$basename_queryfile.tmp";
my $temp_file_aln_all = "$outpath/$basename_queryfile.all";
my $alignment_AFA_oneline_file = "$outpath/$basename_queryfile.afa.oneline";
my $final_topo_file = "$outpath/$basename_queryfile.homo_top";
my $final_alignment_file = "$outpath/$basename_queryfile.total_aligns";
my $file_without_query = "$outpath/$basename_queryfile.aligns_no_query";

multi2oneline ($query_file, $temp_fasta);
open (FINAL_FASTA, ">$query_file_to_use");
open (IN, $temp_fasta);
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

print FINAL_FASTA ">query\n$query_seq\n";
close FINAL_FASTA;

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
system("jackhmmer -E 0.001 -A $alignment_file --tblout $hits_table --noali --acc -o /dev/null $query_file_to_use $DB_fasta");

#esl-reformat to transform the alignment from jackhmmer into AFA format
system ("$rundir/esl-reformat afa $alignment_file > $temp_file_aln");

#keep Evalue for every hit and output it to the alignments file
my %hash_evalues=();
open (TABLE_HITS, $hits_table);
while(<TABLE_HITS>)
{
	chomp $_;
	if($_=~/^3D\_/)
	{
		my @split_hit_line=split(/\s+/, $_);
		my $id_hit = $split_hit_line[0];
		my $evalue_hit = $split_hit_line[4];
		my $score_hit = $split_hit_line[5];
		$hash_evalues{$id_hit} = "Evalue: $evalue_hit | Score: $score_hit";
	}
}
close TABLE_HITS;

#append the evalues to all hits in the alignment file
open (FINAL_AFA, ">$temp_file_aln_all");
open (TMP, $temp_file_aln);
while(<TMP>)
{
	if($_=~/^>(3D.*?)\/.*/)
	{
		chomp $_;
		my $tmp_id=$1;
		print FINAL_AFA $_;
		print FINAL_AFA "\t$hash_evalues{$tmp_id}\n";
	}

	else
	{
		print FINAL_AFA $_;
	}
}
close TMP;
close FINAL_AFA;

#AFA-formatted alignment to oneline
multi2oneline ($temp_file_aln_all, $alignment_AFA_oneline_file);

#grep query id and seq and erase them from the alignment
my $query_id_aln = `tail -2 $alignment_AFA_oneline_file | head -1`;
	chomp $query_id_aln;
my $query_seq_aln = `tail -2 $alignment_AFA_oneline_file | tail -1`;
	chomp $query_seq_aln;
system ("head -n -2 $alignment_AFA_oneline_file > $file_without_query");

#map topologies to hit sequences and erase everything that 
#does not correspond to a letter in the query sequence
my $final_top_3D='';
my $final_query = '';
my $final_3D = '';
my $hit_id_total = '';

open ALN_ALL, ">$final_alignment_file";
open (ALL_SEQ_ALIGN, $file_without_query);
while(<ALL_SEQ_ALIGN>)
{
	if($_=~/^>(3D.*?\|.*?)\/(\d+)\-(\d+)/)
	{
		my $hit_id_wanted = $1;
		$hit_id_total = $_;
		my $hit_seq_start = $2-1;
		my $hit_seq_end = $3;
		my $hit_seq = <ALL_SEQ_ALIGN>;
			chomp $hit_seq;

		my $hit_seq_len = $hit_seq_end - $hit_seq_start;

		#keep only substring of initial hit sequence
		my $respective_hit_seq = substr($hash_top{$hit_id_wanted}, $hit_seq_start, $hit_seq_len );
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
		$final_3D = join('', @final_topo); 	#this is the final version of 3D topology (with "-") | subject to change the terminals
		$final_3D=~s/I/i/g;
		$final_3D=~s/O/o/g;
		$final_3D=~s/X/u/g;
				
		#keep only "significant gaps" and ignore everything which is mapped to a gap in the query seq
		my @split_query = split(//, $query_seq_aln);
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

	while($final_top_3D=~/[iou-](M+\-{0,}M+)[iou-]/g)
	{
		my $TM_seg=$1;
		my $TM_seg_len=length($TM_seg);

		if($TM_seg_len<5)
		{
			my $TM_replacement = '-' x $TM_seg_len;
			$final_top_3D=~s/$TM_seg/$TM_replacement/;
		}
	}
	print ALN_ALL $hit_id_total, $final_top_3D, "\n",$query_id_aln, "\n", "$final_query", "\n\n";
}
close ALL_SEQ_ALIGN;
close ALN_ALL;

#print only the best hit in the final file for web server
my $best_hit = `head -2 $final_alignment_file`;
if($best_hit=~/^>3D_.*?\|(.*?)\/.*?\n(.*)/) 
{
	my $pdb_code = $1;
	my $pdb_top = $2;

	open (HOMO, ">$final_topo_file");
# 	print HOMO $query_id,"\n",">".$pdb_code."\n".$pdb_top."\n";
	print HOMO ">".$pdb_code."\n".$pdb_top."\n";
	close HOMO;
}

else
{
	open (HOMO, ">$final_topo_file");
# 	print HOMO $query_id,"\n","No homologous hits found\n";
	print HOMO "No homologous hits found\n";
	close HOMO;
}

#remove unused files
unlink ($temp_file_aln);
unlink ($temp_file_aln_all);
unlink ($alignment_file);
unlink ($hits_table);
unlink ($alignment_AFA_oneline_file);
unlink ($file_without_query);
unlink ($temp_fasta);
unlink ($query_file_to_use);

=pod
In the end we have:
1) [query_file].aln_all 	: all proteins that had a significant hit and their alignment with the query (with their topologies mapped on it)
2) [query_file].homo_top 	: only the topology inferred from homology (from the top hit)
