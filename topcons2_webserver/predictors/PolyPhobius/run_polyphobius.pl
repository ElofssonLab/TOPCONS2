#we start from when the pipeline creates the blast db and formats it
#we just need to supply the fasta file, check for changing the 5th line
use File::Basename;

#$in_file=$ARGV[0]; #fasta input file
$folder=$ARGV[0];#dirname($in_file);
$outFolder=$ARGV[1];

$in_file="$folder/query.fa";
$blast_db="$folder/query.hits.db"; #uniref.90
$blast_idx="$folder/polyphobius/query.hits.index"; #uniref.90.index
$msa_file="$folder/polyphobius/query.hits.msa";
$bg_file="$folder/polyphobius/query.hits.bg";
$pred_file="$folder/polyphobius/query.poly";

mkdir "$folder/polyphobius/";

print `./blastget -ix $blast_idx -create $blast_db`;
print `perl blastget -db $blast_db -ix $blast_idx $in_file > $bg_file`;
$count_blast_hits = `grep -c "^>" $bg_file`; 
chomp ($count_blast_hits);

if($count_blast_hits>1) #I have a lot of hits, run kalign
{
	print `/usr/bin/kalign -q -f fasta -input $bg_file -output $msa_file`; #kalign run
	print `perl jphobius -poly $msa_file > $pred_file`; #topology prediction
}

else	#No hits, get prediction without homology
{
	print `perl jphobius -poly $in_file > $pred_file`; 
}

#parse the output from PolyPhobius
mkdir("$outFolder/PolyPhobius");
open FINAL, ">$outFolder"."PolyPhobius/query.top";

$final_pred="";
open PREDICTION, $pred_file;
while(<PREDICTION>)
{
	while($_=~/^FT\s+(.*)/mg)
	{
		$topo_part=$1;
		
		if($topo_part=~/^SIGNAL\s+(\d+)\s+(\d+)/)
		{ $len_part = ($2-$1)+1; $sig_part= 'S' x $len_part; $final_pred.=$sig_part;}

		elsif($topo_part=~/^TOPO_DOM\s+(\d+)\s+(\d+)\s+CYTOPLASMIC\./)
	{ $len_part = ($2-$1)+1; $in_part='i' x $len_part; $final_pred.=$in_part;}

		elsif($topo_part=~/^TOPO_DOM\s+(\d+)\s+(\d+)\s+NON CYTOPLASMIC\./)
		{ $len_part = ($2-$1)+1; $out_part='o' x $len_part; $final_pred.=$out_part;}

		elsif($topo_part=~/^TRANSMEM\s+(\d+)\s+(\d+)/)
		{ $len_part = ($2-$1)+1; $tm_part='p' x $len_part; $final_pred.=$tm_part;}
        
	}
}

if($final_pred!~/p/)
{
    $final_pred=~s/i/o/g;
}

print FINAL $final_pred."\n";
close (PREDICTION);
close (FINAL);
#unlink ($blast_idx);
#unlink ($pred_file);
#unlink ($bg_file);
#unlink ($msa_file);
