#transform Philius output to S/IMO format

$IN = $ARGV[0];
$OUT = $ARGV[1];
mkdir($OUT."philius");
$outfile_res = $OUT."philius/query.top";
open INPUT, $IN;
while(<INPUT>)
{

	if ($_=~/#\s+Length\s+(\d+)/)		#get protein length
	{
		$length=$1;
	}

	if ($_=~/# Number of AAs in SP\s+(\d+)\s+[\d\.\w]+/)
	{
		$sp_length=$1;
		if ($sp_length>0)
		{ $signal=$sp_length-1; }
	}

	if($_=~/^ANNOTATION:(.*)/)
	{	
		$seq_philius=$1; chop $seq_philius;
	}
}

open FINAL, ">$outfile_res";
if($seq_philius)		#globular with SP | TM protein | TM protein with SP
{
	if ($signal)
	{
		print FINAL 'S' x $signal;
	}

	@split_seq_phil = split(/\$/, $seq_philius);
	for ($i=0; $i<=$#split_seq_phil; $i++)
	{
		$part_philius = $split_seq_phil[$i];
	
		if($part_philius=~/([I|M|O]):(\d+)-(\d+)-[\w\d\-]+/)
		{
			$lbl=$1;
			$start=$2;
			$end=$3;
				$length=$end-$start+1;

			if($lbl eq 'I') { print FINAL 'i' x $length }
			if($lbl eq 'O') { print FINAL 'o' x $length }
			if($lbl eq 'M') { print FINAL 'p' x $length }
		}
	}
}

else			#globular without SP
{
	print FINAL 'o' x $length;
}
close INPUT;
close FINAL;

