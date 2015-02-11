#!/bin/bash

set -e

if (( $# < 1 )); then
    echo "Usage: $0 <fasta file> <blast installation directory>";
    exit 1;
fi

tmpdir=$1;
infile_path=${tmpdir}/query;
#filename=`basename $infile_path`
blastdir=$2;
#tmpdir=$3;
#db=$3;

if [[ -s ${infile_path}.prf ]]; then
    echo "$0: Found file '${infile_path}.mtx', skipping BLAST-run";
    exit ;
else
    echo "No blast result found for ${infile_path}";
fi

#/scratch/chrisp/tools/blast-2.2.26/bin/rpsblast -i ${infile_path}.fa -d /scratch/chrisp/databases/deltablast/cdd_delta -v 0 -b 500 -m 6 -a 4 > ${infile_path}.blast
#print "$blastdir/bin/blastpgp -i ${infile_path}.fa -d $db -e 1e-3 -v 0 -b 500 -m 6 -a 4 > ${infile_path}.blast"
#exit 1;

#echo "python pfam_scan_to_profile.py ${infile_path}.fa $tmpdir/"

python pfam_scan_to_profile.py ${infile_path}.fa $tmpdir/
#exit;
#$bindir/msa62fasta_oneround.pl ${infile_path}.blast > ${infile_path}.hits.db


#echo "$blastdir/bin/formatdb -i ${infile_path}.hits.db -l /dev/null"
#echo "formatdb"
$blastdir/bin/formatdb -i ${infile_path}.hits.db #-l /dev/null

#echo "$blastdir/bin/blastpgp -j 2 -i ${infile_path}.fa -d ${infile_path}.hits.db -e 10 -v 0 -b 1000000 -a 4 -C ${infile_path}.chk -Q ${infile_path}.psi >/dev/null"

#echo "psi blast"
$blastdir/bin/blastpgp -j 2 -i ${infile_path}.fa -d ${infile_path}.hits.db -e 1.e-5 -v 0 -b 100 -a 4 -C ${infile_path}.chk -Q ${infile_path}.psi -m 6 -o ${infile_path}.blast  >/dev/null
#$blastdir/bin/blastpgp -j 2 -m 6 -F F -e 1.e-5 -i ${infile_path}.fa -d ${infile_path}.hits.db -C ${infile_path}.chk -a4 -Q ${infile_path}.psi -o ${infile_path}.blast  >/dev/null
#echo "Enable for Scampi"
#$bindir/msa62mod_oneround.pl ${infile_path}.blast ${infile_path}.fa > ${infile_path}.raw.prf
#echo "convert"
#python convert_psi_to_prf.py ${infile_path}.psi
mkdir $tmpdir/rawprf/
perl msa62mod_oneround_v2.pl ${infile_path}.blast > ${tmpdir}/rawprf/query.prf
cp ${tmpdir}/rawprf/query.prf ${tmpdir}/rawprf/query.raw.prf
#cp ${infile_path}.raw.prf ${infile_path}.prf
/bin/echo query.chk > ${infile_path}.pn
/bin/echo query.fa > ${infile_path}.sn
$blastdir/bin/makemat -P ${infile_path}
perl mtx2prf.pl ${infile_path}
#rm ${infile_path}.hits.db* ${infile_path}.blast formatdb.log error.log ${infile_path}.pn ${infile_path}.sn ${infile_path}.chk ${infile_path}.mtx ${infile_path}.mn ${infile_path}.aux
