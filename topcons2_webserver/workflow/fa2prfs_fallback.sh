#!/bin/bash

set -e

if (( $# < 2 )); then
    echo "Usage: $0 <fasta file> <blast installation directory> <database>";
    exit 1;
fi

tmpdir=$1;
infile_path=${tmpdir}/query;
blastdir=$2;
db=$3;

$blastdir/bin/blastpgp -i ${infile_path}.fa -d $db -e 1e-3 -v 0 -b 500 -m 6 -a 4 > ${infile_path}.blast
perl msa62fasta_oneround.pl ${infile_path}.blast > ${infile_path}.hits.db
#sed -i 's/U/A/g' ${infile_path}.hits.db
$blastdir/bin/formatdb -i ${infile_path}.hits.db -l /dev/null

rm ${infile_path}.blast
$blastdir/bin/blastpgp -j 2 -i ${infile_path}.fa -d ${infile_path}.hits.db -e 1.e-5 -v 0 -b 100 -a 8 -m 6 -C ${infile_path}.chk -Q ${infile_path}.psi -o ${infile_path}.blast >/dev/null

#python convert_psi_to_prf.py ${infile_path}.psi

#cp ${infile_path}.raw.prf ${infile_path}.prf

mkdir $tmpdir/rawprf/
perl msa62mod_oneround_v2.pl ${infile_path}.blast > ${tmpdir}/rawprf/query.prf
cp ${tmpdir}/rawprf/query.prf ${tmpdir}/rawprf/query.raw.prf
#cp ${infile_path}.raw.prf ${infile_path}.prf
/bin/echo query.chk > ${infile_path}.pn
/bin/echo query.fa > ${infile_path}.sn
$blastdir/bin/makemat -P ${infile_path}
echo mtx
perl mtx2prf.pl ${infile_path}


#perl msa62mod_oneround.pl ${infile_path}.blast ${infile_path}.fa > ${infile_path}.raw.prf

#/bin/echo $filename.chk > ${infile_path}.pn
#/bin/echo $filename.fa > ${infile_path}.sn
#$blastdir/bin/makemat -P ${infile_path}
#perl mtx2prf.pl ${infile_path}
#cp ${infile_path}.raw.prf ${infile_path}.prf
#rm ${infile_path}.hits.db* ${infile_path}.blast formatdb.log error.log ${infile_path}.pn ${infile_path}.sn ${infile_path}.chk ${infile_path}.mtx ${infile_path}.mn ${infile_path}.aux
