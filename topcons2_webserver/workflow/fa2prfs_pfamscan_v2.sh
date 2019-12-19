#!/bin/bash

set -e

if (( $# < 1 )); then
    echo "Usage: $0 <fasta file> <blast installation directory>";
    exit 1;
fi

tmpdir=$1;
infile_path=${tmpdir}/query;  #changed by Nanjiang 2015-02-06, added /
blastdir=$2;

if [[ -s ${infile_path}.prf ]]; then
    echo "$0: Found file '${infile_path}.mtx', skipping BLAST-run";
    exit 0;
else
    echo "No blast result found for ${infile_path}";
fi

python pfam_scan_to_profile.py ${infile_path}.fa $tmpdir/

$blastdir/bin/formatdb -i ${infile_path}.hits.db #-l /dev/null

$blastdir/bin/blastpgp -j 2 -i ${infile_path}.fa -d ${infile_path}.hits.db -e 1.e-5 -v 0 -b 100 -a 8 -C ${infile_path}.chk -Q ${infile_path}.psi -m 6 -o ${infile_path}.blast  >/dev/null

MAKEMAT=$blastdir/bin/makemat

/bin/mkdir -p $tmpdir/PSIBLAST_FILES
/bin/mkdir -p $tmpdir/BLAST_FILES

/bin/mkdir -p $tmpdir/RAW_PRF_FILES
/bin/mkdir -p $tmpdir/PSSM_PRF_FILES

/bin/mkdir -p $tmpdir/CHECK_FILES

/bin/cp ${infile_path}.fa  $tmpdir/CHECK_FILES/query.chd
/bin/cp ${infile_path}.fa  $tmpdir/query.seq
/bin/cp ${infile_path}.blast $tmpdir/BLAST_FILES/query.msa
/bin/cp ${infile_path}.psi $tmpdir/PSSM_PRF_FILES/query.psi
cp ${infile_path}.chk $tmpdir/CHECK_FILES/query.chk

echo query.chk > $tmpdir/CHECK_FILES/DATABASE.pn
echo query.chd > $tmpdir/CHECK_FILES/DATABASE.sn

$MAKEMAT -P $tmpdir/CHECK_FILES/DATABASE

/bin/sh modhmmblast/msa2prf.sh protname $tmpdir $tmpdir/BLAST_FILES $tmpdir/RAW_PRF_FILES
/bin/sh modhmmblast/mtx2prf.sh protname $tmpdir/CHECK_FILES $tmpdir/PSSM_PRF_FILES

python convert_psi_to_prf.py ${infile_path}.psi
