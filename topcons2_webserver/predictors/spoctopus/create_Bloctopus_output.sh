#!/bin/bash

TMPPATH=/tmp
if [ -d /scratch ]; then
    TMPPATH=/scratch
fi

protnamefile=/scratch/KOSTAS/topcons2/topcons2_webserver/predictors/spoctopus/name.list
fastadir=/scratch/KOSTAS/topcons2/topcons2_webserver/input/
workingdir=`/bin/mktemp -d $TMPPATH/BLOCTOPUS_XXXXXXXXXX` || exit 1
BLASTALL=/scratch/KOSTAS/blast-2.2.26/bin/blastall
DATABASE=/scratch/KOSTAS/topcons2/topcons2_webserver/database/blast/uniref90.fasta
BLASTPGP=/scratch/KOSTAS/blast-2.2.26/bin/blastpgp
DATABASE2=/scratch/KOSTAS/topcons2/topcons2_webserver/database/blast/uniref90.fasta
MAKEMAT=/scratch/KOSTAS/blast-2.2.26/bin/makemat

/bin/mkdir $workingdir/PSIBLAST_FILES
/bin/mkdir $workingdir/BLAST_FILES
/bin/mkdir $workingdir/RAW_PRF_FILES
/bin/mkdir $workingdir/PSSM_PRF_FILES
/bin/mkdir $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES
/bin/mkdir $workingdir/SEQNAMEFILES

/bin/sh modhmmblast/run_blast.sh $protnamefile $fastadir $workingdir/BLAST_FILES $BLASTALL $DATABASE

perl msa62fasta_oneround.pl $workingdir/BLAST_FILES/A2RMJ9.msa > $workingdir/A2RMJ9.hits.db
#sed -i 's/U/A/g' ${infile_path}.hits.db
/scratch/KOSTAS/blast-2.2.26/bin/formatdb -i $workingdir/A2RMJ9.hits.db -l /dev/null

/bin/sh modhmmblast/run_psiblast.sh $protnamefile $fastadir $workingdir $workingdir/PSIBLAST_FILES $BLASTPGP $workingdir/A2RMJ9.hits.db $MAKEMAT
/bin/sh modhmmblast/msa2prf.sh $protnamefile $fastadir $workingdir/BLAST_FILES $workingdir/RAW_PRF_FILES
/bin/sh modhmmblast/mtx2prf.sh $protnamefile $workingdir/PSIBLAST_FILES $workingdir/PSSM_PRF_FILES
