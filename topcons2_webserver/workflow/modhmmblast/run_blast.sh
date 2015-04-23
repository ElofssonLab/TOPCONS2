#!/bin/bash


if [ ${#5} -lt 1 ]; then
    /bin/echo "run_blast.sh <protnamefile> <fastadir> <outdir> <blastall-executable> <database>"
    exit
fi

TMPPATH=/tmp
if [ -d /scratch ]; then
    TMPPATH=/scratch
fi

BLAST=$4
DATABASE=$5
################################################
workingdir=`/bin/mktemp -d $TMPPATH/run_blast_XXXXXXXXXX` || exit 1



pid=$$

#args
protnamefile=$1
fastadir=$2
outdir=$3


seqs=`/bin/cat $protnamefile`
for seq in $seqs
do
  seqfile=`/bin/ls ${fastadir}/${seq}*`
  /bin/echo $seqfile
  $BLAST -p blastp -m 6 -F F -e 1.e-5 -i ${seqfile} -d $DATABASE -b 250 > ${outdir}/${seq}.msa
done
/bin/rm -r $workingdir
