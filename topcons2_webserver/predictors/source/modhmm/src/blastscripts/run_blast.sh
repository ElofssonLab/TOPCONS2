#!/bin/bash


if [ ${#4} -lt 1 ]; then
    echo "run_blast.sh <protnamefile> <fastadir> <outdir> <database>"
    exit
fi



DATABASE=$4
################################################
workingdir=`mktemp -d /tmp/run_blast_XXXXXXXXXX` || exit 1



pid=$$

#args
protnamefile=$1
fastadir=$2
outdir=$3

seqs=`cat $protnamefile`
for seq in $seqs
do
  seqfile=`ls ${fastadir}/${seq}*`
  echo $seqfile
  @BLASTALL_EXECUTABLE@ -p blastp -m 6 -F F -e 1.e-5 -i ${seqfile} -d $DATABASE -b 250 > ${outdir}/${seq}.msa
done
rm -r $workingdir