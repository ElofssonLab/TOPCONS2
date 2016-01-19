#!/bin/bash


if [ ${#4} -lt 1 ]; then
    echo "run_psiblast.sh <protnamefile> <fastadir> <outdir> <database>"
    exit
fi

DATABASE=$4
################################################
workingdir=`mktemp -d /tmp/run_psiblast_XXXXXXXXXX` || exit 1

pid=$$

#args
protnamefile=$1
fastadir=$2
outdir=$3

basedir=$workingdir/BLAST_TEMP_$pid
mkdir $basedir
mkdir $basedir/CHECK_FILES


seqs=`cat $protnamefile`
for seq in $seqs
do
  seqfile=`ls ${fastadir}/${seq}*`
@BLASTPGP_EXECUTABLE@ -j 2 -m 6 -F F -e 1.e-5 -i ${seqfile} -d $DATABASE -C ${basedir}/CHECK_FILES/${seq}.chk > /dev/null
  cp ${seqfile} ${basedir}/CHECK_FILES/${seq}.chd
done

cd ${basedir}/CHECK_FILES
ls | grep '.chk' > DATABASE.pn
ls | grep '.chd' > DATABASE.sn

@MAKEMAT_EXECUTABLE@ -P DATABASE;

ls | grep '.mtx' | xargs -I xxx cp xxx ${outdir}/
rm -r $workingdir