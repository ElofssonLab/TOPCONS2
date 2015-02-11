#!/bin/bash


if [ ${#7} -lt 1 ]; then
    /bin/echo "run_psiblast.sh <protnamefile> <fastadir> <outdir> <tmpdir> <blastpgp-executable> <database> <makemat-executable>"
    exit
fi


PSIBLAST=$5
MAKEMAT=$7
DATABASE=$6
workingdir=$3/psiblast_tmp
if [ ! -d $workingdir ]; then
    /bin/mkdir $workingdir;
fi


pid=$$

#args
protnamefile=$1
fastadir=$2
outdir=$4


basedir=$workingdir/BLAST_TEMP_$pid
/bin/mkdir $basedir
/bin/mkdir $basedir/CHECK_FILES


seqs=`/bin/cat $protnamefile`
for seq in $seqs
do
  seqfile=`/bin/ls ${fastadir}/${seq}*`
  $PSIBLAST -j 2 -m 6 -F F -e 1.e-5 -i ${seqfile} -d $DATABASE -C ${basedir}/CHECK_FILES/${seq}.chk > /dev/null
  /bin/cp ${seqfile} ${basedir}/CHECK_FILES/${seq}.chd
done

cd ${basedir}/CHECK_FILES
/bin/ls | /bin/grep '.chk' > DATABASE.pn
/bin/ls | /bin/grep '.chd' > DATABASE.sn

$MAKEMAT -P DATABASE

/bin/ls | /bin/grep '.mtx' | /usr/bin/xargs -I xxx /bin/cp xxx ${outdir}/
