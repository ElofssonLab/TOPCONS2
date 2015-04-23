#!/bin/bash


if [ ${#2} -lt 1 ]; then
    /bin/echo  -e -n "Usage: SPOCTOPUS.sh <protnamefile> <pssm prf dir> <raw prf dir> <outdir> [N/D options]\n\n  -N\t\tsave network predictions\n  -D\tdebug mode"
    exit
fi

TMPPATH=/tmp
if [ -d /scratch ]; then
    TMPPATH=/scratch
fi

#prfFolder=$1
protnamefile=$1
pssmprfdir=$2
rawprfdir=$3
outdir=$4

# /bin/echo query > $prfFolder/query.fa.txt
# protnamefile=$prfFolder/query.fa.txt

#args
#protnamefile=$1
#pssmprfdir=$2
#rawprfdir=$3
#outdir=$4

i=0
for arg in "$@"
do
  if [ $i -gt 3 ] ; then
      args="${args}${arg}"
  fi
i=$(( $i + 1 ))
done

pid=$$

N=`/bin/echo  $args | /bin/grep N`
D=`/bin/echo  $args | /bin/grep D`

octopusdir=.
workingdir=`/bin/mktemp -d $TMPPATH/SPOCTOPUS_XXXXXXXXXX` || exit 1
/bin/mkdir -p $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES
/bin/mkdir -p $workingdir/SEQNAMEFILES

#if [ ${#N} -gt 0 ]; then
#    /bin/mkdir $outdir/NN_PRF_FILES
#fi

/usr/bin/perl $octopusdir/bin/splitseqfile.pl $protnamefile $workingdir/SEQNAMEFILES > /dev/null
seqs=`/bin/ls $workingdir/SEQNAMEFILES | /usr/bin/xargs -I xxx /bin/echo  $workingdir/SEQNAMEFILES/xxx`
for i in $seqs
do
  /bin/sh $octopusdir/bin/run_SPOCTOPUS.sh $i $rawprfdir $pssmprfdir $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES $workingdir $octopusdir > /dev/null
  /bin/cat $i | /usr/bin/xargs -I xxx /bin/cp $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES/xxx.top ${outdir}/SPOCTOPUS.top.tmp 
  mkdir -p ${outdir}/SPOCTOPUS/ > /dev/null
  perl oneline_sp.pl ${outdir}/SPOCTOPUS.top.tmp > ${outdir}/SPOCTOPUS/query.top
  rm ${outdir}/SPOCTOPUS.top.tmp
  #if [ ${#N} -gt 0 ]; then
  #    /bin/cat $i | /usr/bin/xargs -I xxx /bin/cp $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES/xxx.nnprf ${outdir}/NN_PRF_FILES/
  #fi
done

if [ ${#D} -gt 0 ]; then
    /bin/cp -r $workingdir $outdir/TEMP_FILES
fi
/bin/rm -r $workingdir


