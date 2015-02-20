#!/bin/bash


if [ ${#2} -lt 1 ]; then
    /bin/echo  -e -n "Usage: OCTOPUS.sh <protnamefile> <pssm prf dir> <raw prf dir> <outdir> [N/D option]\n\n  -N\t\tsave network predictions\n  -D\tdebug mode\n"
    exit 1
fi

prfFolder=$1
pssmprfdir=$2
rawprfdir=$3
outdir=$4

/bin/echo query > $prfFolder/query.fa.txt
protnamefile=$prfFolder/query.fa.txt

#args
#protnamefile=$1
#pssmprfdir=$prfFolder
#rawprfdir=$prfFolder
#outdir=$4

#timefile="$outdir/exec_times-"`basename $0`".txt";
#/bin/echo echo `/bin/date +"%T %N"` "Start of OCTOPUS.sh" > $timefile;

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
workingdir=`/bin/mktemp -d /tmp/OCTOPUS_XXXXXXXXXX` || exit 1
/bin/mkdir $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES
/bin/mkdir $workingdir/SEQNAMEFILES

#if [ ${#N} -gt 0 ]; then
#    /bin/mkdir $outdir/NN_PRF_FILES
#fi

#/bin/echo `/bin/date +"%T %N"` "Split seqname file" >> $timefile;
/usr/bin/perl $octopusdir/bin/splitseqfile.pl $protnamefile $workingdir/SEQNAMEFILES > /dev/null

#/bin/echo `/bin/date +"%T %N"` "run_OCTOPUS.sh" >> $timefile;
/bin/sh $octopusdir/bin/run_OCTOPUS.sh $protnamefile $rawprfdir $pssmprfdir $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES $workingdir > /dev/null
#cat $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES/exec_times*.txt >> $timefile;

#/bin/echo `/bin/date +"%T %N"` "Back to OCTOPUS: copy result files to outdir" >> $timefile;
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cp $workingdir/DETAILED_TOPOLOGY_FILES_final/xxx.top ${outdir}/OCTOPUS.top.tmp
mkdir ${outdir}/OCTOPUS/ > /dev/null
perl oneline.pl ${outdir}/OCTOPUS.top.tmp > ${outdir}/OCTOPUS/query.top
rm ${outdir}/OCTOPUS.top.tmp

if [ ${#N} -gt 0 ]; then
    /bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cp $workingdir/HMM_PROFILES_2_ALPHA_lonely_revised/xxx.nnprf ${outdir}/NN_PRF_FILES/
fi
if [ ${#D} -gt 0 ]; then
    #/bin/echo `/bin/date +"%T %N"` "Debug: Copy tempdir files to outdir" >> $timefile;
    /bin/cp -r $workingdir $outdir/TEMP_FILES
fi
#/bin/echo `/bin/date +"%T %N"` "Clean up" >> $timefile;
/bin/rm -r $workingdir


