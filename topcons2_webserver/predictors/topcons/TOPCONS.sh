#!/bin/bash

if [ ${#3} -lt 1 ]; then
    /bin/echo -e -n "Usage: TOPCONS.sh <protnamefile> <outdir> <predictionsdir> \n"
    exit
fi

rundir=`dirname $0`
cd $rundir
#args
protnamefile=$1
outdir=$2
predDir=$3

#a=0
#for arg in "$@"
#do
#  if [ $a -gt 1 ] ; then
#      topodirs="${topodirs} ${arg}"
#  fi
#a=$(( $a + 1 ))
#done



##################################################################################
topconsdir=.
##################################################################################

workingdir=`/bin/mktemp -d /tmp/TOPCONS_XXXXXXXXXX` || exit 1
/bin/mkdir $workingdir/SEQNAMEFILES
/bin/mkdir $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES
/bin/mkdir $workingdir/CONSENSUS_PROFILES
seqnamefile=$workingdir/seqnamefile.txt
hmm=topcons

/usr/bin/perl $topconsdir/bin/splitseqfile.pl $protnamefile $workingdir/SEQNAMEFILES
seqs=`/bin/ls $workingdir/SEQNAMEFILES | /usr/bin/xargs -i echo $workingdir/SEQNAMEFILES/{}`
for i in $seqs
do
  seq=`/bin/cat $i`
  #/usr/bin/perl $topconsdir/bin/mk_prf.pl $seq $topodirs > $workingdir/CONSENSUS_PROFILES/${seq}.prf
  python $topconsdir/topcons_make_prf.py $predDir > $workingdir/CONSENSUS_PROFILES/${seq}.prf
  /bin/ls ${workingdir}/CONSENSUS_PROFILES/${seq}* > $seqnamefile

  $topconsdir/../modhmm/bin/modhmms -m $topconsdir/HMM_FILES/${hmm}.txt -s $seqnamefile -f prf -o $workingdir -L --nolabels --viterbi --nopostout > $workingdir/outfile.xml  
  /bin/cat $workingdir/outfile.xml | $topconsdir/../modhmm/bin/modhmmxml2res > $workingdir/${hmm}.hmg.res
  /usr/bin/perl $topconsdir/bin/modhmm_tm2top.pl $workingdir/${hmm}.hmg.res  $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES > /dev/null
  /bin/cp $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES/${seq}.top ${outdir}/topcons.top.tmp
  python reliability_score.py $workingdir/CONSENSUS_PROFILES/${seq}.prf $workingdir/PREDICTED_DETAILED_TOPOLOGY_FILES/${seq}.top ${outdir}/reliability.txt
  perl oneline.pl ${outdir}/topcons.top.tmp > ${outdir}/topcons.top
  rm ${outdir}/topcons.top.tmp
  rm -r ${workingdir}
done

#/bin/mkdir $outdir/CONSENSUS_PROFILES
#/bin/cat $protnamefile | /usr/bin/xargs -i /bin/cp $workingdir/CONSENSUS_PROFILES/{}.prf $outdir/CONSENSUS_PROFILES/
