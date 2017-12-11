#!/bin/bash


if [ ${#3} -lt 1 ]; then
    /bin/echo "mtx2prf.sh <protnamefile> <mtxdir> <outdir>"
    exit
fi


####### hard coded directories and progs###########
modhmmblastdir=modhmmblast
################################################
workingdir=`/bin/mktemp -d $TMPPATH/run_mtx2prf_XXXXXXXXXX` || exit 1



pid=$$

#args
protnamefile=$1
mtxdir=$2
outdir=$3


#cd $outdir
/bin/ls $mtxdir | /bin/grep '.mtx' | /usr/bin/xargs -I xxx /usr/bin/perl $modhmmblastdir/mtx2prf.pl $mtxdir/xxx $outdir

/bin/rm -r $workingdir
