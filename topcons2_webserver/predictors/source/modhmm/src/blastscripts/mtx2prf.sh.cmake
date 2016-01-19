#!/bin/bash


if [ ${#3} -lt 1 ]; then
    echo "mtx2prf.sh <protnamefile> <mtxdir> <outdir>"
    exit
fi


####### hard coded directories and progs###########
modhmmblastdir=@CMAKE_INSTALL_PREFIX@/modhmmblast
################################################
workingdir=`mktemp -d /tmp/run_mtx2prf_XXXXXXXXXX` || exit 1



pid=$$

#args
protnamefile=$1
mtxdir=$2
outdir=$3


#cd $outdir
ls $mtxdir | grep '.mtx' | xargs -I xxx /usr/bin/perl $modhmmblastdir/mtx2prf.pl $mtxdir/xxx $outdir

rm -r $workingdir