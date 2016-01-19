#!/bin/bash


if [ ${#4} -lt 1 ]; then
    echo "msa2prf.sh <protnamefile> <fastadir> <msadir> <outdir>"
    exit
fi


####### hard coded directories and progs###########
modhmmblastdir=@CMAKE_INSTALL_PREFIX@/modhmmblast
################################################
workingdir=`mktemp -d /tmp/run_msa2prf_XXXXXXXXXX` || exit 1



pid=$$

#args
protnamefile=$1
fastadir=$2
msadir=$3
outdir=$4


basedir=$workingdir/BLAST_TEMP_$pid
mkdir $basedir
mkdir $basedir/MOD_FILES
mkdir $basedir/MOD_FILES_QUERY
mkdir $basedir/CHECK_FILES


/usr/bin/perl ${modhmmblastdir}/msa62mod.pl ${protnamefile} 0 ${msadir} $basedir/MOD_FILES ${fastadir}
/usr/bin/perl ${modhmmblastdir}/mod2modquery.pl ${protnamefile} $basedir/MOD_FILES $basedir/MOD_FILES_QUERY
/usr/bin/perl ${modhmmblastdir}/mod_upper_caser.pl $protnamefile $basedir/MOD_FILES_QUERY
/usr/bin/perl ${modhmmblastdir}/mod2prfmod_nolabelfile.pl $protnamefile $fastadir $basedir/MOD_FILES_QUERY $outdir

rm -r $workingdir