#!/bin/bash

# Description: install the uniref90 database for PSIBLAST
# Author: Nanjiang Shu (nanjiang.shu@scilifelab.se)

set -e
dbname=uniref90.fasta
progname=`basename $0`
rundir=`dirname $0`
rundir=`readlink -f $rundir`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname OUTPATH

Examples:
# install the $dbname blastdb at path /data/blastdb
# the old database will be overwritten if exists
    $progname /data/blastdb
"
IsProgExist(){ #{{{
    # usage: IsProgExist prog
    # prog can be both with or without absolute path
    type -P $1 &>/dev/null \
        || { echo The program \'$1\' is required but not installed. \
        Aborting $0 >&2; exit 1; }
    return 0
}
#}}}
exec_cmd(){ #{{{
    local date=`date`
    echo "[$date] $*"
    eval "$*"
}

#}}}
if [ $# -ne 1 ]; then
    echo "$usage"
    exit 1
fi

formatdb=$rundir/tools/blast-2.2.26/bin/formatdb

IsProgExist $formatdb
IsProgExist gzip
IsProgExist readlink

isQuiet=0
outpath=$1

if [ ! -d $outpath ];then
    mkdir -p $outpath
fi

outpath=`readlink -f $outpath`

url=ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/${dbname}.gz
filename=${dbname}.gz

# check whether the database is up-to-date
IS_UP_TO_DATE=false
dbname=${dbname}

if [ -f "$outpath/$dbname" ];then
    date_local_db=$(stat -c  "%y" $outpath/$dbname | awk '{print $1}')
    date_remote_db=$(curl -I $url 2>/dev/null | grep Last-Mo | awk -F: '{print $2}')
    date_remote_db=$(/bin/date -d"$date_remote_db" +%Y-%m-%d)

    if [ "$date_remote_db" == "$date_local_db" -a "$date_remote_db" != "" ];then
        IS_UP_TO_DATE=true
    fi
fi

if [ "$IS_UP_TO_DATE" == "true" ];then
    echo "$dbname is already up to date to $date_remote_db"
    exit 0
fi

tmpdir=$(mktemp -d $outpath/tmpdir.install_${dbname}.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
trap 'rm -rf "$tmpdir"' INT TERM EXIT

cd $tmpdir

exec_cmd "curl -O $url"
exec_cmd "gzip -dN $filename"
exec_cmd "$formatdb -i ${dbname} -p T"

SUCCESS=0
if [ -s ${dbname}.00.phr ] ;then
    SUCCESS=1
fi

if [ $SUCCESS -eq 1 ];then
    cd $outpath
    mv -f $tmpdir/${dbname}* $outpath/
fi

rm -rf $tmpdir
