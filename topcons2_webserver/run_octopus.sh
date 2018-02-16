#!/bin/bash

rundir=`dirname $0`
mainoutdir=


usage="
Description: Run OCTOPUS using this package

usage: $0 -outpath OUTDIR FILE [FILE ...] 
            [-RM] [--debug] [--tmpdir DIR]

Created 2018-02-16, updated 2018-02-16, Nanjiang Shu

    --tmpdir    Set the path for temporary files
    -RM         Delete result file for individual sequences

Example:
    # run OCTOPUS for two sequence files test1.fasta and test2.fasta
    # the result will be output to outdir/test1/ outdir/test2/
    $0 -outpath outdir test1.fasta test2.fasta
"

if [ $# -lt 1 ];then
    echo "$usage"
    exit 1
fi

mainoutdir=
fileList=()
extra_opt=""

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        fileList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -outpath|--outpath) mainoutdir=$2;shift;;
            -debug|--debug) extra_opt="$extra_opt --debug";;
            -RM|--RM) extra_opt="$extra_opt -RM";;
            -tmpdir|--tmpdir) extra_opt="$extra_opt --tmpdir $2";shift;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

if [ "$mainoutdir" == "" ];then
    echo "outdir not set. Exit" >&2 
    exit 1
fi
numFile=${#fileList[@]}
if [ $numFile -eq 0  ]; then
    echo "Fasta seqfile file not set! Exit." >&2
    exit 1
fi


if [ ! -d $mainoutdir ];then
    mkdir -p $mainoutdir
fi


for ((i=0;i<numFile;i++));do
    file=${fileList[$i]}
    if [ ! -s $file ];then
        echo "\"$file\" does not exist. Ignore."
        continue
    fi
    basename=`basename $file`
    rootname=${basename%.*}
    outdir=$mainoutdir/$rootname
    if [ ! -d $outdir ];then
        mkdir -p $outdir
    fi
    echo "$rundir/workflow/pfam_workflow_octopus.py $file $outdir $rundir/tools/blast-2.2.26/ $rundir/database/blast/uniref90.fasta $extra_opt"
    $rundir/workflow/pfam_workflow_octopus.py $file $outdir $rundir/tools/blast-2.2.26/ $rundir/database/blast/uniref90.fasta $extra_opt
done
