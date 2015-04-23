#!/bin/bash

###############################################################################
################### Copyright 2008 University of Washington ###################
###############################################################################

## the first thing we need to do is to translate the input sequence from
## FASTA format to GMTK format ...
##	command-line argument : input FASTA file name
##	outputs : creates two output files in the same directory as the
##		the input FASTA file, one is called <name>.obs (where "name"
##		is taken from the ">" header line in the FASTA file), and
##		the other file is called "input.list"
TMPPATH=/tmp
if [ -d /scratch ]; then
    TMPPATH=/scratch
fi
infile=$1
outFolder=$2
tmpFolder=`/bin/mktemp -d $TMPPATH/runPhilius_XXXXXXXXX` || exit 1
#scratch/
tmpFile=${tmpFolder}/${infile##*/}

cp ${infile} ${tmpFile}
python src/python/faa2gmtk.py ${tmpFile}


## next we need to create two other "list" files
rm -f ${tmpFolder}/curPostpFiles.list
sed "s/\.obs/\.test\.postp/g" ${tmpFolder}/input.list > ${tmpFolder}/curPostpFiles.list
echo $tmpFolder
#rm -f params/testVitFiles.list
sed "s/\.obs/\.testVit\.anyType\.bin/g" ${tmpFolder}/input.list > ${tmpFolder}/testVitFiles.list
#exit
#date
#echo " first pass : gmtkJT "
cp gmtk/params/commonParams.txt $tmpFolder/
cp gmtk/params/model_test2.master.Viterbi $tmpFolder
sed -i 's|scratch|'$tmpFolder'|g' $tmpFolder/model_test2.master.Viterbi
bin/gmtkJT  -strFile gmtk/params/model_test1.str \
            -triFile gmtk/params/model_test1.mytrifile \
            -of1 ${tmpFolder}/input.list -nf1 0 -ni1 2 -fmt1 ascii \
            -inputMasterFile gmtk/params/model_test1.master.Viterbi \
            -inputTrainableParameters gmtk/params/learnedParamsALL.out \
            -pCliquePrintRange 1:2 \
            -cCliquePrintRange 1:1 \
            -eCliquePrintRange 1:2 \
            -doDist \
            -cptNormThreshold 10. \
            -verbosity 1 > ${tmpFolder}/runJT_test.out

#date
#echo " preparing for second pass "

python src/python/prepVEdata10.py ${tmpFolder}/input.list ${tmpFolder}/runJT_test.out test.postp

#date
#echo " second pass : gmtkViterbi "

bin/gmtkViterbiNew \
	-strFile gmtk/params/model_test2_anyType.str \
        -of1 ${tmpFolder}/input.list -nf1 0 -ni1 2 -fmt1 ascii \
        -inputMasterFile ${tmpFolder}/model_test2.master.Viterbi \
        -inputTrainableParameters gmtk/params/learnedParamsALL.out \
        -dumpNames gmtk/params/new.dumpNames.list \
        -ofilelist ${tmpFolder}/testVitFiles.list \
        -cptNormThreshold 10. \
        -verbosity 19 > ${tmpFolder}/runVit_test.out

grep "log(prob(e" ${tmpFolder}/runVit_test.out > ${tmpFolder}/pass2.probs.anyType.out

#date

#echo " writing report file "

python src/python/writeReportOnly.py ${tmpFolder}/testVitFiles.list ${tmpFolder}/runJT_test.out ${tmpFolder}/pass2.probs.anyType.out

#cat scratch/Philius.report.out 
perl philius_convert.pl $tmpFolder/Philius.report.out ${outFolder}
rm -r $tmpFolder
#/*

