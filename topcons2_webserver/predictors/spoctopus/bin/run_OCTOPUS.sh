#!/bin/bash


if [ ${#4} -lt 1 ]; then
    /bin/echo "run_OCTOPUS.sh <protnamefile> <raw_profiledir> <pssm_profile_dir> <outdir> <tmpdir> [R/H usage]"
    exit
fi

#args
protnamefile=$1
raw_prfdir=$2
pssm_prfdir=$3
outdir=$4
RH_usage=$6

#timefile="$outdir/exec_times-"`basename $0`".txt";
#/bin/echo `/bin/date +"%T %N"` "Start of run_OCTOPUS.sh" > $timefile;

pid=$$

################ DIRECTORIES AND PATHS THAT SHOULD BE SET BY USER #################
octopusdir=.
###################################################################################

workingdir=$5;
if [ ! -d $workingdir ]; then
    /bin/mkdir $workingdir
fi

#hmms
hmm=NNZHMM_io_reent_hairpin
H=`/bin/echo $ioRH_usage | /bin/grep H`
R=`/bin/echo $ioRH_usage | /bin/grep R`
if [ ${#RH_usage} -gt 0 ]; then
    hmm=NNZHMM_io
    if [ ${#R} -gt 0 ]; then
	hmm=${hmm}_reent
    fi
    if [ ${#H} -gt 0 ]; then
	hmm=${hmm}_hairpin
    fi
fi

inputdir=$workingdir/INPUT_FILES
/bin/mkdir $inputdir


#global variables
WINSIZE_MLGR_1=29

WINSIZE_M_2=39
WINSIZE_G_2=51

inputdir_MLGR_1=$inputdir/WINSIZE_MLGR_1

/bin/mkdir $inputdir_MLGR_1

#make NN-inputfiles for M,L,G,R
#/bin/echo `/bin/date +"%T %N"` "Make input file for MLGR" >> $timefile;
/usr/bin/perl $octopusdir/bin/prf2NN_input_fast_scoring.pl $protnamefile $pssm_prfdir $inputdir_MLGR_1 $WINSIZE_MLGR_1


#run networks for MLGR
#/bin/echo `/bin/date +"%T %N"` "Run networks for MLGR" >> $timefile;
/bin/mkdir $workingdir/NN_SCORE_FILES
/bin/mkdir $workingdir/NN_RESULT_FILES
for i in M L G R
  do
    /bin/cat $protnamefile | while read line; do $octopusdir/score_nn/pred $octopusdir/NN/${i}.net.ascii $inputdir_MLGR_1/$line.input > $workingdir/NN_RESULT_FILES/$line.${i}.res; done
done

# make NN-profiles
#/bin/echo `/bin/date +"%T %N"` "make NN-profiles" >> $timefile;
for i in M L G R
  do
  /bin/mkdir $workingdir/NN_PRED_PROFILES_${i}
  /bin/cat $protnamefile | /usr/bin/xargs -I xxx /usr/bin/perl $octopusdir/bin/NN_pred2prf.pl $workingdir/SEQNAMEFILES/xxx.seq $workingdir/NN_RESULT_FILES/xxx.${i}.res $pssm_prfdir $workingdir/NN_PRED_PROFILES_${i}
  /bin/mkdir $workingdir/NN_PRED_PROFILES_OLD_${i}
  /bin/cp  $workingdir/NN_PRED_PROFILES_${i}/* $workingdir/NN_PRED_PROFILES_OLD_${i}/
done


#make io inputfiles
#/bin/echo `/bin/date +"%T %N"` "make io inputfiles" >> $timefile;
/bin/mkdir $workingdir/INPUT_FILES/I_vs_O
/usr/bin/perl $octopusdir/bin/prf2KR_WY_NN_input_scoring.pl $protnamefile $raw_prfdir $workingdir/INPUT_FILES/I_vs_O 25 31

#run io net
# c-code gives wrong values in 7th decimal when using two outputs
# < 0.890549017680243 0.108731644327106
# > 0.890549019819532 0.108731642192651
# => using perl-code instead
#/bin/echo `/bin/date +"%T %N"` "run io net" >> $timefile;
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /usr/bin/perl $octopusdir/bin/score_net.pl $octopusdir/NN/io.net.ascii $workingdir/INPUT_FILES/I_vs_O/xxx.input $workingdir/NN_RESULT_FILES/xxx.io.res



#make NN-profiles
#/bin/echo `/bin/date +"%T %N"` "make NN-profiles" >> $timefile;
/bin/mkdir $workingdir/NN_PRED_PROFILES_io
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /usr/bin/perl $octopusdir/bin/NN_pred2prf_2.pl $workingdir/SEQNAMEFILES/xxx.seq $workingdir/NN_RESULT_FILES/xxx.io.res $raw_prfdir $workingdir/NN_PRED_PROFILES_io


#make second layer input for MG
#/bin/echo `/bin/date +"%T %N"` "make second layer input for MG" >> $timefile;
/bin/mkdir $workingdir/INPUT_FILES/WIN_2_39
/usr/bin/perl $octopusdir/bin/NN_profiles2NN_input_scoring.pl $protnamefile $workingdir/NN_PRED_PROFILES_M $workingdir/NN_PRED_PROFILES_G $workingdir/NN_PRED_PROFILES_L $workingdir/NN_PRED_PROFILES_io $workingdir/INPUT_FILES/WIN_2_39 $WINSIZE_M_2
/bin/mkdir $workingdir/INPUT_FILES/WIN_2_51
/usr/bin/perl $octopusdir/bin/NN_profiles2NN_input_scoring.pl $protnamefile $workingdir/NN_PRED_PROFILES_M $workingdir/NN_PRED_PROFILES_G $workingdir/NN_PRED_PROFILES_L $workingdir/NN_PRED_PROFILES_io $workingdir/INPUT_FILES/WIN_2_51 $WINSIZE_G_2

#run second layer nets
#/bin/echo `/bin/date +"%T %N"` "run second layer nets" >> $timefile;
/bin/cat $protnamefile | while read line; do $octopusdir/score_nn/pred $octopusdir/NN2/G.net.ascii $workingdir/INPUT_FILES/WIN_2_51/$line.input > $workingdir/NN_RESULT_FILES/$line.G2.res; done
/bin/cat $protnamefile | while read line; do $octopusdir/score_nn/pred $octopusdir/NN2/M.net.ascii $workingdir/INPUT_FILES/WIN_2_39/$line.input > $workingdir/NN_RESULT_FILES/$line.M2.res; done

# make NN-profiles
#/bin/echo `/bin/date +"%T %N"` "make NN-profiles" >> $timefile;
for i in M G
  do
  /bin/cat $protnamefile | /usr/bin/xargs -I xxx /usr/bin/perl $octopusdir/bin/NN_pred2prf.pl $workingdir/SEQNAMEFILES/xxx.seq $workingdir/NN_RESULT_FILES/xxx.${i}2.res $pssm_prfdir $workingdir/NN_PRED_PROFILES_${i}
done

# make HMM-profiles
#/bin/echo `/bin/date +"%T %N"` "make HMM-profiles" >> $timefile;
/bin/mkdir $workingdir/HMM_PROFILES_MIOGR
#/usr/bin/perl $octopusdir/bin/NN_pred2prf_MIOGR.pl $protnamefile $workingdir/NN_PRED_PROFILES_M $workingdir/NN_PRED_PROFILES_L $workingdir/NN_PRED_PROFILES_G $workingdir/NN_PRED_PROFILES_R $workingdir/HMM_PROFILES_MIOGR
/usr/bin/perl $octopusdir/bin/NN_pred2prf_MIOGR_with_normalization.pl $protnamefile $workingdir/NN_PRED_PROFILES_M $workingdir/NN_PRED_PROFILES_L $workingdir/NN_PRED_PROFILES_G $workingdir/NN_PRED_PROFILES_R $workingdir/HMM_PROFILES_MIOGR



#make 2-alpha-files
#/bin/echo `/bin/date +"%T %N"` "make 2-alpha-files" >> $timefile;
/bin/mkdir $workingdir/HMM_PROFILES_2_ALPHA
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $workingdir/HMM_PROFILES_MIOGR $workingdir/NN_PRED_PROFILES_io $workingdir/HMM_PROFILES_2_ALPHA


#run HMM
#/bin/echo `/bin/date +"%T %N"` "run HMM" >> $timefile;
/bin/ls $workingdir/HMM_PROFILES_2_ALPHA/* > $workingdir/HMM_seqfiles.txt
../modhmm/bin/modhmms_octopus -m $octopusdir/HMM_FILES/${hmm}.txt -s $workingdir/HMM_seqfiles.txt -f prf -o $workingdir -r $octopusdir/util/replacement_letter_multi.rpl -L -M DP -v --nopostout -u --nolabels --viterbi > $workingdir/outfile.xml
#/bin/echo `/bin/date +"%T %N"` "convert modhmmxml to top-files" >> $timefile;
/bin/cat $workingdir/outfile.xml | ../modhmm/bin/modhmmxml2res > $workingdir/${hmm}.hmg.res
/bin/mkdir $workingdir/DETAILED_TOPOLOGY_FILES
/usr/bin/perl $octopusdir/bin/modhmm_tm2top.pl $workingdir/${hmm}.hmg.res $workingdir/DETAILED_TOPOLOGY_FILES/ $workingdir/DETAILED_TOPOLOGY_FILES/



# lonely revise
#/bin/echo `/bin/date +"%T %N"` "lonely revise" >> $timefile;
/bin/mkdir $workingdir/HMM_PROFILES_MIOGR_lonely_revised
/usr/bin/perl $octopusdir/bin/lonely_revise.pl $protnamefile $workingdir/HMM_PROFILES_MIOGR $workingdir/DETAILED_TOPOLOGY_FILES/ $workingdir/HMM_PROFILES_MIOGR_lonely_revised


/bin/mkdir $workingdir/HMM_PROFILES_2_ALPHA_lonely_revised
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $workingdir/HMM_PROFILES_MIOGR_lonely_revised $workingdir/NN_PRED_PROFILES_io $workingdir/HMM_PROFILES_2_ALPHA_lonely_revised



#rerun HMM
#/bin/echo `/bin/date +"%T %N"` "rerun HMM" >> $timefile;
/bin/ls $workingdir/HMM_PROFILES_2_ALPHA_lonely_revised/* > $workingdir/HMM_seqfiles.txt
../modhmm/bin/modhmms_octopus -m $octopusdir/HMM_FILES/${hmm}.txt -s $workingdir/HMM_seqfiles.txt -f prf -o $workingdir -r $octopusdir/util/replacement_letter_multi.rpl -L -M DP -v --nopostout --viterbi -u --nolabels > $workingdir/outfile2.xml
#/bin/echo `/bin/date +"%T %N"` "convert modhmmxml to top-files" >> $timefile;
/bin/cat $workingdir/outfile2.xml | ../modhmm/bin/modhmmxml2res > $workingdir/${hmm}.hmg.res
/bin/mkdir $workingdir/DETAILED_TOPOLOGY_FILES_lr
/usr/bin/perl $octopusdir/bin/modhmm_tm2top.pl $workingdir/${hmm}.hmg.res $workingdir/DETAILED_TOPOLOGY_FILES_lr/ $workingdir/DETAILED_TOPOLOGY_FILES_lr/

#/bin/echo `/bin/date +"%T %N"` "post process preds" >> $timefile;
/bin/mkdir $workingdir/DETAILED_TOPOLOGY_FILES_final
/usr/bin/perl $octopusdir/bin/post_process_preds.pl $protnamefile $workingdir/DETAILED_TOPOLOGY_FILES_lr $workingdir/DETAILED_TOPOLOGY_FILES_final
