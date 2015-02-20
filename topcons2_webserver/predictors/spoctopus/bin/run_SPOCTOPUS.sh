#!/bin/bash


if [ ${#6} -lt 1 ]; then
    /bin/echo "run_SPOCTOPUS.sh <protnamefile> <raw_profiledir> <pssm_profile_dir> <outdir> <tmpdir> <octopusdir>"
    exit
fi

#args
protnamefile=$1
raw_prfdir=$2
pssm_prfdir=$3
outdir=$4
tmpdir=$5
octopusdir=$6


pid=$$

################ DIRECTORIES AND PATHS THAT SHOULD BE SET BY USER #################
octopusdir=.
###################################################################################

workingdir=$tmpdir/spoctopus_tmp_`basename $protnamefile .seq`;
if [ ! -d $workingdir ]; then
    /bin/mkdir $workingdir;
fi

modhmmdir=${octopusdir}/modhmmts_spoctopus

basedir=${workingdir}/TEMP

#hmms
hmm=NNZHMM_io_reent_hairpin_sigpep

/bin/mkdir $basedir
inputdir=$basedir/INPUT_FILES
/bin/mkdir $inputdir


#global variables
WINSIZE_MLGR_1=29
WINSIZE_S=29

WINSIZE_M_2=39
WINSIZE_G_2=51

inputdir_MLGR_1=$inputdir/WINSIZE_MLGR_1
inputdir_S=$inputdir/WINSIZE_S

/bin/mkdir $inputdir_MLGR_1
/bin/mkdir $inputdir_S

#make NN-inputfiles for M,L,G,R
/bin/mkdir $basedir/PSSMDIR
/bin/mkdir $basedir/RAWDIR

/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cp $pssm_prfdir/xxx.prf  $basedir/PSSMDIR/
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cp $raw_prfdir/xxx.prf  $basedir/RAWDIR/
raw_prfdir=$basedir/RAWDIR
pssm_prfdir=$basedir/PSSMDIR


/usr/bin/perl $octopusdir/bin/prf2NN_input_fast_scoring.pl $protnamefile $pssm_prfdir $inputdir_MLGR_1 $WINSIZE_MLGR_1
/usr/bin/perl $octopusdir/bin/prf2NN_input_fast_scoring.pl $protnamefile $raw_prfdir $inputdir_S $WINSIZE_S

/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cat $inputdir_MLGR_1/xxx.input > $basedir/input_MLGR.dat
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cat $inputdir_S/xxx.input > $basedir/input_S.dat

#run networks for MLGR
/bin/mkdir $basedir/NN_SCORE_FILES
/bin/mkdir $basedir/NN_RESULT_FILES
for i in M L G R
  do
  $octopusdir/score_nn/pred $octopusdir/NN/${i}.net.ascii $basedir/input_MLGR.dat > $basedir/NN_RESULT_FILES/${i}.res
done
$octopusdir/score_nn/pred $octopusdir/NN/S.net.ascii $basedir/input_S.dat > $basedir/NN_RESULT_FILES/S.res


# make NN-profiles
for i in M L G R
  do
  /bin/mkdir $basedir/NN_PRED_PROFILES_${i}
  /usr/bin/perl $octopusdir/bin/NN_pred2prf.pl $protnamefile $basedir/NN_RESULT_FILES/${i}.res $pssm_prfdir $basedir/NN_PRED_PROFILES_${i}
  /bin/mkdir $basedir/NN_PRED_PROFILES_OLD_${i}
  /bin/cp  $basedir/NN_PRED_PROFILES_${i}/* $basedir/NN_PRED_PROFILES_OLD_${i}/
done
/bin/mkdir $basedir/NN_PRED_PROFILES_S
/usr/bin/perl $octopusdir/bin/NN_pred2prf.pl $protnamefile $basedir/NN_RESULT_FILES/S.res $pssm_prfdir $basedir/NN_PRED_PROFILES_S
/bin/mkdir $basedir/NN_PRED_PROFILES_OLD_S
/bin/cp  $basedir/NN_PRED_PROFILES_S/* $basedir/NN_PRED_PROFILES_OLD_S/



#make io inputfiles
/bin/mkdir $basedir/INPUT_FILES/I_vs_O
/usr/bin/perl $octopusdir/bin/prf2KR_WY_NN_input_scoring.pl $protnamefile $raw_prfdir $basedir/INPUT_FILES/I_vs_O 25 31

/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cat  $basedir/INPUT_FILES/I_vs_O/xxx.input > $basedir/input_i_vs_o.dat

#run io net
# c-code dysfunctional, see run_OCTOPUS.sh
$octopusdir/bin/score_net.pl $octopusdir/NN/io.net.ascii $basedir/input_i_vs_o.dat $basedir/NN_RESULT_FILES/io.res

#make NN-profiles
 /bin/mkdir $basedir/NN_PRED_PROFILES_io
/usr/bin/perl $octopusdir/bin/NN_pred2prf_2.pl $protnamefile $basedir/NN_RESULT_FILES/io.res $raw_prfdir $basedir/NN_PRED_PROFILES_io


#make second layer input for MG
/bin/mkdir $basedir/INPUT_FILES/WIN_2_39
/usr/bin/perl $octopusdir/bin/NN_profiles2NN_input_scoring.pl $protnamefile $basedir/NN_PRED_PROFILES_M $basedir/NN_PRED_PROFILES_G $basedir/NN_PRED_PROFILES_L $basedir/NN_PRED_PROFILES_io $basedir/INPUT_FILES/WIN_2_39 $WINSIZE_M_2
/bin/mkdir $basedir/INPUT_FILES/WIN_2_51
/usr/bin/perl $octopusdir/bin/NN_profiles2NN_input_scoring.pl $protnamefile $basedir/NN_PRED_PROFILES_M $basedir/NN_PRED_PROFILES_G $basedir/NN_PRED_PROFILES_L $basedir/NN_PRED_PROFILES_io $basedir/INPUT_FILES/WIN_2_51 $WINSIZE_G_2

/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cat  $basedir/INPUT_FILES/WIN_2_51/xxx.input > $basedir/input_G.dat
/bin/cat $protnamefile | /usr/bin/xargs -I xxx /bin/cat  $basedir/INPUT_FILES/WIN_2_39/xxx.input > $basedir/input_M.dat

#run second layer nets
$octopusdir/score_nn/pred $octopusdir/NN2/G.net.ascii $basedir/input_G.dat > $basedir/NN_RESULT_FILES/G2.res
$octopusdir/score_nn/pred $octopusdir/NN2/M.net.ascii $basedir/input_M.dat > $basedir/NN_RESULT_FILES/M2.res


# make NN-profiles
for i in M G
  do
  /usr/bin/perl $octopusdir/bin/NN_pred2prf.pl $protnamefile $basedir/NN_RESULT_FILES/${i}2.res $pssm_prfdir $basedir/NN_PRED_PROFILES_${i}
done

/bin/mkdir  $basedir/S_PREDS
# make sigpep prediction
/usr/bin/perl $octopusdir/bin/make_simple_sigpep_pred.pl $protnamefile $basedir/NN_PRED_PROFILES_S/ 13 0.56 $basedir/S_PREDS
 
# run sigpep hmm
/usr/bin/perl $octopusdir/bin/NN_pred2prf_SCL.pl $protnamefile $basedir/NN_PRED_PROFILES_S $basedir/NN_PRED_PROFILES_S $basedir/S_PREDS
 
/bin/ls $basedir/S_PREDS/ | /bin/grep 'sigpeppa' | /usr/bin/xargs -I xxx /bin/echo $basedir/S_PREDS/xxx > $basedir/SIGPEP_HMM_SEQFILES.txt
/bin/sed -i "s|sigpeppa|prf|g" $basedir/SIGPEP_HMM_SEQFILES.txt
../modhmm/bin/modhmms_spoctopus -m $octopusdir/HMM_FILES/spoctopus.txt -s $basedir/SIGPEP_HMM_SEQFILES.txt -f prf -o $basedir -L -M DP -v --nopostout --viterbi > $basedir/outfile0.xml
/bin/cat $basedir/outfile0.xml | ../modhmm/bin/modhmmxml2res > $basedir/spoctopus.hmg.res
/usr/bin/perl $octopusdir/bin/modhmm_tm2top.pl $basedir/spoctopus.hmg.res $basedir/S_PREDS/ $basedir/S_PREDS/

# make HMM-profiles with adjusted S values
/bin/mkdir $basedir/HMM_PROFILES_MIOGR
/usr/bin/perl $octopusdir/bin/adjust_SMLGR_profiles.pl $protnamefile $basedir/S_PREDS $basedir/NN_PRED_PROFILES_S $basedir/NN_PRED_PROFILES_M $basedir/NN_PRED_PROFILES_L  $basedir/NN_PRED_PROFILES_G  $basedir/NN_PRED_PROFILES_R $basedir/HMM_PROFILES_MIOGR

/bin/mkdir $basedir/HMM_PROFILES_MIOGR_CLEAN
/usr/bin/perl $octopusdir/bin/nonadjust_SMLGR_profiles.pl $protnamefile $basedir/S_PREDS $basedir/NN_PRED_PROFILES_S $basedir/NN_PRED_PROFILES_M $basedir/NN_PRED_PROFILES_L  $basedir/NN_PRED_PROFILES_G  $basedir/NN_PRED_PROFILES_R $basedir/HMM_PROFILES_MIOGR_CLEAN


#make 2-alpha-files
/bin/mkdir $basedir/HMM_PROFILES_2_ALPHA
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $basedir/HMM_PROFILES_MIOGR $basedir/NN_PRED_PROFILES_io $basedir/HMM_PROFILES_2_ALPHA

#make 2-alpha-files
/bin/mkdir $basedir/HMM_PROFILES_2_ALPHA_CLEAN
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $basedir/HMM_PROFILES_MIOGR_CLEAN $basedir/NN_PRED_PROFILES_io $basedir/HMM_PROFILES_2_ALPHA_CLEAN




#run HMM
/bin/ls $basedir/HMM_PROFILES_2_ALPHA/* > $basedir/HMM_seqfiles.txt
../modhmm/bin/modhmms_spoctopus -m $octopusdir/HMM_FILES/${hmm}.txt -s $basedir/HMM_seqfiles.txt -f prf -o $basedir -L -M DP -v --nopostout -u --nolabels --viterbi > $basedir/outfile.xml
/bin/cat $basedir/outfile.xml | ../modhmm/bin/modhmmxml2res > $basedir/${hmm}.hmg.res
/bin/mkdir $basedir/DETAILED_TOPOLOGY_FILES
/usr/bin/perl $octopusdir/bin/modhmm_tm2top.pl $basedir/${hmm}.hmg.res $basedir/DETAILED_TOPOLOGY_FILES/ $basedir/DETAILED_TOPOLOGY_FILES/



# lonely revise
/bin/mkdir $basedir/HMM_PROFILES_MIOGR_lonely_revised
/usr/bin/perl $octopusdir/bin/lonely_revise_sigpep.pl $protnamefile $basedir/HMM_PROFILES_MIOGR $basedir/DETAILED_TOPOLOGY_FILES/ $basedir/HMM_PROFILES_MIOGR_lonely_revised
/bin/mkdir $basedir/HMM_PROFILES_MIOGR_lonely_revised_CLEAN
/usr/bin/perl $octopusdir/bin/lonely_revise_sigpep.pl $protnamefile $basedir/HMM_PROFILES_MIOGR_CLEAN $basedir/DETAILED_TOPOLOGY_FILES/ $basedir/HMM_PROFILES_MIOGR_lonely_revised_CLEAN


/bin/mkdir $basedir/HMM_PROFILES_2_ALPHA_lonely_revised
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $basedir/HMM_PROFILES_MIOGR_lonely_revised $basedir/NN_PRED_PROFILES_io $basedir/HMM_PROFILES_2_ALPHA_lonely_revised
/bin/mkdir $basedir/HMM_PROFILES_2_ALPHA_lonely_revised_CLEAN
/usr/bin/perl $octopusdir/bin/merge_HMM_profile_files.pl $protnamefile $basedir/HMM_PROFILES_MIOGR_lonely_revised $basedir/NN_PRED_PROFILES_io $basedir/HMM_PROFILES_2_ALPHA_lonely_revised_CLEAN


#rerun HMM
/bin/ls $basedir/HMM_PROFILES_2_ALPHA_lonely_revised/* > $basedir/HMM_seqfiles.txt
../modhmm/bin/modhmms_spoctopus -m $octopusdir/HMM_FILES/${hmm}.txt -s $basedir/HMM_seqfiles.txt -f prf -o $basedir -L -M DP -v --nopostout --viterbi -u --nolabels > $basedir/outfile2.xml
/bin/cat $basedir/outfile2.xml | ../modhmm/bin/modhmmxml2res > $basedir/${hmm}.hmg.res

/bin/mkdir $basedir/DETAILED_TOPOLOGY_FILES_lr
/usr/bin/perl $octopusdir/bin/modhmm_tm2top.pl $basedir/${hmm}.hmg.res $basedir/DETAILED_TOPOLOGY_FILES_lr/ $basedir/DETAILED_TOPOLOGY_FILES_lr/

/bin/mkdir $basedir/DETAILED_TOPOLOGY_FILES_final
/usr/bin/perl $octopusdir/bin/post_process_preds.pl $protnamefile $basedir/DETAILED_TOPOLOGY_FILES_lr $basedir/DETAILED_TOPOLOGY_FILES_final

/bin/ls $basedir/DETAILED_TOPOLOGY_FILES_final | /usr/bin/xargs -I xxx /bin/cp $basedir/DETAILED_TOPOLOGY_FILES_final/xxx $outdir/
/bin/mkdir $basedir/NN_PRED_OUTFILES
/usr/bin/perl $octopusdir/bin/make_pred_output_sigpep.pl $protnamefile $basedir/HMM_PROFILES_2_ALPHA_CLEAN $basedir/NN_PRED_OUTFILES
/bin/ls $basedir/NN_PRED_OUTFILES/ | /usr/bin/xargs -I xxx /bin/cp $basedir/NN_PRED_OUTFILES/xxx $outdir/


