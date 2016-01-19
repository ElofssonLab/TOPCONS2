package "modhmms@suffix@"
version "@PACKAGE_VERSION@"
purpose "score sequences on hidden markov models"

#compulsory options
option "hmmnamefile" m "model namefile for models in hmg format" string typestr="filename" yes
option "seqnamefile" s "sequence namefile (for seuences in fasta, smod, msamod or prfmod format)" string typestr="filename" yes 
option "seqformat" f "format of input sequences (fa=fasta, s=smod, msa=msamod, prf=prfmod)" string yes


#non compulsory options
option "outpath" o "output directory" string typestr="dir" no
option "freqfile" q "background frequency file" string typestr="filename" no
option "smxfile" x "substitution matrix file" string typestr="filename" no	
option "replfile" r "replacement letter file" string typestr="filename" no
option "priorfile" p "sequence prior file (for msa input files)" string typestr="filename" no
option "nullfile" n "null model file" string typestr="filename" no

#output options
option "anchor" a "hmm=results are hmm-ancored (default), seq=results are sequence anchored" string no
option "labeloutput" L "output will print predicted labeling and posterior label probabilities" flag off
option "alignmentoutput" A "output will print log likelihood, log odds and reversi scores" flag off



#non compulsory msa only options
option "msascoring" M "scoring method for alignment and profile data options = DP/DPPI/GM/GMR/DPPI/PI/PIS default=GM" string no
option "usecolumns" c "specify which columns to use for alignment input data, options = all/nr, where all means use all columns
and nr specifies a sequence in the alignment and the columns where this sequence have non-gap symbls are used
default = all" string no

#flags
option "nolabels" - "do not use labels even though the input sequences are labeled" flag off
option "verbose" v "print some information about what is going on" flag off
@commentout_userscore_ggo@option "userscore" u "use special emission score" flag off

#algorithm usage flags
defgroup "score_algs"
groupoption "viterbi" - "Use viterbi algorithm for alignment and/or label scoring (default no)" group="score_algs" 
groupoption "nbest" - "Use n-best (=1-best) algorithm for label scoring (default yes)" group="score_algs"
groupoption "forward" - "Use forward algorithm for alignment scoring (default yes)" group="score_algs"

option "max_d" - "Retrain model on each sequence using Baum-Welch before scoring" flag off
option "savehmm" - "Save retrained HMM to file" flag off

@commentout_globmemscan_ggo@option "globmemscan" g "Run prefilter to filter out globular proteins" flag off
@commentout_xml_on_stdin_ggo@option "xml-on-stdin" X "input data as xml on stdin" flag off

section "options for specific output control"
option "path" - "Print most likely statepath" flag off
option "nopostout" - "no posterior probability information for label scoring" flag off
option "nolabelout" - "no predicted labeling for label scoring" flag off
option "nollout" - "no log likelihood score for alignment scoring" flag off
option "nooddsout" - "no log odds score for alignment scoring" flag off
option "norevout" - "no reversi score for alignment scoring" flag off

option "alignpostout" - "print posterior probability information for alignment scoring" flag off
option "alignlabelout" - "print predicted labeling for alignment scoring" flag off
option "labelllout" - "print log likelihood score for label scoring" flag off
option "labeloddsout" - "print log odds score for label scoring" flag off
option "labelrevout" - "print reversi score for label scoring" flag off



#option <long> <short> <desc> <argtype> {typestr="<type descr>"} {default="<default value>"} <required> {multiple}
#option <long> <short> <desc> flag      <onoff>
#option <long> <short> <desc> no


