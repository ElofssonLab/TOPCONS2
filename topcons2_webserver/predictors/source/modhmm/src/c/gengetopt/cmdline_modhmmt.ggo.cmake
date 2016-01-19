package "modhmmt@suffix@"
version "@PACKAGE_VERSION@"
purpose "train a hidden markov model"

#compulsory options
option "hmminfile" i "modelfile (in .hmg format)" string typestr="filename" yes
option "seqnamefile" s "sequence namefile (for sequences in fasta, smod, msamod or prfmod format)" string typestr="filename" yes 
option "seqformat" f "format of input sequences (fa=fasta, s=smod, msa=msamod, prf=prfmod)" string yes
option "outfile" o "model outfile" string typestr="filename" yes

#non compulsory options
option "freqfile" q "background frequency file" string typestr="filename" no
option "smxfile" x "substitution matrix file" string typestr="filename" no	
option "replfile" r "replacement letter file" string typestr="filename" no
option "alg" a "training algorithm (cml=conditional maximum likelihood, bw=baum-welch (default), disc=discriminative training)" string no
option "negseqnamefile" n "sequence namefile for negative training sequences (for sequences in fasta, smod, msamod
or prfmod format)" string typestr="filename" no 
option "optalpha" z "alphabet to optimize (parameters for transitions and all other alphabets will not be changed" string no

#non compulsory msa only options
option "msascoring" M "scoring method for alignment and profile data options = DP/DPPI/GM/GMR/DPPI/PI/PIS default=GM" string no
option "usecolumns" c "specify which columns to use for alignment input data, options = all/nr, where all means use all columns
and nr specifies a sequence in the alignment and the columns where this sequence have non-gap symbls are used
default = all" string no

#flags
option "nolabels" - "do not use labels even though the input sequences are labeled" flag off
option "noprior" - "do not use priors when training even though the the model file has prior files specified" flag off
option "tpcounts" - "use pseudocounts for transition parameter updates" flag off
option "epcounts" - "use pseudocounts for emission parameter updates" flag off
option "transonly" - "only update transition parameters" flag off
option "emissonly" - "only update emission parameters" flag off
option "verbose" v "print some information about what is going on" flag off


#option <long> <short> <desc> <argtype> {typestr="<type descr>"} {default="<default value>"} <required> {multiple}
#option <long> <short> <desc> flag      <onoff>
#option <long> <short> <desc> no




