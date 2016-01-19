package "modseqalign@suffix@"
version "@PACKAGE_VERSION@"
purpose "align 2 sequences to each other using their most likely state path through a given HMM"

#compulsory options
option "hmmfile" m "model namefile for models in hmg format" string typestr="filename" yes
option "target" s "sequence namefile (for seuences in fasta, smod, msamod or prfmod format)" string typestr="filename" yes
option "template" t "sequence namefile (for seuences in fasta, smod, msamod or prfmod format)" string typestr="filename" yes 
option "seqformat" f "format of input sequences (fa=fasta, s=smod, msa=msamod, prf=prfmod)" string yes
option "outfile" o "model outfile" string typestr="filename" yes

#non compulsory options
option "freqfile" q "background frequency file" string typestr="filename" no
option "smxfile" x "substitution matrix file" string typestr="filename" no	
option "replfile" r "replacement letter file" string typestr="filename" no
option "priorfile" p "sequence prior file (for msa input files)" string typestr="filename" no



#non compulsory msa only options
option "msascoring" M "scoring method for alignment and profile data options = DP/DPPI/GM/GMR/DPPI/PI/PIS default=GM" string no
option "usecolumns" c "specify which columns to use for alignment input data, options = all/nr, where all means use all columns
and nr specifies a sequence in the alignment and the columns where this sequence have non-gap symbls are used
default = all" string no

#flags
option "nolabels" - "do not use labels even though the input sequences are labeled" flag off
option "verbose" v "print some information about what is going on" flag off




#option <long> <short> <desc> <argtype> {typestr="<type descr>"} {default="<default value>"} <required> {multiple}
#option <long> <short> <desc> flag      <onoff>
#option <long> <short> <desc> no


