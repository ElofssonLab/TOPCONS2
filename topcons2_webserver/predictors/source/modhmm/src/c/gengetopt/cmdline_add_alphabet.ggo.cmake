package "add_alphabet@suffix@"
version "@PACKAGE_VERSION@"
purpose "add alphabet with predefined emission probabilities to a given model file"

option "hmminfile" i "modelfile (in .hmg format)" string typestr="filename" yes
option "outfile" o "model outfile" string typestr="filename" yes
option "alphafile" a "alphabet file" string typestr="filename" yes

option "verbose" v "print some information about what is going on" flag off

