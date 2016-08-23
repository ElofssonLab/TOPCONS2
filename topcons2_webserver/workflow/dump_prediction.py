#!/usr/bin/python
# Description:
# Dump TOPCONS2 prediction in one text file
# and also output the concensus prediction in FASTA format
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s SEQFILE PATH_RESULT OUTFILE [-odg] [-orel]
"""%(progname)

usage_ext="""
Description:
    Dump TOPCONS2 prediction in one text file and also output 
    the concensus prediction in FASTA format

    SEQFILE     The original sequence file in FASTA format used in TOPCONS2 prediction
    PATH_RESULT The output path containing the prediction from TOPCONS2
    OUTFILE     The file the result will be output. In addition, the concensus
                prediction in FASTA format will be output to ${OUTFILE}.fa
OPTIONS:
  -h, --help    Print this help message and exit
  -odg, --odg   Output also the deltaG scores
  -orel, --orel Output also the reliability scores

Created 2016-08-18, updated 2016-08-18, Nanjiang Shu
"""
usage_exp="""
Examples:

    %s input.fa out/  input.pred_topcons2.txt

"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def DumpPredictionTOPCONS2(seqfile, path_result, outfile, isWriteDG, isWriteRel):#{{{
    (seqidlist, seqannolist, seqlist) = myfunc.ReadFasta(seqfile)
    outfile_fa = "%s.fa"%(outfile)

    fpout = None
    try:
        fpout = open(outfile, "w")
    except IOError:
        print >> sys.stderr, "Failed to write to file \"%s\""%(outfile)
        return 1

    fpout_fa = None
    try:
        fpout_fa = open(outfile_fa, "w")
    except IOError:
        print >> sys.stderr, "Failed to write to file \"%s\""%(outfile_fa)
        return 1

    methodlist = ['TOPCONS', 'OCTOPUS', 'Philius', 'PolyPhobius', 'SCAMPI',
            'SPOCTOPUS', 'Homology']

    for i in xrange(len(seqidlist)):
        subdirname = "seq_%d"%(i)
        subdir = "%s/%s"%(path_result, subdirname)
        seq = seqlist[i]
        length = len(seq)
        desp = seqannolist[i]
        print >> fpout, "Sequence number: %d"%(i+1)
        print >> fpout, "Sequence name: %s"%(desp)
        print >> fpout, "Sequence length: %d aa."%(length)
        print >> fpout, "Sequence:\n%s\n\n"%(seq)
        topo_consensus = ""
        for i in xrange(len(methodlist)):
            method = methodlist[i]
            seqid = ""
            seqanno = ""
            top = ""
            if method == "TOPCONS":
                topfile = "%s/%s/topcons.top"%(subdir, "Topcons")
            elif method == "Philius":
                topfile = "%s/%s/query.top"%(subdir, "philius")
            elif method == "SCAMPI":
                topfile = "%s/%s/query.top"%(subdir, method+"_MSA")
            else:
                topfile = "%s/%s/query.top"%(subdir, method)
            if os.path.exists(topfile):
                (seqid, seqanno, top) = myfunc.ReadSingleFasta(topfile)
            else:
                top = ""
            if top == "":
                #top = "***No topology could be produced with this method topfile=%s***"%(topfile)
                top = "***No topology could be produced with this method***"

            if method == "TOPCONS":
                topo_consensus = top

            if method == "Homology":
                showtext_homo = method
                if seqid != "":
                    showtext_homo = seqid
                print >> fpout, "%s:\n%s\n\n"%(showtext_homo, top)
            else:
                print >> fpout, "%s predicted topology:\n%s\n\n"%(method, top)

            if isWriteDG:
                dgfile = "%s/dg.txt"%(subdir)
                dg_content = ""
                if os.path.exists(dgfile):
                    dg_content = myfunc.ReadFile(dgfile)
                lines = dg_content.split("\n")
                dglines = []
                for line in lines:
                    if line and line[0].isdigit():
                        dglines.append(line)
                if len(dglines)>0:
                    print >> fpout,  "\nPredicted Delta-G-values (kcal/mol) "\
                            "(left column=sequence position; right column=Delta-G)\n"
                    print >> fpout, "\n".join(dglines)

            if isWriteRel:
                reliability_file = "%s/Topcons/reliability.txt"%(subdir)
                reliability = ""
                if os.path.exists(reliability_file):
                    reliability = myfunc.ReadFile(reliability_file)
                if reliability != "":
                    print >> fpout, "\nPredicted TOPCONS reliability (left "\
                            "column=sequence position; right column=reliability)\n"
                    print >> fpout, reliability

        print >> fpout, "##############################################################################"



        # write the concensus prediction in FASTA format
        print >> fpout_fa, ">%s"%(desp)
        print >> fpout_fa, topo_consensus

    if fpout:
        try:
            fpout.close()
        except IOError:
            pass
    if fpout_fa:
        try:
            fpout_fa.close()
        except IOError:
            pass

    return 0

#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1


    isWriteDG = False
    isWriteRel = False

    i = 1
    posArgList = [] #list of positional arguments
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-odg", "--odg"]:
                isWriteDG = True
                i += 1
            elif argv[i] in ["-orel", "--orel"]:
                isWriteRel = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            posArgList.append(argv[i])
            if len(posArgList) > 3:
                print >> sys.stderr, "Error! Two many positional arguments. Exit!"
                print >> sys.stdout, usage_short
                return 1
            i += 1

    if len(posArgList) < 3:
        print >> sys.stderr, "Error! Two few positional arguments. Exit!"
        print >> fpout, usage_short
        return 1
    seqfile = posArgList[0]
    path_result = posArgList[1]
    outfile = posArgList[2]
    outfile_fa = "%s.fa"%(outfile)

    if not os.path.exists(seqfile):
        print >> sys.stderr, "seqfile \"%s\" does not exist. Exit!"%(seqfile)
        return 1
    if not os.path.exists(path_result):
        print >> sys.stderr, "path_result \"%s\" does not exist. Exit!"%(path_result)
        return 1

    DumpPredictionTOPCONS2(seqfile, path_result, outfile, isWriteDG, isWriteRel)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
