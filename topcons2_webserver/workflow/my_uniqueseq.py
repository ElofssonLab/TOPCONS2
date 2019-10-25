#!/usr/bin/env python
# Description:
# remove duplicated sequences from fasta files
import os
import sys
import myfunc
import md5
import subprocess
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [FILE ...] [-outpath DIR]
"""%(progname)

usage_ext="""
Description:
    remove duplicated sequences from fasta files
    the output file will be $outpath/$rootname.nr.fa

OPTIONS:
  -outpath DIR   Output the result to outpath, (default: the same as input file)
  -m id|seq      Set method, by unique seqid or unique seq, (default: id)
  -l LISTFILE    Set the listfile
  -md5 yes|no    Whether use md5 encoding, (default: yes)
                 When enabled, it may use fewer memory if -m = seq
  -h, --help     Print this help message and exit

Created 2014-12-12, updated 2014-12-12, Nanjiang Shu
"""
usage_exp="""\
Examples:
    %s test.fa -m id    # remove sequences in test.fa has duplicated
    %s                  # seqIDs and output the result to test.nr.fa
""" %(progname, wspace)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def RemoveDupSeq(infile, g_outpath, method, isUseMD5):#{{{
    if g_outpath == "":
        outpath = myfunc.my_dirname(infile)
    else:
        outpath = g_outpath
    rootname = os.path.basename(os.path.splitext(infile)[0])

    outfile = "%s%s%s"%(outpath, os.sep, rootname)

    fpout = myfunc.myopen(outfile, None, "w", False)
    if fpout == None:
        return 1

    hdl = myfunc.ReadFastaByBlock(infile)
    if hdl.failure:
        return -1

    myset = set([])

    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            if method == "id":
                key = rd.seqid
            elif method == "seq":
                if isUseMD5:
                    key = md5.new(rd.seq).digest()
                else:
                    key = rd.seq
            if not key in myset:
                myset.add(key)
                fpout.write(">%s\n%s\n"%(rd.description, rd.seq))
        recordList = hdl.readseq()

    hdl.close()
    myfunc.myclose(fpout)
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    fileListFile = ""
    fileList = []

    i = 1
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
            elif argv[i] in ["-outpath", "--outpath",]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-m", "--m"]:
                (g_params['method'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-md5", "--md5"]:
                (tmpstr, i) = myfunc.my_getopt_str(argv, i)
                if tmpstr.lower() == "yes":
                    g_params['isUseMD5'] = True
                elif tmpstr.lower() == "no":
                    g_params['isUseMD5'] = False
                else:
                    print >> sys.stderr, "Bad syntax. option -md5 must be followed by yes or no"
                    return 1
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if outpath != "" and not os.path.exists(outpath):
        cmd = ["mkdir","-p",outpath]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e
            raise

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    if len(fileList) < 1:
        print >> sys.stderr, "%s: no input file is set. exit"%(sys.argv[0])

    if not g_params['method'] in ["id","seq"]:
        print >> sys.stderr, "%s: bad method \"%s\""%(sys.argv[0], g_params['method'])


    for i in xrange(len(fileList)):
        RemoveDupSeq(fileList[i], outpath, g_params['method'], g_params['isUseMD5'])

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method'] = "id"
    g_params['isUseMD5'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
