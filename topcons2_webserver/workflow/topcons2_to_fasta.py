#!/usr/bin/env python
# Description:
# convert topcons result file to fasta format


# output agreement list
# SeqID maxNumIDT numpredictor RLTY
# method 1:
# output files with different agreed numbers
# compare the topology of topcons to subpredictors, 
# can be 0 identical, 1 identical, up to 5
# 
import os
import sys
import myfunc
import libtopologycmp as lcmp
progname = os.path.basename(sys.argv[0])

min_TM_overlap = 5

usage = """
Usage:  %s  topcons_result_file [-outpath DIR]

Description: 
    Extract topologies from the Topcons result file and 
    output predicted topologies of each method in FASTA format

OPTIONS:
  -outpath DIR    Set ouput path, default: the same as the result file
  -q              Quiet mode
  -m INT          Set extracting method, (default: 0)
                  method 0: all predicted topologies are output
                  method 1: only those predicted the same (same numTM and same
                            NtermStatus) for all sub-predictors are extracted.
  -h, --help      Print this help message and exit

Created 2016-08-23, updated 2016-08-23, Nanjiang Shu

Examples:
    %s query.result.txt -outpath outdir
"""%(progname, progname)
BLOCK_SIZE=100000

def PrintHelp():
    print usage

def ExtractFromTopconsResult(recordContent):#{{{
#     print
#     print "recordContent:"
#     print "=============================="
#     print recordContent
#     print "=============================="
#     print
    record = {}
    record['anno'] = ""
    record['rlty'] = -1
    topoNameList =  ['predtopo_TOPCONS', 'predtopo_OCTOPUS',
            'predtopo_Philius','predtopo_PolyPhobius','predtopo_SCAMPI_msa',
            'predtopo_SPOCTOPUS']
    for name in topoNameList:
        record[name] = ""

    lines = recordContent.split("\n")
    numLine = len(lines)
    i = 0
    while i < numLine:
#         print i
        if lines[i].find("Sequence name") == 0:
            record['anno'] = lines[i][15:]
            i += 1
        elif lines[i].find('SCAMPI predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_SCAMPI_msa'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find('Philius predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_Philius'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find('PolyPhobius predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_PolyPhobius'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find('OCTOPUS predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_OCTOPUS'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find('SPOCTOPUS predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_SPOCTOPUS'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find('TOPCONS predicted topology') == 0:
            j = 1
            while lines[i+j] != "":
                record['predtopo_TOPCONS'] += lines[i+j]
                j+=1
            i += j
        elif lines[i].find("Predicted TOPCONS reliability") == 0:
            j = 1
            sumrlty = 0.0
            cnt = 0
            while lines[i+j][0:1].isdigit():
                ss = lines[i+j].split()
                if len(ss) == 2:
                    sumrlty += float(ss[1])
                    cnt += 1
                j+=1
            if cnt > 0:
                record['rlty'] = sumrlty / (cnt)
            i += j
        else:
            i += 1
    record['rlty'] *= 100.0
    record['seqid'] = myfunc.GetSeqIDFromAnnotation(record['anno'])

    if record['predtopo_TOPCONS'] != "":
        record['seqlength'] = len(record['predtopo_TOPCONS'])
        for name in topoNameList:
            if record[name].find("No TM-regions predicted") != -1:
                record[name] = ""
        #print record
        return record
    else:
        return {}
#}}}
def IsAllIdenticalTopology(topoList):#{{{
    numSeq = len(topoList)
    if numSeq <= 1:
        return True
    else:
        posTMList = [myfunc.GetTMPosition(topo) for topo in topoList]
        NtermStateList = [ lcmp.GetNtermState(topo) for topo in topoList]
        numTMList = [len(posTM) for posTM in posTMList]

        for i in xrange(numSeq-1):
            for j in xrange(i+1, numSeq):
                if not lcmp.IsIdenticalTopology(NtermStateList[i],
                        NtermStateList[j], numTMList[i], numTMList[j],
                        posTMList[i], posTMList[j], topoList[i], topoList[j],
                        min_TM_overlap):
                    return False
        return True
#}}}


def Read_topcons_result_from_buffer(buff, recordList, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer="";
    beg=0;
    end=0;
    while 1:
        beg=buff.find("Sequence number",beg);
        if beg >= 0:
            end=buff.find("\nSequence number",beg+1);
            if end >=0:
                recordContent = buff[beg:end];
                record = ExtractFromTopconsResult(recordContent)
                if record != {}:
                    recordList.append(record)
                beg = end;
            else:
                unprocessedBuffer = buff[beg:];
                break;
        else:
            unprocessedBuffer = buff[end:];
            break;
    if isEOFreached and unprocessedBuffer:
        recordContent = unprocessedBuffer
        record = ExtractFromTopconsResult(recordContent)
        if record != {}:
            recordList.append(record)
        unprocessedBuffer = ""
    return unprocessedBuffer;
    #}}}

def Topcons2Fasta(infile, outpath):#{{{
    try:
        rootname=os.path.basename(os.path.splitext(infile)[0])
        outfile_TOPCONS = outpath + os.sep + rootname + "_TOPCONS.topo"
        outfile_OCTOPUS = outpath + os.sep + rootname + "_OCTOPUS.topo"
        outfile_SPOCTOPUS = outpath + os.sep + rootname + "_SPOCTOPUS.topo"
        outfile_SCAMPI_msa = outpath + os.sep + rootname + "_SCAMPI_msa.topo"
        outfile_Philius = outpath + os.sep + rootname + "_Philius.topo"
        outfile_PolyPhobius = outpath + os.sep + rootname + "_PolyPhobius.topo"
        outRLTYFile = outpath + os.sep + rootname + "_TOPCONS.rlty"
        fpout_TOPCONS = open(outfile_TOPCONS, "w")
        fpout_OCTOPUS = open(outfile_OCTOPUS, "w")
        fpout_SPOCTOPUS = open(outfile_SPOCTOPUS, "w")
        fpout_SCAMPI_msa = open(outfile_SCAMPI_msa, "w")
        fpout_Philius = open(outfile_Philius, "w")
        fpout_PolyPhobius = open(outfile_PolyPhobius, "w")
        fpout_rlty = open(outRLTYFile, "w")
        fpin = open(infile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            unprocessedBuffer = Read_topcons_result_from_buffer(
                    buff, recordList, isEOFreached);
            if len(recordList) > 0: 
                for record in recordList:
                    seqid = myfunc.GetSeqIDFromAnnotation(record['anno'])

                    if record['predtopo_TOPCONS'] != "":
                        fpout_TOPCONS.write(">%s predtopo_TOPCONS2 rlty=%.2f\n"%(
                            record['anno'], record['rlty']))
                        fpout_TOPCONS.write("%s\n"%record['predtopo_TOPCONS'])
                    if record['predtopo_OCTOPUS'] != "":
                        fpout_OCTOPUS.write(">%s predtopo_OCTOPUS\n"%(
                            record['anno']))
                        fpout_OCTOPUS.write("%s\n"%record['predtopo_OCTOPUS'])
                    if record['predtopo_SPOCTOPUS'] != "":
                        fpout_SPOCTOPUS.write(">%s predtopo_SPOCTOPUS\n"%(
                            record['anno']))
                        fpout_SPOCTOPUS.write("%s\n"%record['predtopo_SPOCTOPUS'])
                    if record['predtopo_SCAMPI_msa'] != "":
                        fpout_SCAMPI_msa.write(">%s predtopo_SCAMPI_msa\n"%(
                            record['anno']))
                        fpout_SCAMPI_msa.write("%s\n"%record['predtopo_SCAMPI_msa'])
                    if record['predtopo_Philius'] != "":
                        fpout_Philius.write(">%s predtopo_Philius\n"%(
                            record['anno']))
                        fpout_Philius.write("%s\n"%record['predtopo_Philius'])
                    if record['predtopo_PolyPhobius'] != "":
                        fpout_PolyPhobius.write(">%s predtopo_PolyPhobius\n"%(
                            record['anno']))
                        fpout_PolyPhobius.write("%s\n"%record['predtopo_PolyPhobius'])

                    if record['rlty'] != -100.0:
                        fpout_rlty.write("%s %.2f\n"%(seqid, record['rlty']))
            if isEOFreached == True:
                break;
        fpin.close()
        fpout_TOPCONS.close()
        fpout_PolyPhobius.close()
        fpout_OCTOPUS.close()
        fpout_Philius.close()
        fpout_SCAMPI_msa.close()
        fpout_rlty.close()

        print "Results have been output to"
        print "\t%s"%outfile_TOPCONS
        print "\t%s"%outfile_OCTOPUS
        print "\t%s"%outfile_SPOCTOPUS
        print "\t%s"%outfile_PolyPhobius
        print "\t%s"%outfile_Philius
        print "\t%s"%outfile_SCAMPI_msa
        print "\t%s"%outRLTYFile

    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        raise
#}}}
def Topcons2Fasta_method1(infile, outpath):#{{{
    try:
        rootname=os.path.basename(os.path.splitext(infile)[0])
        outfile_SPlist = outpath + os.sep + rootname + "_TOPCONS.sp_list"
        outfile_TOPCONS = outpath + os.sep + rootname + "_TOPCONS.topo"
        outfile_TOPCONS_m1 = outpath + os.sep + rootname + "_TOPCONS.m1.topo"
        outfile_TOPCONS_filterSP = outpath + os.sep + rootname + "_TOPCONS_filterSP.topo"
        outfile_agreement = outpath + os.sep + rootname + ".agreement.stat.txt"
        logfile1 = outpath + os.sep + rootname + ".m1.idt.log"
        logfile2 = outpath + os.sep + rootname + ".m1.nonidt.log"
        logfile3 = outpath + os.sep + rootname + ".m1.idt.but.numpred.lt.4.log"
        fpout_SPlist = open(outfile_SPlist, "w")
        fpout_TOPCONS = open(outfile_TOPCONS, "w")
        fpout_TOPCONS_m1 = open(outfile_TOPCONS_m1, "w")
        fpout_TOPCONS_filterSP = open(outfile_TOPCONS_filterSP, "w")
        fplog1 = open(logfile1, "w")
        fplog2 = open(logfile2, "w")
        fplog3 = open(logfile3, "w")
        fpout_agree = open(outfile_agreement, "w")
        fpin = open(infile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            unprocessedBuffer = Read_topcons_result_from_buffer(
                    buff, recordList, isEOFreached);
            if len(recordList) > 0: 
                for record in recordList:
                    if record['predtopo_TOPCONS'] == "":
                        continue
                    seqid = myfunc.GetSeqIDFromAnnotation(record['anno'])
                    topoList = []
                    topoList.append(record['predtopo_OCTOPUS'])
                    topoList.append(record['predtopo_SPOCTOPUS'])
                    topoList.append(record['predtopo_SCAMPI_msa'])
                    topoList.append(record['predtopo_PolyPhobius'])
                    topoList.append(record['predtopo_Philius'])
#                     print "======================="
#                     print seqid
#                     print topoList
#                     print "======================="
# Annotation: matchList, matching the target topology to ordered topology list
# 1 for identical, 0 for non identical and -1 for empty topology
                    (matchList, numIDTtopo, numPredictor) = lcmp.MatchTopology(
                            record['predtopo_TOPCONS'], topoList,
                            min_TM_overlap, seqid)

                    fpout_agree.write("%s\t%d\t%d"%(seqid, numIDTtopo,
                        numPredictor))
                    for tt in matchList:
                        fpout_agree.write("\t%d"%(tt))
                    fpout_agree.write("\t%6.2f"%(record['rlty']))
                    fpout_agree.write("\n")

                    msg =  ">%s TOPCONS RLTY=%.2f" \
                            "numIDTtopo=%d numPredictor=%d\n"
                    if record['predtopo_TOPCONS'].find('S')>=0:
                        pp = record['predtopo_TOPCONS'].rfind('S')
                        fpout_SPlist.write("%s %d %s\n"%(record['seqid'], pp, 'Y'))

                    if record['predtopo_TOPCONS'].find('M')>=0:
                        fpout_TOPCONS.write(msg%(record['anno'] ,
                            record['rlty'], numIDTtopo, numPredictor))
                        fpout_TOPCONS.write("%s\n"%record['predtopo_TOPCONS'])
                        if record['predtopo_TOPCONS'].find('S')>=0:
                            pp = record['predtopo_TOPCONS'].rfind('S')
                            iostat = record['predtopo_TOPCONS'][pp+1]
                            top = record['predtopo_TOPCONS'].replace('S', iostat)
                        else:
                            top = record['predtopo_TOPCONS']
                        fpout_TOPCONS_filterSP.write(msg%(record['anno'] ,
                            record['rlty'], numIDTtopo, numPredictor))
                        fpout_TOPCONS_filterSP.write("%s\n"%top)


                    if numIDTtopo == numPredictor:
                        if numPredictor >= 4:
                            msg =  ">%s TOPCONS RLTY=%.2f" \
                                    "numIDTtopo=%d numPredictor=%d\n"
                            fpout_TOPCONS_m1.write(msg%(record['anno'] ,
                                record['rlty'], numIDTtopo, numPredictor))
                            fpout_TOPCONS_m1.write("%s\n"%record['predtopo_TOPCONS'])
                            msg = "%s RLTY= %.2f numIDTtopo= %d numPredictor= %d"
                            print >> fplog1, msg%(seqid, record['rlty'],
                                numIDTtopo, numPredictor)
                        else:
                            msg = "%s RLTY= %.2f numIDTtopo= %d numPredictor= %d"
                            print >> fplog3,  msg%(seqid, record['rlty'],
                                numIDTtopo, numPredictor)
                    else:
                        msg = "%s RLTY= %.2f numIDTtopo= %d numPredictor= %d"
                        print >> fplog2, msg%(seqid, record['rlty'],
                            numIDTtopo, numPredictor)
            if isEOFreached == True:
                break;
        fpin.close()
        fpout_SPlist.close()
        fpout_TOPCONS.close()
        fpout_TOPCONS_m1.close()
        fpout_TOPCONS_filterSP.close()
        fplog1.close()
        fplog2.close()
        fplog3.close()
        fpout_agree.close()
        print "Result have been output to"
        print "\t%s"%outfile_TOPCONS
        print "\t%s"%outfile_agreement
        print "\t%s"%logfile1
        print "\t%s"%logfile2
        print "\t%s"%logfile3

    except IOError:
        msg = "Failed to read file {} in function {}"
        print >> sys.stderr, msg.format(infile, sys._getframe().f_code.co_name)
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    topcons_result_file = ""
    method  = 0

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            topcons_result_file = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-m", "--m", "-method", "--method"]:
                try:
                    method = int(argv[i+1])
                    i += 2
                except (TypeError, ValueError, IndexError):
                    print >> sys.stderr, "argument error: -m INT"
                    return 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            topcons_result_file = argv[i]
            i += 1

    if topcons_result_file == "":
        print >> sys.stderr, "Infile not set. Exit."
        return 1
    elif not os.path.exists(topcons_result_file):
        print >> sys.stderr, "topcons_result_file %s is empty.  Exit"%(
                topcons_result_file)
        return 1

    if outpath == "":
        outpath = os.path.dirname(topcons_result_file)
        if outpath == "":
            outpath = "."
    elif not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    if method == 0:
        Topcons2Fasta(topcons_result_file, outpath)
    elif method == 1:
        Topcons2Fasta_method1(topcons_result_file, outpath)
    else:
        print >> sys.stderr, "Wrong method, must be 0 or 1, but %d was set." \
                % (method)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}

if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
