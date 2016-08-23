# Created 2011-11-08, updated 2013-02-21, Nanjiang Shu

import os
import sys
import myfunc
GAP = '-'

def FilterSignalPeptideInTopology(topo, sp_pos):#{{{
    """
    Filter signal peptide in topology
    sp_pos: location of signal peptide
    """
    gapless_topo = topo.replace(GAP, '')
    #get position of the N-terminal TM  
    (fb, fe) = myfunc.GetFirstTMPosition(gapless_topo)
    if fb != -1:
        cov = myfunc.coverage(0, sp_pos, fb, fe)
#        print "cov = %d, LTM = %d"%(cov, fe-fb)
        if cov/float(fe-fb) > 0.5: # the first TM helix is signal peptide
# get the position of the first TM helix in aligned form
            (fb_aln, fe_aln) = myfunc.GetFirstTMPosition(topo)
            Nterm = GetNtermState(topo)
            if Nterm == 'i':
                newNterm = 'o'
            else:
                newNterm = 'i'
            newtopo = topo[:fe_aln].replace(Nterm, newNterm).replace('M',
                    newNterm) + topo[fe_aln:]
            return newtopo
        else:
            return topo
    else:
        return topo
#}}}
def RemoveUnnecessaryGap_old(seqList): #{{{
    """Remove unnecessary gaps in the alignment, i.e. colums with all gaps
    """
    numSeq = len(seqList)
    if numSeq < 1:
        return seqList

    lengthList = [len(seq) for seq in seqList]
    if max(lengthList) != min(lengthList):
        print >> sys.stderr, "Error! seqList with un-equal length."
        return seqList
    lengthAlignment = lengthList[0]
    cnt_GAP = [0]*lengthAlignment

    for i in xrange(numSeq):
        seq = seqList[i]
        for j in xrange(lengthAlignment):
            if seq[j] == GAP:
                cnt_GAP[j] += 1

    newSeqList = []
    for i in xrange(numSeq):
        seq = seqList[i]
        newseqli = []
        for j in xrange(lengthAlignment):
            if cnt_GAP[j] < numSeq:
                newseqli.append(seq[j])
        newseq = "".join(newseqli)
        newSeqList.append(newseq)
    return newSeqList
#}}}
def RemoveUnnecessaryGap(seqList): #{{{
    """Remove unnecessary gaps in the alignment, i.e. colums with all gaps
    this method is 10x faster than RemoveUnnecessaryGap_old
    There might be bugs with this function, 2014-10-02, do not use it
    """
    numSeq = len(seqList)
    if numSeq < 1:
        return seqList

    lengthList = [len(seq) for seq in seqList]
    if max(lengthList) != min(lengthList):
        print >> sys.stderr, "Error! seqList with un-equal length."
        return seqList
    lengthAlignment = lengthList[0]

    #cnt_GAP = [0]*lengthAlignment
    allGapPosSet = set(range(lengthAlignment))

    for i in xrange(numSeq):
        seq = seqList[i]
        nonGapPosList = []
        for j in allGapPosSet:
            if seq[j] != GAP:
                nonGapPosList.append(j)
        allGapPosSet = allGapPosSet - set(nonGapPosList)

    isAllGapList = []
    for j in xrange(lengthAlignment):
        if j in allGapPosSet:
            isAllGapList.append(True)
        else:
            isAllGapList.append(False)
#    isAllGapList = [True]*lengthAlignment
    trimmedSeqPosList = []
    j = 0
    while j < lengthAlignment:
        if not isAllGapList[j]:
            j1 = 0
            while j+j1 < lengthAlignment and not isAllGapList[j+j1]:
                j1 += 1
            trimmedSeqPosList.append((j,j+j1))
            j += j1
        else:
            j += 1

    newSeqList = []
    for i in xrange(numSeq):
        seq = seqList[i]
        slist = [seq[p[0]:p[1]] for p in trimmedSeqPosList]
        trimmedseq = "".join(slist)
        newSeqList.append(trimmedseq)
    return newSeqList
#}}}

def Get_IOState_upstream(topo, begin_TM):#{{{
    """
    Get inside/outside state for the loop before the current TM helix
    Input:
        topo        topology sequence of the protein
        begin_TM    sequence position at the beginning of the TM helix
                    (begin_TM, end_TM) defines the location of the TM helix 
                    in the sequence
    Output:
        state       'i' or 'o', if all gaps, return empty string ""
    """
    i = begin_TM
    while i >= 0:
        if topo[i] in ['i','o']:
            return topo[i]
        i -= 1
    return ''
#}}}
def Get_IOState_downstream(topo, end_TM):#{{{
    """
    Get inside/outside state for the loop after the current TM helix
    Input:
        topo        topology sequence of the protein
        end_TM      sequence position at the end of the TM helix
                    (begin_TM, end_TM) defines the location of the TM helix 
                    in the sequence
    Output:
        state       'i' or 'o', if all gaps, return empty string ""
    """
    i = end_TM
    length = len(topo)
    while i < length:
        if topo[i] in ['i','o']:
            return topo[i]
        i += 1
    return ''
#}}}

def WriteUnmappedRecord(record, prefix, text, fpout):#{{{
    fpout.write("%-12s : %2d "%("%s%s"%(prefix, text), record["numTMunmapped"]))
    fpout.write(" | index ")
    for j in xrange(record["numTMunmapped"]):
        fpout.write(" %2d " %(record["index"][j]))
    fpout.write(" | gap ")
    for j in xrange(record["numTMunmapped"]):
        fpout.write(" %6.3f " %(record["gapFraction"][j]))
    fpout.write(" | DG ")
    for j in range(record["numTMunmapped"]):
        fpout.write(" %6.3f " %(record["DGvalue"][j]))
    fpout.write("\n")
#}}}
def WriteAna(ana, fpout, text):#{{{
    if 'Nterm' in ana and ana["Nterm"]["numTMunmapped"] > 0:
        WriteUnmappedRecord(ana['Nterm'], 'Nterm', text, fpout)
    if 'Cterm' in ana and ana["Cterm"]["numTMunmapped"] > 0:
        WriteUnmappedRecord(ana['Cterm'], 'Cterm', text, fpout)
    if 'internal' in ana:
        for m in range(len(ana["internal"])):
            newtext = "%s-%d"%(text,m)
            WriteUnmappedRecord(ana['internal'][m], 'Inter', newtext, fpout)
#}}}
def WriteOverallInfo_pairwise(id1, id2, seqIdentity, class_global, numTM1, #{{{
        numTM2, seqLength1, seqLength2,  fpout = None, 
        uniprot2pdbMap = {}, swissprotAcSet = set([])):
    """
    Write the general information of the pairwise topology comparison
    Input:
        id1
        id2
        seqIdentity
        class_global        class of the topology comparison
        numTM1              number of TM helices for protein 1
        numTM2
        seqLength1
        seqLength2
        uniprot2pdbMap      uniprotac: {pdbid1:{'type':X-ray, 'resolution':1.5}, pdbid2:{}}
        swissprotAcSet      set([uniprotac, uniprotac])

    """
    if fpout != None:
        fpout.write ("%s %s %s %.2f %s %3d %3d %4d %4d"%(
            "PairwiseComparison:", id1, id2, seqIdentity, class_global, numTM1,
            numTM2, seqLength1, seqLength2) )

        if swissprotAcSet != set([]):
            isSwissprot1 = False
            isSwissprot2 = False
            if id1 in swissprotAcSet:
                isSwissprot1 = True
            if id2 in swissprotAcSet:
                isSwissprot2 = True
#             print "numPDBID1=%d, numPDBID2=%d"%(numPDBID1, numPDBID2)
            fpout.write("  | ")
            if isSwissprot1 ==  False and isSwissprot2 == 0:
                fpout.write("NO_SPROT")
            else:
                if isSwissprot1  and isSwissprot2:
                    fpout.write("BOTH_SPROT")
                elif isSwissprot1:
                    fpout.write("ONE_SPROT_1")
                elif isSwissprot2:
                    fpout.write("ONE_SPROT_2")

        if uniprot2pdbMap != {}:
            try:
                dt1 = uniprot2pdbMap[id1]
                numPDBID1 = len(dt1)
            except KeyError:
                numPDBID1 = 0

            try:
                dt2 = uniprot2pdbMap[id2]
                numPDBID2 = len(dt2)
            except KeyError:
                numPDBID2 = 0
            fpout.write("  | ")
            if numPDBID1 ==  0 and numPDBID2 == 0:
                fpout.write("NO_PDB")
            else:
                if numPDBID1 > 0 and numPDBID2 > 0:
                    fpout.write("BOTH_PDB")
                else:
                    fpout.write("ONE_PDB")
                if numPDBID1 > 0:
                    fpout.write(" Num_PDB_1=%d 1st=%s"%(numPDBID1, dt1.keys()[0]))
                if numPDBID2 > 0:
                    fpout.write(" Num_PDB_2=%d 1st=%s"%(numPDBID2, dt2.keys()[0]))

        fpout.write("\n")
#}}}
def WriteOverallInfo_msa(comparisonClassNameList, entryID, numSeq, numTM_IDT, #{{{
        numTMcons, numIDTTopo, indexClass, fpout):
    if fpout != None:
        tag = "MultipleComparison:"
        sizeEntryID = len(entryID)
        fpout.write("%-*s " % (len(tag), "#"))
        fpout.write("%*s %5s %6s %7s %5s" %(max(7, sizeEntryID),
            "EntryID", "Nseq", "nTMIDT", "nTMcons", "nIDT"))
        for i in range(len(comparisonClassNameList)):
            fpout.write(" %5s"% ("N%d"%(i+1)))
        fpout.write(" %6s"% ("NIDT%"))
        for i in range(len(comparisonClassNameList)):
            fpout.write(" %6s"% ("N%d%%"%(i+1)))
        fpout.write("\n")

        fpout.write("%s "%(tag))
        fpout.write("%*s %5d %6d %7d %5d" %(max(7, sizeEntryID),
            entryID, numSeq, numTM_IDT, numTMcons, numIDTTopo))
        for i in range(len(comparisonClassNameList)):
            fpout.write(" %5d"% len(indexClass[i]))
        fpout.write(" %6.2f"% (numIDTTopo/float(numSeq)*100))
        for i in range(len(comparisonClassNameList)):
            fpout.write(" %6.2f"% (len(indexClass[i])/float(numSeq)*100))
        fpout.write("\n")
#}}}
def ScanfMemberInfo(line): #{{{
    memberrecord = {}
    if not line:
        return memberrecord
    else:
        strs = line.split(':')
        ss = strs[1].split()
        memberrecord['id'] = ss[0]
        memberrecord['seqLength'] = int(ss[1])
        ss = strs[2].split()
        memberrecord['dgscore'] = [float(x) for x in ss]
        return memberrecord
#}}}
def WritePairCmpRecord(recordList,cntTotalOutputRecord, fpout):#{{{
    """ Return (status, cntTotalOutputRecord)"""
    if recordList == None or recordList == []:
        return (1, cntTotalOutputRecord)
    numRecord = len(recordList)
    for i in xrange(numRecord):
        record=recordList[i]
        if record == {}:
            continue
        print >> fpout, "//Begin record", cntTotalOutputRecord+1
        if 'general_info_line' in record and record['general_info_line'] != "":
            print >> fpout, record['general_info_line']
        else:
            WriteOverallInfo_pairwise(record['id1'], record['id2'],
                    record['seqidt1'], record['cmpclass'], record['numTM1'],
                    record['numTM2'], record['seqLength1'],
                    record['seqLength2'], fpout)
        print >> fpout, record['mapTMline'][0]
        print >> fpout, record['mapTMline'][1]
        if  record['ana1'] != {}:
            #print record['ana1']
            WriteAna(record['ana1'],fpout, "1") 
        if  record['ana2'] != {}:
            #print record['ana2']
            WriteAna(record['ana2'],fpout, "2") 
        print >> fpout, "//End record", cntTotalOutputRecord+1
        cntTotalOutputRecord += 1
    return (0, cntTotalOutputRecord)
#}}}

def ReadDupPairList(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return []

    li = []
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split()
            if len(strs) >= 2:
                if strs[1] == 'y':
                    ss = strs[0].split("-")
                    if len(ss) == 2:
                        li.append((ss[0], ss[1]))
        lines = hdl.readlines()
    hdl.close()
    return (li)
#}}}
def ParseDupHit(text):#{{{
    """
    Parse a duplication hit 
    e.g. 
    1-575(nTM=2) 2-571(nTM=2) 
    """
    li = []
    strs = text.split()
    for ss in strs:
        strs1 = ss.rstrip(')').split('(nTM=') # ['35-345', '7']
        li1 = [int(x) for x in strs1[0].split('-')]
        li1.append(int(strs1[1])) # li1 is a list of three values (begin, end, numTM)
        li.append(li1)
    return li
#}}}

def ReadDupPairDict(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}

    dt = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split()
            if line == "" or line[0] == "#":
                continue
            if len(strs) >= 2:
                if strs[1] == 'y': # it is a duplicated pair
                    ss = strs[0].split("-")
                    if len(ss) == 2:
                        key = (ss[0], ss[1])
                        dt[key] = {}
                        dt[key]['isDup'] = 'y'
                        li = []
                        strs1 = line.split('|')
                        for j in xrange(1, len(strs1)):
                            hit = ParseDupHit(strs1[j].strip()) # hit is a list of 
                                                                # two segments
                                                                # from query
                                                                # and template
                            li.append(hit)
                        dt[key]['hit'] = li
        lines = hdl.readlines()
    hdl.close()
    return (dt)
#}}}

def ReadPfamDefFile(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return ({}, {})

    dtClan = {}
    dtPfam = {}
    lines = hdl.readlines()

    while lines != None:
        for line in lines:
            strs = line.split("\t")
            try:
                pfamid = strs[0]
                pfamDefShort = strs[3]
                dtPfam[pfamid] = pfamDefShort

                clanid = strs[1]
                clanDefShort = strs[2]
                if clanid != "\N":
                    dtClan[clanid] = clanDefShort
                else:
                    dtClan[pfamid] = pfamDefShort
            except IndexError:
                pass
        lines = hdl.readlines()
    hdl.close()
    return (dtPfam, dtClan)
#}}}
def ReadSignalPDict(infile):#{{{
# format of signalp file
# SeqID location Y
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}

    signalpDict = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) >= 2:
                try:
                    signalpDict[strs[0]] = int(strs[1])
                except (ValueError):
                    pass
        lines = hdl.readlines()
    hdl.close()
    return signalpDict
#}}} 
def GetSeqIDT(record, seqidttype):#{{{
    seqidt = -1.0
    if seqidttype == 0:
        if 'seqidt' in record:
            seqidt = record['seqidt']
        else:
            print >> sys.stderr, "%s %s, no seqidt info" % (
                    record['id1'], record['id2'])
    elif seqidttype == 1:
        if 'seqidt1' in record:
            seqidt = record['seqidt1']
        else:
            print >> sys.stderr, "%s %s, no seqidt1 info" % (
                    record['id1'], record['id2'])
    elif seqidttype == 2:
        if 'seqidt2' in record:
            seqidt = record['seqidt2']
        else:
            print >> sys.stderr, "%s %s, no seqidt2 info" % (
                    record['id1'], record['id2'])
    return seqidt
#}}}

# topology related functions
def GetNtermState(topo):#{{{
    if topo[0] != GAP:
        return topo[0]
    else:
        topo=topo.lstrip(GAP)
        if topo != "":
            return topo[0]
        else:
            return None
#}}}
def GetCtermState(topo):#{{{
    length = len(topo)
    if topo[length-1] != GAP:
        return topo[length-1]
    else:
        topo = topo.rstrip(GAP)
        if topo != "":
            return topo[len(topo)-1]
        else:
            return None
#}}}
def IsIdenticalTopology(Nterm1, Nterm2, numTM1, numTM2, posTM1, posTM2, #{{{
        topo1, topo2, min_TM_overlap = 5):

    """Check whether topo1 and topo2 are identical"""
# Created 2011-11-15, updated 2011-11-15
# Two topologies are considered identical (Krogh et al. 2001) if
# 1. numTM1 == numTM2
# 2. Each helix of the compared topology should overlap by at least N (e.g. 5)
#    residues
# 3. Each helix is oriented in the same way
    if numTM1 != numTM2:
        return False
    else:
        if Nterm1 != Nterm2:
            return False
        else:
            for i in range (numTM1): 
                (b1,e1) = posTM1[i]
                (b2,e2) = posTM2[i]
                (common_b, common_e) = (max(b1,b2), min(e1,e2))
                overlap = common_e - common_b
                if overlap <= 0:
                    return False
                else:
                    cntCommonM = 0
                    for j in xrange(common_b, common_e):
                        if topo1[j] == 'M' and topo2[j] == 'M':
                            cntCommonM += 1
                        if cntCommonM >= min_TM_overlap:
                            break
#                     print ("cntCommonM=", cntCommonM, "min_TM_overlap=",
#                             min_TM_overlap)
                    if cntCommonM < min_TM_overlap:
                        return False
    return True
#}}}
def IsIdenticalTopology_simple( topo1, topo2, min_TM_overlap = 5):#{{{

    """Check whether topo1 and topo2 are identical"""
# Created 2011-11-15, updated 2011-11-15
# Two topologies are considered identical (Krogh et al. 2001) if
# 1. numTM1 == numTM2
# 2. Each helix of the compared topology should overlap by at least N (e.g. 5)
#    residues
# 3. Each helix is oriented in the same way
    numTM1 = myfunc.CountTM(topo1)
    numTM2 = myfunc.CountTM(topo2)
    Nterm1 = GetNtermState(topo1)
    Nterm2 = GetNtermState(topo2)
    posTM1 = myfunc.GetTMPosition(topo1)
    posTM2 = myfunc.GetTMPosition(topo2)

    if numTM1 != numTM2:
        return False
    else:
        if Nterm1 != Nterm2:
            return False
        else:
            for i in range (numTM1): 
                (b1,e1) = posTM1[i]
                (b2,e2) = posTM2[i]
                (common_b, common_e) = (max(b1,b2), min(e1,e2))
                overlap = common_e - common_b
                if overlap <= 0:
                    return False
                else:
                    cntCommonM = 0
                    for j in xrange(common_b, common_e):
                        if topo1[j] == 'M' and topo2[j] == 'M':
                            cntCommonM += 1
                        if cntCommonM >= min_TM_overlap:
                            break
#                     print ("cntCommonM=", cntCommonM, "min_TM_overlap=",
#                             min_TM_overlap)
                    if cntCommonM < min_TM_overlap:
                        return False
    return True
#}}}
def IsInvertedTopology(Nterm1, Nterm2, numTM1, numTM2, posTM1, posTM2, #{{{
        topo1, topo2, min_TM_overlap = 5):
    """Check whether topo1 and topo2 are inverted"""
# Created 2013-07-11, updated 2013-07-11, Nanjiang Shu 
# Two topologies are considered inverted (Krogh et al. 2001) if
# 1. numTM1 == numTM2
# 2. Each helix of the compared topology should overlap by at least N (e.g. 5)
#    residues
# 3. Each helix is oriented in different way
    if numTM1 != numTM2 or Nterm1 == Nterm2:
        return False
    else:
        for i in range (numTM1): 
            (b1,e1) = posTM1[i]
            (b2,e2) = posTM2[i]
            (common_b, common_e) = (max(b1,b2), min(e1,e2))
            overlap = common_e - common_b
            if overlap <= 0:
                return False
            else:
                cntCommonM = 0
                for j in xrange(common_b, common_e):
                    if topo1[j] == 'M' and topo2[j] == 'M':
                        cntCommonM += 1
                    if cntCommonM >= min_TM_overlap:
                        break
                if cntCommonM < min_TM_overlap:
                    return False
    return True
#}}}
def MatchTopology(targetTopo, topoList, min_TM_overlap = 5, seqid = ""):#{{{
## compare targetTopo to all topologies in the topoList
# return (matchList, numIDTtopo, numPredictor
    numList = len(topoList)
    matchList = []
# 0 for different topology
# 1 for identical topology
# -1 for empty topology

# debug
#     print "SeqID: %s"%(seqid)
#     print GetNtermState(targetTopo), myfunc.GetTMPosition(targetTopo)
#     print
#     for tt in topoList:
#         if tt:
#             print GetNtermState(tt), myfunc.GetTMPosition(tt)
#         else:
#             print "Null"
#     print

    NtermStateTarget = GetNtermState(targetTopo)
    posTMtarget = myfunc.GetTMPosition(targetTopo)
    numTMtarget = len(posTMtarget)
    for i in xrange(numList):
        if topoList[i] == "":
            matchList.append(-1)
        else:
            NtermState = GetNtermState(topoList[i])
            posTM = myfunc.GetTMPosition(topoList[i])
            numTM = len(posTM)
            if IsIdenticalTopology(NtermStateTarget, NtermState,
                    numTMtarget, numTM, posTMtarget, posTM, targetTopo,
                    topoList[i], min_TM_overlap):
                matchList.append(1)
            else:
                matchList.append(0)

    numIDTtopo = matchList.count(1)
    numPredictor = matchList.count(1) + matchList.count(0)
    return (matchList, numIDTtopo, numPredictor)
#}}}
def GetSeq2AlignMap(seq, alignedseq):#{{{
    """
    map index in the original sequence to the aligned sequence
    e.g. map[0] = 0
         map[1] = 5
    """
    mp = {}
    j = 0
    for i in xrange(len(alignedseq)):
        if alignedseq[i] != GAP:
            mp[j] = i
            j += 1
    return mp
#}}}
def GetAlign2SeqMap(alignedseq, seq):#{{{
    """
    map index in the aligned sequence to the original sequence
    e.g. map[0] = 0
         map[5] = 1
    """
    mp = {}
    j = 0
    for i in xrange(len(alignedseq)):
        if alignedseq[i] != GAP:
            mp[i] = j
            j += 1
        else:
            mp[i] = j  # assign to the neast residue to the left
    return mp
#}}}
    

def MatchToAlignedSeq(unalignedseq, alignedseq, seqID): #{{{
    """match the unaligned seq to the aligned seq, gaps are added
    return alignedseq at failure"""
    newseq = ""
    j = 0
    for i in xrange(len(alignedseq)):
        if alignedseq[i] != GAP:
            newseq += unalignedseq[j]
            j += 1
        else:
            newseq += GAP
    if len(newseq) != len(alignedseq):
        print >> sys.stderr, "failed to match sequence for ID %s" %seqID
        return alignedseq
    else:
        return newseq
#}}}

# alignment related functions
def GetAlignmentFactorFromPairAlignment(seq1,seq2, isLocalAlignment):#{{{
    """
    Return alignment factor as a dictionary
    """
    alignFactor = {}
    alnLength = len(seq1)
    if isLocalAlignment is True:
        cntLocalAlnLength = 0
        cntLocalIDT = 0
        cntLocalGap = 0
        cntLocalLen1 = 0
        cntLocalLen2 = 0
        cntLocalUnAligned = 0
        for i in xrange(alnLength):
            if ((seq1[i].isalpha() and seq1[i].islower()) or (seq2[i].isalpha()
                and seq2[i].islower())):
                cntLocalUnAligned += 1
            else:
                cntLocalAlnLength += 1
                if seq1[i] == seq2[i]:
                    cntLocalIDT += 1
                elif seq1[i] == "-" or seq2[i] == "-":
                    cntLocalGap += 1
                if seq1[i] != "-":
                    cntLocalLen1 += 1
                if seq2[i] != "-":
                    cntLocalLen2 += 1
        alignFactor['numIDT'] = cntLocalIDT
        alignFactor['numGap'] = cntLocalGap
        alignFactor['alnLength'] = cntLocalAlnLength
        alignFactor['seqidt0'] = myfunc.FloatDivision(cntLocalIDT,
                cntLocalAlnLength) * 100
        alignFactor['seqidt1'] = myfunc.FloatDivision(cntLocalIDT,
                min(cntLocalLen1,cntLocalLen2))*100
        alignFactor['seqidt2'] = myfunc.FloatDivision(cntLocalIDT,
                cntLocalAlnLength - cntLocalGap)*100
        alignFactor['seqLength1'] = cntLocalLen1
        alignFactor['seqLength2'] = cntLocalLen2
        alignFactor['numUnaligned'] = cntLocalUnAligned

    else:
        len1 = len(seq1.replace("-", ""))
        len2 = len(seq2.replace("-", ""))
        cntIDT = 0
        cntGap = 0
        for i in xrange(alnLength):
            if seq1[i] == seq2[i]:
                cntIDT += 1
            elif seq1[i] == "-" or seq2[i] == "-":
                cntGap += 1
        alignFactor['numIDT'] = cntIDT
        alignFactor['numGap'] = cntGap
        alignFactor['alnLength'] = alnLength
        alignFactor['seqidt0'] = myfunc.FloatDivision(cntIDT, alnLength) * 100
        alignFactor['seqidt1'] = myfunc.FloatDivision(cntIDT, min(len1,len2)) * 100
        alignFactor['seqidt2'] = myfunc.FloatDivision(cntIDT, alnLength - cntGap) * 100
        alignFactor['seqLength1'] = len1
        alignFactor['seqLength2'] = len2
    return alignFactor
#}}}

def ExtractFromPairCmpRecordContent(recordContent):#{{{
    """
    Extract pairwise topology comparison from the record content in the file
    *.paircmp
    updated 2011-11-21
    """
    record = {}
    lines = recordContent.split('\n')
    if len(lines) <= 1: # record is empty
        print >> sys.stderr, "record is empty\n", recordContent
        return {}

    record['mapTMline']=[]
    record['general_info_line']= ""
    record['mapArray'] = []
    record['ana1'] = {}
    record['ana2'] = {}
    record['member'] = []

    for line in lines:
        tag = myfunc.GetFirstWord(line)
        if tag == "PairwiseComparison:":
            ScanfOverallInfo_pairwise(line, record)
            record['general_info_line'] = line
        elif tag == "SeqID" or tag == "TMMap":
            record['mapTMline'].append(line)
            record['mapArray'].append([int(x) for x in
                line.split(':')[1].replace('-','').split()])
        elif tag[0:6] == "Member":
            record['member'].append(ScanfMemberInfo(line))
        elif tag == "NtermTopo1":
            record['NterTopo1'] = line.split()[1]
        elif tag == "NtermTopo2":
            record['NterTopo2'] = line.split()[1]
        elif tag == "Nterm1":
            record['ana1']['Nterm'] = ScanfUnmappedRecord(line)
        elif tag == "Nterm2":
            record['ana2']['Nterm'] = ScanfUnmappedRecord(line)
        elif tag == "Cterm1":
            record['ana1']['Cterm'] = ScanfUnmappedRecord(line)
        elif tag == "Cterm2":
            record['ana2']['Cterm'] = ScanfUnmappedRecord(line)
        elif tag[0:6] == "Inter1":
            if 'internal' not in record['ana1']:
                record['ana1']['internal'] = []
            record['ana1']['internal'].append(ScanfUnmappedRecord(line))
        elif tag[0:6] == "Inter2":
            if 'internal' not in record['ana2']:
                record['ana2']['internal'] = []
            record['ana2']['internal'].append(ScanfUnmappedRecord(line))
    return record

#}}}
def ReadPairCmpResultFromBuffer(buff,recordList):#{{{
    """ 
    Return unprocessedBuffer 
    """
    if not buff:
        return ""
    unprocessedBuffer=""
    beg=0
    end=0
    while 1:
        beg=buff.find("//Begin",beg)
        if beg >= 0:
            end=buff.find("\n//End",beg+1)
            if end >=0:
                recordContent=buff[beg:end]
                recordList.append(ExtractFromPairCmpRecordContent(recordContent))
                beg=end
            else:
                unprocessedBuffer=buff[beg:]
                break
        else:
            unprocessedBuffer=buff[end:]
            break
    return unprocessedBuffer
#}}}

def ScanfUnmappedRecord(line):#{{{
# Example:
# CtermCons    :  1  | index   1 | gap   1.000  | DG  -1.924
    subAna={}
    subAna['numTMunmapped'] = 0
    subAna['index'] = []
    subAna['gapFraction'] = []
    subAna['DGvalue'] = []

    strs=line.split('|')
    for i in range(len(strs)):
        if i == 0:
            substrs=strs[i].split(':')
            subAna['numTMunmapped'] = int(substrs[1])
        else:
            substrs=strs[i].split()
            if substrs[0] == 'index':
                for j in range(1,len(substrs)):
                    subAna['index'].append(int(substrs[j]))
            elif substrs[0] == 'gap':
                for j in range(1,len(substrs)):
                    subAna['gapFraction'].append(float(substrs[j]))
            elif substrs[0] == 'DG':
                for j in range(1,len(substrs)):
                    subAna['DGvalue'].append(float(substrs[j]))
    return subAna
#}}}
def ScanfOverallInfo_pairwise(line, cmpRecord):#{{{
    try:
        strs = line.split()
        cmpRecord['id1']=strs[1]
        cmpRecord['id2']=strs[2]
        cmpRecord['seqidt']=float(strs[3])
        if strs[4].find(";") == -1:
            cmpRecord['cmpclass'] = strs[4]
            cmpRecord['isLocalAlignment'] = False
        else:
            strs1 = strs[4].split(";")
            cmpRecord['isLocalAlignment'] = True
            cmpRecord['cmpclass'] = strs1[0]
            cmpRecord['alignrange'] = strs1[1]
            cmpRecord['aligned_numTM1'] = int(strs1[2])
            cmpRecord['aligned_numTM2'] = int(strs1[3])
            cmpRecord['numTM_unaligned_Nterm1'] = int(strs1[4])
            cmpRecord['numTM_unaligned_Nterm2'] = int(strs1[5])
            cmpRecord['numTM_unaligned_Cterm1'] = int(strs1[6])
            cmpRecord['numTM_unaligned_Cterm2'] = int(strs1[7])
        cmpRecord['numTM1'] = int(strs[5])
        cmpRecord['numTM2'] = int(strs[6])
        cmpRecord['seqLength1'] = int(strs[7])
        cmpRecord['seqLength2'] = int(strs[8])
        return 0
    except (IndexError, ValueError):
        print >> sys.stderr, "Bad paircmp line, \"%s\""%(line)
        return 1
#}}}

def CopyGeneralInfo_pairwise(newRecord, record):#{{{
    newRecord['id1'] = record['id1']
    newRecord['id2'] = record['id2']
    newRecord['seqidt'] = record['seqidt']
    if 'seqidt1' in record:
        newRecord['seqidt1'] = record['seqidt1']
    if 'seqidt2' in record:
        newRecord['seqidt2'] = record['seqidt2']
    if 'seqdef1' in record:
        newRecord['seqdef1'] = record['seqdef1']
    if 'seqdef2' in record:
        newRecord['seqdef2'] = record['seqdef2']
    if 'pfamid1' in record:
        newRecord['pfamid1'] = record['pfamid1']
    if 'pfamid2' in record:
        newRecord['pfamid2'] = record['pfamid2']
    if 'pfamid-inter' in record:
        newRecord['pfamid-inter'] = record['pfamid-inter']
    if 'pfamid-union' in record:
        newRecord['pfamid-union'] = record['pfamid-union']
    newRecord['cmpclass'] = record['cmpclass']
    newRecord['numTM1'] = record['numTM1']
    newRecord['numTM2'] = record['numTM2']
    newRecord['seqLength1'] = record['seqLength1']
    newRecord['seqLength2'] = record['seqLength2']
    newRecord['mapTMline'] = record['mapTMline']
    if 'mapArray' in record:
        newRecord['mapArray'] = record['mapArray']
    if 'dgscore' in record:
        newRecord['dgscore'] = record['dgscore']

#}}}

def IsUnammpedTMSatisfied(record, parameters):#{{{
    gapFracThreshold  = parameters['minGapFraction']
    DGvalueThreshold = parameters['maxDGvalue']
    if record == {}:
        return False
    elif record['numTMunmapped'] <= 0:
        return False
    else:
        for j in range(record['numTMunmapped']):
            if  (record['gapFraction'][j] >= gapFracThreshold and
                    record['DGvalue'][j] <= DGvalueThreshold): 
                return True
        return False
#}}}
def SelectAnaDIFF(ana, parameters):#{{{
# ChangeLog 2012-04-03
#   "> gapFracThreshold" changed to ">= gapFracThreshold"
#   this solved the problem that many pairs with gap fraction == 0.0 are
#   ignored.
#   "< DGvalueThreshold" changed to "<= DGvalueThreshold"
    if 'minGapFraction' in parameters:
        gapFracThreshold  = parameters['minGapFraction']
    else:
        print >> sys.stderr, "minGapFraction not set in parameters, set to 0.0"
        gapFracThreshold = 0.0
    if 'maxDGvalue' in parameters:
        DGvalueThreshold = parameters['maxDGvalue']
    else:
        print >> sys.stderr, "DGvalueThreshold not set in parameters, set to 10.0"
        DGvalueThreshold = 10.0
    selecttype = 'all'
    if 'selecttype' in parameters:
        selecttype = parameters['selecttype']

    newana={}
    newana["Nterm"] ={}
    newana["Nterm"]["numTMunmapped"] = 0
    newana["Nterm"]["index"] = []
    newana["Cterm"] ={}
    newana["Cterm"]["numTMunmapped"] = 0
    newana["Cterm"]["index"] = []
    newana["internal"] = []
    for tp in ['Nterm', 'Cterm']:
        if tp in ana and ana[tp]["numTMunmapped"] > 0:
            if IsUnammpedTMSatisfied(ana[tp], parameters):
                newana[tp] = ana[tp]
    if 'internal' in ana:
        for m in range(len(ana['internal'])):
            if  ana['internal'][m]['numTMunmapped'] > 0:
                for j in range(ana['internal'][m]["numTMunmapped"] ):
                    if IsUnammpedTMSatisfied(ana['internal'][m], parameters):
                        newana['internal'].append(ana['internal'][m])

    if (selecttype == 'all' and ( newana["Nterm"]["numTMunmapped"] > 0 or
        newana["Cterm"]["numTMunmapped"]>0 or len(newana["internal"]) > 0)):
            return newana
    elif (selecttype == 'nterm' and newana["Nterm"]["numTMunmapped"] > 0):
            return newana
    elif (selecttype == 'cterm' and newana["Cterm"]["numTMunmapped"] > 0):
            return newana
    elif (selecttype == 'internal' and len(newana["internal"]) > 0):
            return newana
    else:
        return {} 
#}}}

def GetTopoStateFraction(topoSeqList):#{{{
    """return (cnt_i, cnt_o, cntM, cnt_GAP, per_i, per_o, per_M, per_GAP)"""
    lengthAlignment=len(topoSeqList[0])
    numSeq = len(topoSeqList)
    cnt_i = [0]*lengthAlignment
    cnt_o = [0]*lengthAlignment
    cnt_M = [0]*lengthAlignment
    cnt_GAP = [0]*lengthAlignment
    for i in xrange(numSeq):
        topo = topoSeqList[i]
        for j in xrange(lengthAlignment):
            s = topo[j]
            if s == GAP:
                cnt_GAP[j] += 1
            elif s == 'i':
                cnt_i[j] += 1
            elif s == 'o':
                cnt_o[j] += 1
            elif s == 'M':
                cnt_M[j] += 1

    per_GAP = [0.0]*lengthAlignment
    per_o = [0.0]*lengthAlignment
    per_M = [0.0]*lengthAlignment
    per_i = [0.0]*lengthAlignment
    numSeq_float = float(numSeq)
    for j in xrange(lengthAlignment):
        per_o[j] = cnt_o[j]/(numSeq_float)
        per_M[j] = cnt_M[j]/(numSeq_float)
        per_i[j] = cnt_i[j]/(numSeq_float)
        per_GAP[j] = cnt_GAP[j]/(numSeq_float)
    return (cnt_i, cnt_o, cnt_M, cnt_GAP, per_i, per_o, per_M, per_GAP)

#}}}


def ReadDGScore(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        msg = "Failed to read file %s in function %s"
        print >> sys.stderr, msg%(infile, sys._getframe().f_code.co_name)
        return {}
    dgScoreDict = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line and line[0] != "#":
                strs = line.split()
                numStr = len(strs)
                if numStr >= 2:
                    try:
                        seqid = strs[0]
                        if numStr == 2:
                            dgscore = float(strs[1])
                        elif numStr == 3:
                            dgscore = float(strs[2])
                        if not seqid in dgScoreDict:
                            dgScoreDict[seqid] = []
                        dgScoreDict[seqid].append(dgscore)
                    except (ValueError, TypeError):
                        pass
        lines = hdl.readlines()
    hdl.close()
    return dgScoreDict
#}}}
def ReadRLTYInfo(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        msg = "Failed to read file %s in function %s"
        print >> sys.stderr, msg%(infile, sys._getframe().f_code.co_name)
        return {}
    rltyDict = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split()
            if len(strs) == 2:
                try:
                    seqid = strs[0]
                    rlty = float(strs[1])
                    rltyDict[seqid] = rlty
                except (ValueError, TypeError, KeyError):
                    pass
        lines = hdl.readlines()
    hdl.close()
    return rltyDict
#}}}

def ReadPairAlnTableInfo(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        msg = "Failed to read file %s in function %s"
        print >> sys.stderr, msg%(infile, sys._getframe().f_code.co_name)
        return {}
    pairalnStat = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#":
                strs = line.split()
                if len(strs) == 13:
                    try:
                        id1 = strs[0]
                        id2 = strs[1]
                        seqidt = float(strs[2])
                        alignLen = float(strs[4])
                        seqlen1 = int(strs[5])
                        seqlen2 = int(strs[6])
                        seqidt1 = float(strs[11])
                        seqidt2 = float(strs[12])
                        pairid = id1+'-'+id2
                        pairalnStat[pairid] = {}
                        tmpdict = pairalnStat[pairid] 
                        tmpdict['seqidt'] = seqidt
                        tmpdict['seqidt1'] = seqidt1
                        tmpdict['seqidt2'] = seqidt1
                        tmpdict['seqLength1'] = seqlen1
                        tmpdict['seqLength2'] = seqlen2
                    except (IndexError, ValueError, TypeError, KeyError):
                        pass
        lines = hdl.readlines()
    hdl.close()
    return pairalnStat
#}}}
