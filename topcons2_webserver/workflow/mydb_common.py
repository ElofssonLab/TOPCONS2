#!/usr/bin/env python

# Description:
#   A collection of classes and functions used by the class MyDB
#
# Author: Nanjiang Shu (nanjiang.shu@scilifelab.se)
#
# Address: Science for Life Laboratory Stockholm, Box 1031, 17121 Solna, Sweden

import sys
import os
from array import array
import mybase

FORMAT_BINARY = 0
FORMAT_TEXT = 1
TYPE_DICT = 0
TYPE_LIST = 1
LargeFileThresholdSize = 1.5*1024*1024*1024
version = "1.4"

def GetIndexFileHeaderText(headerinfo):#{{{
    """
    Get the header information of the index file in ASCII format
    """
    (dbname, version, ext, prefix) = headerinfo
    indexFileHeaderText = []
    indexFileHeaderText.append("DEF_VERSION %s"%(version))
    indexFileHeaderText.append("DEF_DBNAME %s"%(dbname))
    indexFileHeaderText.append("DEF_EXTENSION %s"%(ext))
    indexFileHeaderText.append("DEF_PREFIX %s"%(prefix))
    return indexFileHeaderText
#}}}
def GetIndexFile(dbname, formatindex):#{{{
    """
    Get the name of the index file
    """
# return (indexfile, formatindex)
    indexfile = ""
    if formatindex == FORMAT_BINARY:
        indexfile = dbname + ".indexbin"
        if not os.path.exists(indexfile):
            formatindex = FORMAT_TEXT
            msg = "Binary index file {} does not exist. "\
                   "Try looking for text index file"
            print >> sys.stderr, msg.format(indexfile)
            indexfile = dbname + ".index"
            if not os.path.exists(indexfile):
                msg = "Text index file {} does not exist"
                print >> sys.stderr, msg.format(indexfile)
                indexfile = ""
    else:
        indexfile = dbname+".index"
        if not os.path.exists(indexfile):
            formatindex = FORMAT_BINARY
            msg = "Text index file {} does not exist."\
                    "Try looking for binary index file"
            print >> sys.stderr, msg.format(indexfile)
            indexfile=dbname+".indexbin"
            if not os.path.exists(indexfile):
                msg = "Binary index file {} does not exist"
                print >> sys.stderr, msg.format(indexfile)
                indexfile = ""
    return (indexfile, formatindex)
#}}}
def WriteIndexHeader(indexFileHeaderText, formatindex, fpindex):#{{{
    """
    Write the header information to the index file of the database
    """
    if formatindex == FORMAT_TEXT:
        for s in indexFileHeaderText:
            print >> fpindex, s
    else:
        dumpedtext='\n'.join(s for s in indexFileHeaderText)
        vI = array('I')
        vI.append(len(dumpedtext))
        vI.tofile(fpindex)
        fpindex.write(dumpedtext)
#}}}
def WriteIndexContent(indexList, formatindex, fpindex):#{{{
    """
    Write the content to the index file of the database
    """
    if formatindex == FORMAT_TEXT:
        numRecord = len(indexList[0])
        idList = indexList[0]
        v1 = indexList[1]
        v2 = indexList[2]
        v3 = indexList[3]
        for i in xrange(numRecord):
            print >> fpindex, idList[i], v1[i], v2[i], v3[i]
    else:
        maxOffset = max(indexList[2])

        numRecord = len(indexList[0])

        idList = indexList[0]
        v1 = indexList[1]
        v3 = indexList[3]
        if maxOffset > LargeFileThresholdSize:
            v2 = indexList[2]
        else: #'I'
            v2 = array('I', [x for x in indexList[2]])

        dumpedliststr = '\n'.join(s for s in idList)

        vI=array('I')
        vI.append(len(dumpedliststr))
        vI.tofile(fpindex)
        fpindex.write(dumpedliststr)

        vI=array('I')
        vI.append(numRecord)
        vI.tofile(fpindex)

        v1.tofile(fpindex)
        v2.tofile(fpindex)
        v3.tofile(fpindex)
#}}}
def ReadIndex_binary(indexfile, isPrintWarning = False):#{{{
    """
    Read the index file of the binary format
    """
# return (indexList, headerinfo, dbfileindexList)
    indexList = []
    indexFileHeaderText = []
    size_indexfile = os.path.getsize(indexfile)
    cntReadByte = 0
    try:
        fpin=open(indexfile, "rb")
        vI = array('I')
        vI.fromfile(fpin,1)
        cntReadByte += vI.itemsize
        dumpedtext = fpin.read(vI[0])
        cntReadByte += vI[0]

        strs = dumpedtext.split("\n")
        origdbname = ""
        origversion = ""
        origext = ""
        origprefix = ""
        for line in strs:
            if not line or line[0] == "#":
                continue
            ss=line.split()
            if ss[0] == "DEF_DBNAME":
                if len(ss)>=2:
                    origdbname=ss[1]
            elif ss[0] == "DEF_VERSION":
                if len(ss)>=2:
                    origversion=ss[1]
            elif ss[0] == "DEF_EXTENSION":
                if len(ss)>=2:
                    origext=ss[1]
            elif ss[0] == "DEF_PREFIX":
                if len(ss)>=2:
                    origprefix=ss[1]
        if isPrintWarning:
            if origversion == "": 
                msg = "{}: Warning! No version info in the index file {}"
                print >> sys.stderr, msg.format(sys.argv[0],indexfile)
            elif origversion != version:
                msg = "{}: Warning! Version conflicts. "\
                        "Version of the index file {} ({}) "\
                        "!= version of the program ({})"
                print >> sys.stderr, msg.format(sys.argv[0], indexfile,
                        origversion, version)

        headerinfo = (origdbname, origversion, origext, origprefix)
        #read in other information
        vI = array('I')
        vI.fromfile(fpin,1)
        cntReadByte += vI.itemsize

        dumpedidlist=fpin.read(vI[0])
        cntReadByte += vI[0]

        idlist = dumpedidlist.split("\n")
        vI=array('I')
        vI.fromfile(fpin,1)
        cntReadByte += vI.itemsize

        numRecord = vI[0]
        if numRecord != len(idlist):
            msg = "{}: numID ({}) != numRecord ({}) for indexfile {} "
            print >> sys.stderr, msg.format(sys.argv[0], len(idlist),
                    numRecord, indexfile)

        sizeRecord_I = (array('B').itemsize + array('I').itemsize +
                array('I').itemsize)
        sizeRecord_L = (array('B').itemsize + array('L').itemsize +
                array('I').itemsize)
        sizeRecord = int(mybase.FloatDivision(size_indexfile - cntReadByte, numRecord))
        if abs(sizeRecord - sizeRecord_I)  < abs(sizeRecord - sizeRecord_L):
            vIarray=[array('B'), array('I'), array('I')]
        else:
            vIarray=[array('B'), array('L'), array('I')]
        for i in range(3):
            vIarray[i].fromfile(fpin,numRecord)

        lastDBFileIndex = vIarray[0][numRecord-1]
        dbfileindexList = range(lastDBFileIndex+1)

        indexList.append(idlist)
        for i in range(3):
            indexList.append(vIarray[i])
        fpin.close()
        return (indexList, headerinfo, dbfileindexList)
    except IOError:
        msg = "Failed to read index file {} in function {}"
        print >> sys.stderr, msg.format(indexfile, sys._getframe().f_code.co_name)
        return (None, None, None)
#}}}
def ReadIndex_text(indexfile, isPrintWarning = False):#{{{
    """
    Read the index file of the ASCII format
    """
# return (indexList, headerinfo, dbfileindexList)
    indexList = []
    idList = []
    v1 = array('B') # dbfile index
    v2 = array('L') # offset
    v3 = array('I') # block size
    apd1 = idList.append
    apd2 = v1.append
    apd3 = v2.append
    apd4 = v3.append
    indexFileHeaderText = []
    origdbname=""
    origversion=""
    origext=""
    origprefix=""
    try:

        hdl = mybase.ReadLineByBlock(indexfile)
        lines = hdl.readlines()
        while lines != None:
            for line in lines:
                if not line or line[0] == "#":
                    continue
                strs = line.split()
                if strs[0] == "DEF_DBNAME":
                    if len(strs)>=2:
                        origdbname=strs[1]
                elif strs[0] == "DEF_VERSION":
                    if len(strs)>=2:
                        origversion=strs[1]
                elif strs[0] == "DEF_EXTENSION":
                    if len(strs)>=2:
                        origext=strs[1]
                elif strs[0] == "DEF_PREFIX":
                    if len(strs)>=2:
                        origprefix = strs[1]
                else:
                    apd1(strs[0])
                    apd2(int(strs[1]))
                    apd3(int(strs[2]))
                    apd4(int(strs[3]))
            lines = hdl.readlines()

        indexList.append(idList)
        indexList.append(v1)
        indexList.append(v2)
        indexList.append(v3)

        headerinfo = (origdbname, origversion, origext, origprefix)

        numRecord = len(idList)
        lastDBFileIndex = v1[numRecord-1]
        dbfileindexList = range(lastDBFileIndex+1)

        if isPrintWarning:
            if origversion == "":
                msg = "{}: Warning! No version info in the index file {}"
                print >> sys.stderr,msg.format(sys.argv[0],indexfile)
            elif origversion != version:
                msg = "{}: Warning! Version conflicts. "\
                        "Version of the index file {} ({}) "\
                        "!= version of the program ({})"
                print >> sys.stderr, msg.format(sys.argv[0],indexfile,
                        origversion, version)
        return (indexList, headerinfo, dbfileindexList)
    except IOError:
        msg = "Failed to read index file {} in function {}"
        print >> sys.stderr, msg.format(indexfile, sys._getframe().f_code.co_name)
        return (None, None, None)
#}}}
