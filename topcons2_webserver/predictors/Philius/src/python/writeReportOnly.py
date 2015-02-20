###############################################################################
################### Copyright 2008 University of Washington ###################
###############################################################################
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

"""
This program reads a binary gmtkViterbi output file and compares it to the
labeled observation file ...

label			allowed states 					topoClass
0  "."			anything goes					0
1  "i"			1,15,16,17		'inside'		1
2  "M"			2,3,4,12,13,14		'membrane'		2
3  "o"			5,6,7,8			'outside' (short)	3
4  "O"	 		9,10,11,25,27,28	'outside' (long)	3
5  "n"			18,26			SP n-region		4
6  "h"			19			SP h-region		4
7  "c"			20,21,22,23		SP c-region		4
8  "C"			24			SP CUT site		4

"""

__author__ = "Sheila M. Reynolds (sheila@ee.washington.edu)"

#------------------------------------------------------------------------------

numTypes = 4

pString = [0] * numTypes
pString[0] = "GLOBULAR"
pString[1] = "GLOBULAR with SIGNAL PEPTIDE"
pString[2] = "TM PROTEIN"
pString[3] = "TM PROTEIN with SIGNAL PEPTIDE"

topoString = [0] * 6
topoString[0] = "BAD TOPOLOGY INDEX !!!"
topoString[1] = "inside         "
topoString[2] = "TMhelix        "
topoString[3] = "outside        "
topoString[4] = "BAD TOPOLOGY INDEX ???"
topoString[5] = "signal-peptide "

topoString[1] = "I"
topoString[2] = "M"
topoString[3] = "O"

#------------------------------------------------------------------------------
# here we just want to return a string that indicates the predicted protein
# type :	GLOBULAR				type0
#		GLOBULAR with SIGNAL PEPTIDE		type1
#		TM PROTEIN				type2			
#		TM PROTEIN with SIGNAL PEPTIDE		type3

def getPredictedType ( labelVec ):

    # first check to see if there is a signal-peptide ...
    if ( 5 in labelVec ):
	spBit = 1
    else:
	spBit = 0

    # then check to see if there is a TM segment ...
    if ( 2 in labelVec ):
	tmBit = 1
    else:
	tmBit = 0

    return ( spBit + 2*tmBit )

#------------------------------------------------------------------------------

def getSPlength ( labelVec ):

    if ( 5 not in labelVec ): return ( 0 )

    nn = len(labelVec)

    numAAs = 0
    for ii in range(nn):
	if ( labelVec[ii] == 5 ): numAAs += 1

    return ( numAAs )

#------------------------------------------------------------------------------

def getNumberTMHelices ( labelVec ):

    # print labelVec

    if ( 2 not in labelVec ): return ( 0, 0 )

    nn = len(labelVec)
    tmpVec = labelVec + [99]

    numTMHs = 0
    numAAs = 0
    ii = 0
    while ( ii < nn ):
	while ( tmpVec[ii] != 2  and  ii < nn ): ii += 1
	if ( tmpVec[ii] == 2 ):
	    numTMHs += 1
	    while ( tmpVec[ii] == 2 ): 
		ii += 1
		numAAs += 1

    return ( numTMHs, numAAs )

#------------------------------------------------------------------------------

def topoClass2label ( iTopoClass ):

    if ( iTopoClass == 1 ): return 1
    if ( iTopoClass == 2 ): return 2
    if ( iTopoClass == 3 ): return 3
    if ( iTopoClass == 4 ): return 5

    print ' ERROR in topoClass2label ??? ', iTopoClass
    sys.exit(-1)

#------------------------------------------------------------------------------

def label2topoClass ( iLabel ):

    if ( iLabel == 1 ): return 1
    if ( iLabel == 2 ): return 2
    if ( iLabel == 3 ): return 3
    if ( iLabel == 5 ): return 4

    print ' ERROR in label2topoClass ??? ', iLabel
    sys.exit(-1)

#------------------------------------------------------------------------------

def getViterbiSeqs ( vitFilename ):

    # print ' '
    # print ' '
    # print ' in getViterbiSeqs ... ', vitFilename

    # number of variables dumped out by the Viterbi pass ...
    numDump = 4

    try:
	fhVit = file ( vitFilename )
    except:
	print ' ERROR in getViterbiSeqs : error opening file ', vitFilename
	sys.exit(-1)

    wholeFile = fhVit.read()
    fhVit.close()

    numSamplesPerRec = numDump
    numBytesPerInt = 4
    numBytesPerRec = numSamplesPerRec * numBytesPerInt

    numBytes = len(wholeFile)
    numSamples = numBytes / numBytesPerRec

    vitStates = [0] * numSamples
    vitLabels = [0] * numSamples
    iByte = 0

    for ii in range(numSamples):

	# the hidden variables currently being written out by gmtkViterbi
	# are : iState, iToppoClass, spBit, and tmBit
	( iState, iTopoClass, spBit, tmBit ) = \
	    unpack ( '<IIII', wholeFile[iByte:iByte+numBytesPerRec] )
	iByte += numBytesPerRec
	vitStates[ii] = iState
	vitLabels[ii] = topoClass2label ( iTopoClass )


    return ( vitStates, vitLabels )

#------------------------------------------------------------------------------

def getTopology ( inLabels ):

    # start with the first label ...
    topoList = [ inLabels[0] ]

    # then walk along looking for label changes ...
    for ii in range(1,len(inLabels)):
	# if the label is 0, skip it ...
	if ( inLabels[ii] == 0 ): continue
	# otherwise, if the label has changed, add it to the topology vector
	if ( inLabels[ii] != topoList[-1] ):
	    topoList += [ inLabels[ii] ]

    # minor hack to the topology : labels '3' and '4' are both 'outside'
    # loops so we'll just replace all 4's with 3's ...
    numChanges = 0
    copyList = [0] * len(topoList)
    for ii in range(len(topoList)):
	copyList[ii] = topoList[ii]

    for ii in range(len(topoList)):
	if ( topoList[ii] == 4 ): 
	    topoList[ii] = 3
	    numChanges += 1

    # if ( numChanges > 0 ):
    #	print '         input topology  : ', copyList
    #	print '         output topology : ', topoList

    # another minor hack : if the topology starts with a SP, then it will
    # pass through [5,6,7,8] ... we'll simplify this down to just [5] ...

    if ( topoList[0] == 5 ):
	if ( topoList[1] == 6 ):
	    newList = [5] + topoList[4:]
	    return ( newList )

    return ( topoList )

#------------------------------------------------------------------------------

def getProteinName ( aFilename ):

    nn = len(aFilename) - 1

    i1 = nn
    while ( aFilename[i1] != '/'  and  i1>0 ): i1 -= 1
    i1 += 1

    i2 = i1 + 1
    while ( aFilename[i2] != '.'  and  i2<nn ): i2 += 1

    return ( aFilename[i1:i2] )

#------------------------------------------------------------------------------

def writeDetailedTopology ( fhOut, vitStates, vitLabels, segScores ):

    nn = len(vitLabels)

    tmpLabels = vitLabels + [99]
    # print tmpLabels

    lastSegLabel = -1

    i1 = 0
    iSeg = 0
    while ( i1 < nn ):
	i2 = i1+1
	while ( tmpLabels[i2] == tmpLabels[i1]  and  i2 < nn ): i2 += 1
	# print i1, i2, nn
	if ( i1 < nn ):

	    # if the previous segment was a signal peptide, we need to decrement
	    # the start position by one ...
	    if ( lastSegLabel == 5 ):
	        fhOut.write ( '%s:%d-%d-%4.2f$' % \
		 ( topoString[tmpLabels[i1]], i1, i2, segScores[iSeg][1] ) )
	    # if the current segment is a signal peptide, we need to decrement 
	    # the stop position by one ...
	    elif ( tmpLabels[i1] == 5 ):
	        fhOut.write ( '%s:%d-%d-%4.2f$' % \
		 ( topoString[tmpLabels[i1]], i1, i2-1, segScores[iSeg][1] ) )
	    # otherwise everything is ok
	    else:
	        fhOut.write ( '%s:%d-%d-%4.2f$' % \
		 ( topoString[tmpLabels[i1]], i1+1, i2, segScores[iSeg][1] ) )
	    iSeg += 1

	    # if this portion represents a signal-peptide, then we should also
	    # write out a more detailed description :
	    if ( tmpLabels[i1] == 5 ):
		if ( vitStates[i1] == 18  or  vitStates[i1] == 26  or  vitStates[i1] == 29 ):
		    n1 = i1
		    n2 = n1 + 1
		    while ( vitStates[n2] == 18  or  vitStates[n2] == 26  or  vitStates[n2] == 29 ): n2 += 1
		    h1 = n2
		    h2 = h1 + 1
		    while ( vitStates[h2] == 19 ): h2 += 1
		    c1 = h2
		    c2 = c1 + 1
		    while ( vitStates[c2] < 24 ): c2 += 1
		    if ( vitStates[c2] != 24 ):
			print ' no CUT site ??? ' 
			sys.exit(-1)
		else:
		    print ' UNEXPECTED VALUE ??? ', vitStates[i1]
		    sys.exit(-1)
		fhOut.write ( '                n-region    %6d  %6d\n' % ( n1+1, n2 ) )
		fhOut.write ( '                h-region    %6d  %6d\n' % ( h1+1, h2 ) )
		fhOut.write ( '                c-region    %6d  %6d\n' % ( c1+1, c2 ) )
		fhOut.write ( '           cleavage between %d and %d \n\n\nANNOTATION:' % ( c2, (c2+1) ) )   #ktsirig

	    lastSegLabel = tmpLabels[i1]

	i1 = i2

#------------------------------------------------------------------------------

def getPass2Probs ( fhP2P ):

    aLine = fhP2P.readline()
    tokenList = aLine.split()
    # print len(tokenList), tokenList

    totalProb = float ( tokenList[ 7][:-1] )
    perFrameP = float ( tokenList[14] )

    return ( totalProb, perFrameP )

#------------------------------------------------------------------------------

def getSegmentScores ( vitLabels, fhJT ):

    # print ' in getSegmentScores : ', len(vitLabels), fhJT.tell()
    numLabels = len(vitLabels)
    # print vitLabels

    tmpLabels = vitLabels + [99]
    tmpProbs = [-1] * numLabels

    # we should be at the "Segment" line ...
    aLine = fhJT.readline()
    # print aLine

    # print ' looping over iLabel = 0 to %d ' % ( numLabels-1 )
    for iLabel in range(numLabels):

	curClass = label2topoClass ( vitLabels[iLabel] )
	# print ' curClass : ', curClass

	# skip the '----' line ...
	aLine = fhJT.readline()

	# next read the first Partition ...
	aLine = fhJT.readline()
	# print aLine
	tokenList = aLine.split()
	# print len(tokenList), tokenList
	numEntries = int ( tokenList[10] )
	# print ' reading %d entries ... ' % numEntries

	for ii in range(numEntries):
	    aLine = fhJT.readline()
	    tokenList = aLine.split()
	    prob = float ( tokenList[1] )
	    iTenClass = int ( tokenList[2][-1] )
	    iclass = mapDownFromTen ( iTenClass )
	    # print tokenList, prob, iTenClass, iclass
	    if ( iclass == curClass ):
		if ( abs(tmpProbs[iLabel]+1.) < 0.001 ):
		    tmpProbs[iLabel] = prob
		else:
		    tmpProbs[iLabel] += prob

	if ( tmpProbs[iLabel] < 0 ):
	    print ' ERROR ??? did not get a probability for this label ??? '
	    print iLabel, tmpProbs[iLabel-10:iLabel+1], tmpProbs
	    sys.exit(-1)

    # now we want to figure out the individual segment scores ...

    # print ' now have filled up the tmpProbs vector  ... '
    # print tmpProbs

    lastLabel = tmpLabels[0]
    sumLog = log ( tmpProbs[0] )
    sumLin = tmpProbs[0]
    sumN = 1
    segScores = []

    for iLabel in range(1,numLabels+1):

	if ( tmpLabels[iLabel] != lastLabel ):

	    logScore = sumLog / float(sumN)
	    logScore = exp ( logScore )

	    linScore = sumLin / float(sumN)

	    # to be a touch conservative, we will clip the scores at 0.99 since
	    # they will be written out with 2 digits of precision and we don't
	    # want it looking like "1.00" ...
	    if ( logScore > 0.99 ): logScore = 0.99
	    if ( linScore > 0.99 ): linScore = 0.99

	    # print ' computing new segment score : ', sumLog, sumN, logScore, linScore
	    segScores += [ ( logScore, linScore ) ]

	    if ( iLabel < numLabels ):
		lastLabel = tmpLabels[iLabel]
	        sumLog = log ( tmpProbs[iLabel] )
		sumLin = tmpProbs[iLabel]
	        sumN = 1

	else:

	    sumLog += log ( tmpProbs[iLabel] )
	    sumLin += tmpProbs[iLabel]
	    sumN += 1
	    

    # print ' segScores : ', segScores
    # print ' range     : ', min(segScores), max(segScores)

    # print fhJT.tell()
	    
    return ( segScores )

#------------------------------------------------------------------------------

def mapDownFromTen ( iTenClass ):

    if ( iTenClass == 0 ): 
	print ' ERROR ??? how can this be ??? '
	sys.exit(-1)

    if ( iTenClass < 5 ): return ( iTenClass )

    if ( iTenClass == 5 ): 
	print ' ERROR ??? how can this be ??? '
	sys.exit(-1)

    return ( iTenClass-5 )

#------------------------------------------------------------------------------

def getJTposteriors ( fhJT ):

    # print ' in getJTposteriors ... ', fhJT.tell()
    # print fhJT

    expNumTMaa = 0.
    expNumSPaa = 0.

    numTopoClasses = 5
    topoProbs0 = [0.] * numTopoClasses
    topoProbsN = [0.] * numTopoClasses

    # we should be at the "Segment" line ...
    aLine = fhJT.readline()
    tokenList = aLine.split()
    # print len(tokenList), tokenList
    totalProb = float ( tokenList[6][:-1] )
    perFrameP = float ( tokenList[13] )
    # print totalProb, perFrameP

    # skip the '----' line ...
    aLine = fhJT.readline()

    # next read the first Partition ...
    aLine = fhJT.readline()
    # print aLine
    tokenList = aLine.split()
    # print len(tokenList), tokenList
    numEntries0 = int ( tokenList[10] )
    cliqueH0 = float ( tokenList[12][2:] )
    # print ' reading %d entries ... ', numEntries0
    for ii in range(numEntries0):
	aLine = fhJT.readline()
	tokenList = aLine.split()
	prob = float ( tokenList[1] )
	iTenClass = int ( tokenList[2][-1] )
	iclass = mapDownFromTen ( iTenClass )
	# print tokenList, prob, iTenClass, iclass
	topoProbs0[iclass] = prob

	if ( iclass == 2 ): expNumTMaa += prob
	if ( iclass == 4 ): expNumSPaa += prob

    # print ' got the first partition ... '

    # and then we need to skip forward until we get to the last partition ...
    # but in the meanwhile we want to look at each partition and find the one
    # with the highest probability of having topoClass=2

    # 18-Jan-08 : also want to count up the total probability of being in 
    # the membrane ... ( topoClass=2 )

    maxTM_p = -1.
    maxTM_h = -1.
    
    done = 0
    while not done:
	aLine = fhJT.readline()
	# print ' A : ', aLine
	if ( aLine.startswith("---") ):
	    aLine = fhJT.readline()
	    # print ' B : ', aLine
	if ( aLine.startswith("Partition ") ):
	    if ( aLine.find(" (E), Clique") > 0 ): 
		done = 1
	    else:
	        tokenList = aLine.split()
	        # print ' C : ', tokenList
	        tmpNum = int ( tokenList[10] )
	        tmpH   = float ( tokenList[12][2:] )
	        for ii in range(tmpNum):
		    aLine = fhJT.readline()
		    tokenList = aLine.split()
		    prob = float ( tokenList[1] )
		    iclass = int ( tokenList[2][-1] )
		    if ( iclass == 2 ):
		        if ( prob > maxTM_p ):
			    maxTM_p = prob
			    maxTM_h = tmpH
			expNumTMaa += prob
		    if ( iclass == 4 ): expNumSPaa += prob

    # print ' maxTM : ', maxTM_p, maxTM_h

    # and then we read the last Partition ...
    # print aLine
    tokenList = aLine.split()
    # print len(tokenList), tokenList
    numEntriesN = int ( tokenList[10] )
    cliqueHN = float ( tokenList[12][2:] )
    for ii in range(numEntriesN):
	aLine = fhJT.readline()
	tokenList = aLine.split()
        # print len(tokenList), tokenList
	prob = float ( tokenList[1] )
	iTenClass = int ( tokenList[2][-1] )
	iclass = mapDownFromTen ( iTenClass )
	topoProbsN[iclass] = prob
	if ( iclass == 2 ): expNumTMaa += prob
	if ( iclass == 4 ): expNumSPaa += prob

    # and we also need to take a look at the information in the new
    # extra clique in the last Partition ...
    aLine = fhJT.readline()
    # print aLine
    aLine = fhJT.readline()
    # print aLine
    tokenList = aLine.split()
    numEntriesN = int ( tokenList[10] )
    typeHN = float ( tokenList[12][2:] )
    typeProbs = [0] * 4
    for ii in range(numEntriesN):
	aLine = fhJT.readline()
	# print aLine
	tokenList = aLine.split()
	prob = float ( tokenList[1] )
	itype = int ( tokenList[2][-1] )
	typeProbs[itype] = prob
	if ( iclass == 2 ): expNumTMaa += prob
	if ( iclass == 4 ): expNumSPaa += prob


    # print topoProbs0
    # print topoProbsN
    # print ' typeHN    : ', typeHN
    # print ' typeProbs : ', typeProbs

    # sys.exit(-1)

    return ( totalProb,  perFrameP, \
	     cliqueH0,   topoProbs0, \
	     cliqueHN,   topoProbsN, \
	     maxTM_p,    maxTM_h, \
	     typeHN,     typeProbs, \
	     expNumSPaa, expNumTMaa )


#------------------------------------------------------------------------------

def writeOneReport ( fhOut, vitFilename, fhJT, fhP2P ):

    # print ' '
    # print ' '
    # print ' >>>>>>>>>>>>>>>>>>>>> '
    # print ' in writeOneReport ... '
    # print vitFilename

    pName = getProteinName ( vitFilename )
    # print pName
    fhOut.write ( '\n#\n' )
    fhOut.write ( '# Name                       %s\n' % pName )

    ( vitStates, vitLabels ) = getViterbiSeqs ( vitFilename )
    # print vitStates
    # print vitLabels

    if ( 1 ):

	jtPos0 = fhJT.tell()

        ( jtTotalProb, jtPerFrameP, \
          cliqueH0,    topoProbs0, \
          cliqueHN,    topoProbsN, \
	  maxTM_p,     maxTM_h, \
	  typeHN,      typeProbs, \
	  expNumSPaa,  expNumTMaa ) = getJTposteriors ( fhJT )

	jtPos1 = fhJT.tell()
	# print ' fhJT file positions : ', jtPos0, jtPos1

        # print jtTotalProb, jtPerFrameP, cliqueH0, cliqueHN, maxTM_p, maxTM_h, typeHN, typeProbs
	# print topoProbs0
	# print topoProbsN
        # sys.exit(-1)

	( vitTotalProb, vitPerFrameP ) = getPass2Probs ( fhP2P )
	# print vitTotalProb, vitPerFrameP
	# sys.exit(-1)

	# we need to get something from fhJT regarding this protein, so we need
	# to rewind ...
	fhJT.seek(jtPos0)
	segScores = getSegmentScores ( vitLabels, fhJT )
	# and then go back to where getJTposteriors had left us before ...
	fhJT.seek(jtPos1)

    pType = getPredictedType ( vitLabels )
    fhOut.write ( '# Annotation                 %s\n' % pString[pType] )

    fhOut.write ( '# Length                     %6d\n' % len(vitLabels) )

    numAAs = getSPlength ( vitLabels )
    fhOut.write ( '# Number of AAs in SP           %3d   %5.1f\n' % ( numAAs, expNumSPaa ) )

    ( numTMHs, numAAs ) = getNumberTMHelices ( vitLabels )
    fhOut.write ( '# Number of predicted TMHs      %3d\n' % numTMHs )
    if ( numTMHs > 0 ):
        fhOut.write ( '# Number of AAs in TMHs         %3d   %5.1f   %4.1f\n' % \
		 ( numAAs, expNumTMaa, (float(numAAs)/float(numTMHs)) ) )
    else:
        fhOut.write ( '# Number of AAs in TMHs         %3d   %5.1f\n' % \
		 ( numAAs, expNumTMaa ) )


    # first we write out the overall protein-type confidence score
    if ( pType < 2 ):
	fhOut.write ( '# Score                           %4.2f\n' % min(0.99,typeProbs[pType]) )
	
    else:
	fhOut.write ( '# Scores                          %4.2f' % min(0.99,typeProbs[pType]) )
	
	# a second topology confidence score is written out ONLY for TM proteins
	# ( including SP+TM proteins )
	if ( pType == 2 ):
	    # for TM proteins, we take the MIN() over the first and last IN/OUT
	    # segments and over the individual TM helices
	    newScore = min ( segScores[0][1], segScores[-1][1] )
	    for iSeg in range(1,len(segScores),2):
		if ( newScore > segScores[iSeg][1] ): newScore = segScores[iSeg][1]

	elif ( pType == 3 ):
	    # for SP+TM proteins, we do the same thing, but which segments are the
	    # TM helices is shifted over by 1 ...
	    newScore = min ( segScores[0][1], segScores[-1][1] )
	    for iSeg in range(2,len(segScores),2):
		if ( newScore > segScores[iSeg][1] ): newScore = segScores[iSeg][1]

	fhOut.write ( '  %4.2f\n' % newScore )
	
    fhOut.write ( '#\nANNOTATION:' )	#ktsirig

    if ( pType > 0 ):
	writeDetailedTopology ( fhOut, vitStates, vitLabels, segScores )

    fhOut.write ( '\n\n\n\n' )	#ktsirig
    fhOut.write ( '----------------------------------------------------------------------\n' )

    # sys.exit(-1)

#------------------------------------------------------------------------------

def writeAllReports ( fhOut, vitFileList, jtFileName, p2pFileName ):

    # print ' in writeAllReports : '
    # print '     vitFileList : ', vitFileList
    # print '     jtFileName  : ', jtFileName
    # print '     p2pFileName : ', p2pFileName

    fhVit = file ( vitFileList )
    fhJT  = file ( jtFileName )
    fhP2P = file ( p2pFileName )

    for aLine in fhVit:
	vitFilename = aLine[:-1]
	writeOneReport ( fhOut, vitFilename, fhJT, fhP2P )

#------------------------------------------------------------------------------

def getDirName ( faaFilename ):

    if ( faaFilename.find('/') < 0 ): return ( '' )

    ii = faaFilename.find('/')
    jj = ii
    while ( jj >= 0 ):
        jj = faaFilename.find('/',jj+1)
        if ( jj >= 0 ): ii = jj

    dirName = faaFilename[:ii]

    return ( dirName )

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# test-driver is here

import sys

from   struct  import  *
from   math    import  *

if __name__ == "__main__":

    if ( len(sys.argv)!=4 ):
	print ' Usage : %s <vitFileList> <jtFileName> <p2pFileName> ' % sys.argv[0]
	sys.exit(-1)

    if ( len(sys.argv) == 4 ):

	vitFileList = sys.argv[1]
	jtFileName  = sys.argv[2]
	p2pFileName = sys.argv[3]

	i1 = vitFileList.find('vit') + 3
	i2 = vitFileList.find('.list') + 1
	if ( i2 <= i1 ):
	    print ' ERROR ??? unexpected file name ??? '
	    print vitFileList
	    print i1, i2
	    sys.exit(-1)

	dirName = getDirName ( vitFileList )
	reportName = dirName + '/Philius.report.out'
	# print ' output report name : ', reportName

        fhOut = file ( reportName, 'w' )

        writeAllReports ( fhOut, vitFileList, jtFileName, p2pFileName )

        
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
