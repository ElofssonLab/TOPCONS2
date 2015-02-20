###############################################################################
################### Copyright 2008 University of Washington ###################
###############################################################################
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
"""
"""

__author__ = "Sheila M. Reynolds (sheila@ee.washington.edu)"

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

from math import *

def myLN ( x ):

    if ( x < 1.e-43 ):
	return ( -99. )
    else:
	return ( log(x) )

#------------------------------------------------------------------------------
# this function reads the output of gmtkJT (which has been piped to a file)
# and returns a matrix of topology-class probabilities (ignore class-0 which
# is not currently being used
#
# the gmtkJT output looks like this :
#
#    Segment 0, after CE, log(prob(evidence)) = -1497.470041, per frame =-2.830756, per numUFrams = -2.830756
#    --------
#    Partition 0 (P), Clique 1: Printing Clique with 1 variables, 2 entries
#    0: 5.83349709e-02 topoClass(0)=3
#    1: 9.41665029e-01 topoClass(0)=1
#    --------
#    Partition 1 (C), Clique 1: Printing Clique with 1 variables, 3 entries
#    0: 5.83349709e-02 topoClass(1)=3
#    1: 9.41665029e-01 topoClass(1)=1
#    2: 9.77297704e-16 topoClass(1)=2
#    --------
#    Partition 2 (C), Clique 1: Printing Clique with 1 variables, 3 entries
#    0: 5.83349709e-02 topoClass(2)=3
#    1: 9.41665029e-01 topoClass(2)=1
#    2: 8.13225946e-12 topoClass(2)=2
#    --------
#	...
#    --------
#    Partition 527 (E), Clique 1: Printing Clique with 1 variables, 2 entries, H=2.923754e-02
#    0: 9.97027021e-01 topoClass(528)=3
#    1: 2.97297927e-03 topoClass(528)=1
#    --------
#    Partition 527 (E), Clique 2: Printing Clique with 1 variables, 4 entries, H=6.064373e-01
#    0: 4.97263717e-03 pType(528)=2
#    1: 2.50773065e-03 pType(528)=3
#    2: 1.22812216e-01 pType(528)=0
#    3: 8.69707416e-01 pType(528)=1
# 
# Note that if the topoClass probability for a class is 0, it will not be
# written out (as in the first frame, above, where topoClass(0)=2 is absent),
# and the ordering of the topoClasses is not fixed
#
# before any of the clique probabilities are written out, a total log-probability
# of the evidence is written out which could maybe be used to give a sense of
# how well that model appears to 'fit' the data
#

def getTopologyProbs ( fhJT ):

    # the topology probability matrix will have, for each topology class
    # the probability of being in that class in this particular frame
    #		class 0 : unused
    #		class 1 : inside
    #		class 2 : membrane
    #		class 3 : outside
    #		class 4 : signal-peptide

    # NEW : 01Feb08 : now instead of this being a posterior just on topoClass,
    # it is a posterior on topoChange x topoClass :
    #		0 : impossible
    #		1 : topoChange=0 AND topoClass=1 (inside)
    #		2 : topoChange=0 AND topoClass=2 (membrane)
    #		3 : topoChange=0 AND topoClass=3 (outside)
    #		4 : topoChange=0 AND topoClass=4 (signal-peptide)
    #		5 : impossible
    #		6 : topoChange=1 AND topoClass=1 
    #		7 : topoChange=1 AND topoClass=2 
    #		8 : topoChange=1 AND topoClass=3 
    #		9 : topoChange=1 AND topoClass=4 

    topoProbs = [ [], [], [], [], [], [], [], [], [], [] ]
    topoMax = []

    done1 = 0
    while not done1:

	# read a line from the JT output file ...
	aLine = fhJT.readline()

	# this should not happen !!!
	if ( len(aLine) == 0 ):
	    print fhJT
	    print aLine
	    print " ERROR !!! empty input line from JT output file ??? "
	    sys.exit(-1)

	# print ' (0) aLine : ', aLine, len(aLine)

	if ( aLine.startswith('---') ): continue

	if ( aLine.startswith('Partition ') ):
	    # print ' aLine starts with Partition ... '

	    # take a look at the partition information ...
	    tokenList = aLine.split()
	    iFrame = int ( tokenList[1] )
	    frameLabel = tokenList[2][:-1]
	    numEntries = int ( tokenList[10] )
	    # print tokenList, iFrame, frameLabel, numEntries
	    if ( frameLabel == "(E)" ): 
		done1 = 1
		# print ' >>>>> frameLabel is (E) ... setting done1 to TRUE ... '

	    # now read all of the probabilities for this partition ...
	    curProbs = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
	    for ii in range(numEntries):
	        curPos = fhJT.tell()
		bLine = fhJT.readline()
		# print ' (1) bLine : ', bLine

		if ( 1 ):
		    tokenList = bLine.split()
		    if ( len(tokenList) != 3 ):
			print ' ERROR ??? more than 3 tokens ??? ', len(tokenList)
			print bLine, tokenList
			sys.exit(-1)
		    pVal = float ( tokenList[1] )
		    try:
		        iClass = int ( tokenList[2][-1] )
		    except:
			print bLine
			print tokenList
			print pVal
			sys.exit(-1)
		    curProbs[iClass] = pVal
		    # print bLine, iClass, pVal

	    # print ' got : ', curProbs
	    maxVal = -1.
	    maxClass = -1
	    for iClass in range(len(curProbs)):
		topoProbs[iClass] += [ curProbs[iClass] ]
		if ( maxVal < curProbs[iClass] ):
		    maxVal = curProbs[iClass]
		    maxClass = iClass
	    topoMax += [ maxClass ]

	    # now skip the '-----' line
	    aLine = fhJT.readline()

	    # if we have gotten to the "E" partition, there is one more clique
	    # that we need to read (the one that has the pType posteriors) ...

	    if ( done1 ):

		aLine = fhJT.readline()
		if not ( aLine.startswith('Partition ') ):
		    print ' ERROR ??? this line does not start with Partition ??? '
		    print aLine
		    sys.exit(-1)

	        # take a look at the partition information ...
	        tokenList = aLine.split()
	        iFrame = int ( tokenList[1] )
	        frameLabel = tokenList[2][:-1]
	        numEntries = int ( tokenList[10] )
	        # print tokenList, iFrame, frameLabel, numEntries
	        if ( frameLabel != "(E)" ): 
		    print ' ERROR ??? this is not an epilogue clique ??? '
		    print aLine
		    sys.exit(-1)

	        # now read all of the probabilities for this partition ...
	        pTypeProbs = [ 0, 0, 0, 0 ]
		for ii in range(numEntries):
	            curPos = fhJT.tell()
		    bLine = fhJT.readline()
		    # print ' (1) bLine : ', bLine

		    tokenList = bLine.split()
		    if ( len(tokenList) != 3 ):
			print ' ERROR ??? more than 3 tokens ??? ', len(tokenList)
			print bLine, tokenList
			sys.exit(-1)
		    pVal = float ( tokenList[1] )
		    try:
		        iClass = int ( tokenList[2][-1] )
		    except:
			print bLine
			print tokenList
			print pVal
			sys.exit(-1)
		    pTypeProbs[iClass] = pVal
		    # print bLine, iClass, pVal

	        # print ' got : ', pTypeProbs
		

	
    # print len(topoProbs[0])
    # print ' leaving getTopologyProbs ... '
    return ( topoProbs, topoMax )

#------------------------------------------------------------------------------
# this is the top-level function ...

def writeVEdata ( obsFileList, jtFilename, fileSuffix ):

    # print ' opening input files ... '
    # print obsFileList
    # print jtFilename

    try:
        fhObs = file ( obsFileList )
        fhJT  = file ( jtFilename )
    except:
	print '     ERROR opening input file ??? '
	sys.exit(-1)

    numProbs = 10

    jthSegment = 0
    done = 0
    while not done:

	bLine = fhObs.readline()

	if ( len(bLine) < 5 ): done = 1

	else:
	    if ( 1 ):

		obsFilename = bLine[:-1]
		postpFilename = obsFilename[:obsFilename.find(".obs")] + fileSuffix

		# print ' obsFilename   : ', obsFilename
		# print ' postpFilename : ', postpFilename

		( jtTopoProbs, jtTopoMax ) = getTopologyProbs ( fhJT )
		# print jtTopoMax

		fhTmp = file ( postpFilename, 'w' )
		for jj in range(len(jtTopoProbs[0])):
		    for kk in range(numProbs):
			lnProb = myLN ( jtTopoProbs[kk][jj] )
			fhTmp.write ( ' %12.8f ' % lnProb )
		    fhTmp.write ( '\n' )
		fhTmp.write ( '\n' )

		fhTmp.close()

		jthSegment += 1
		continue


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# test-driver is here

import sys

from  struct import *
from  math   import *

if __name__ == "__main__":

    # this program expects the following two command-line arguments:
    # the name of the 'obs' file list (eg ./params/testFiles0.list)
    # and the name of the JT output file (eg ./params/runJT.0.out)

    if ( len(sys.argv) < 3 ):
	print ' Usage : %s <obs-file-list> <jt-out-file> [file-suffix] ' % sys.argv[0]

    obsFileList = sys.argv[1]
    jtFilename  = sys.argv[2]
    if ( len(sys.argv) == 4 ):
	fileSuffix = sys.argv[3]

    if ( fileSuffix[0] != '.' ):
	fileSuffix = '.' + fileSuffix

    writeVEdata ( obsFileList, jtFilename, fileSuffix )


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
