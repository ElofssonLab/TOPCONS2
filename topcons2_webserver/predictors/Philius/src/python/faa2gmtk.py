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

def stripBlanks ( inString ):

    outString = ''
    for ii in range(len(inString)):
	if ( inString[ii] != ' ' ): outString += inString[ii]

    return ( outString )

#------------------------------------------------------------------------------

def appendSegment ( inSeq, newLine ):

    if ( newLine[-1] == '\n' ):
	outSeq = inSeq + stripBlanks ( newLine[:-1] )
    else:
	outSeq = inSeq + stripBlanks ( newLine )

    return ( outSeq )

#------------------------------------------------------------------------------

def readNextSequenceFromFASTA ( fh ):

    hdrLine = fh.readline()
    if ( len(hdrLine) == 0 ): return ( '', '', '', '' )

    if ( hdrLine[0] != '>' ):
	print ' FATAL ERROR in readNextSequenceFromFASTA : file not properly positioned '
	print fh
	print hdrLine
	sys.exit(-1)

    hdrLine = hdrLine[1:]
    hdrTokens = hdrLine.split()

    done = 0
    eof = 0
    rewind = 0

    newSeq = ''
    lblSeq = ''

    # when we start reading, we assume that we are reading sequence ...
    inSeq = 1
    inLbl = 0

    while not done:
	curPos = fh.tell()
	aLine = fh.readline()
	if ( len(aLine) == 0 ): 
	    eof = 1
	    done = 1
	elif ( aLine[0] == '>' ):
	    rewind = 1
	    done = 1
	elif ( aLine[0] == '#' ):
	    # there are some labels for this sequence ...
	    inSeq = 0
	    inLbl = 1
	    lblSeq = appendSegment ( lblSeq, aLine[1:] )
	else:
	    if ( inSeq ): newSeq = appendSegment ( newSeq, aLine )
	    if ( inLbl ): lblSeq = appendSegment ( lblSeq, aLine )

    if ( rewind ): fh.seek(curPos)

    #### print ' returning ', hdrTokens[0], hdrTokens[1]
    #### print newSeq
    #### print lblSeq

    if ( len(hdrTokens) == 1 ):
	return ( hdrTokens[0], '', newSeq, lblSeq )
    else:
        return ( hdrTokens[0], hdrTokens[1], newSeq, lblSeq )

#------------------------------------------------------------------------------

# the 20 canonical amino acids ...
CanonicalAminos = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
		    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]

# we will actually use a list of length 24, where:
#	the unknown "X" occupies the "special" position 0
# 	"U"(21) is the 21st amino acid, selenocysteine
#	the ambiguous "B"(22) can represent either 'D'(3) or 'N'(12)
# and 	the ambiguous "Z"(23) can represent either 'E'(4) or 'Q'(14)

aaList = [ 'X' ] + CanonicalAminos + [ 'U', 'B', 'Z' ]

def aaId ( aaInitial ):
    try:
	return ( aaList.index ( aaInitial ) )
    except:
	return ( 0 )

#------------------------------------------------------------------------------

labelList = [ '.', 'i', 'M', 'o', 'O', 'n', 'h', 'c', 'C' ]

def labelId ( aLabel ):
    try:
	return ( labelList.index ( aLabel ) )
    except:
	return ( 0 )

#------------------------------------------------------------------------------

def makeGMTKobsFilename ( hdr1 ):

    filename = hdr1 + '.obs'

    return ( filename )

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

def writeGMTKobsFile ( aSeq, lSeq, gmtkFilename ):

    try:
	fh = file ( gmtkFilename, 'w' )
    except:
	print ' ERROR in writeGMTKobsFile : failed to open output file ', gmtkFilename
	sys.exit(-1)

    seqLength = len(aSeq)

    if ( len(lSeq) == seqLength ):
	useLabels = 1
    else:
	useLabels = 0

    # write out the sequence 
    for ii in range(seqLength):

	tmpId = aaId ( aSeq[ii] )
	if ( useLabels ): tmpLabel = labelId ( lSeq[ii] )

	# if any of the ambiguous codes are used, map them to something
	# in the range [0,20] ...
	if ( tmpId == 21 ): tmpId = 0
	if ( tmpId == 22 ):
	    r = random.random()
	    if ( r < 0.5 ): 
		tmpId =  3
	    else:
		tmpId = 12
	if ( tmpId == 23 ):
	    r = random.random()
	    if ( r < 0.5 ):
		tmpId =  4
	    else:
		tmpId = 14

	# and finally write it out
	if ( useLabels ):
            fh.write ( '%2d  %1d \n' % ( tmpId, tmpLabel ) )
	else:
            fh.write ( '%2d  0 \n' % tmpId )

    fh.write('\n')
    fh.close()

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# main program is here -- assumes command-line argument giving name of .faa
# file to be read

import os
import commands
import sys
import random

if __name__ == "__main__":

    if not ( len(sys.argv) == 2 ):
	print ' Usage : %s <faa filename> ' % sys.argv[0] 

    faaFilename = sys.argv[1]

    # we need to assume that this FASTA file might be HUGE ... so we only
    # want to handle ONE sequence at a time ...

    try:
        fh = file ( faaFilename )
    except:
	print ' FATAL ERROR : failed to open %s ' % faaFilename
	sys.exit(-1)

    listFileName = ''

    done = 0
    numSeqs = 0
    while not done:

	( hdr1, hdr2, aSeq, lSeq ) = readNextSequenceFromFASTA ( fh )
	if ( hdr1==''  and  hdr2==''  and aSeq=='' ):
	    done = 1
	else:

	    if ( len(aSeq) < 3 ): 
		print ' '
		print ' WARNING : skipping too short sequence ', hdr1, hdr2, aSeq
		print ' '
		continue

	    dirName = getDirName ( faaFilename )
	    gmtkFilename = makeGMTKobsFilename ( hdr1 ) 
	    fullName = dirName + '/' + gmtkFilename

            print listFileName

	    if ( listFileName == '' ):
	        listFileName = dirName + '/' + 'input.list'
	        fhList = file ( listFileName, 'w' )

	    writeGMTKobsFile ( aSeq, lSeq, fullName )

	    fhList.write ( '%s/%s\n' % ( dirName, gmtkFilename ) )

	    numSeqs += 1
    
    fhList.close()

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
