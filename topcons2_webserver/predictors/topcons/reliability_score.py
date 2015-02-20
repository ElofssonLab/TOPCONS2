#!/usr/bin/python

import sys,re,os;

def prf2reliability(prffile, topfile, outfile) :
    intWindowSize = 10; # no positions to each side of aa
    
    # Parse prffile
    strSeqName = '';
    intSeqLength = 0;
    dctProfile = dict();
    flhPrf = file(prffile);
    for line in flhPrf :
        if line.startswith('Sequence') :
            strSeqName = line.strip().split()[-1];
        elif line.startswith('Length of query sequence') :
            intSeqLength = int(re.match('Length of query sequence:\s+(\d+)\s+$', line).group(1))
        # ALPHABET:   i       M       o       -     <SPACE> <LABEL> <QUERY>
        # COL    1:   0.00    0.00    100.00  0.00    0.00    0.00    .
        # ...
        elif line.startswith('COL') :
            lstCols = line.strip().split();
            dctProfile[lstCols[1].rstrip(':')] = {'i': float(lstCols[2]), 
                                                  'M': float(lstCols[3]),
                                                  'o': float(lstCols[4]),
                                                  'S': float(lstCols[5]),
						  'p': float(lstCols[6])};
    flhPrf.close();

    if len(dctProfile) != intSeqLength :
        print >>sys.stderr, "Profile '%s' doesn't have the right number of 'COL' elements: %d (header says %d)" % (prffile, len(Profile), intSeqLength);
        sys.exit(1);
    
    # Parse topfile
    flhTop = file(topfile);
    strTitle = flhTop.readline();
    strTop = flhTop.read().strip().replace('\n','');
    flhTop.close();
    
#     print >>sys.stderr, strTop;

    # compute reliability scores
    lstRelScores = [];
    intAAPos = intWindowSize + 1;
    # loop over AA positions
    while intAAPos < intSeqLength - intWindowSize :
        # loop over the window around the AA pos and
        # sum the 'certainty' for the positions' particular topologies
        fltSum = 0.0;
        for intWindowPos in range(-intWindowSize, intWindowSize + 1) :
            i = intAAPos + intWindowPos;
            chrTop = strTop[i-1]; # working with one-based index, strings use zero-based
            fltSum += dctProfile[str(i)][chrTop];
#             print >>sys.stderr, intAAPos, i, chrTop, dctProfile[str(i)][chrTop], fltSum;
        # rel score is the avg of the scores in the window
        lstRelScores.append( fltSum/(intWindowSize * 2 + 1) );
#         print >>sys.stderr, intAAPos, fltSum, intWindowSize * 2 + 1,  fltSum/(intWindowSize * 2 + 1);
        intAAPos += 1;
        

    # write to outfile
    flhOut = file(outfile, 'w');
    intCounter = intWindowSize + 1;
    for fltScore in lstRelScores :
        print >>flhOut, "%d\t%.2f" % (intCounter, fltScore);
        intCounter += 1;
    flhOut.close();
    
    return;


if __name__ == '__main__' :
    # Check, parse $argv
    strUsage = "Usage: %s <topology profile file> <topology file> <reliability score (out) file>";
    if ( len(sys.argv) < len(strUsage.split('<')) ) :
        print >>sys.stderr, strUsage % os.path.basename(sys.argv[0]);
        sys.exit(1);

    prf2reliability(sys.argv[1], sys.argv[2], sys.argv[3]);
    

