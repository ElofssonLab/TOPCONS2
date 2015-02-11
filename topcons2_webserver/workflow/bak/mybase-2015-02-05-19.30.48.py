#!/usr/bin/env python

# Description:
#   this is a base module, it will not import any self created modules
#
# Author: Nanjiang Shu (nanjiang.shu@scilifelab.se)
#
# Address: Science for Life Laboratory Stockholm, Box 1031, 17121 Solna, Sweden

def FloatDivision(x1, x2):#{{{
    """
    Return the division of two values
    """
    try:
        return float(x1)/x2
    except ZeroDivisionError:
        return 0.0
#}}}
class ReadLineByBlock:#{{{
# Description: readlines by BLOCK reading, end of line is not included in each
#              line. Empty lines are not ignored 
#              reading lines at about 3 times faster than normal readline()
# Function:
#   readlines()
#   close()
#
# Usage:
# handel = ReadLineByBlock(infile)
# if handel.failure:
#   print "Failed to init ReadLineByBlock for file", infile
#   return 1
# lines = handel.readlines()
# while lines != None:
#       do_something
#       lines = handel.readlines()
    def __init__(self, infile, BLOCK_SIZE=100000):#{{{
        self.failure = False
        self.filename = infile
        self.BLOCK_SIZE = BLOCK_SIZE
        self.isEOFreached = False
        try: 
            self.fpin = open(infile, "rb")
        except IOError:
            print >> sys.stderr, "Failed to read file %s"%(self.filename)
            self.failure = True
            return None
        self.unprocessedBuff = ""
#}}}
    def __del__(self):#{{{
        try:
            self.fpin.close()
        except IOError:
            print >> sys.stderr, "Failed to close file %s"%(self.filename)
            return 1
#}}}
    def close(self):#{{{
        try:
            self.fpin.close()
        except IOError:
            print >> sys.stderr, "Failed to close file %s"%(self.filename)
            return 1
#}}}
    def readlines(self):#{{{
        if self.isEOFreached and not self.unprocessedBuff:
            return None
        else:
            buff = self.fpin.read(self.BLOCK_SIZE)
            if not buff:
                self.isEOFreached = True
            if self.unprocessedBuff:
                buff = self.unprocessedBuff + buff

            strs = buff.split('\n')
            numStrs = len(strs)
            if not self.isEOFreached:
                self.unprocessedBuff = strs[numStrs-1]
                return strs[:numStrs-1]
            else:
                self.unprocessedBuff = ""
                return strs
#}}}
#}}}
