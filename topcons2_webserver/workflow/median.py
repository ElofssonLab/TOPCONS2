import numpy
import os
import sys

def median(lst):
    return numpy.median(numpy.array(lst))

def main(filepath):
    num_list = []
    more30 = 0
    more60 = 0
    with open(filepath) as inFile:
        for line in inFile:
            num_list.append(float(line.split(";")[1]))
            if float(line.split(";")[1]) > 30:
                more30 += 1
            if float(line.split(";")[1]) > 60:
                more60 += 1

    print median(num_list)
    print numpy.average(num_list)
    print max(num_list)
    print min(num_list)
    print more30
    print more60

if __name__=="__main__":
    main(sys.argv[1])
