import os
import sys
from Bio import SeqIO

def readFiles(inPath):
    seqTemp = ""
    with open(inPath + "query.top") as inFile:
        for line in inFile:
            seqTemp += line.strip()
    return seqTemp

def printPrf(prf_List):
    outString = ""
    outString += "Sequence: query\n"
    outString += "NR of aligned sequences: 100\n"
    outString += "Length of query sequence: " + str(len(prf_List)) + "\n\n"
    outString += "START 1\n"
    outString += "ALPHABET:   i       M       o       S       p       -     <SPACE> <LABEL> <QUERY>\n"
    for index, entry in enumerate(prf_List):
        outString += "COL" + str.rjust(str(index + 1), 5) + ":   " + entry + "0.00    0.00    0.00    0.00    .\n"

    outString += "END 1"

    print outString


def main(argvs):
    inTopoDir = argvs[1]

    topSeqs = []
    method_list = []
    number_sp_methods = 0

    try:
        topSeqs.append(readFiles(inTopoDir + "SPOCTOPUS/"))
        method_list.append("spoctopus")
        number_sp_methods += 1
    except:
        # file not generated
        pass
    try:
        topSeqs.append(readFiles(inTopoDir + "philius/"))
        method_list.append("philius")
        number_sp_methods += 1
    except:
        # file not generated
        pass
    try:
        topSeqs.append(readFiles(inTopoDir + "PolyPhobius/"))
        method_list.append("polyPhobius")
        number_sp_methods += 1
    except:
        # file not generated
        pass
    try:
        topSeqs.append(readFiles(inTopoDir + "OCTOPUS/"))
        method_list.append("octopus")
    except:
        # file not generated
        pass
    try:
        topSeqs.append(readFiles(inTopoDir + "SCAMPI_MSA/"))
        method_list.append("scampimsa")
    except:
        # file not generated
        pass

    if len(topSeqs) == 0:
        print "Error: All predictors failed!"
        exit()

    prf_list = []
    foundTM = False
    for i in range(0, len(topSeqs)):
        if topSeqs[i].count("M") > 0:
            foundTM = True
            break

    for i in range(0, len(topSeqs[0])):
        dict_pos = {"i":0, "M":0, "o":0, "S":0, "p":0}
        for entry in topSeqs:
            dict_pos[entry[i]] += 1
        if dict_pos["S"] > 0 and foundTM is False:
            dict_pos = {"i":0, "M":0, "o":0, "S":0, "p":0}
            for entry in topSeqs[0:number_sp_methods]:
                    dict_pos[entry[i]] += 1
        predictors_used = sum(dict_pos.values())
        prf_list.append(str.ljust(str("%.2f" % (dict_pos["i"]/float(predictors_used)*100)), 8) + str.ljust(str("%.2f" % (dict_pos["M"]/float(predictors_used)*100)), 8) + str.ljust(str("%.2f" % (dict_pos["o"]/float(predictors_used)*100)), 8) + str.ljust(str("%.2f" % (dict_pos["S"]/float(predictors_used)*100)), 8) + str.ljust(str("%.2f" % (dict_pos["p"]/float(predictors_used)*100)), 8))

    printPrf(prf_list)

if __name__=="__main__":
    main(sys.argv)
