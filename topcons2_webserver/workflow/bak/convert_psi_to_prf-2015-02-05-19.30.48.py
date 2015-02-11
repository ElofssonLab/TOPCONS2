import os
import sys

def fixAAorder(count_array):
    dic_aa = {}
    aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    int_count_array = [int(entry) for entry in count_array ]
    if sum(int_count_array) == 0:
        int_count_array = [entry + 1 for entry in int_count_array]
    freq_array_temp = [entry / float(sum(int_count_array)) for entry in int_count_array]
    freq_array = []

    for i in range(0, len(aa_list)):
        entry = aa_list[i]
        dic_aa[entry] = freq_array_temp[i]

    for key in sorted(dic_aa.iterkeys()):
        freq_array.append(str("%.2f" % dic_aa[key]))

    return freq_array



def main(file_name):
    freqList = []
    seqList = []
    with open(file_name) as inFile:
        bFoundStart = False
        for line in inFile:
            if line.strip() == "" and bFoundStart is True:
                break

            if bFoundStart is True:
                freqList.append(line.strip().split()[22:42])
                seqList.append(line.strip().split()[1])

            if line.replace(" ","").find("ARN") != -1:
                bFoundStart = True

    with open(file_name.replace(".psi", ".raw.prf"), "w") as outFile:#file_name.replace(".psi", ".raw.prf"), "w") as outFile:
        outFile.write("Sequence: query" + "\n") #+ file_name.replace(".psi", "") + "\n")
        outFile.write("NR of aligned sequences: 100"  + "\n")
        outFile.write("Length of query sequence: " + str(len(freqList)) + "\n\n\n")
        outFile.write("START 1" + "\n")
        outFile.write("ALPHABET:       A       C       D       E       F       G       H       I       K       L       M       N       P       Q       R       S       T       V       W       Y       -     <SPACE> <LABEL> <QUERY>" + "\n")
        for i in range(0, len(seqList)):
            freq_str = fixAAorder(freqList[i])
            #freq_str = [str.ljust(entry, 8) for entry in freq_str]
            outFile.write("COL" + str.rjust(str(i + 1), 5) + ":       " + "".join([str.ljust(entry, 8) for entry in freq_str]) + "0.00    0.00    .       " + seqList[i] + "       \n")
        outFile.write("END 1")




if __name__=="__main__":
    main(sys.argv[1])
