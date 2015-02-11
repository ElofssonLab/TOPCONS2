import sys
import os
import linecache
from random import randrange
import myfunc

cddseqdb = "../database/cddfamseqdb_nr90/uniref100.cddfamseq.nr90"
cdd_hmm = "../database/cdd_hmm/cdd.hmm"

def createHitDB(pfamList, prot_name, work_dir):
    hdl = myfunc.MyDB(cddseqdb)
    if hdl.failure:
        print "Error"
        return 1
    with open(work_dir + prot_name + ".hits.db.temp", "w") as outFile:
        for pfamid in pfamList:
            record = hdl.GetRecord(pfamid)
            if record:
                outFile.write(record)
        hdl.close()

    os.system("python my_uniqueseq.py " + work_dir + prot_name + ".hits.db.temp")


def main(argvs):
    input_file = argvs[1]
    work_dir = argvs[2]

    name_temp = (input_file[input_file.rfind("/")+1:])
    name = name_temp[:name_temp.rfind(".")]
    startDir = os.getcwd()
    #os.chdir(os.path.abspath("../database/pfam_seq/PfamScan/"))
    sCmd = "hmmscan -E 0.1 --tblout " + work_dir + name + ".txt " + cdd_hmm + " " + input_file
    os.system(sCmd)
    #os.chdir(startDir)
    pfamList = []
    bFoundStart = False
    with open(work_dir + name + ".txt") as inFile:
        for line in inFile:
            if line.find("#") == -1:
                pfamList.append(line.split(" ")[0])

    createHitDB(list(set(pfamList)), name, work_dir)

if __name__ == "__main__":
    main(sys.argv)
