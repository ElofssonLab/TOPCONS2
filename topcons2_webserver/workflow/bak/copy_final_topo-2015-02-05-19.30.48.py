import os

inPath = "../out/"
for entry in os.listdir(inPath):
    prot_name = entry.split("_")[0]
    os.system("cp " + inPath + entry + "/Topcons/topcons.top " + "../out_topo/" + prot_name + ".top")
