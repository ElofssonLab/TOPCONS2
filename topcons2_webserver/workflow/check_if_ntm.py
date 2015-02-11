import os
import sys

def main(args):
    scampi_res_file = args[1]
    octopus_res_file = args[2]
    spoctopus_res_file = args[3]

    predTopo = ""
    # Check if Scampi predicts a nTM protein
    with open(scampi_res_file, "rU") as predictedHandle:
        for line in predictedHandle:
            if line.find(">") == -1:
                predTopo = line.strip()

    # If there is a membrane we do nothing
    if predTopo.count("M") > 0:
        exit()

    # Else we make sure that (SP)Octopus and Scampi have only "S" and "o"
    predTopo.replace("i", "o")

    with open(scampi_res_file, "w") as outFile:
        outFile.write(predTopo)

    with open(octopus_res_file, "w") as outFile:
        outFile.write(predTopo)

    predTopo = ""
    with open(spoctopus_res_file, "rU") as predictedHandle:
        for line in predictedHandle:
            if line.find(">") == -1:
                predTopo = line.strip()

    with open(spoctopus_res_file, "w") as outFile:
        outFile.write(predTopo.replace("i", "o").replace("M", "o"))

if __name__=="__main__":
    main(sys.argv)
