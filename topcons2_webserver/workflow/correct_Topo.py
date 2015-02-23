import sys
import os

def main(args):
    with open(args[1] + "/SPOCTOPUS/query.top", "r+") as inFile:
        temp_line = ""
        for line in inFile:
            temp_line += line.strip()
        inFile.seek(0)
        inFile.write(temp_line.replace("p", "M"))

    with open(args[1] + "/PolyPhobius/query.top", "r+") as inFile:
        temp_line = ""
        for line in inFile:
            temp_line += line.strip()
        inFile.seek(0)
        inFile.write(temp_line.replace("p", "M"))

    with open(args[1] + "/philius/query.top","r+") as inFile:
        temp_line = ""
        for line in inFile:
            temp_line += line.strip()
        inFile.seek(0)
        inFile.write(temp_line.replace("p", "M"))


if __name__=="__main__":
    sys.exit(main(sys.argv))
