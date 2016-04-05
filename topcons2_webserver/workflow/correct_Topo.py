import sys
import os

def main(args):
    filename =  args[1] + "/SPOCTOPUS/query.top"
    if os.path.exists(filename):
        with open(filename, "r+") as inFile:
            temp_line = ""
            for line in inFile:
                temp_line += line.strip()
            inFile.seek(0)
            inFile.write(temp_line.replace("p", "M"))

    filename =  args[1] + "/PolyPhobius/query.top"
    if os.path.exists(filename):
        with open(filename, "r+") as inFile:
            temp_line = ""
            for line in inFile:
                temp_line += line.strip()
            inFile.seek(0)
            inFile.write(temp_line.replace("p", "M"))

    filename =  args[1] + "/philius/query.top"
    if os.path.exists(filename):
        with open(filename,"r+") as inFile:
            temp_line = ""
            for line in inFile:
                temp_line += line.strip()
            inFile.seek(0)
            inFile.write(temp_line.replace("p", "M"))


if __name__=="__main__":
    sys.exit(main(sys.argv))
