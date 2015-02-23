#!/usr/bin/env python
import sys
import os
import time
from Bio import SeqIO
import tempfile
import subprocess
import module_locator
import ntpath

DEBUG=False
# it seems apache on centos does not like /usr/local/bin
# /usr/local/bin can not be added to path
usage="""
Usage: %s inFile out_path blastDir blastDB
"""%(sys.argv[0])

def main(args):
    try:
        inFile = os.path.realpath(args[1])
        out_path= os.path.realpath(args[2]) + os.sep
        blastDir = args[3]
        blastDB = args[4]
    except IndexError:
        print >> sys.stderr, "Bad syntax"
        print usage
        sys.exit(1)

    # Set the working dir to the script location
    my_path = module_locator.module_path()
    os.chdir(my_path)

    # Timing remove from final version
    #print "Timing remove from final version"
    timingfile = "%s/%s"%(out_path, "time.txt")
    with open(timingfile, "w") as timingFileOut:
        with open(inFile, "rU") as seqFile:
            for index, entry in enumerate(list(SeqIO.parse(seqFile, "fasta"))):
                # Timing remove from final version
#                 print "Timing remove from final version"
                start = time.time()


                #Create folders for tmp data and output
                used_pfam = "pfam"
                tmpDir = tempfile.mkdtemp(prefix="seq_" + str(index) + "_") + "/"
#                 os.chmod(tmpDir, 0755)
                tmpDir_pfam = tmpDir
                tmpDir_cdd = ""
                tmpDir_uniref = ""

# Changed by Nanjiang at 2015-02-05 23:07:44, no random suffix in the folder,
# since the specified out_path should be exclusively for this query 
        #outDir = tempfile.mkdtemp(dir=out_path, prefix="seq_" + str(index) + "_") + "/"
                protnamefile = "%s/query.fa.txt"%(tmpDir)
                try:
                    fpout = open(protnamefile, "w")
                    print >> fpout, "query"
                    fpout.close()
                except IOError:
                    print >> sys.stderr, "Failed to write to protnamefile %s. "\
                            "Exit."%(protnamefile)
                    sys.exit(1)


                outDir = "%s%s%s/"%(out_path, os.sep, "seq_%d"%(index))
                if os.path.exists(tmpDir) is False:
                    os.mkdir(tmpDir)

                if os.path.exists(outDir) is False:
                    os.mkdir(outDir)

                if os.path.exists(outDir + "Topcons/") is False:
                    os.mkdir(outDir + "Topcons/")

#                 outfile = "%s/%s"%(tmpDir, "query.fa")
                with open(tmpDir + "query.fa", "w") as outFile:
                    outFile.write(">query" + "\n" + str(entry.seq))

                with open(outDir + "seq.fa", "w") as outFile:
                    outFile.write(">query" + "\n" + str(entry.seq))

                startDir = os.getcwd()

                # Run Philius because it does not need the profile
                os.chdir(os.path.abspath("../predictors/Philius/"))
                cmd = ["./runPhilius.sh", tmpDir + "query.fa", outDir]
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline
                p_philius = subprocess.Popen(cmd)
                os.chdir(startDir)

                # At the same time the profiles can be created
                cmd = ["./fa2prfs_pfamscan_v2.sh", tmpDir, blastDir]
                rmsg = ""
                try:
                    rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError, e:
                    print e
                    print rmsg
                # In case we do not find a hit, we have to run hmmscan on the cdd database
                query_seqdbfile = "%s/%s"%(tmpDir, "query.hits.db")
#                 if os.stat(tmpDir + "query.hits.db").st_size == 0:
                filesize = 0
                try:
                    filesize = os.stat(query_seqdbfile).st_size
                except OSError:
                    filesize = -1

                if filesize <= 0:
                    tmpDir_cdd = tempfile.mkdtemp(prefix="seq_" + str(index) + "_") + "/"
                    with open(tmpDir_cdd + "query.fa", "w") as outFile:
                        outFile.write(">query" + "\n" + str(entry.seq))
                    used_pfam = "cdd"
                    cmd = ["./fa2prfs_hmmscan.sh", tmpDir_cdd, blastDir]
                    try:
                        rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError, e:
                        print e
                        print rmsg

                    tmpDir = tmpDir_cdd
                # In case we do not find a hit, we have to run the old script

                try:
                    filesize = os.stat(query_seqdbfile).st_size
                except OSError:
                    filesize = -1

#                 if os.stat(tmpDir + "query.hits.db").st_size == 0:
                if filesize <= 0:
                    tmpDir_uniref = tempfile.mkdtemp(prefix="seq_" + str(index) + "_") + "/"
                    with open(tmpDir_uniref + "query.fa", "w") as outFile:
                        outFile.write(">query" + "\n" + str(entry.seq))
                    used_pfam = "uniref"
                    cmd =  ["./fa2prfs_fallback_v2.sh", tmpDir_uniref, blastDir, blastDB]
                    cmdline = " ".join(cmd)
                    if DEBUG:
                        print "cmdline:", cmdline
                    try:
                        rmsg = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError, e:
                        print e
                        print rmsg
                    tmpDir = tmpDir_uniref

                # Once the profile is created start all other predictors
                os.chdir(os.path.abspath("../predictors/scampi-msa/"))
                cmd = ["perl", "run_SCAMPI_MSA.pl", tmpDir , outDir]
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline
                p_scampi = subprocess.Popen(cmd)
                os.chdir(startDir)

                os.chdir(os.path.abspath("../predictors/spoctopus/"))
                cmd = ["./SPOCTOPUS.sh", protnamefile ,tmpDir + "PSSM_PRF_FILES/", tmpDir + "RAW_PRF_FILES/", outDir]
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline
                p_spoctopus = subprocess.Popen(cmd)
                os.chdir(startDir)

                os.chdir(os.path.abspath("../predictors/spoctopus/"))
                cmd =  ["./OCTOPUS.sh", protnamefile,tmpDir + "PSSM_PRF_FILES/", tmpDir + "RAW_PRF_FILES/", outDir]
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline

                p_octopus = subprocess.Popen(cmd)
                os.chdir(startDir)

                os.chdir(os.path.abspath("../predictors/PolyPhobius/"))
                cmd =  ["perl", "run_polyphobius.pl", tmpDir, outDir]
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline
                p_polyphobius = subprocess.Popen(cmd)
                os.chdir(startDir)

                p_philius.communicate() #now wait for philius
                p_scampi.communicate() #now wait for SCAMPI MSA
                p_spoctopus.communicate() #now wait for OCTOPUS
                p_octopus.communicate() #now wait for SPOCTOPUS
                p_polyphobius.communicate() #now wait for PolyPhobius

                # Kostas wants (SP)Octopus to be "o"*len when Scampi does not predict a membrane
                subprocess.call(["python","./check_if_ntm.py", outDir + "SCAMPI_MSA/query.top", outDir + "OCTOPUS/query.top", outDir + "SPOCTOPUS/query.top"])

                # When the predictors are done we can run topcons
                os.chdir(os.path.abspath("../predictors/topcons/"))
                results = []
                count_pred = 0
                try:
                    if os.stat(outDir + "OCTOPUS/query.top").st_size != 0:
                        results.append("%s/%s"%(outDir, "OCTOPUS/"))
                        count_pred += 1
                except OSError:
                    # In case of a random crash we predict without spoctopus 
                    pass

                try:
                    if os.stat(outDir + "SPOCTOPUS/query.top").st_size != 0:
                        results.append("%s/%s"%(outDir , "SPOCTOPUS/"))
                        count_pred += 1
                except OSError:
                    # In case of a random crash we predict without octopus
                    pass

                try:
                    if os.stat(outDir + "philius/query.top").st_size != 0:
                        results.append("%s/%s"%( outDir, "philius/"))
                        count_pred += 1
                except OSError:
                    pass

                try:
                    if os.stat(outDir + "PolyPhobius/query.top").st_size != 0:
                        results.append("%s/%s"%(outDir, "PolyPhobius/"))
                        count_pred += 1
                except OSError:
                    pass


                if os.stat(outDir + "SCAMPI_MSA/query.top").st_size != 0:
                    results.append("%s/%s"%( outDir, "SCAMPI_MSA/"))
                    count_pred += 1

                outdir_topcons = "%s/Topcons/"%(outDir)
                #cmd = ["./TOPCONS.sh", protnamefile, outdir_topcons] + results
                # 2015-02-23, the syntax of TOPCONS.sh has been changed
                cmd = ["./TOPCONS.sh", protnamefile, outdir_topcons, outDir]
                cmdline = " ".join(cmd)
                #os.system("./TOPCONS.sh " + tmpDir + "query.fa.txt " + outDir + "Topcons/ " + results)
                #print "./TOPCONS.sh " + tmpDir + "query.fa.txt " + outDir + "Topcons/ " + results 
                #subprocess.call(["./TOPCONS.sh", tmpDir + "query.fa.txt", outDir + "Topcons/", outDir + "OCTOPUS/", outDir + "philius/", outDir + "PolyPhobius/", outDir + "SPOCTOPUS/", outDir + "SCAMPI_MSA/"])

                if DEBUG:
                    print "cmdline: cmdline"
                try:
                    rmsg = subprocess.check_output(cmd)
                except subprocess.CalledProcessError, e:
                    print >> sys.stderr, str(e)
                os.chdir(startDir)
                os.chdir(os.path.abspath("../tools/dgpred_standalone/"))
                os.system("perl myscanDG.pl -o " + outDir + "dg.txt " + tmpDir + "query.fa")
                os.chdir(startDir)

                cmd = ["perl", "create_topcons_plot.pl", outDir + "/"]
                try:
                    rmsg = subprocess.check_output(cmd)
                except subprocess.CalledProcessError, e:
                    print e
                    print rmsg

                end = time.time()
                lines = 0
                with open(tmpDir + "query.hits.db") as inFile:
                    for line in inFile:
                        if line.find(">") == -1:
                            lines += 1
                timingFileOut.write(str(entry.id) + ";" + str(end - start) + ";" + used_pfam + ";" + str(lines) + ";" + str(count_pred) + "\n")
                #Remove the tmpFolder

                if not DEBUG: #debugging
                    if os.path.exists(tmpDir) is True:
                        p = subprocess.call(["rm", "-r", tmpDir])
                    if os.path.exists(tmpDir_cdd) is True:
                        p = subprocess.call(["rm", "-r", tmpDir_cdd])
                    if os.path.exists(tmpDir_uniref) is True:
                        p = subprocess.call(["rm", "-r", tmpDir_uniref])
                    if os.path.exists(tmpDir_pfam) is True:
                        p = subprocess.call(["rm", "-r", tmpDir_pfam])
                else:
                    print "tmpDir=%s"%(tmpDir)
                #exit()


if __name__=="__main__":
    sys.exit(main(sys.argv))
