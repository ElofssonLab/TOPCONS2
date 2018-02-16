#!/usr/bin/env python
# Description:
#   This script is a variant of pfam_workflow.py which use the same pipeline as
#   TOPCONS2 but runs only one of the sub-predictors, OCTOPUS

# Created 2016-04-05, updated 2016-04-05, Nanjiang Shu

import sys
import os
import time
from Bio import SeqIO
import tempfile
import subprocess
import module_locator
import ntpath
import myfunc
import shutil
import argparse


rundir = os.path.dirname(os.path.realpath(__file__))
basedir = os.path.realpath("%s/../"%(rundir))  #path to topcons2_webserver
os.environ['TOPCONS2'] = basedir

def main(args, g_params):
    parser = argparse.ArgumentParser(
            description='TOPCONS2_OCTOPUS workflow master script',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Created 2015-05-05, updated 2018-02-16, Peters Christoph and Nanjiang Shu

Examples:
''')
    parser.add_argument('inFile', metavar='inFile',
            help='Specify the input amino acid sequence file in FASTA format')
    parser.add_argument('out_path', metavar='out_path',
            help='Specify the outpath for result')
    parser.add_argument('blastDir', metavar='blastDir',
            help='Specify the path for psiblast, which contains bin/blastpgp')
    parser.add_argument('blastDB', metavar='blastDB',
            help='Specify the name of the blastdb, including the path')
    parser.add_argument('-tmpdir', '--tmpdir', metavar='DIR', dest='TMPPATH', 
            help='Specify the directory where the temporary files will be written to')
    parser.add_argument('-debug', '--debug', action='store_true', default=False,  dest='isDEBUG', 
            help='Output debug info')
    parser.add_argument('-RM', '--remove-individual-files', action='store_true', default=False,  dest='isRemoveFile', 
            help='Delete result for individual sequences')

    args = parser.parse_args()

    g_params['DEBUG'] = args.isDEBUG
    g_params['REMOVE_IND_FILES'] = args.isRemoveFile
    inFile = os.path.abspath(args.inFile)
    out_path = os.path.abspath(args.out_path)
    blastDir = os.path.abspath(args.blastDir)
    blastDB = os.path.abspath(args.blastDB)
    if args.TMPPATH != None:
        g_params['TMPPATH'] = os.path.abspath(args.TMPPATH)

    if not os.access(g_params['TMPPATH'], os.W_OK):
        print >> sys.stderr, "Error. TMPPATH '%s' not writable. Exit."%(g_params['TMPPATH'])
        return 1
    if not os.access(out_path, os.W_OK):
        print >> sys.stderr, "Error. out_path '%s' not writable. Exit."%(out_path)
        return 1

    os.environ['TMPPATH'] = g_params['TMPPATH']

    DEBUG = g_params['DEBUG']
    TMPPATH = g_params['TMPPATH']
    if not os.path.exists(inFile):
        print >> sys.stderr, "inFile %s does not exist. Exit."%(inFile)
        sys.exit(1)
    if not os.path.exists(out_path):
        try:
            os.makedirs(out_path)
        except OSError:
            print >> sys.stderr, "Failed to create out_path %s. Exit."%(out_path)
            sys.exit(1)

    if not "BLASTDB" in os.environ: # this fixed the warning message of unset 'BLASTDB'
        try:
            blastdbpath = os.path.realpath(os.path.dirname(blastDB))
            os.environ['BLASTDB'] = blastdbpath
        except:
            pass

    # Set the working dir to the script location
    my_path = module_locator.module_path()
    os.chdir(my_path)
    inFile_rootname = os.path.basename(os.path.splitext(inFile)[0])

    # Timing remove from final version
    #print "Timing remove from final version"
    timingfile = "%s/%s"%(out_path, "time.txt")
    topfile_OCTOPUS = "%s/%s.OCTOPUS.topfa"%(out_path, inFile_rootname)
    topfile_SPOCTOPUS = "%s/%s.SPOCTOPUS.topfa"%(out_path, inFile_rootname)
    fpout_OCTOPUS = open(topfile_OCTOPUS, "w")
    fpout_SPOCTOPUS = open(topfile_SPOCTOPUS, "w")
    with open(timingfile, "w") as timingFileOut:
        with open(inFile, "rU") as seqFile:
            for index, entry in enumerate(list(SeqIO.parse(seqFile, "fasta"))):
                # Timing remove from final version
#                 print "Timing remove from final version"
                start = time.time()


                #Create folders for tmp data and output
                used_pfam = "pfam"
                tmpDir = tempfile.mkdtemp(prefix="%s/seq_"%(TMPPATH) + str(index) + "_") + "/"
                os.chmod(tmpDir, 0755)
                tmpDir_pfam = tmpDir
                tmpDir_cdd = ""
                tmpDir_uniref = ""

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

#                 if os.path.exists(outDir + "Topcons/") is False:
#                     os.mkdir(outDir + "Topcons/")

#                 outfile = "%s/%s"%(tmpDir, "query.fa")
                with open(tmpDir + "query.fa", "w") as outFile:
                    outFile.write(">query" + "\n" + str(entry.seq))

                with open(outDir + "seq.fa", "w") as outFile:
                    outFile.write(">query" + "\n" + str(entry.seq))

                startDir = os.getcwd()


                # At the same time the profiles can be created
                cmd = ["./fa2prfs_pfamscan_v2.sh", tmpDir_pfam, blastDir]
                cmdline = " ".join(cmd)
                rmsg = ""
                try:
                    print "cmdline: ", cmdline
                    rmsg = subprocess.check_call(cmd, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError, e:
                    print "errmsg:", e
                    print "rmsg:", rmsg
                    pass
                query_seqdbfile = "%s/%s"%(tmpDir_pfam, "query.hits.db")
                filesize = 0
                try:
                    filesize = os.path.getsize(query_seqdbfile)
                except OSError:
                    filesize = -1
                    pass
                if DEBUG:
                    print "After fa2prfs_pfamscan_v2.sh filesize(%s)=%d"%(query_seqdbfile, filesize)

                # In case we do not find a hit, we have to run hmmscan on the cdd database
                if filesize <= 0:
                    tmpDir_cdd = tempfile.mkdtemp(prefix="%s/seq_cdd_"%(TMPPATH) + str(index) + "_") + "/"
                    os.chmod(tmpDir_cdd, 0755)
                    with open(tmpDir_cdd + "query.fa", "w") as outFile:
                        outFile.write(">query" + "\n" + str(entry.seq))
                    used_pfam = "cdd"
                    cmd = ["./fa2prfs_hmmscan.sh", tmpDir_cdd, blastDir]
                    cmdline = " ".join(cmd)
                    try:
                        print "\ncmdline:",cmdline
                        rmsg = subprocess.check_call(cmd, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError, e:
                        print "errmsg:", e
                        print "rmsg:", rmsg
                        pass

                    tmpDir = tmpDir_cdd

                    query_seqdbfile = "%s/%s"%(tmpDir_cdd, "query.hits.db")
                    try:
                        filesize = os.path.getsize(query_seqdbfile)
                    except OSError:
                        filesize = -1
                        pass

                    if DEBUG:
                        print "After fa2prfs_hmmscan.sh filesize(%s)=%d"%(query_seqdbfile, filesize)
                # In case we do not find a hit, we have to run the old script
                if filesize <= 0:
                    tmpDir_uniref = tempfile.mkdtemp(prefix="%s/seq_uniref_"%(TMPPATH) + str(index) + "_") + "/"
                    os.chmod(tmpDir_uniref, 0755)
                    with open(tmpDir_uniref + "query.fa", "w") as outFile:
                        outFile.write(">query" + "\n" + str(entry.seq))
                    used_pfam = "uniref"
                    cmd =  ["./fa2prfs_fallback_v2.sh", tmpDir_uniref, blastDir, blastDB]
                    cmdline = " ".join(cmd)
                    try:
                        print "\ncmdline:", cmdline
                        rmsg = subprocess.check_call(cmd, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError, e:
                        print e
                        print rmsg
                        pass
                    tmpDir = tmpDir_uniref

                    query_seqdbfile = "%s/%s"%(tmpDir_uniref,"query.hits.db")
                    try:
                        filesize = os.path.getsize(query_seqdbfile)
                    except OSError:
                        filesize = -1
                        pass

                    if DEBUG:
                        print "After fa2prfs_fallback_v2.sh filesize(%s)=%d"%(query_seqdbfile, filesize)

                # Once the profile is created start all other predictors

                os.chdir(os.path.abspath("../predictors/spoctopus/"))
                outDir_SPOCTOPUS = outDir + os.sep + "SPOCTOPUS"
                if not os.path.exists(outDir_SPOCTOPUS):
                    os.makedirs(outDir_SPOCTOPUS)
                cmd = ["./SPOCTOPUS.sh", protnamefile ,tmpDir + "PSSM_PRF_FILES/", tmpDir + "RAW_PRF_FILES/", outDir_SPOCTOPUS, "-N"] #output also the ANN result for SPOCTOPUS, changed 2016-01-26
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline
                p_spoctopus = subprocess.Popen(cmd)
                os.chdir(startDir)

                os.chdir(os.path.abspath("../predictors/spoctopus/"))
                outDir_OCTOPUS = outDir + os.sep + "OCTOPUS"
                if not os.path.exists(outDir_OCTOPUS):
                    os.makedirs(outDir_OCTOPUS)
                cmd =  ["./OCTOPUS.sh", protnamefile,tmpDir + "PSSM_PRF_FILES/", tmpDir + "RAW_PRF_FILES/", outDir_OCTOPUS, "-N"] #output also the ANN result for OCTOPUS, changed 2016-01-26
                cmdline = " ".join(cmd)
                if DEBUG:
                    print "cmdline:", cmdline

                p_octopus = subprocess.Popen(cmd)
                os.chdir(startDir)


                p_spoctopus.communicate() #now wait for OCTOPUS
                p_octopus.communicate() #now wait for SPOCTOPUS
                count_pred = 2


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
                        p = subprocess.call(["rm", "-rf", tmpDir])
                    if os.path.exists(tmpDir_cdd) is True:
                        p = subprocess.call(["rm", "-rf", tmpDir_cdd])
                    if os.path.exists(tmpDir_uniref) is True:
                        p = subprocess.call(["rm", "-rf", tmpDir_uniref])
                    if os.path.exists(tmpDir_pfam) is True:
                        p = subprocess.call(["rm", "-rf", tmpDir_pfam])
                else:
                    print "tmpDir=%s"%(tmpDir)

                p = subprocess.call(["python","correct_Topo.py", outDir])

                topfile = "%s/%s/%s"%(outDir, "OCTOPUS", "query.top")
                if os.path.exists(topfile):
                    top = myfunc.ReadFile(topfile).strip()
                    if top:
                        fpout_OCTOPUS.write(">%s\n"%(entry.description))
                        fpout_OCTOPUS.write("%s\n"%(top))

                topfile = "%s/%s/%s"%(outDir, "SPOCTOPUS", "query.top")
                if os.path.exists(topfile):
                    top = myfunc.ReadFile(topfile).strip()
                    if top:
                        fpout_SPOCTOPUS.write(">%s\n"%(entry.description))
                        fpout_SPOCTOPUS.write("%s\n"%(top))

                if g_params['REMOVE_IND_FILES']:
                    shutil.rmtree(outDir)

    fpout_SPOCTOPUS.close()
    fpout_OCTOPUS.close()

if __name__=="__main__":
    g_params = {}
    g_params['DEBUG'] = False
    g_params['REMOVE_IND_FILES'] = False
    g_params['TMPPATH']="/tmp"
    if os.path.exists("/scratch") and os.access("/scratch", os.W_OK):
        g_params['TMPPATH']="/scratch"
    sys.exit(main(sys.argv, g_params))
