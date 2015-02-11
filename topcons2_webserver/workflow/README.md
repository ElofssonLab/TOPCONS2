# FastPSSM

Before running:

There are three pathes in pfam_scan_to_profile.py that need to be adjusted:

+ pfam_Dir: pfam database for pfam_scan
+ pfamseqdb: database containing the full length sequences for each family
+ pfamScan: location of pfam_scan


Running the workflow:

./pfam_workflow.py inFile outDir blastDir

inFile can contain one or many sequences in fasta format.
