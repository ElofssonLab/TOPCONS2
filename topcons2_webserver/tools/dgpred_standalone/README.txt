==LOG: 2010-09-01 15:58:38 Wednesday  Week 35 <nanjiang@illergard>
Three script files

analyze.pl
calc_dG.pl
dG_predTM.pl

have been rewritten into the standalone version, so that the programs can be
run by e.g.

./analyze.pl  YIYLGGAILAEVIGTTLMKF -o test/rst.html
./calc_dG.pl test/fragfile.txt -o test/rst.html
./dG_predTM.pl test/test.fa -o test/rst.html

==LOG: 2010-09-26 14:21:55 Sunday  Week 39 <nanjiang@illergard>
a script myscanDG.pl has been modified from DGpred_multi_scan.pl
for calculating DG values from either single sequence or multiple sequence
alignment

# 1. calculate DG values for 'test.fa' based on single sequence
     ./myscanDG.pl test/test.fa
# 2. calculate DG values for 'test.fa' based on multiple sequence alignment
     ./myscanDG.pl -multi test/test.fa
# 3. calculate DG values from the prf file 
     ./myscanDG.pl test/test.prf

