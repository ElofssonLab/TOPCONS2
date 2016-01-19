#!/bin/sh

set -e
/data3/modhmm-projects/trunk/modhmm/src/c/modhmms_octopus -m /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/hmg.list -s /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/query.list -f prf -o /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1 -r /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/replacement_letter_multi.rpl -L -M DP -v --nopostout -u --nolabels --viterbi  | /usr/bin/xsltproc /data3/modhmm-projects/trunk/modhmm/src/xslt/xml2res.xsl - > /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/NNZHMM_io_reent_hairpin.hmg.res 

diff --brief /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/NNZHMM_io_reent_hairpin.hmg.res /data3/modhmm-projects/trunk/modhmm/src/c/test/octopus/test1/NNZHMM_io_reent_hairpin.hmg.res
