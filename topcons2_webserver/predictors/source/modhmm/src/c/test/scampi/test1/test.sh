#!/bin/sh

set -e
/data3/modhmm-projects/trunk/modhmm/src/c/modhmms_scampi  -f fa -s /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/snf.1 -m /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/hmg.list -o /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1 -r /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/replacement_letter_multi.rpl --nopostout --nolabels --viterbi -u -L -g | /usr/bin/xsltproc /data3/modhmm-projects/trunk/modhmm/src/xslt/xml2res.xsl - > /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/DGHMM_KR_21_single.1.hmg.res

diff --brief /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/DGHMM_KR_21_single.1.hmg.res /data3/modhmm-projects/trunk/modhmm/src/c/test/scampi/test1/DGHMM_KR_21_single.1.hmg.res
