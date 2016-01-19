#!/bin/bash

set -e
/data3/modhmm-projects/trunk/modhmm/src/c/modhmms_spoctopus -m /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1/hmg.list -s /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1/query.list -f prf -o /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1 -L -M DP -v --nopostout -u --nolabels --viterbi  | /usr/bin/xsltproc /data3/modhmm-projects/trunk/modhmm/src/xslt/xml2res.xsl - > /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1/NNZHMM_io_reent_hairpin_sigpep.hmg.res

diff --brief /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1/NNZHMM_io_reent_hairpin_sigpep.hmg.res /data3/modhmm-projects/trunk/modhmm/src/c/test/spoctopus/test1/NNZHMM_io_reent_hairpin_sigpep.hmg.res
