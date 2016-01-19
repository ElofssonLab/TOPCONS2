#!/bin/sh
set -e
gengetoptversion=$1
installDir="$2"
tar xfz "../download/gengetopt-${gengetoptversion}.tar.gz"
cd gengetopt-${gengetoptversion}
./configure "--prefix=$installDir" && make install
