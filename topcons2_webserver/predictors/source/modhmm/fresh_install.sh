#!/bin/bash

prefix=$1

if [ "$prefix" == "" ]; then 
    echo "usage: $0 ROOT_DIR"
    exit 1
fi

if [ ! -d "$prefix" ]; then 
    mkdir -p $prefix
fi

tmpdir=$(mktemp -d $PWD/tmpdir.install.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }  

cd $tmpdir


cmake  -DCMAKE_INSTALL_PREFIX=$prefix/modhmm ../

make
make install

rm -rf $tmpdir
#echo $tmpdir
