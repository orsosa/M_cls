#!/bin/bash
if [ -z "$1" ];then
    echo "you must supply a directory"
    exit 1
fi

cd $1
ls */p*/*.tar.gz | xargs -I V bash -c 'a=V; d=`dirname $a`; fn=`basename $a`; tar -C $d -xvzf $a'
cd ..
