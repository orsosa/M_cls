#!/bin/bash
outdir=$1
source ~/.bashrc

echo $LD_LIBRARY_PATH
ps

shift
opt="$@"

echo "options in bash: "$opt

./run_M "*.root" $opt quiet
echo "outdir: "$outdir

dir="`ls -d */ | xargs -I V bash -c 'a=V;echo ${a/\//}'`"
echo "dir to be compressed: "$dir
tar -czf $dir.tar.gz $dir
ls

swif outfile $dir.tar.gz file:$outdir/$dir.tar.gz
