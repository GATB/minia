#!/bin/bash

# look for minia binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/minia" ]
then
 bindir="../bin"
elif [ -f "../build/bin/minia" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled Minia binary"
 exit 1
fi

verbose="0"
if [[ "$1" != "" ]]
then
    verbose="1"
fi

function test()
{
    file=$1
    k=$2
    traversal=$3
    args="-in "$file".fa -kmer-size "$k" -abundance-min 1 -traversal $traversal"
    if [[ $verbose == "0" ]]
    then
        eval $bindir/minia $args > /dev/null 2>/dev/null
    else
        eval $bindir/minia $args
    fi
    python compare_fasta.py "$file".solution.fa "$file".contigs.fa
}

#test X 5 unitig
#echo "has a known bug where the AAAA node gets repeated a few times.., so it should print NOT EQUAL. but should print 6 below:"
#grep ">" X.contigs.fa|wc -l

test tip 21 contig

test bubble 21 contig

#test ec 21 contig
