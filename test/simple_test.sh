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
    args="-in "$file".fa -kmer-size "$k" -abundance-min 1 -traversal $traversal -minimizer-size 3"
    if [[ $verbose == "0" ]]
    then
        eval $bindir/minia $args > /dev/null 2>/dev/null
    else
        eval $bindir/minia $args
    fi
    python compare_fasta.py "$file".solution.fa "$file".contigs.fa
}

# Test One
test X 5 unitig
echo "has a known issue where the AAAA node gets repeated a few times, so it should print NOT EQUAL. But should print '6', now:"
lines=`grep ">" X.contigs.fa|wc -l`
lines="$(sed -e 's/[[:space:]]*$//' <<<${lines})"
echo "  $lines"
if [ $lines -eq 6 ]; then
  echo  PASSED
else
  echo  FAILED
  exit 1
fi

# Test Two
test tip 21 contig
if [ $? -eq 0 ]; then
  echo  PASSED
else
  echo  FAILED
  exit 1
fi

# Test Three
test bubble 21 contig
if [ $? -eq 0 ]; then
  echo  PASSED
else
  echo  FAILED
  exit 1
fi

# Test Four
test ec 21 contig
if [ $? -eq 0 ]; then
  echo  PASSED
else
  echo  FAILED
  exit 1
fi

# do some cleanup
rm -f *.h5 *.contigs.fa
