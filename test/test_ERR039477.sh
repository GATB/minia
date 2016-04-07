#! /bin/bash

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

################################################################################
# we download some banks
################################################################################
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"

################################################################################
# we launch minia; note that we use only one thread (no real time issues with
# potential different results)
################################################################################
$bindir/minia -nb-cores 1 -in ERR039477.fastq.gz

################################################################################
# we check the result
################################################################################
md5sum ERR039477.contigs.fa > ERR039477.check

diff ./ERR039477.md5 ./ERR039477.check

if [ $? -eq 0 ]; then
   echo "TEST OK"
else
   echo "TEST KO"
fi

################################################################################
# clean up
################################################################################
rm -f  ERR039477.fastq.gz  ERR039477.check
