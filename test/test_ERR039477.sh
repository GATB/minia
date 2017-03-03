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
# we download a sample bank from EBI
################################################################################
# if wget is not installed, you may use "curl -O ..."
DATA_SAMPLE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR039/ERR039477/ERR039477.fastq.gz"
WGET_PATH=`which wget`
echo ">>> Retrieving data sample: ${DATA_SAMPLE}"
if [ ! -z "$WGET_PATH" ] ; then
  echo "    using '$WGET_PATH'..."
  if [ $SILENT_MODE=="true"  ] ; then
    wget --quiet ${DATA_SAMPLE}
  else
    wget ${DATA_SAMPLE}
  fi
else
   CURL_PATH=`which curl`
  if [ ! -z "$CURL_PATH" ] ; then
    echo "    using '$CURL_PATH'..."
    if [ $SILENT_MODE=="true"  ] ; then
      curl --silent -O ${DATA_SAMPLE}
    else
      curl -O ${DATA_SAMPLE}
    fi
  else
    echo "    /!\ error: unable to find 'wget' or 'curl'"
    exit 1
  fi
fi

################################################################################
# we launch minia; note that we use only one thread (no real time issues with
# potential different results)
################################################################################
$bindir/minia -nb-cores 1 -in ERR039477.fastq.gz

################################################################################
# we check the result
################################################################################
MD5SUM_PATH=`which md5sum`
if [ ! -z "$MD5SUM_PATH" ] ; then
  # Linux
  md5sum ERR039477.fastq.contigs.fa > ERR039477.check
else
  # OSX: 'md5 -r' is equivalent to Linux md5sum
  md5 -r ERR039477.fastq.contigs.fa > ERR039477.check
fi

# 'diff' cannot be used: md5sum produces a single space between value and file
# name whereas md5 produces two spaces...
#   diff ./ERR039477.md5 ./ERR039477.check

REF_CHKSUM=`cut -d ' ' -f 1 ERR039477.md5`
CHKSUM=`cut -d ' ' -f 1 ERR039477.check`

if [ "$REF_CHKSUM" == "$CHKSUM" ]; then
   echo "TEST OK"
else
   echo "some debug: $REF_CHKSUM $CHKSUM"
   head -n 4 ERR039477.fastq.contigs.fa
   wc -l ERR039477.fastq.contigs.fa
   echo $bindir/minia -nb-cores 1 -in ERR039477.fastq.gz
   ls -l ERR039477.fastq.gz
   ls -l $bindir/minia
   echo "TEST KO"
   exit 1
fi

################################################################################
# clean up
################################################################################
rm -f  ERR039477.fastq.contigs.fa ERR039477.fastq.gz ERR039477.fastq.h5 ERR039477.check
