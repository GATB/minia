rm -f buchnera.reads*.fq
~/tools/samtools-1.2/misc/wgsim  -N 150000 -r 0.0005 buchnera.fasta buchnera.reads1.fq /dev/null > /dev/null
~/tools/samtools-1.2/misc/wgsim  -N 150000 -r 0.0005 buchnera.fasta buchnera.reads2.fq /dev/null > /dev/null
#~/tools/samtools-1.2/misc/wgsim  -N 150000 -r 0.0005 buchnera.fasta buchnera.reads3.fq /dev/null > /dev/null
ls -1 buchnera.reads*.fq > buchnera
