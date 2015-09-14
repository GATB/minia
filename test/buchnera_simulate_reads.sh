rm -f buchnera.reads*.fq
N=15000
for i in $(seq 1 15)
do
    ~/tools/samtools-1.2/misc/wgsim  -N $N -r 0.008 buchnera.fasta buchnera.reads$i.fq /dev/null > /dev/null
done
ls -1 buchnera.reads*.fq > buchnera
