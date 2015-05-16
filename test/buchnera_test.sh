file=buchnera
k=22

../build/minia -in buchnera -abundance-min 2 -kmer-size $k

echo "tigops, bandage:"
echo "../../tigops/build/tigops fasta2fastg -tigs "$file".contigs.fa -out "$file".contigs.fastg -kmer-size $k -rename"

echo "~/tools/bandage/Bandage/Bandage $file.contigs.fastg"

