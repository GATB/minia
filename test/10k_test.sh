file=read50x_ref10K_e001
k=21
t=2
traversal=contig

../build/minia -in "$file".fa -kmer-size "$k" -no-length-cutoff -abundance-min $t -traversal $traversal

echo "legacy: ~/gatb-pipeline/minia/minia -in "$file".fa -kmer-size "$k" -no-length-cutoff -abundance-min $t -traversal $traversal"

# running tigops
echo "tigops: ../../tigops/build/tigops fasta2fastg -tigs "$file".contigs.fa -out "$file".contigs.fastg -kmer-size $k"

echo "Bandage: ~/tools/bandage/Bandage/Bandage $file.contigs.fastg"

