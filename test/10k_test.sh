#!/bin/bash
file=read50x_ref10K_e001
k=21
t=2
traversal=contig

set -x
../build/minia -in "$file".fa -kmer-size "$k" -abundance-min $t -traversal $traversal
set +x

echo "legacy: ~/gatb-pipeline/minia/minia-legacy -in "$file".fa -kmer-size "$k" -no-length-cutoff -abundance-min $t -traversal $traversal"

# running tigops
echo "tigops, bandage:"
echo "../../tigops/build/tigops fasta2fastg -tigs "$file".contigs.fa -out "$file".contigs.fastg -kmer-size $k -rename"

echo "~/tools/bandage/Bandage/Bandage $file.contigs.fastg"

