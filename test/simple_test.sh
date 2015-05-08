
function test()
{
    file=$1
    k=$2
    traversal=$3
    ../build/minia -in "$file".fa -kmer-size "$k" -abundance-min 1 -traversal $traversal
    python compare_fasta.py "$file".fa "$file".contigs.fa
}

test X 5 unitig

test tip 21 contig

test bubble 21 contig 
