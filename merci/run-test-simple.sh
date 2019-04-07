cd test-simple; 
../../buildk32/merci -kmer-size 11 -reads reads.fa -assembly assembly.fa -verbose  1 ;
cd ..

echo "wc, should have 4 lines:"
wc -l test-simple/assembly.fa
echo "wc, should have 2 lines:"
wc -l test-simple/assembly.fa.merci
