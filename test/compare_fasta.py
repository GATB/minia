#very basic, one seq per line only, no pysam for max portability
import sys

fasta1 = sys.argv[1]
fasta2 = sys.argv[2]

def read_seqs(fasta):
    seqs = set()
    for  line in open(fasta):
        if line[0] == '>': continue
        seqs.add(line.strip())
    return seqs

s1 = read_seqs(fasta1)
s2 = read_seqs(fasta2)
if s1 == s2:
    exit("OK")

exit("NOT EQUAL: %d sequence(s) in %s not in %s" % (len(s1.difference(s2)), fasta1, fasta2))
