#very basic, one seq per line only, no pysam for max portability
#doesn't care if sequences are in revcomp or forward order
import sys

fasta1 = sys.argv[1]
fasta2 = sys.argv[2]

revcomp = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
def normalize(seq):
    rev = revcomp(seq)
    return min(rev,seq)

def read_seqs(fasta):
    seqs = set()
    for  line in open(fasta):
        if line[0] == '>': continue
        seqs.add(normalize(line.strip()))
    return seqs

s1 = read_seqs(fasta1)
s2 = read_seqs(fasta2)

if s1 == s2:
    sys.exit(0)

print("NOT EQUAL: %d sequence(s) in %s not in %s" % (len(s1.difference(s2)), fasta1, fasta2))
sys.exit(1)
