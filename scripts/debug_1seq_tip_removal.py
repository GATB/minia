# take an assembly and a reference genome.
# determines whether each unitig is either fully erroneous or fully genomic, but that'd be an assembly error if it was both
import sys
fasta = open(sys.argv[1])
genome = [
"CTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGC",
"GCCGAACGGCCCTGGTCACCGTCTGGCAGGTTAAAGGAAACGGTAGCTCGTTAGCCAGAGCGTGGTCCGCCCCACAATTCCCCGCGACAG"
]
id=""
k=31
for line in fasta:
    if line.startswith('>'):
        #print abundance
        id=line
    else:
        seq = line.strip()
        is_error = False
        is_genome = False
        for i in xrange(len(seq)-k+1):
            kmer = seq[i:i+k]
            if kmer not in genome[0] and kmer not in genome[1]:
                is_error = True
            else:
                is_genome = True

        if is_error and is_genome:
            print "mix of erroneous kmers and true genomic kmers inside the same assembled contig!!"
            print id
            print seq
            exit(1)
        else:
            print id,is_error,is_genome
