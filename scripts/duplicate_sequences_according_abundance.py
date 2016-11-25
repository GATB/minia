import sys
fasta = open(sys.argv[1])
id=""
for line in fasta:
    if line.startswith('>'):
        ma = line.find("MA=")
        abundance = float(line[ma+3:])
        #print abundance
        id=line
    else:
        for i in xrange(int(abundance)):
            print id,
            print line,

