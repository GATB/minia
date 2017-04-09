import sys
file=sys.argv[1]
for line in open(file):
     line=line.strip()
     if not line.startswith(">"):
        print(line.strip())
        continue
     if "L:" in line:
         s = line.split("L:")
     else:
         s = [line.strip()]
     print(s[0])

