#!/usr/bin/env python
import sys
max_mem = 0
for line in sys.stdin: 
    s = line.replace(']','').split()
    if 'MB' not in s:
        continue
    i=s.index('MB')
    if s[i-1].isdigit(): 
        mem=int(s[i-1])
        max_mem=max(mem,max_mem)

print max_mem
