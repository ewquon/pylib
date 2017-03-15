#!/usr/bin/env python
import sys
import os

dirs = os.walk('.').next()[1]

dirlist = []
numlist = []
for d in dirs:
    try: 
        step = float(d)
	numlist.append(step)
        dirlist.append(d)
    except ValueError: pass

# sort list of floats
indices = [i[0] for i in sorted(enumerate(numlist), key=lambda x:x[1])]

if len(sys.argv) > 1:
    sep = sys.argv[1]
else: sep = ' '

#print ' '.join(dirlist)
#print ' '.join([dirlist[i] for i in indices])
#print ' '.join([dirlist[i] for i in indices]).strip()
print sep.join([dirlist[i] for i in indices]).strip()
