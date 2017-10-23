#!/usr/bin/env python
import sys

if len(sys.argv) <= 1:
    sys.exit('Specify solver log file')

nodes = []
with open(sys.argv[1]) as f:
    line = f.readline()
    while not line=='':
        if line.startswith('Host'):
            hostname = line.split(':')[1].strip().strip('"')
            nodes.append(hostname)
        elif line.startswith('Slaves'):
            Nslaves = int(f.readline())
            f.readline()
            for i in range(Nslaves):
                line = f.readline().strip('"')
                node = line.split('.')[0]
                if node not in nodes:
                    nodes.append(node)
            break
        line = f.readline()
print ' '.join([ hostname+'*' if n==hostname else n for n in sorted(nodes)])
