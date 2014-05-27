#coding:utf-8

import os, sys, re

filename = sys.argv[1]
counts = {}
with open(filename) as fi:
    for line in fi:
        if line.startswith('#'): continue
        items = line.split('\t', 1)
        if items[0] in counts:
            counts[items[0]] += 1
        else:
            counts[items[0]] = 1
            pass
for chrm in sorted(counts.keys()):
    print('{}\t{}'.format(chrm, counts[chrm]))
