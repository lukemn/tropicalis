#!/usr/bin/env python

import sys, signal
from sys import stderr, stdout, stdin
signal.signal(signal.SIGPIPE, signal.SIG_DFL) #silent on head

# quick strip of genotype calls per sample

samples = 0
for l in stdin:
    col = l.strip().split()
    if col[0] == '#CHROM':
        samples = col[9:]
    if samples:
        gts = col[:2] + col[3:7] + [i.split(':')[0] for i in col[9:]]
        stdout.write('\t'.join(gts)+'\n') 

    
