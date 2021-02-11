#!/usr/bin/env python3

'''
Print consensus sequence from a multiple alignment
First argument is a alignment fasta
Second argument is either a float determining the proportion a base should
be in to be called, or it is 'r' which will choose the one with highest proportion
and if there are ties it will pick a random one among the ties
'''

import sys
import numpy as np

from Bio import AlignIO
from Bio.Align import AlignInfo

ali = AlignIO.read(sys.argv[1], 'fasta')

if sys.argv[2] == 'r':

    def randargmax(b, **kw):
      return np.argmax(np.random.random(b.shape) * (b==b.max()), **kw)

    arr = np.transpose(np.asarray([list(str(x.seq).upper()) for x in ali]))
    cc = [np.unique(x, return_counts=True) for x in arr]
    
    gg = 0
    while gg < 10:
        cons = ''.join([x[0][randargmax(x[1])] for x in cc])
        stripped = cons.replace('-', '')
        gg += 1
        if '-' not in stripped:
            break

    print(stripped)
else:
    print(AlignInfo.SummaryInfo(ali).dumb_consensus(float(sys.argv[2])).upper())

