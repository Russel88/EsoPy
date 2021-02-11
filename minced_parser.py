#!/usr/bin/env python

import sys
import re

# Define the CRISPR class
class CRISPR:
    count = 0
    def __init__(self, sequence):
        self.sequence = sequence.rstrip()
        CRISPR.count += 1
        self.crispr = 'CRISPR_{}'.format(CRISPR.count)
        self.repeats = []
        self.spacers = []
    def setPos(self, start, end):
        self.start = start.rstrip()
        self.end = end.rstrip()
    def addRepeat(self, repeat):
        self.repeats.append(repeat.rstrip())
    def addSpacer(self, spacer):
        self.spacers.append(spacer.rstrip())
    def getConsensus(self):
        self.cons = max(set(self.repeats), key = self.repeats.count) 

# Read file line by line
file = open(sys.argv[1], 'r')

crisprs = []
for ll in file:
    # Record sequence accession
    if ll.startswith('Sequence'):
        sequence_current = re.sub('\' \(.*', '', re.sub('Sequence \'', '', ll))
    # Create instance of CRISPR and add positions
    if ll.startswith('CRISPR'):
        crisp_tmp = CRISPR(sequence_current)
        pos = re.sub('.*Range: ', '', ll)
        start = re.sub(' - .*', '', pos)
        end = re.sub('.* - ', '', pos)
        crisp_tmp.setPos(start, end)
    # Add Repeats and Spacers to the current instance
    if ll[:1].isdigit():
        lll = ll.split()
        if len(lll) == 7:
            crisp_tmp.addRepeat(lll[1])
            crisp_tmp.addSpacer(lll[2])
        if len(lll) == 2:
            crisp_tmp.addRepeat(lll[1])
    # Save the instance
    if ll.startswith('Repeats'):
        crisp_tmp.getConsensus()
        crisprs.append(crisp_tmp)

file.close()

# Output
# Fasta output
def getSeq(what):
    for crisp in crisprs:
        n = 0
        for sq in getattr(crisp, what):
            n += 1
            print('>' + crisp.crispr + '@' + crisp.sequence + '_' + str(n) + '\n' + sq)

# Tab output
def getSeqT(what):
    for crisp in crisprs:
        for sq in getattr(crisp, what):
            print(crisp.crispr + '@' + crisp.sequence + '\t' + sq)

# Tab output trimmed for flanking repeats
def getSeqTT():
    for crisp in crisprs:
        for sq in crisp.repeats[1:len(crisp.repeats)-1]:
            print(crisp.crispr + '@' + crisp.sequence + '\t' + sq)


if sys.argv[2] == 'tab':
    for crisp in crisprs:
        print(crisp.crispr + '\t' + crisp.sequence + '\t' + crisp.start + '\t' + crisp.end + '\t' + crisp.cons)
if sys.argv[2] == 'repeats':
    getSeq('repeats')
if sys.argv[2] == 'spacers':
    getSeq('spacers')
if sys.argv[2] == 'trepeats':
    getSeqT('repeats')
if sys.argv[2] == 'tspacers':
    getSeqT('spacers')
if sys.argv[2] == 'ttrepeats':
    getSeqTT()

