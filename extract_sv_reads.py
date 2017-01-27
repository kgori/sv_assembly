#!/usr/bin/env python

"""
Extract just the 'SV' reads --- and their mates --- from the reads picked out by Tigra
"""

from Bio import SeqIO
from collections import defaultdict
import re
rgx = re.compile(',SV$')



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    args = parser.parse_args()
    INFILE = args.infile
    OUTFILE = args.outfile

    seqs = [(seq.name, seq.seq) for seq in SeqIO.parse(INFILE, 'fasta')]
    d = defaultdict(list)
    for (name, seq) in seqs:
        d[rgx.sub('', name)].append((name, seq))

    with open(OUTFILE, 'w') as outfile:
        for key, seqlist in d.items():
            if any([rgx.search(name) for (name, seq) in seqlist]):
                for (_, seq) in seqlist:
                    outfile.write('>{}\n{}\n'.format(key, seq))

