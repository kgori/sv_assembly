import skbio
from skbio import DNA
import os
import sys

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str)
    parser.add_argument('outfile', type=str)
    parser.add_argument('translation_file', type=str)
    args = parser.parse_args()
    if not os.path.isfile(args.infile):
        sys.stderr.write("File not found: {}".format(args.infile))
        sys.exit(1)

    fasta = [DNA(str(seq).upper(), seq.metadata) for seq in skbio.read(args.infile, format='fasta')]
    names = []
    new_names = []

    with open(args.translation_file, 'w') as fl:
        for (i, seq) in enumerate(fasta, start=1):
            name = seq.metadata['id'] + seq.metadata['description']
            short_name = '^'.join(name.split('^')[:4]) + '^{}'.format(i)
            fl.write('{}\t{}\n'.format(short_name, name))
            seq.metadata['id'] = short_name
            del seq.metadata['description']

    skbio.io.write((seq for seq in fasta), format='fasta', into=args.outfile)

