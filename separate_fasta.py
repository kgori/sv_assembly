#!/usr/bin/env python

"""
Script to split the contigs fasta file produced by tigra into separate files for each SV
"""

from Bio import SeqIO
import os, sys

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    parser.add_argument('-o', '--outdir', type=str, default='.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Inform the user of file closure')
    return parser.parse_args()


class FileHandleDictionary(dict):
    """
    Like a regular dict, but will close any open files it contains
    when it's finished with.

    Usage:
        with FileHandleDictionary() as fhd:
            make some files...

        ... now all files are closed

    """
    def __init__(self, inform_of_closures=False):
        super(FileHandleDictionary, self).__init__()
        self._inform = inform_of_closures

    def __enter__(self):
        return self

    def __exit__(self, *args):
        for file_handle in self.values():
            if isinstance(file_handle, file) and not file_handle.closed:
                if self._inform:
                    sys.stdout.write('INFO: Closing {}\n'.format(file_handle.name))
                file_handle.close()


if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.filename):
        sys.stderr.write("ERROR: {} not found\nExiting - no results have been written.\n".format(args.filename))
        sys.exit()

    if args.outdir is not None and not os.path.exists(args.outdir):
        sys.stderr.write("ERROR: {} not found\nExiting - no results have been written.\n".format(args.outdir))
        sys.exit()

    seqs = SeqIO.parse(args.filename, "fasta")
    with FileHandleDictionary(args.verbose) as file_handles:
        for seq in seqs:
            # tigra bundles info fields together separated by '^'
            # all but the last field is SV specific - i.e.
            # all contigs built for the same SV will have
            # the same leading fields. The last field is contig specific.
            id_ = seq.description.rpartition("^")
            if not id_[0] in file_handles:
                file_handles[id_[0]] = open(os.path.join(args.outdir, '{}.fa'.format(id_[0])), 'w')
            handle = file_handles[id_[0]]
            SeqIO.write(seq, handle, 'fasta')
