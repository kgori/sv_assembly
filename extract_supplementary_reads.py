#!/usr/bin/env python
import argparse
import pysam
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help='Input BAM file to extract reads from')
    parser.add_argument('outfile', type=str, help='Full path and filename for new BAM file to write reads into')
    return parser.parse_args()

def args_ok(arguments):
    """
    Return True if arguments are OK, otherwise return False
    """
    if not os.path.exists(arguments.infile):
        sys.stderr.write("ERROR: Can't find input file {}\n".format(arguments.infile))
        return False
    dirname = os.path.dirname(arguments.outfile)
    if dirname == '':
        return True
    if not os.path.isdir(dirname):
        sys.stderr.write("ERROR: Directory {} is not writable\n".format(dirname))
        return False
    return True

if __name__ == '__main__':
    args = parse_args()
    if args_ok(args):
        with pysam.Samfile(args.infile, 'rb') as inbam:
            with pysam.Samfile(args.outfile, 'wb', template=inbam) as outbam:
                for read in inbam:
                    if read.is_supplementary:
                        outbam.write(read)
        pysam.index(args.outfile)
