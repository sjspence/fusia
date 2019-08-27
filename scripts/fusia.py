#!/usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO

def main(options):
    
    test_seq = 'CATCGCCGCGCAAACGTTTAAGCAG'

    for f in os.listdir(options.i):
        if f.endswith('.fa') or f.endswith('.fasta') or \
                                f.endswith('.fsa'):

            file_name = '%s/%s' % (options.i, f)
            
            # START HERE
            SeqIO.read()

            for line in file_handle:
                if test_seq in line:
                    print(line)

            file_handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", dest="i", required=True,
                        help="Input directory of assemblies in fasta format.",
                        metavar="INPUT")
    parser.add_argument("-s", "--sequence", dest="s",
                        help="Sequence to check for uniqueness within a set \
                        of assemblies.")
    parser.add_argument("-o", "--output", dest="o",
                        help="Output TBD.",
                        metavar="OUTPUT")

    options = parser.parse_args()
    main(options)
