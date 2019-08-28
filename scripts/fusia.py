#!/usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


def main(options):
    
    # Format input sequence to search both forward and reverse complements
    search_seq = (options.s).upper()
    search_seq_rc = str(Seq(search_seq).reverse_complement())

    # Prepare out_file handle for recording search results
    try:
        assert options.o.endswith('.csv')
    except:
        raise ValueError('Output file must be .CSV')
    out_handle = open(options.o, 'w')
    out_handle.write('# search_seq: %s\n' % search_seq)
    out_handle.write('# search_seq_rc: %s\n' % search_seq_rc)
    out_handle.write('file_name,record_id,record_description\n')

    # Launch search by iterating through available fasta files
    for f in os.listdir(options.i):
        if f.endswith('.fa') or f.endswith('.fasta') or \
                                f.endswith('.fna'):

            file_name = '%s/%s' % (options.i, f)
            
            records = SeqIO.parse(file_name, "fasta")

            for record in records:
                if (search_seq in record.seq) or (search_seq_rc in record.seq):
                    out_string = ','.join([file_name, record.id, 
                                           record.description])
                    out_handle.write('%s\n' % out_string)

    out_handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", dest="i", required=True, type=str,
                        help="Input directory of assemblies in fasta format.",
                        metavar="INPUT")
    parser.add_argument("-s", "--sequence", dest="s", required=True, type=str,
                        help="Sequence to check for uniqueness within a set \
                        of assemblies.")
    parser.add_argument("-o", "--out_file", dest="o", required=True, type=str,
                        help="Output file summarizing locations of the \
                        substring.", metavar="OUTPUT")

    options = parser.parse_args()
    main(options)
