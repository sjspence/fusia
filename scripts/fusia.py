#!/usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq

ACCEPTED_IUPAC = set('GATCN')
EXTENDED_IUPAC = set('URYSWKMBDHVN.-')

def check_seq(sequence):
    """
    Ensure sequence is uppercase for comparison, and with only canonical
    IUPAC DNA characters
    """

    sequence = sequence.upper()

    try:
        assert set(sequence).issubset(ACCEPTED_IUPAC)
    except:
        raise ValueError('Unexpected character\n%s' % sequence)

    return sequence


def main(options):
    
    # Format input sequence to search both forward and reverse complements
    search_seq = check_seq(options.s)
    search_seq_rc = str(Seq(search_seq).reverse_complement())

    # Prepare out_file handle for recording search results
    try:
        assert options.o.endswith('.csv')
    except:
        raise ValueError('Output file must be .CSV')

    out_handle = open(options.o, 'w')
    out_handle.write('# search_seq: %s\n' % search_seq)
    out_handle.write('# search_seq_rc: %s\n' % search_seq_rc)
    out_handle.write('file_name,record_id,record_description, count\n')

    # Launch search by iterating through available fasta files
    for f in os.listdir(options.i):
        if f.endswith('.fa') or f.endswith('.fasta') or \
                                f.endswith('.fna'):

            print("Searching for seq in %s" % f)
            file_name = '%s/%s' % (options.i, f)
            records = SeqIO.parse(file_name, "fasta")

            for record in records:

                # Ensure sequence is uppercase for comparison, only canonical
                # IUPAC DNA characters
                record_seq = str(record.seq).upper()
                try:
                    assert set(record_seq).issubset(ACCEPTED_IUPAC)
                except:
                    raise ValueError('Unexpected character\n%s' % record_seq)

                if (search_seq in record_seq) or (search_seq_rc in record_seq):

                    count = record_seq.count(search_seq)
                    count_rc = record_seq.count(search_seq_rc)
                    count_total = count + count_rc

                    out_string = ','.join([file_name, record.id, 
                                           record.description, 
                                           str(count_total)])
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
