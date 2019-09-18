#!/usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import config
import util_funcs
import alignment_funcs

ACCEPTED_IUPAC = set('GATCN')
EXTENDED_IUPAC = set('URYSWKMBDHVN.-')


def launch_substring(options):

    # Format input sequence to search both forward and reverse complements
    search_seq = util_funcs.check_seq(options.s)
    search_seq_rc = str(Seq(search_seq).reverse_complement())

    # Prepare out_file handle for recording search results
    try:
        assert options.o.endswith('.csv')
    except ValueError:
        print('Output file must be .CSV')

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
                except ValueError:
                    print('Unexpected character\n%s' % record_seq)

                if (search_seq in record_seq) or (search_seq_rc in record_seq):

                    count = record_seq.count(search_seq)
                    count_rc = record_seq.count(search_seq_rc)
                    count_total = count + count_rc

                    out_string = ','.join([file_name, record.id,
                                           record.description,
                                           str(count_total)])
                    out_handle.write('%s\n' % out_string)

    out_handle.close()


def launch_unique(options):

    print("##", "-" * 50)
    print("## Launching \'fusia.py unique\' to identify unique substrings")
    print("## among genome assemblies")
    print("##", "-" * 50)

    # Prepare output directories
    if not os.path.exists(options.o):
        os.makedirs(options.o)

    # Prepare output files
    out_xmfa = "%s/output.xmfa" % options.o
    out_tree = "%s/output.tree" % options.o
    out_backbone = "%s/output.backbone" % options.o
    out_subaln_csv = "%s/output_subalignments.csv" % options.o

    print("##")
    print("## Iterate through input fasta files and convert to .gbk")
    print("## -----------------------")

    gbk_files = []
    for f in os.listdir(options.i):
        if f.endswith('.fa') or f.endswith('.fasta') or \
                                f.endswith('.fna'):

            input_fasta = "%s/%s" % (options.i, f)
            gbk_filename = '.'.join(f.split('.')[0:-1]) + '.gbk'
            output_gbk = "%s/%s" % (options.o, gbk_filename)

            assert output_gbk not in gbk_files
            gbk_files.append(output_gbk)

            if not os.path.exists(output_gbk):
                print("Converting %s to .gbk" % f)
                input_handle = open(input_fasta, 'r')
                sequences = list(SeqIO.parse(input_handle, "fasta"))
                for s in sequences:
                    full_description = s.id
                    s.name = '_'.join(full_description.split('_')[0:2])
                    s.seq.alphabet = generic_dna

                output_handle = open(output_gbk, 'w')
                count = SeqIO.write(sequences, output_handle, "genbank")

                input_handle.close()
                output_handle.close()
                print("Converted %i records" % count)
            else:
                print("Detected existing gbk file: %s" % output_gbk)

    print("##")
    print("## Launch progressiveMauve alignment on the list of .gbk files")
    print("## -----------------------")

    if not os.path.exists(out_xmfa):
        alignment_funcs.launch_progressivemauve(out_xmfa, out_tree,
                                                out_backbone, gbk_files)
    else:
        print("Detected existing .xmfa file: %s" % out_xmfa)

    print("##")
    print("## Parse .xmfa output to identify regions unique to one genome")
    print("## -----------------------")

    if not os.path.exists(out_subaln_csv):
        alignment_funcs.parse_xmfa(out_xmfa, out_subaln_csv)

    else:
        print("Detected existing subalignments file: %s" % out_subaln_csv)

    alignment_funcs.find_unique_regions(out_subaln_csv)


def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest='command')

    substring = subparsers.add_parser('substring', help="Find occurences of \
                                      input substring.")
    substring.add_argument("-i", "--in_dir", dest="i", required=True,
                           type=str, metavar="INPUT",
                           help="Input directory of assemblies in fasta \
                           format.")
    substring.add_argument("-s", "--sequence", dest="s", required=True,
                           type=str,
                           help="Sequence to check for uniqueness within a \
                           set of assemblies.")
    substring.add_argument("-o", "--out_file", dest="o", required=True,
                           type=str, metavar="OUTPUT",
                           help="Output file summarizing locations of the \
                           substring.")

    unique = subparsers.add_parser('unique', help="Find unique region across \
                                   input fasta files.")
    unique.add_argument("-i", "--in_dir", dest="i", required=True,
                        type=str, metavar="INPUT",
                        help="Input directory of assemblies in fasta \
                        format.")
    unique.add_argument("-o", "--out_dir", dest="o", required=True,
                        type=str, metavar="OUTPUT",
                        help="Output directory of assemblies in gbk \
                        format.")

    options = parser.parse_args()

    if options.command == 'substring':
        launch_substring(options)

    elif options.command == 'unique':
        launch_unique(options)

    else:
        raise ValueError("Unrecognized command: %s" % options.command)


if __name__ == '__main__':
    main()
