#!/usr/bin/env python

import subprocess
import pandas as pd


def parse_xmfa(out_xmfa, out_csv):
    """Convert a progressiveMauve xmfa file into a parsable .CSV.

    Args:
        out_xmfa (str): Path to input file with xmfa subalignments.
        out_csv (str): Desired output .CSV file parsed from .XMFA.
    """

    assert out_csv.endswith('.csv')
    xmfa_handle = open(out_xmfa, 'r')

    df = pd.DataFrame(columns = ['subalignment', 'file', 'seq_num', 
                                 'seq_strand', 'seq_start', 'seq_end', 
                                 'seq_aln'])

    curr_subalignment = 0
    curr_seq = ''

    for line in xmfa_handle:
        if line.startswith("#"):
            continue

        elif line.startswith("> "):

            # Save information for current aligned sequence
            if curr_seq != '':
                curr_row = pd.Series([curr_subalignment, file_id, seq_num, 
                                      seq_strand, seq_start, seq_end, 
                                      curr_seq], index=df.columns)
                df = df.append(curr_row, ignore_index=True)

            # Parse information for next aligned sequence
            curr_seq = ''
            seq_id = line.strip().replace('> ', '')
            seq_coords = seq_id.split()[0]
            seq_strand = seq_id.split()[1]
            file_id = seq_id.split()[2]

            seq_num = seq_coords.split(':')[0]
            seq_start = int(seq_coords.split(':')[1].split('-')[0])
            seq_end = int(seq_coords.split(':')[1].split('-')[1])
            
            assert seq_start <= seq_end

        elif line.startswith("="):

            # Save information for current aligned sequence
            if curr_seq != '':
                curr_row = pd.Series([curr_subalignment, file_id, seq_num, 
                                      seq_strand, seq_start, seq_end, 
                                      curr_seq], index=df.columns)
                df = df.append(curr_row, ignore_index=True)

            # Iterate to next set of subalignments
            if curr_subalignment % 100 == 0:
                print("Recorded subalignment %s" % curr_subalignment)
            curr_seq = ''
            curr_subalignment += 1

        else:
            curr_seq += line.strip()

    df.to_csv(out_csv)

def launch_progressivemauve(out_xmfa, out_tree, out_backbone, gbk_files):
    """Launch progressiveMauve for multiple genome alignment.

    Args:
        out_xmfa (str): Output filename for .xmfa alignment file.
        out_tree (str): Output filename for .tree alignment guide tree.
        out_backbone (str): Output filename for .backbone file.
        gbk_files (list(str)): A list of all input .gbk files to be aligned.

    Todo:
        * Save printed shell output into a log file
    """

    cmd = 'progressiveMauve'
    cmd += ' --output=%s' % out_xmfa
    cmd += ' --output-guide-tree=%s' % out_tree
    cmd += ' --backbone-output=%s' % out_backbone
    cmd += ' ' + ' '.join(gbk_files)

    print(cmd)
    subprocess.check_call(cmd, shell=True)
