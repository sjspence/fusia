#!/usr/bin/env python

import subprocess


def launch_progressivemauve(out_xmfa, out_tree, out_backbone, gbk_files):
    """Launch progressiveMauve for multiple genome alignment.

    Args:
        out_xmfa (str): Output filename for .xmfa alignment file.
        out_tree (str): Output filename for .tree alignment guide tree.
        out_backbone (str): Output filename for .backbone file.
        gbk_files (list(str)): A list of all input .gbk files to be aligned.
    """

    cmd = 'progressiveMauve'
    cmd += ' --output=%s' % out_xmfa
    cmd += ' --output-guide-tree=%s' % out_tree
    cmd += ' --backbone-output=%s' % out_backbone
    cmd += ' ' + ' '.join(gbk_files)

    print(cmd)
    subprocess.check_call(cmd, shell=True)
