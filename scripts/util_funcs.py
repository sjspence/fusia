#!/usr/bin/env python

ACCEPTED_IUPAC = set('GATCN')
EXTENDED_IUPAC = set('URYSWKMBDHVN.-')


def check_seq(sequence):
    """Check input DNA sequence and modify for uppercase and IUPAC characters.

    Args:
        sequence (str): Input DNA sequence.
    """

    sequence = sequence.upper()

    try:
        assert set(sequence).issubset(ACCEPTED_IUPAC)
    except:
        raise ValueError('Unexpected character\n%s' % sequence)

    return sequence
