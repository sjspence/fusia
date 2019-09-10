#!/usr/bin/python

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