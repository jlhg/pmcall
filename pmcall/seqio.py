from table import DNA_CODON_TABLE
from itertools import izip


class Seq(object):

    def __init__(self, name, seq, line):
        self.name = name
        self.seq = seq
        self.line = line


def parse_fasta(path):
    """A generator that iterates two elements (header
    and sequence) each time"""
    fi = open(path, 'r')

    header = ''
    sequence = []
    while 1:
        line = fi.readline()
        if line == '':
            if header and sequence:
                yield header, ''.join(sequence)
            break

        line = line.strip()
        if line == '':
            continue
        elif line[0] == '>':
            if header and sequence:
                yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                header = line[1:]
        else:
            sequence.append(line)

    fi.close()


def reverse_complement(seq):
    table = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    rvseq = []
    for i in seq.upper()[::-1]:
        rvseq.append(table.get(i))

    return ''.join(rvseq)


def translate(seq, frame):
    if frame < 0:
        aa = []
        for i in izip(*[iter(reverse_complement(seq)[-frame - 1:])] * 3):
            aa.append(DNA_CODON_TABLE.get(''.join(i)))
    else:
        for i in izip(*[iter(seq[frame - 1:])] * 3):
            aa.append(DNA_CODON_TABLE.get(''.join(i)))

    return ''.join(aa)
