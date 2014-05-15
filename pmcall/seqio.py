from table import DNA_CODON_TABLE
from itertools import izip


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
        if i in table:
            rvseq.append(table.get(i))
        else:
            rvseq.append('N')

    return ''.join(rvseq)


def translate(seq, frame):
    seq = seq.upper()
    aa = []

    if frame < 0:
        for i in izip(*[iter(reverse_complement(seq)[-frame - 1:])] * 3):
            codon = ''.join(i)
            if codon in DNA_CODON_TABLE:
                aa.append(DNA_CODON_TABLE.get(''.join(i)))
            else:
                aa.append('X')
    else:
        for i in izip(*[iter(seq[frame - 1:])] * 3):
            codon = ''.join(i)
            if codon in DNA_CODON_TABLE:
                aa.append(DNA_CODON_TABLE.get(''.join(i)))
            else:
                aa.append('X')

    return ''.join(aa)
