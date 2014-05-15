import re
from itertools import izip
from seqio import reverse_complement

CLU_CONSV_NAME = '_consv_'


class Clustal(object):

    def __init__(self):
        self.title = ''
        self.records = {}


class ClustalParser(object):

    def __init__(self):
        self._pattern = re.compile('(\S+)(\s+)(\S+)')

    def parse(self, fclustal):
        clustal = Clustal()
        slen = 0
        with open(fclustal, 'r') as fi:
            clustal.title = fi.readline().rstrip('\n')
            for line in fi:
                if line.startswith(' '):
                    line = line.rstrip('\n')
                    if CLU_CONSV_NAME in clustal.records:
                        clustal.records.get(CLU_CONSV_NAME).append(line[-slen:])
                    else:
                        clustal.records.update({CLU_CONSV_NAME: [line[-slen:]]})
                else:
                    match = self._pattern.match(line)
                    if match:
                        name = match.group(1)
                        if name == CLU_CONSV_NAME:
                            raise Exception('{0} is the built-in variable for Clustal.'.format(CLU_CONSV_NAME))
                        seq = match.group(3)
                        slen = len(seq)
                        if name in clustal.records:
                            clustal.records.get(name).append(seq)
                        else:
                            clustal.records.update({name: [seq]})

        for i, j in clustal.records.items():

            clustal.records.update({i: ''.join(j)})

        return clustal


def aa_to_nt(fclustal_aa, path_out, ntseq, chunk_size=60):
    clustalparser = ClustalParser()
    clustal_aa = clustalparser.parse(fclustal_aa)
    clustal_nt = Clustal()
    clustal_nt.title = 'CLUSTAL format alignment converted from aa clustal'
    namepattern = re.compile('(\S+)\((ss|rs|rc),([+-]*\d+)\)')
    namelength = 0

    clustal_aa.records.pop(CLU_CONSV_NAME)
    for name, aaseq in clustal_aa.records.items():
        if len(name) > namelength:
            namelength = len(name)

        match = namepattern.match(name)
        seq = ntseq.get(match.group(1))
        frame = int(match.group(3))
        orf = []
        if frame < 0:
            for i in izip(*[iter(reverse_complement(seq)[-frame - 1:])] * 3):
                orf.append(''.join(i))
        else:
            for i in izip(*[iter(seq[frame - 1:])] * 3):
                orf.append(''.join(i))

        ncodon = 0
        nt = []
        for aa in aaseq:
            if aa == '-':
                nt.append('---')
            else:
                nt.append(orf[ncodon])
                ncodon += 1
        clustal_nt.records.update({name: chunk_string(''.join(nt), chunk_size)})

    consv = []

    for i in izip(*[''.join(x) for x in clustal_nt.records.values()]):
        if all(x == i[0] for x in i):
            consv.append('*')
        else:
            consv.append(' ')

    consvseq = chunk_string(''.join(consv), chunk_size)

    nchunk = len(consvseq)
    namelength += 4
    with open(path_out, 'w') as fo:
        fo.write(clustal_nt.title)
        fo.write('\n\n\n')
        fo.flush()
        for i in range(nchunk):
            for name, chunkseq in clustal_nt.records.items():
                fo.write('{0: <{1}}{2}\n'.format(name, namelength, chunkseq[i]))
            fo.write('{0: <{1}}{2}\n\n'.format('', namelength, consvseq[i]))
            fo.flush()


def chunk_string(text, length):
    return [text[0 + x:length + x] for x in range(0, len(text), length)]
