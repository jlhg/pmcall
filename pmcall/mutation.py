import re
from itertools import izip
from clustal import ClustalParser, CLU_CONSV_NAME


class MutationParser(object):

    def __init__(self, blocklen=10, perfect_match_percent=0.9, sidelen=5, side_match_percent=0.8):
        assert isinstance(blocklen, int)
        assert isinstance(perfect_match_percent, float)
        assert isinstance(sidelen, int)
        assert isinstance(side_match_percent, float)

        self.blocklen = blocklen
        self.perfect_match_percent = perfect_match_percent
        self.sidelen = sidelen
        self.side_match_percent = side_match_percent
        self._unknown_code = {
            'nt': 'N',
            'aa': 'X',
        }
        self._headerpattern = re.compile('\S+\((ss|rs|rc)(,[+-]*\d+)*\)')

    def parse(self, fclustal, seqtype):
        """
        fclustal: clustal file
        seqtype: 'nt' for nucleotide; 'aa' for amino acid
        """
        if seqtype not in ('nt', 'aa'):
            raise Exception("Argument 'seqtype' must be 'nt' or 'aa'.")

        clustal_mutation_profile = ClustalMutationProfile()
        clustalparser = ClustalParser()
        clustal = clustalparser.parse(fclustal)

        clustalseq_consv = clustal.records.pop(CLU_CONSV_NAME)
        clustalseqs_ss = []
        clustalseqs_rs = []
        clustalseqs_rc = []
        for header, seq in clustal.records.items():
            match = self._headerpattern.match(header)
            if match:
                line = match.group(1)
            else:
                raise Exception('Unkonwn line: {0}'.format(header))
            if line == 'ss':
                clustalseqs_ss.append(seq)
            elif line == 'rs':
                clustalseqs_rs.append(seq)
            elif line == 'rc':
                clustalseqs_rc.append(seq)
            else:
                raise Exception('Unkonwn line: {0}'.format(header))

        if clustalseqs_rc:
            apply_rc = True
        else:
            apply_rc = False

        nogap_indexes = []
        for i, j in enumerate(izip(*clustal.records.values())):
            if all(x != '-' for x in j):
                nogap_indexes.append(i)

        block = group_continuous_number(nogap_indexes)
        mutation_indexes = []
        for i in block:
            clustal_mutation_profile.nblock += 1
            blocklen = i[1] - i[0] + 1
            nast = 0
            for i in range(i[0], i[1] + 1):
                if clustalseq_consv[i] == '*':
                    nast += 1
            perfect_match_percent = nast / float(blocklen)

            if self.blocklen > blocklen or self.perfect_match_percent > perfect_match_percent:
                continue

            block_profile = BlockProfile()
            block_profile.start = i[0] + 1
            block_profile.end = i[1] + 1
            block_profile.length = i[1] - i[0] + 1
            block_profile.perfect_match_percent = round(perfect_match_percent, 2)
            clustal_mutation_profile.block_profile.append(block_profile)

            for j in range(i[0], i[1] + 1):
                bases_ss = [x[j] for x in clustalseqs_ss]
                bases_rs = [x[j] for x in clustalseqs_rs]
                bases_rc = [x[j] for x in clustalseqs_rc]
                bases_ss_set = set(bases_ss)
                bases_rs_set = set(bases_rs)
                bases_rc_set = set(bases_rc)

                if len(bases_ss_set) != 1 or len(bases_rs_set) != 1:
                    continue

                if apply_rc and len(bases_rc_set) != 1:
                    continue

                if self._unknown_code.get(seqtype) in bases_ss_set & bases_rs_set:
                    continue

                if apply_rc and self._unknown_code.get(seqtype) in bases_rc_set:
                    continue

                retrieve_mutation = False
                if self.sidecheck(seq=clustalseq_consv, position=j, symbol='*'):
                    if apply_rc:
                        if bases_ss_set == bases_rc_set and bases_ss_set != bases_rs_set:
                            retrieve_mutation = True
                    elif bases_ss_set != bases_rs_set:
                        retrieve_mutation = True

                if retrieve_mutation:
                    clustal_mutation_profile.nmutation += 1
                    mutation_indexes.append(j)

        if mutation_indexes:
            for header, seq in clustal.records.items():
                clustal_mutation_profile.mutation_profiles.update({header: []})
                for i in mutation_indexes:
                    mutation_profile = ''.join([
                        bases_ss[0],
                        len(seq[0: i + 1].replace('-', '')),
                        bases_rs[0],
                    ])
                    clustal_mutation_profile.mutation_profiles.get(header).append(mutation_profile)

        return clustal_mutation_profile

    def sidecheck(self, seq, position, symbol='*'):
        if position < self.sidelen or position + self.sidelen > len(seq) - 1:
            # Mutation position is on sides, return False.
            return False

        # Check two sides
        nsymbol = 0
        for i in range(position + 1, position + 1 + self.sidelen):
            if seq[i] == symbol:
                nsymbol += 1

        if nsymbol < self.sidelen * self.side_match_percent:
            return False

        nsymbol = 0
        for i in range(position - self.sidelen, position):
            if seq[i] == symbol:
                nsymbol += 1

        if nsymbol >= self.sidelen * self.side_match_percent:
            return True
        else:
            return False


def group_continuous_number(numbers):
    if not numbers:
        yield 0, 0
    else:
        # Yield a numbers of tuple of continuous numbers
        first = last = numbers[0]
        for n in numbers[1:]:
            if n - 1 == last:
                # Part of the group
                last = n
            else:
                # Not part of the group
                yield first, last
                first = last = n
        # Yield the last group
        yield first, last


class BlockProfile(object):

    def __init__(self):
        self.start = None
        self.end = None
        self.length = None
        self.perfect_match_percent = None


class ClustalMutationProfile(object):

    def __init__(self):
        self.nblock = 0
        self.nmutation = 0
        self.block_profiles = []
        self.mutation_profiles = {}
