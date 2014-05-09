class Parser(object):

    def __init__(self, block_len=10, perfect_matches=0.9, neighbor_len=5, neighbor_matches=0.8):
        assert isinstance(block_len, int)
        assert isinstance(perfect_matches, float)
        assert isinstance(neighbor_len, int)
        assert isinstance(neighbor_matches, float)

        self.block_len = block_len
        self.perfect_matches = perfect_matches
        self.neighbor_len = neighbor_len
        self.neighbor_matches = neighbor_matches

    def parse(self, path_in, seqtype):
        """
        path_in: clustal file
        alignment order: ss, rs, rc
        seqtype: 'nt' for nucleotide; 'aa' for amino acid
        """
        if seqtype not in ('nt', 'aa'):
            raise Exception("Argument 'seqtype' must be 'nt' or 'aa'.")

        mutdisc = _MutDisc(block_len=self.block_len,
                           perfect_matches=self.perfect_matches,
                           neighbor_len=self.neighbor_len,
                           neighbor_matches=self.neighbor_matches)
        mutdisc.parse(path_in, seqtype)

    def get_clustal_html(self):
        if not self._parse:
            return None

        html = []
        html.append('<div class="%s">' % (self.htmltag.get('clustal')))
        html.append(self._title.get_html())
        html.append('<br><br>')

        for i in range(1, len(self._star.count)):
            html.append(self._ss.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._rs.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._rc.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._star.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append('<br>')

        html.append('</div>')

        return ''.join(html)


class _Clutitle(object):

    def __init__(self, tag='alignment'):
        self.line = None
        self.tag = tag

    def parse(self, line):
        self.line = line

    def get_html(self):
        return '<span class="%s">%s</span></br>' % (self.tag, self.line)


class _CluSequence(object):

    def __init__(self, tag_name, tag_gap, tag_seqtype, tag_var):
        import re
        self.name = None
        self.blank = None
        self.linelength = []
        self._base = []
        self.pattern = re.compile('(\S+)\((\S+)\)(\s+)(\S+)')
        self.tag_name = tag_name
        self.tag_gap = tag_gap
        self.tag_seqtype = tag_seqtype
        self.tag_var = tag_var
        self._mutposlist = []

    def parse(self, line):
        match = self.pattern.match(line)

        if match is not None:
            self.name = match.group(1)
            extra = match.group(2).split(',')
            if len(extra) == 2:
                self.line = extra[0]
                self.frame = int(extra[1])
            else:
                self.line = extra[0]
                self.frame = None
            assert self.line in {'ss', 'rs', 'rc'}

            self.blank = match.group(3)
            self.linelength.append(len(match.group(4)))
            self._base = self._base + list(match.group(4))
        else:
            raise Exception('Format error.')

    def get_sequence_length(self):
        return len(self._base)

    def get_base(self, position):
        return self._base[position]

    def get_bases(self, start_position, stop_position):
        return ''.join(self._base[start_position:stop_position + 1])

    def get_raw_position(self, position):
        return len(''.join(self._base[0:position + 1]).replace('-', ''))

    def add_mutposition(self, position):
        self._mutposlist.append(position)

    def get_html(self, start_position, end_position):
        values = []
        values.append('<span class="%s">%s</span>' % (self.tag_name, self.name))
        values.append('&nbsp;' * len(self.blank))

        for i in range(start_position, end_position):
            if i in self._mutposlist:
                values.append('<span class="%s %s">%s</span>' % (self.tag_var, self._base[i], self._base[i]))
            elif i == '-':
                values.append('<span class="%s">%s</span>' % (self.tag_gap, self._base[i]))
            else:
                values.append('<span class="%s %s">%s</span>' % (self.tag_seqtype, self._base[i], self._base[i]))

        values.append('<br>')

        return ''.join(values)


class _CluAst(object):

    def __init__(self, tag_ast, tag_dot, tag_col, tag_count):
        self.blank = None
        self._base = []
        self.tag_ast = tag_ast
        self.tag_dot = tag_dot
        self.tag_col = tag_col
        self.tag_count = tag_count
        self.count = [0]

    def parse(self, line, linelength):
        line = line.rstrip('\n')
        self.blank = line[0: len(line) - linelength]
        self._base = self._base + list(line[len(line) - linelength:len(line)])

        if self.count:
            self.count.append(self.count[-1] + linelength)
        else:
            self.count.append(linelength)

    def get_nasterisk(self, start, end):
        star_num = 0

        for i in range(start, end + 1):
            if self._base[i] == '*':
                star_num += 1

        return star_num

    def neighbor_star_check(self, position, star_check_num):
        if position < star_check_num or position + star_check_num > len(self._base) - 1:
            # Mutation position is on sides, return False.
            return False

        # Check two sides
        star_num = 0

        for i in range(position + 1, position + 1 + star_check_num):
            if self._base[i] == '*':
                star_num += 1

        if star_num < star_check_num * 0.8:
            return False

        star_num = 0

        for i in range(position - star_check_num, position):
            if self._base[i] == '*':
                star_num += 1

        if star_num >= star_check_num * 0.8:
            return True
        else:
            return False

    def get_html(self, start_position, end_position):
        values = []
        values.append('&nbsp;' * len(self.blank))

        for i in range(start_position, end_position):
            if self._base[i] == '*':
                values.append('<span class="%s">%s</span>' % (self.tag_ast, self._base[i]))
            elif self._base[i] == ':':
                values.append('<span class="%s">%s</span>' % (self.tag_col, self._base[i]))
            elif self._base[i] == '.':
                values.append('<span class="%s">%s</span>' % (self.tag_dot, self._base[i]))
            elif self._base[i] == ' ':
                values.append('&nbsp;')
            else:
                raise Exception('Extra symbol was found.')

        values.append('<span class="%s">%s</span>' % (self.tag_count, end_position))
        values.append('<br>')

        return ''.join(values)


def group_continuous_number(numbers):
    if not numbers:
        yield 0, 0
    else:
        # yield a numbers of tuple of continuous numbers
        first = last = numbers[0]
        for n in numbers[1:]:
            if n - 1 == last:
                # part of the group
                last = n
            else:
                # not part of the group
                yield first, last
                first = last = n
        # yield the last group
        yield first, last


class MutationProfile(object):

    def __init__(self):
        self.query_id = None
        self.profile = None
        self.line = None
        self.frame = None


class _MutDisc(object):

    def __init__(self, block_len=10, perfect_matches=0.9, neighbor_len=5, neighbor_matches=0.8):
        assert isinstance(block_len, int)
        assert isinstance(perfect_matches, float)
        assert isinstance(neighbor_len, int)
        assert isinstance(neighbor_matches, float)

        self.block_len = block_len
        self.perfect_matches = perfect_matches
        self.neighbor_len = neighbor_len
        self.neighbor_matches = neighbor_matches

    def parse(self, path_in, seqtype, htmltag=None):
        self.unknown_base = {
            'nt': 'N',
            'aa': 'X',
        }
        self.htmltag = {
            'clustal': 'alignment',
            'clutitle': 'title',
            'ss_name': 'ss',
            'rs_name': 'rs',
            'rc_name': 'rc',
            'ast': 'ast',
            'gap': 'gap',
            'dot': 'dot',
            'col': 'col',
            'count': 'count',
            'origin': 'o',
            'mutation': 'm',
            'seqtype': seqtype,
        }
        if htmltag:
            self.htmltag.update(htmltag)

        self._title = _Clutitle(tag=self.htmltag.get('clutitle'))
        self._seq = {}
        # self._ss = _CluSequence(tag_name=self.htmltag.get('ss_name'),
        #                         tag_gap=self.htmltag.get('gap'),
        #                         tag_seqtype=self.htmltag.get('seqtype'),
        #                         tag_var=self.htmltag.get('origin'))
        # self._rs = _CluSequence(tag_name=self.htmltag.get('rs_name'),
        #                         tag_gap=self.htmltag.get('gap'),
        #                         tag_seqtype=self.htmltag.get('seqtype'),
        #                         tag_var=self.htmltag.get('mutation'))
        # self._rc = _CluSequence(tag_name=self.htmltag.get('rc_name'),
        #                         tag_gap=self.htmltag.get('gap'),
        #                         tag_seqtype=self.htmltag.get('seqtype'),
        #                         tag_var=self.htmltag.get('origin'))
        self._ast = _CluAst(tag_ast=self.htmltag.get('ast'),
                            tag_dot=self.htmltag.get('dot'),
                            tag_col=self.htmltag.get('col'),
                            tag_count=self.htmltag.get('count'))
        self.nmut = 0
        self.position = []
        self.profile = []
        self.block = []

        with open(path_in, 'r') as fi:
            title = fi.readline().rstrip('\n')
            if title[0:7] == 'CLUSTAL':
                self._title.parse(title)
            else:
                raise Exception('Need a clustal file.')

            while 1:
                line = fi.readline()
                if line == '':
                    # end of file
                    break
                elif line == '\n':
                    nseq = 0
                    continue
                elif line.startswith(' '):
                    self._ast.parse(line, self._seq.get(1).linelength[-1])
                else:
                    nseq += 1
                    if nseq not in self._seq:
                        self._seq.update({nseq: _CluSequence()})
                    self._seq.get(nseq).parse(line)

        assert all(self._seq.get(1).get_sequence_length() == self._seq.get(x).get_sequence_length() for x in self._seq)
        nogap_indexs = []

        for i in range(self._seq.get(1).get_sequence_length()):
            # if self._ss._base[i] != '-' and self._rs._base[i] != '-' and self._rc._base[i] != '-':
            if all(x._base[i] != '-' for x in self._seq.items()):
                nogap_indexs.append(i)

        block = group_continuous_number(nogap_indexs)

        for i in block:
            block_len = i[1] - i[0] + 1
            perfect_matches = self._ast.get_nasterisk(i[0], i[1]) / float(block_len)

            if self.block_len > block_len or self.perfect_matches > perfect_matches:
                # does not pass
                continue

            self.block.append('{0}..{1} L={2} SP={3}'.format(i[0] + 1,
                                                             i[1] + 1,
                                                             i[1] - i[0] + 1,
                                                             round(perfect_matches, 2)))

            for j in range(i[0], i[1] + 1):
                base_ss = [x.get_base(j) for x in self._seq if x.line == 'ss']
                base_rs = [x.get_base(j) for x in self._seq if x.line == 'rs']
                base_rc = [x.get_base(j) for x in self._seq if x.line == 'rc']
                if base_rc:
                    apply_rc = True
                else:
                    apply_rc = False

                # base_ss = self._ss.get_base(j)
                # base_rs = self._rs.get_base(j)
                # base_rc = self._rc.get_base(j)

                base_ss_set = set(base_ss)
                base_rs_set = set(base_rs)
                base_rc_set = set(base_rc)

                if len(base_ss_set) != 1 or len(base_rs_set) != 1:
                    continue

                if apply_rc and len(base_rc_set) != 1:
                    continue

                if self._unknown_base.get(seqtype) in base_ss_set & base_rs_set:
                    continue

                if apply_rc and self._unknown_base.get(seqtype) in base_rc_set:
                    continue

                # if any([base_ss, base_rs, base_rc]) == self._unknown_base.get(seqtype):
                #     continue

                ret_profile = False
                if self._ast.neighbor_star_check(j, self.neighbor_len):
                    if apply_rc:
                        if base_ss_set == base_rc_set and base_ss_set != base_rs_set:
                            ret_profile = True
                    elif base_ss_set != base_rs_set:
                        ret_profile = True

                if ret_profile:
                    self.nmut += 1
                    for _seq in self._seq:
                        mutation_profile = MutationProfile()
                        mutation_profile.query_id = _seq.name
                        mutation_profile.line = _seq.line
                        mutation_profile.frame = _seq.frame
                        mutation_profile.profile = '{0}{1}{2}'.format(
                            base_ss[0], _seq.get_raw_position(j), base_rs[0])
                        self.profile.append(mutation_profile)
                        # for _seq in self._seq:
                        #     _seq.add_mutposition(j)
