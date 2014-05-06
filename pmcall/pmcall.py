#!/usr/bin/env python
#
# Copyright (C) 2014 Jian-Long Huang <jianlong@ntu.edu.tw>

__version__ = '1.0'

import sys
import logging
import argparse
from operator import and_
from os import makedirs
from os.path import join, exists, abspath, basename
from logging import StreamHandler, Formatter
from shutil import rmtree
from blast import blastx, makeblastdb
from blast import besthit, parse_blastext_frame
from seqio import parse_fasta, translate, Seq
# from obejct import Sequence


def prompt(msg):
    while 1:
        resp = raw_input('{0} [yes/no] '.format(msg))
        resp = resp.strip().lower()
        if resp in {'yes', 'no'}:
            return resp
        else:
            sys.stdout.write('Your response (\'{0}\') was not one of the expected '
                             'responses: yes, no\n'.format(resp))
            sys.stdout.flush()


def run(args):
    if exists(args.outdir):
        resp = prompt('The directory {0} existed, remove it?'.format(args.outdir))
        if resp == 'yes':
            rmtree(args.outdir)
            makedirs(args.outdir)
        else:
            sys.stdout.write('Cancel to perform the pmcall.\n')
            sys.stdout.flush()
            return

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=join(args.outdir, 'log'),
        filemode='w',
    )
    stream_handler = StreamHandler(sys.stdout)
    formatter = Formatter('%(asctime)s [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(stream_handler)

    if args.rc:
        apply_rc = True
    else:
        apply_rc = False

    # Change to absolute path
    for i, j in enumerate(args.ss):
        args.ss[i] = abspath(j)

    for i, j in enumerate(args.rs):
        args.rs[i] = abspath(j)

    if apply_rc:
        for i, j in enumerate(args.rc):
            args.rc[i] = abspath(j)

    if args.blastx:
        for i in enumerate(args.blastx):
            args.blastx[i] = abspath(j)

    args.ref = abspath(args.ref)

    makedirs(join(args.out, 'blast'))
    makedirs(join(args.out, 'msain_aa'))
    makedirs(join(args.out, 'msain_nt'))
    makedirs(join(args.out, 'clu_aa'))
    makedirs(join(args.out, 'clu_nt'))

    # Perform blastx and besthit
    if not args.blastx:
        blastdbout = join(args.out, 'blast', basename(args.ref))
        makeblastdb(path_in=args.ref, dbtype='prot', path_out=blastdbout)
        args.blast = []
        for i in args.ss:
            blastout = join(args.out, 'blast', '_'.join(basename(i), basename(args.ref)))
            blastx(query=i, db=args.ref, out=blastout)
            args.blast.append(blastout)

    for i in args.blast:
        besthit(path_in=i, path_out=join(args.out, 'blast', basename(i) + '.besthit'))

    # Parse the blastx results
    sidsets = []
    idmap_qs = {}
    sid_qidset = {}
    qframe = {}
    for i in args.blast:
        sidset = set()
        for qid, sid, qframe, sframe in parse_blastext_frame(i):
            sidset.add(sid)
            idmap_qs.update({qid: sid})
            if sid in sid_qidset:
                sid_qidset.get(sid).add(qid)
            else:
                sid_qidset.update({sid: {qid}})
            qframe.update({''.join(qid, sid): int(qframe)})
        sidsets.append(sidset)

    common_sid = reduce(and_, sidsets)

    # Parse FASTA files
    seq = {}
    for i in args.ss:
        for header, sequence in parse_fasta(i):
            if idmap_qs.get(header) in common_sid:
                seq.update({header: Seq(header, sequence, 'ss')})

    for i in args.rs:
        for header, sequence in parse_fasta(i):
            if idmap_qs.get(header) in common_sid:
                seq.update({header: Seq(header, sequence, 'rs')})

    if apply_rc:
        for i in args.rc:
            for header, sequence in parse_fasta(i):
                if idmap_qs.get(header) in common_sid:
                    seq.update({header: Seq(header, sequence, 'rc')})

    for sid in common_sid:
        with open(join(args.out, 'msain_aa', sid), 'w') as fo:
            for qid in sid_qidset.get(sid):
                fo.write('>')
                fo.write(seq.get(qid).header)
                fo.write('\n')
                fo.write(translate(seq.get(qid).sequence), qframe.get(''.join(qid, sid)))
                fo.write('\n')

    #


def main():
    parser = argparse.ArgumentParser(prog='pmcall',
                                     description='Point mutation discovery',
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=70))

    parser.add_argument('-ss', dest='ss',
                        metavar='<file>',
                        nargs='+',
                        required=True,
                        help='susceptible line sequences in FASTA format')
    parser.add_argument('-rs', dest='rs',
                        metavar='<file>',
                        nargs='+',
                        required=True,
                        help='resistant line sequences in FASTA format')
    parser.add_argument('-rc', dest='rc',
                        metavar='<file>',
                        nargs='+',
                        help='recovered line sequences in FASTA format')
    parser.add_argument('-blastx', dest='blastx',
                        metavar='<file>',
                        nargs='+',
                        required=True,
                        help='blastx supports in extended tabular format')
    parser.add_argument('-ref', dest='refseq',
                        metavar='<file>',
                        required=True,
                        help='reference sequences in FASTA format')
    parser.add_argument('-out', dest='outdir',
                        metavar='<dir>',
                        required=True,
                        help='output directory')
    parser.add_argument('-bl', dest='block_len',
                        metavar='<int>',
                        type=int,
                        default=10,
                        help='block length (default: 10)')

    parser.add_argument('-pp', dest='perfact_matches_percent',
                        metavar='<float>',
                        type=float,
                        default=0.9,
                        help='percentage of the number of perfect matches '
                        'in a clustal block (default: 0.9)')

    parser.add_argument('-ln', dest='neighbor_length',
                        metavar='<int>',
                        type=int,
                        default=5,
                        help='length of each neighboring side for a '
                        'point mutation (default: 5)')

    parser.add_argument('-np', dest='neighbor_matches_percent',
                        metavar='<float>',
                        type=float,
                        default=0.8,
                        help='percentage of the number of perfect matches in '
                        'each neighboring side (default: 0.8)')

    parser.add_argument('-cpu', dest='cpu',
                        metavar='<int>',
                        type=int,
                        default=1,
                        help='number of threads to perform this program (default: 1)')
    args = parser.parse_args()
    print(args)
    # run(args)


if __name__ == '__main__':
    main()
