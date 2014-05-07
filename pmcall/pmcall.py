#!/usr/bin/env python
#
# Copyright (C) 2014 Jian-Long Huang <jianlong@ntu.edu.tw>

__version__ = '1.0'

import sys
import logging
import argparse
import signal
from operator import and_
from os import makedirs, listdir
from os.path import join, exists, abspath, basename, isfile
from logging import StreamHandler, Formatter
from shutil import rmtree
from multiprocessing import Pool
from multiprocessing.managers import SyncManager
from blast import blastx, makeblastdb
from blast import besthit, parse_blastext_frame
from seqio import parse_fasta, translate, Seq
from msa import mafft
import msap
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
    pool_timeout = 65535
    task = []

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

    blastx_dir = join(args.out, 'blastx')
    msain_aa_dir = join(args.out, 'msain_aa')
    msain_nt_dir = join(args.out, 'msain_nt')
    clu_aa_dir = join(args.out, 'clu_aa')
    clu_nt_dir = join(args.out, 'clu_nt')
    makedirs(blastx_dir)
    makedirs(msain_aa_dir)
    makedirs(msain_nt_dir)
    makedirs(clu_aa_dir)
    makedirs(clu_nt_dir)

    # Perform blastx and besthit
    if not args.blastx:
        blastdbout = join(blastx_dir, basename(args.ref))
        makeblastdb(path_in=args.ref, dbtype='prot', path_out=blastdbout)

        task[:] = []
        args.blast = []
        for i in args.ss:
            blastout = join(blastx_dir, '{0}_to_{1}'.format(basename(i), basename(args.ref)))
            task.append((i, args.ref, blastout))
            args.blast.append(blastout)
        Pool(args.nthread).map_async(blastx, task).get(pool_timeout)

    task[:] = []
    for i in args.blast:
        task.append((i, join(blastx_dir, '{0}.besthit'.format(basename(i)))))
    Pool(args.nthread).map_async(besthit, task).get(pool_timeout)

    # Parse blastx files
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

    # Multiple sequence alignment
    task[:] = []
    for sid in common_sid:
        namelength = 0
        with open(join(msain_aa_dir, sid), 'w') as fo:
            for qid in sid_qidset.get(sid):
                if len(qid) > namelength:
                    namelength = len(qid)
                fo.write('>')
                fo.write(seq.get(qid).header)
                fo.write('\n')
                fo.write(translate(seq.get(qid).sequence), qframe.get(''.join(qid, sid)))
                fo.write('\n')
        path_fa = join(msain_aa_dir, '{0}.fa'.format(sid))
        path_clu = join(clu_aa_dir, '{0}.clu'.format(sid))
        task.append((path_fa, path_clu, namelength + 4))

    Pool(args.nthread).map_async(mafft, task).get(pool_timeout)

    # Parse clustal files
    manager = SyncManager()
    manager.start(lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
    shared_list = manager.list()
    parser = msap.Parse(perfact_matches_percent=args.perfect_matches_percent,
                        block_len=args.block_len,
                        neighbor_length=args.neighbor_length,
                        neighbor_matches_percent=args.neighbor_matches_percent)
    task[:] = [(join(msain_aa_dir, f), parser, shared_list) for f in listdir(msain_aa_dir) if isfile(f)]
    Pool(args.nthread).map_async(parse_clustal, task).get(pool_timeout)

    with open(join(args.out, 'result.txt'), 'w') as fo:
        for i in shared_list:
            fo.write(i)
            fo.write('\n')
            fo.flush()


def parse_clustal(path_in, parser, shared_list):
    result = parser.parse(path_in)
    shared_list.append(result.data)


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

    parser.add_argument('-pp', dest='perfect_matches_percent',
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

    parser.add_argument('-nthread', dest='nthread',
                        metavar='<int>',
                        type=int,
                        default=1,
                        help='number of threads to perform this program (default: 1)')
    args = parser.parse_args()
    print(args)
    # run(args)


if __name__ == '__main__':
    main()
