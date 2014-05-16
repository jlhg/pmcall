#!/usr/bin/env python
#
# Copyright (C) 2014 Jian-Long Huang <jianlong@ntu.edu.tw>

__version__ = '1.0'

import re
import sys
import logging
import argparse
import signal
from operator import and_
from os import makedirs, listdir
from os.path import join, exists, basename, isfile, splitext
from logging import StreamHandler, Formatter
from shutil import rmtree
from multiprocessing import Pool
from multiprocessing.managers import SyncManager
from pmcall.thirdparty.shutil33 import which
from pmcall.blast import blastx, makeblastdb
from pmcall.blast import besthit, parse_blastext_frame
from pmcall.seqio import parse_fasta, translate
from pmcall.msa import mafft
from pmcall.mutation import MutationParser
from pmcall.clustal import aa_to_nt


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


def check_required():
    commands = [
        'mafft',
        'blastx',
    ]
    for i in commands:
        if not which(i):
            return i


def replace_nonalphanum(text, symbol='_'):
    return re.sub('[^0-9a-zA-Z_.-]+', symbol, text)


def fwrite(f, contents):
    with open(f, 'w') as fo:
        fo.write(contents)
        fo.flush()


def blastx_wrapper(args):
    return blastx(*args)


def besthit_wrapper(args):
    return besthit(*args)


def mafft_wrapper(args):
    return mafft(*args)


def aa_to_nt_wrapper(args):
    return aa_to_nt(*args)


def parse_clustal_wrapper(args):
    return parse_clustal(*args)


def run(args):
    # Check the required programs
    required_program = check_required()
    if required_program:
        sys.exit('{0} is required to perform the mutation discovery.'.format(required_program))

    flog = join(args.dout, 'log')
    fstep = join(args.dout, 'step')
    fprofile_aa = join(args.dout, 'mutation_profiles_aa.txt')
    fprofile_nt = join(args.dout, 'mutation_profiles_nt.txt')
    dmakeblastdb = join(args.dout, 'blastdb')
    dblastx = join(args.dout, 'blastx')
    dbesthit = join(args.dout, 'besthit')
    dfa_aa = join(args.dout, 'fa_aa')
    dfa_nt = join(args.dout, 'fa_nt')
    dclu_aa = join(args.dout, 'clu_aa')
    dclu_nt = join(args.dout, 'clu_nt')

    if exists(args.dout) and not exists(fstep):
        resp = prompt('The directory {0} existed, remove it?'.format(args.dout))
        if resp == 'yes':
            rmtree(args.dout)
            makedirs(args.dout)
        else:
            sys.stdout.write('Cancel to perform the mutation discovery.\n')
            sys.stdout.flush()
            return

    if not exists(args.dout):
        makedirs(args.dout)

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=flog,
        filemode='w',
    )
    stream_handler = StreamHandler(sys.stdout)
    formatter = Formatter('%(asctime)s [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(stream_handler)

    logger = logging.getLogger('main')
    logger.info('Starting to perform mutation discovery.')
    logger.info('Block length (-bl): {0}'.format(args.blocklen))
    logger.info('Perfect matching percent (-pp): {0}'.format(args.perfect_match_percent))
    logger.info('Side length (-sl): {0}'.format(args.sidelen))
    logger.info('Side perfect matching percent (-sp): {0}'.format(args.side_match_percent))

    if exists(fstep):
        step = int(open(fstep, 'r').read())
    elif args.blastx:
        step = 2
        fwrite(fstep, '2')
    else:
        step = 1
        fwrite(fstep, '1')

    pool_timeout = 65535
    task = []

    if args.rc:
        apply_rc = True
    else:
        apply_rc = False

    if step == 1:
        # Perform blastx
        logger.info('Performing blastx.')
        if exists(dblastx):
            rmtree(dblastx)
        makedirs(dblastx)

        fmakeblastdbout = join(dblastx, basename(args.ref))
        makeblastdb(path_in=args.ref, dbtype='prot', path_out=fmakeblastdbout)

        task[:] = []
        args.blastx = []

        if apply_rc:
            fqueries = args.ss + args.rs + args.rc
        else:
            fqueries = args.ss + args.rs

        for i in fqueries:
            fblastxout = join(dblastx, '{0}_to_{1}.blastx'.format(basename(i), basename(args.ref)))
            task.append((i, fmakeblastdbout, fblastxout))
            args.blastx.append(fblastxout)
        Pool(args.nthread).map_async(blastx_wrapper, task).get(pool_timeout)

        step = 2
        fwrite(fstep, '2')

    if step == 2:
        # Perform best hit filter
        logger.info('Performing best hit filter.')
        if exists(dbesthit):
            rmtree(dbesthit)
        makedirs(dbesthit)

        task[:] = []

        if args.blastx:
            for i in args.blastx:
                task.append((i, join(dbesthit, '{0}.besthit'.format(basename(i)))))
        else:
            for i in (join(dblastx, x) for x in listdir(dblastx) if isfile(join(dblastx, x)) and splitext(x)[-1] == '.blastx'):
                task.append((i, join(dbesthit, '{0}.besthit'.format(basename(i)))))
        Pool(args.nthread).map_async(besthit_wrapper, task).get(pool_timeout)

        step = 3
        fwrite(fstep, '3')

    if step == 3:
        # Generate protein fasta files for multiple sequence alignment
        logger.info('Generating protein fasta files for multiple sequence alignment.')
        if exists(dfa_aa):
            rmtree(dfa_aa)
        makedirs(dfa_aa)

        # Parse best hit blastx
        logger.info('Reading blastx files.')
        sid_sets = []
        idmap_qs = {}
        sid_qidset = {}
        qframe = {}
        for i in (join(dbesthit, x) for x in listdir(dbesthit) if isfile(join(dbesthit, x))):
            sid_set = set()
            for qid, sid, qf, sf in parse_blastext_frame(i):
                sid_set.add(sid)
                idmap_qs.update({qid: sid})
                if sid in sid_qidset:
                    sid_qidset.get(sid).add(qid)
                else:
                    sid_qidset.update({sid: {qid}})
                qframe.update({''.join([qid, sid]): int(qf)})
            sid_sets.append(sid_set)

        common_sid_set = reduce(and_, sid_sets)

        # Parse FASTA files
        logger.info('Reading fasta files.')
        qseq = {}
        qline = {}
        for i in args.ss:
            for header, sequence in parse_fasta(i):
                if idmap_qs.get(header) in common_sid_set:
                    qseq.update({header: sequence})
                    qline.update({header: 'ss'})
        for i in args.rs:
            for header, sequence in parse_fasta(i):
                if idmap_qs.get(header) in common_sid_set:
                    qseq.update({header: sequence})
                    qline.update({header: 'rs'})
        if apply_rc:
            for i in args.rc:
                for header, sequence in parse_fasta(i):
                    if idmap_qs.get(header) in common_sid_set:
                        qseq.update({header: sequence})
                        qline.update({header: 'rc'})

        logger.info('Generating protein fasta files.')
        for sid in common_sid_set:
            namelength = 0
            data = []
            for qid in sid_qidset.get(sid):
                qf = qframe.get(''.join([qid, sid]))
                name = '{0}({1},{2})'.format(qid, qline.get(qid), qf)
                if len(name) > namelength:
                    namelength = len(name)
                data.append('>{0}\n'.format(name))
                data.append('{0}\n'.format(translate(qseq.get(qid), qf)))
            namelength += 4
            with open(join(dfa_aa, '{0}.{1}.fa'.format(replace_nonalphanum(sid), namelength)), 'w') as fo:
                fo.write(''.join(data))
                fo.flush()

        step = 4
        fwrite(fstep, '4')

    if step == 4:
        # Perform MAFFT (multiple sequence alignment)
        logger.info('Performing MAFFT.')
        if exists(dclu_aa):
            rmtree(dclu_aa)
        makedirs(dclu_aa)

        task[:] = []
        for i in (x for x in listdir(dfa_aa) if isfile(join(dfa_aa, x))):
            seg = i.split('.')
            namelength = int(seg[-2])
            task.append((join(dfa_aa, i), join(dclu_aa, '{0}.clu'.format('.'.join(seg[:-1]))), namelength))
        Pool(args.nthread).map_async(mafft_wrapper, task).get(pool_timeout)

        step = 5
        fwrite(fstep, '5')

    if step == 5:
        # Convert aa clustal files to nt clustal files
        logger.info('Converting protein clustal files to nucleotide clustal files.')
        if exists(dclu_nt):
            rmtree(dclu_nt)
        makedirs(dclu_nt)

        if 'qseq' not in locals():
            sid_sets = []
            idmap_qs = {}
            sid_qidset = {}
            for i in (join(dbesthit, x) for x in listdir(dbesthit) if isfile(join(dbesthit, x))):
                sid_set = set()
                for qid, sid, qf, sf in parse_blastext_frame(i):
                    sid_set.add(sid)
                    idmap_qs.update({qid: sid})
                    if sid in sid_qidset:
                        sid_qidset.get(sid).add(qid)
                    else:
                        sid_qidset.update({sid: {qid}})
                sid_sets.append(sid_set)

            common_sid_set = reduce(and_, sid_sets)
            qseq = {}
            for i in args.ss:
                for header, sequence in parse_fasta(i):
                    if idmap_qs.get(header) in common_sid_set:
                        qseq.update({header: sequence})
            for i in args.rs:
                for header, sequence in parse_fasta(i):
                    if idmap_qs.get(header) in common_sid_set:
                        qseq.update({header: sequence})
            if apply_rc:
                for i in args.rc:
                    for header, sequence in parse_fasta(i):
                        if idmap_qs.get(header) in common_sid_set:
                            qseq.update({header: sequence})

        task[:] = []
        for i in (x for x in listdir(dclu_aa) if isfile(join(dclu_aa, x))):
            task.append((join(dclu_aa, i), join(dclu_nt, i), qseq, 60))
        Pool(args.nthread).map_async(aa_to_nt_wrapper, task).get(pool_timeout)

        step = 6
        fwrite(fstep, '6')

    if step == 6:
        # Parse clustal files and report mutation profiles
        logger.info('Begining to parse clustal files and report mutation profiles.')
        manager = SyncManager()
        manager.start(lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))

        profile_aa = manager.dict()
        mutationparser = MutationParser(
            blocklen=args.blocklen,
            perfect_match_percent=args.perfect_match_percent,
            sidelen=args.sidelen,
            side_match_percent=args.side_match_percent)

        task[:] = []
        for i in (join(dclu_aa, x) for x in listdir(dclu_aa) if isfile(join(dclu_aa, x))):
            task.append((i, mutationparser, 'aa', profile_aa))
        Pool(args.nthread).map_async(parse_clustal_wrapper, task).get(pool_timeout)

        write_profile(fprofile_aa, 'aa', profile_aa)

        profile_nt = manager.dict()

        task[:] = []
        for i in (join(dclu_nt, x) for x in listdir(dclu_nt) if isfile(join(dclu_nt, x))):
            task.append((i, mutationparser, 'nt', profile_nt))
        Pool(args.nthread).map_async(parse_clustal_wrapper, task).get(pool_timeout)

        write_profile(fprofile_nt, 'nt', profile_nt)

        step = 7
        fwrite(fstep, '7')

    logger.info('Mutation discovery finished.')


def parse_clustal(path_in, mutationparser, seqtype, clustal_mutation_profiles):
    profile = mutationparser.parse(path_in, seqtype)
    clustal_mutation_profiles.update({basename(path_in).rstrip('.fa.clu'): profile})


def write_profile(fprofile, seqtype, clustal_mutation_profiles):
    namepattern = re.compile('(\S+)\((ss|rs|rc),([+-]*\d+)\)')

    with open(fprofile, 'w') as fo:
        fo.write('# Mutation Profiles\n')
        fo.write('\t'.join([
            'n_block',
            'n_mutation',
            'subject_id',
            'query_id',
            'mutation_profile',
        ]))
        fo.write('\n')

        for sid, cm_profile in clustal_mutation_profiles.items():
            for qid, m_profile in cm_profile.mutation_profiles.items():
                if seqtype == 'nt':
                    frame = int(namepattern.match(qid).group(3))
                    posfix_m_profile = []
                    for i in m_profile:
                        posfix = re.sub('(\d+)', lambda x: str(int(x.group(1)) + (abs(frame) - 1)), i)
                        posfix_m_profile.append(posfix)
                    fo.write('\t'.join([
                        str(cm_profile.nblock),
                        str(cm_profile.nmutation),
                        sid,
                        qid,
                        ', '.join('{0}({1})'.format(m_profile[x], posfix_m_profile[x]) for x in range(len(m_profile))),
                    ]))
                else:
                    fo.write('\t'.join([
                        str(cm_profile.nblock),
                        str(cm_profile.nmutation),
                        sid,
                        qid,
                        ', '.join(m_profile),
                    ]))
                fo.write('\n')
                fo.flush()


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
    parser.add_argument('-ref', dest='ref',
                        metavar='<file>',
                        required=True,
                        help='reference sequences in FASTA format')
    parser.add_argument('-out', dest='dout',
                        metavar='<dir>',
                        required=True,
                        help='output directory')
    parser.add_argument('-blastx', dest='blastx',
                        metavar='<file>',
                        nargs='+',
                        help='blastx supports in extended tabular format')
    parser.add_argument('-bl', dest='blocklen',
                        metavar='<int>',
                        type=int,
                        default=10,
                        help='block length (default: 10)')
    parser.add_argument('-pp', dest='perfect_match_percent',
                        metavar='<float>',
                        type=float,
                        default=0.9,
                        help='percentage of the number of perfect matches '
                        'in a clustal block (default: 0.9)')
    parser.add_argument('-sl', dest='sidelen',
                        metavar='<int>',
                        type=int,
                        default=5,
                        help='length of each neighboring side for a '
                        'point mutation (default: 5)')
    parser.add_argument('-sp', dest='side_match_percent',
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

    run(args)


if __name__ == '__main__':
    main()
