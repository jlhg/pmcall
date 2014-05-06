from subprocess import call


def makeblastdb(path_in, dbtype, path_out):
    cmd = """
    set -e
    makeblastdb -in {0} -dbtype {1} -out {2} -parse_seqids
    """.format(path_in, dbtype, path_out)

    return call(cmd, shell=True)


def blastx(query, db, out):
    cmd = """
    set -e
    blastx -query {0} -db {1} -out {2} -outfmt "6 std qlen slen qframe sframe \
    gaps nident positive ppos qcovhsp staxids salltitles"
    """.format(query, db, out)

    return call(cmd, shell=True)


def besthit(path_in, path_out):
    cmd = """
    set -e
    cat {0} |
    awk -F'\t' '$1 !~ /^#/ {{OFS="\t"; print $0\}}' | \
    sort -t$'\t' -k1d,1 -k11g,11 -k3gr,3 -k4gr,4 -k21gr,21 -k13gr | \
    sort -t$'\t' -k1,1 -u >{1}
    """.format(path_in, path_out)

    return call(cmd, shell=True)


def parse_blastext_frame(blastext_path):
    with open(blastext_path, 'r') as fi:
        for line in fi:
            if line.startswith('#'):
                continue
            data = line.split('\t')
            yield data[0], data[1], data[14], data[15]
