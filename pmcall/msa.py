from subprocess import call


def mafft(path_in, path_out, namelength=30):
    cmd = """
    mafft --genafpair --maxiterate 1000 \
    --preservecase --anysymbol \
    --clustalout --namelength {0} \
    {1} >{2} 2>/dev/null
    """.format(namelength, path_in, path_out)

    return call(cmd, shell=True, executable='/bin/bash')
