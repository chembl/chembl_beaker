__author__ = 'mnowotka'

import tempfile
import os.path
from subprocess import call
import requests

FTP_ROOT = "ftp.ebi.ac.uk/pub/databases/chembl"
EBI_FTP = FTP_ROOT + "/ChEMBLdb/releases/"
CHEMBL_API_ROOT = "https://www.ebi.ac.uk/chembl/api/data"

# ----------------------------------------------------------------------------------------------------------------------


def _get_latest_chembl_version():
    ret = requests.get('{0}/{1}'.format(CHEMBL_API_ROOT, 'status.json'))
    return ret.json()['chembl_db_version'].lower()

# ----------------------------------------------------------------------------------------------------------------------


def _parse_results(file_name):
    ret = {}
    if not file_name or not os.path.exists(file_name):
        return ret
    with open(file_name, 'r') as f:
        for line in f.readlines():
            data = line.strip().split('\t')
            ids = data[1].split(',')
            seq = data[0]
            bit_score = data[11]
            if seq not in ret:
                ret[seq] = {}
            for id in ids:
                if id not in ret[seq]:
                    ret[seq][id] = [float(bit_score)]
                else:
                    ret[seq][id].append(float(bit_score))
    for seq in ret:
        for id in ret[seq]:
            ret[seq][id] = sum(ret[seq][id]) / len(ret[seq][id])
    try:
        os.remove(file_name)
    except:
        pass
    return ret

# ----------------------------------------------------------------------------------------------------------------------


def _blast(sequence, type, blast_bindir, chembl_db_release=None, seg=False, blast_datadir=None, blast_tmpdir=None):
    if type not in ('COMPOUNDS', 'TARGETS'):
        return
    if not chembl_db_release:
        chembl_db_release = _get_latest_chembl_version()
    if not blast_datadir:
        blast_datadir = tempfile.gettempdir()
    if not blast_tmpdir:
        blast_tmpdir = tempfile.gettempdir()
    fd, fpath = tempfile.mkstemp(prefix='blast_', suffix='.fa', dir=blast_tmpdir)

    fname = os.path.splitext(os.path.basename(fpath))[0]
    os.write(fd, sequence)
    os.close(fd)

    arguments = [os.path.join(blast_bindir, 'blastall')]
    flags = ['-m', '8', '-Y', '1000000000', '-p', 'blastp']
    additional_flags = ["-F", 'T']
    if seg:
        additional_flags = ["-F", 'T']

    blast_database_path = os.path.join(blast_datadir, chembl_db_release)
    if type == 'COMPOUNDS':
        blast_database_path += '_bio'
    out_file_name_tabbed = os.path.join(blast_tmpdir, "{0}.tabbed".format(fname))
    db = ['-d', blast_database_path]
    io = ['-i', fpath, '-o', out_file_name_tabbed]

    call(arguments + flags + additional_flags + db + io)

    return _parse_results(out_file_name_tabbed)


# ----------------------------------------------------------------------------------------------------------------------

def _blast_all(sequence, blast_bindir, chembl_db_release=None, seg=False, blast_datadir=None, blast_tmpdir=None):
    return {
        'compounds': _blast(sequence, 'COMPOUNDS', blast_bindir, chembl_db_release, seg, blast_datadir, blast_tmpdir),
        'targets': _blast(sequence, 'TARGETS', blast_bindir, chembl_db_release, seg, blast_datadir, blast_tmpdir),
    }

# ----------------------------------------------------------------------------------------------------------------------
