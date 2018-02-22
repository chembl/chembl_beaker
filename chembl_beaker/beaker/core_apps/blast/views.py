__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from bottle import request
import base64
from chembl_beaker.beaker import app, config
from chembl_beaker.beaker.core_apps.blast.impl import _blast_all

# ----------------------------------------------------------------------------------------------------------------------


def blastView(seq, params):

    kwargs = dict()
    kwargs['seg'] = params.get('seg', False)
    kwargs['blast_bindir'] = config.get('blast_bindir')
    kwargs['chembl_db_release'] = config.get('chembl_db_release')
    kwargs['blast_datadir'] = config.get('blast_datadir')
    kwargs['blast_tmpdir'] = config.get('blast_tmpdir')

    return _blast_all(seq, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/blast/<sequence>', method=['OPTIONS', 'GET'], name="blast")
def blast(sequence):
    """
Uses BLAST to search in sequences stored in ChEMBL:
(Please note that this method uses old BLAST API and is intended for the internal use at ChEMBL).

    curl -X GET ${BEAKER_ROOT_URL}blast/$(cat seq.fa | base64 -w 0 | tr "+/" "-_")
    """

    seq = base64.urlsafe_b64decode(sequence)
    return blastView(seq, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/blast', method=['OPTIONS', 'POST'], name="blast")
def image2ctab():
    """
Uses BLAST to search in sequences stored in ChEMBL:
(Please note that this method uses old BLAST API and is intended for the internal use at ChEMBL).

    curl -X POST --data-binary @seq.fa ${BEAKER_ROOT_URL}blast
    curl -X POST -F "file=@seq.fa" ${BEAKER_ROOT_URL}blast
    """
    seq = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return blastView(seq, request.params)

# ----------------------------------------------------------------------------------------------------------------------

