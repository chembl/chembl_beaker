__author__ = 'efelix'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import AllChem
from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.D2Coords.impl import _ctab22D, _smiles22D
from chembl_beaker.beaker.utils.io import _parseFlag
import base64

# ----------------------------------------------------------------------------------------------------------------------


def ctab22DView(data, params):
    kwargs = dict()
    kwargs['multi'] = int(params.get('multi', False))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return _ctab22D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab22D/<ctab>', method=['OPTIONS', 'GET'], name="ctab22D")
def ctab22D(ctab):
    """
Generate 2D coordinates from molfile using Schrodinger's coordgen.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cuRL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab22D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab22DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab22D', method=['OPTIONS', 'POST'], name="ctab22D")
def ctab22D():
    """
Generate 2D coordinates from molfile using Schrodinger's coordgen.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}ctab22D
    curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}ctab22D
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab22DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles22DView(data, params):
    kwargs = dict()
    kwargs['multi'] = int(params.get('multi', False))
    kwargs['computeCoords'] = False
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles22D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles22D/<smiles>', method=['OPTIONS', 'GET'], name="smiles22D")
def smiles22D(smiles):
    """

Generate 2D coordinates from SMILES using Schrodinger's coordgen.
CTAB is urlsafe_base64 encoded string containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles22D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles22D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles22DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles22D', method=['OPTIONS', 'POST'], name="smiles22D")
def smiles22D():
    """
Generate 2D coordinates from SMILES using Schrodinger's coordgen.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles22D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles22D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles22D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles22D
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles22DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
