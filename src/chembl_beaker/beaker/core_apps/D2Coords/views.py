__author__ = 'efelix'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import AllChem
from beaker import app
from bottle import request
from beaker.core_apps.D2Coords.impl import _ctab22D, _smiles22D, _is3D
from beaker.utils.io import _parseFlag

# ----------------------------------------------------------------------------------------------------------------------


def ctab22DView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))

    return _ctab22D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


def is3DView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', True))

    return _is3D(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab22D', method=['OPTIONS', 'POST'], name="ctab22D")
def ctab22D():
    """
Generate 2D coordinates for a molecule using Schrodinger's coordgen.
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


@app.route('/is3D', method=['OPTIONS', 'POST'], name="is3D")
def is3D():
    """
Check if molecule has any 3D coordinate.
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}is3D
    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}is3D
    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}is3D
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}is3D
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return is3DView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
