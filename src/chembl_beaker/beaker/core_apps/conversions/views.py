__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from beaker import app
from bottle import request
from beaker.core_apps.conversions.impl import _ctab2smiles, _smiles2ctab, _inchi2ctab, _ctab2smarts
from beaker.core_apps.conversions.impl import _ctab2inchi, _inchi2inchiKey
from beaker.core_apps.conversions.impl import _canonicalize_smiles, _ctab2inchiKey
from beaker.core_apps.conversions.impl import _smiles2inchi, _smiles2inchiKey
from beaker.core_apps.conversions.impl import _smarts2ctab
from beaker.utils.io import _parseFlag

# ----------------------------------------------------------------------------------------------------------------------

def ctab2smilesView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['nameHeader'] = params.get('nameHeader', 'Name')
    kwargs['includeHeader'] = _parseFlag(params.get('includeHeader', True))
    return _ctab2smiles(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------

@app.route('/ctab2smiles', method=['OPTIONS', 'POST'], name="ctab2smiles")
def ctab2smiles():
    """
Converts CTAB to SMILES format. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST -F "file=@isomeric.mol" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@isomeric.mol" -F "isomericSmiles=1" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=0" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=1" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@explicitHs.mol" -F "removeHs=0" ${BEAKER_ROOT_URL}ctab2smiles
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab2smilesView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2smartsView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', True))
    return _ctab2smarts(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2smarts', method=['OPTIONS', 'POST'], name="ctab2smarts")
def ctab2smarts():
    """
Converts CTAB to SMILES format. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST -F "file=@isomeric.mol" ${BEAKER_ROOT_URL}ctab2smiles
    curl -X POST -F "file=@isomeric.mol" -F "isomericSmiles=1" ${BEAKER_ROOT_URL}ctab2smarts
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smarts
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=0" ${BEAKER_ROOT_URL}ctab2smarts
    curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=1" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smarts
    curl -X POST -F "file=@explicitHs.mol" -F "removeHs=0" ${BEAKER_ROOT_URL}ctab2smarts
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab2smartsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2ctabView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2ctab(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2ctab', method=['OPTIONS', 'POST'], name="smiles2ctab")
def smiles2ctab():
    """
Converts SMILES to CTAB. This method accepts single or multiple SMILES or *.smi file.
cURL examples:

    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2ctab
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2ctab
    curl -X POST -F "file=@rules.smi" -F "computeCoords=0"  ${BEAKER_ROOT_URL}smiles2ctab
    curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2ctab
    curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2ctab
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles2ctabView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smarts2ctabView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    return _smarts2ctab(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smarts2ctab', method=['OPTIONS', 'POST'], name="smarts2ctab")
def smarts2ctab():
    """
Converts SMARTS to CTAB. This method accepts single or multiple SMARTS.
cURL examples:

    curl -X POST -F "file=@amino.sma" ${BEAKER_ROOT_URL}smarts2ctab
    curl -X POST -F "file=@amino.sma" -F "computeCoords=0"  ${BEAKER_ROOT_URL}smarts2ctab
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smarts2ctabView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2inchiView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2inchi(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2inchi', method=['OPTIONS', 'POST'], name="smiles2inchi")
def smiles2inchi():
    """
Converts SMILES to InChi. This method accepts single or multiple SMILES or *.smi file.
cURL examples:

    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
    curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2inchi
    curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles2inchiView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2inchiKeyView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2inchiKey(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2inchiKey', method=['OPTIONS', 'POST'], name="smiles2inchiKey")
def smiles2inchiKey():
    """
Converts SMILES to InChi Key. This method accepts single or multiple SMILES or *.smi file.
cURL examples:

    curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2inchiKey
    curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
    curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2inchiKey
    curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchiKey
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return smiles2inchiKeyView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def canonicalizeSmilesView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['in_delimiter'] = params.get('in_delimiter', ' ')
    kwargs['out_delimiter'] = params.get('out_delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['nameHeader'] = params.get('nameHeader', 'Name')
    kwargs['includeHeader'] = _parseFlag(params.get('includeHeader', True))

    if params.get('titleLine') is None and not data.startswith(b'SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _canonicalize_smiles(data, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------


@app.route('/canonicalizeSmiles', method=['OPTIONS', 'POST'], name="canonicalizeSmiles")
def canonicalizeSmiles():
    """
Converts SMILES to canonical SMILES. This method accepts single or multiple SMILES or *.smi file.
cURL examples:

    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@aspirin_with_header.smi" -F "out_delimiter=|" -F "nameHeader=foo" ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=0" -F "sanitize=0" ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=0" -F "sanitize=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=1" -F "sanitize=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@isomeric.smi" ${BEAKER_ROOT_URL}canonicalizeSmiles
    curl -X POST -F "file=@isomeric.smi" -F "isomericSmiles=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return canonicalizeSmilesView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2ctab', method=['OPTIONS', 'POST'], name="inchi2ctab")
def inchi2ctab():
    """
Converts InChi to CTAB. This method accepts one or multiple InChis.
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2ctab
    curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2ctab
    """

    inchis = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return _inchi2ctab(inchis)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2inchiView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _ctab2inchi(data, **kwargs) 

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchi', method=['OPTIONS', 'POST'], name="ctab2inchi")
def ctab2inchi():
    """
Converts CTAB to InChis. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchi
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchi
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab2inchiView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2inchiKeyView(data, params):
    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _ctab2inchiKey(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchiKey', method=['OPTIONS', 'POST'], name="ctab2inchiKey")
def ctab2inchiKey():
    """
Converts CTAB to InChi Keys. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchiKey
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchiKey
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return ctab2inchiKeyView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2inchiKey', method=['OPTIONS', 'POST'], name="inchi2inchiKey")
def inchi2inchiKey():
    """
Converts InChis to InChiKeys. This method accepts one or multiple InChis.
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2inchiKey
    curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2inchiKey
    """

    inchis = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return _inchi2inchiKey(inchis)

# ----------------------------------------------------------------------------------------------------------------------
