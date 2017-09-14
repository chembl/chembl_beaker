__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.conversions.impl import _ctab2smiles, _smiles2ctab, _inchi2ctab, _ctab2smarts
from chembl_beaker.beaker.core_apps.conversions.impl import _ctab2inchi, _inchi2inchiKey, _ctab2xyz
from chembl_beaker.beaker.core_apps.conversions.impl import _canonicalize_smiles, _ctab2inchiKey
from chembl_beaker.beaker.core_apps.conversions.impl import _smiles2inchi, _smiles2inchiKey
from chembl_beaker.beaker.core_apps.conversions.impl import _smarts2ctab
from chembl_beaker.beaker.utils.io import _parseFlag
import base64

# ----------------------------------------------------------------------------------------------------------------------


def ctab2smilesView(data, params):
    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['nameHeader'] = params.get('nameHeader', 'Name')
    kwargs['includeHeader'] = _parseFlag(params.get('includeHeader', True))
    kwargs['isomericSmiles'] = _parseFlag(params.get('isomericSmiles', False))
    kwargs['kekuleSmiles'] = _parseFlag(params.get('kekuleSmiles', False))
    return _ctab2smiles(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2smiles/<ctab>', method=['OPTIONS', 'GET'], name="ctab2smiles")
def ctab2smiles(ctab):
    """
Converts CTAB to SMILES format. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2smiles/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}ctab2smiles/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")?isomericSmiles=1
    curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=0"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=1&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat explicitHs.mol | base64 -w 0 | tr "+/" "-_")"?removeHs=0"
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2smilesView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2smilesView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2smartsView(data, params):
    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['isomericSmiles'] = _parseFlag(params.get('isomericSmiles', False))
    return _ctab2smarts(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2smarts/<ctab>', method=['OPTIONS', 'GET'], name="ctab2smarts")
def ctab2smarts(ctab):
    """
Converts CTAB to SMARTS format. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2smarts/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}ctab2smarts/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")?isomericSmiles=1
    curl -X GET "${BEAKER_ROOT_URL}ctab2smarts/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smarts/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=0"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smarts/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=1&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}ctab2smarts/"$(cat explicitHs.mol | base64 -w 0 | tr "+/" "-_")"?removeHs=0"
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2smartsView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2smartsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2xyzView(data, params):
    computeCoords = _parseFlag(params.get('computeCoords', True))
    return _ctab2xyz(data, computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2xyz/<ctab>', method=['OPTIONS', 'GET'], name="ctab2xyz")
def ctab2xyz(ctab):
    """
Converts the molecules in the CTAB to xyz format. CTAB is a single molfile
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2xyz/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET "${BEAKER_ROOT_URL}ctab2xyz/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}ctab2xyz/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?computeCoords=0"
    """
    data = base64.urlsafe_b64decode(ctab)
    return ctab2xyzView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2xyz', method=['OPTIONS', 'POST'], name="ctab2xyz")
def ctab2xyz():
    """
Converts the molecules in the CTAB to xyz format. CTAB is a single molfile
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2xyz
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2xyz
    curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}ctab2xyz
    curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0"  ${BEAKER_ROOT_URL}ctab2xyz
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2xyzView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2ctabView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2ctab(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2ctab/<smiles>', method=['OPTIONS', 'GET'], name="smiles2ctab")
def smiles2ctab(smiles):
    """
Converts SMILES to CTAB. This method accepts urlsafe_base64 encoded string containing single or multiple SMILES
optionally containing header line, specific to *.smi format.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET "${BEAKER_ROOT_URL}smiles2ctab/"$(cat rules.smi | base64 -w 0 | tr "+/" "-_")"?computeCoords=0"
    curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2ctabView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2ctabView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smarts2ctabView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', True))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    return _smarts2ctab(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smarts2ctab/<smarts>', method=['OPTIONS', 'GET'], name="smarts2ctab")
def smarts2ctab(smarts):
    """
Converts SMARTS to CTAB. This method accepts urlsafe_base64 encoded string containing single or multiple SMARTS.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smarts2ctab/$(echo '[CH2X4][CX3](=[OX1])[NX3H2]' | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smarts2ctab/$(cat amino.sma | base64 -w 0 | tr "+/" "-_")
    curl -X GET "${BEAKER_ROOT_URL}smarts2ctab/"$(cat amino.sma | base64 -w 0 | tr "+/" "-_")"?computeCoords=0"
    """

    data = base64.urlsafe_b64decode(smarts)
    return smarts2ctabView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smarts2ctab', method=['OPTIONS', 'POST'], name="smarts2ctab")
def smarts2ctab():
    """
Converts SMARTS to CTAB. This method accepts single or multiple SMARTS.
cURL examples:

    curl -X POST -F "file=@amino.sma" ${BEAKER_ROOT_URL}smarts2ctab
    curl -X POST -F "file=@amino.sma" -F "computeCoords=0"  ${BEAKER_ROOT_URL}smarts2ctab
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smarts2ctabView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2inchiView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2inchi(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2inchi/<smiles>', method=['OPTIONS', 'GET'], name="smiles2inchi")
def smiles2inchi(smiles):
    """
Converts SMILES to InChi. This method accepts urlsafe_base64 encoded string containing single or multiple SMILES
optionally containing header line, specific to *.smi format.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2inchiView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2inchiView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2inchiKeyView(data, params):
    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _smiles2inchiKey(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2inchiKey/<smiles>', method=['OPTIONS', 'GET'], name="smiles2inchiKey")
def smiles2inchiKey(smiles):
    """
Converts SMILES to InChi Key. This method accepts urlsafe_base64 encoded string containing single or multiple SMILES
optionally containing header line, specific to *.smi format.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2inchiKeyView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
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
    kwargs['isomericSmiles'] = _parseFlag(params.get('isomericSmiles', False))
    kwargs['kekuleSmiles'] = _parseFlag(params.get('kekuleSmiles', False))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    return _canonicalize_smiles(data, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------


@app.route('/canonicalizeSmiles/<smiles>', method=['OPTIONS', 'GET'], name="canonicalizeSmiles")
def canonicalizeSmiles(smiles):
    """
Converts SMILES to canonical SMILES. This method accepts urlsafe_base64 encoded string containing single or multiple
SMILES optionally containing header line, specific to *.smi format.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}canonicalizeSmiles/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}canonicalizeSmiles/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")"?out_delimiter=|&nameHeader=foo"
    curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=0"
    curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=1&sanitize=1"
    curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat isomeric.smi | base64 -w 0 | tr "+/" "-_")"?isomericSmiles=1"
    """

    data = base64.urlsafe_b64decode(smiles)
    return canonicalizeSmilesView(data, request.params)

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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return canonicalizeSmilesView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2ctab/<inchi>', method=['OPTIONS', 'GET'], name="inchi2ctab")
def inchi2ctab(inchi):
    """
Converts InChi to CTAB. This method accepts urlsafe_base64 encoded string containing one or multiple InChis.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}inchi2ctab/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")
    """

    inchis = base64.urlsafe_b64decode(inchi)
    return _inchi2ctab(inchis)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2ctab', method=['OPTIONS', 'POST'], name="inchi2ctab")
def inchi2ctab():
    """
Converts InChi to CTAB. This method accepts one or multiple InChis.
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2ctab
    curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2ctab
    """

    inchis = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _inchi2ctab(inchis)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2inchiView(data, params):
    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _ctab2inchi(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchi/<ctab>', method=['OPTIONS', 'GET'], name="ctab2inchi")
def ctab2inchi(ctab):
    """
Converts CTAB to InChis. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2inchi/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2inchiView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchi', method=['OPTIONS', 'POST'], name="ctab2inchi")
def ctab2inchi():
    """
Converts CTAB to InChis. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchi
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchi
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2inchiView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def ctab2inchiKeyView(data, params):
    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    return _ctab2inchiKey(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchiKey/<ctab>', method=['OPTIONS', 'GET'], name="ctab2inchiKey")
def ctab2inchiKey(ctab):
    """
Converts CTAB to InChi Keys. CTAB is urlsafe_base64 encoded string containing single molfile or concatenation
of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2inchiKey/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2inchiKeyView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2inchiKey', method=['OPTIONS', 'POST'], name="ctab2inchiKey")
def ctab2inchiKey():
    """
Converts CTAB to InChi Keys. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchiKey
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchiKey
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2inchiKeyView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2inchiKey/<inchi>', method=['OPTIONS', 'GET'], name="inchi2inchiKey")
def inchi2inchiKey(inchi):
    """
Converts InChis to InChiKeys. This method accepts urlsafe_base64 encoded string containing one or multiple InChis.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}inchi2inchiKey/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")
    """

    inchis = base64.urlsafe_b64decode(inchi)
    return _inchi2inchiKey(inchis)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/inchi2inchiKey', method=['OPTIONS', 'POST'], name="inchi2inchiKey")
def inchi2inchiKey():
    """
Converts InChis to InChiKeys. This method accepts one or multiple InChis.
cURL examples:

    curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2inchiKey
    curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2inchiKey
    """

    inchis = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return _inchi2inchiKey(inchis)

# ----------------------------------------------------------------------------------------------------------------------

