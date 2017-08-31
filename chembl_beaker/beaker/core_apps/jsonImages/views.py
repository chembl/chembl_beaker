__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.jsonImages.impl import _ctab2json, _smiles2json
from bottle import request, response
from chembl_beaker.beaker.utils.io import _parseFlag
import base64

# ----------------------------------------------------------------------------------------------------------------------


def ctab2jsonView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 200))
    kwargs['legend'] = params.get('legend','')
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    response.content_type = 'application/json'

    return _ctab2json(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2json/<ctab>', method=['OPTIONS', 'GET'], name="ctab2json")
def ctab2json(ctab):
    """
Converts CTAB to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. CTAB is urlsafe_base64 encoded string containing
single molfile or concatenation of multiple molfiles. Size is the optional size of image in pixels (default value
is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}ctab2json/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}ctab2json/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0
    curl -X GET ${BEAKER_ROOT_URL}ctab2json/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1
    curl -X GET ${BEAKER_ROOT_URL}ctab2json/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin
    curl -X GET ${BEAKER_ROOT_URL}ctab2json/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400
    """

    data = base64.urlsafe_b64decode(ctab)
    return ctab2jsonView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/ctab2json', method=['OPTIONS', 'POST'], name="ctab2json")
def ctab2json():
    """
Converts CTAB to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. CTAB is either single molfile or SDF file. Size
is the optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2json
    curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2json
    curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2json
    curl -X POST -F "file=@aspirin.mol" -F "legend=aspirin" ${BEAKER_ROOT_URL}ctab2json
    curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2json
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return ctab2jsonView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def smiles2jsonView(data, params):

    kwargs = dict()
    kwargs['size'] = int(params.get('size', 200))
    kwargs['legend'] = params.get('legend','')
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'application/json'
    return _smiles2json(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2json/<smiles>', method=['OPTIONS', 'GET'], name="smiles2json")
def smiles2json(smiles):
    """
Converts SMILES to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. This method accepts urlsafe_base64 encoded string
containing single or multiple SMILES optionally containing header line, specific to *.smi format. Size is the
optional size of image in pixels (default value is 200 px). Legend is optional label in the bottom of image.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin
    curl -X GET ${BEAKER_ROOT_URL}smiles2json/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400
    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2jsonView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2json', method=['OPTIONS', 'POST'], name="smiles2json")
def smiles2json():
    """
Converts SMILES to JSON. Resulting JSON contains data about 2D graphic representation of compounds and can be used
directly by Raphael.js library to draw compound image in browser. This method accepts single or multiple SMILES or
*.smi file. Size is the optional size of image in pixels (default valueis 200 px). Legend is optional label in
the bottom of image.
cURL examples:

    curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2json
    curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2json
    curl -X POST -F "file=@aspirin_no_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2json
    curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2json
    curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2json
    curl -X POST -F "file=@aspirin_with_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2json
    curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2json
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2jsonView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

