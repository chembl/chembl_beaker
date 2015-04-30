__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.fingerprints.impl import _sdf2fps
from bottle import request
import base64
import json

#-----------------------------------------------------------------------------------------------------------------------

def sdf2fpsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['type'] = params.get('type', 'morgan')
    kwargs['radius'] = int(params.get('radius', 2))
    kwargs['n_bits'] = int(params.get('n_bits', 2048))

    return _sdf2fps(data, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2fps/<ctab>', method=['OPTIONS', 'GET'], name="sdf2fps")
def sdf2fps(ctab):
    """
Computes fingerprints for given compounds and writes them as fps format. Type is a optional type of fingerprints,
can be 'morgan', 'pair' or 'maccs' (default to 'morgan'). CTAB is urlsafe_base64 encoded string containing single
molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}sdf2fps/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?n_bits=1024&radius=5"
    curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?type=pair"
    curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?type=maccs"

    """

    data = base64.urlsafe_b64decode(ctab)
    return sdf2fpsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2fps', method=['OPTIONS', 'POST'], name="sdf2fps")
def sdf2fps():
    """
Computes fingerprints for given compounds and writes them as fps format. Type is a optional type of fingerprints,
can be 'morgan', 'pair' or 'maccs' (default to 'morgan'). CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}sdf2fps
    curl -X POST -F "file=@aspirin.mol" -F "n_bits=1024" -F "radius=5" ${BEAKER_ROOT_URL}sdf2fps
    curl -X POST -F "file=@aspirin.mol" -F "type=maccs" ${BEAKER_ROOT_URL}sdf2fps

    """
    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return sdf2fpsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------
