__author__ = 'mnowotka'

from chembl_beaker.beaker import app
from bottle import request
from chembl_beaker.beaker.core_apps.ringInfo.impl import _isAtomInRing, _atomRings
from chembl_beaker.beaker.core_apps.ringInfo.impl import _bondRings, _isBondInRing
from chembl_beaker.beaker.core_apps.ringInfo.impl import _numAtomRings, _numBondRings, _numRings
from chembl_beaker.beaker.utils.io import _parseFlag
import base64
import json

#-----------------------------------------------------------------------------------------------------------------------

def atomRingsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_atomRings(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomRings/<ctab>', method=['OPTIONS', 'GET'], name="atomRings")
def atomRings(ctab):
    """
Returns atom indexes group by SSSR rings. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}atomRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return atomRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomRings', method=['OPTIONS', 'POST'], name="atomRings")
def atomRings():
    """
Returns atom indexes group by SSSR rings. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}atomRings
    curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}atomRings
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return atomRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def bondRingsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_bondRings(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/bondRings/<ctab>', method=['OPTIONS', 'GET'], name="bondRings")
def bondRings(ctab):
    """
Returns bond indexes group by SSSR rings. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}bondRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return bondRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/bondRings', method=['OPTIONS', 'POST'], name="bondRings")
def bondRings():
    """
Returns bond indexes group by SSSR rings. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}bondRings
    curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}bondRings
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return bondRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def atomIsInRingView(data, params):

    kwargs = dict()
    kwargs['index'] = int(params.get('index'))
    kwargs['size'] = int(params.get('size'))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_isAtomInRing(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomIsInRing/<ctab>/<index>/<size>', method=['OPTIONS', 'GET'], name="atomIsInRing")
def atomIsInRing(ctab, index, size):
    """
A predicate, checking if atom is in ring. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}atomIsInRing/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")/1/6

    """

    data = base64.urlsafe_b64decode(ctab)

    params = request.params
    params['index'] = index
    params['size'] = size

    return atomIsInRingView(data, params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomIsInRing', method=['OPTIONS', 'POST'], name="atomIsInRing")
def atomIsInRing():
    """
A predicate, checking if atom is in ring. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST -F "file=@addHs.mol" -F "index=1" -F "size=6" ${BEAKER_ROOT_URL}atomIsInRing
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return atomIsInRingView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def bondIsInRingView(data, params):

    kwargs = dict()
    kwargs['index'] = int(params.get('index'))
    kwargs['size'] = int(params.get('size'))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_isBondInRing(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/bondIsInRing/<ctab>/<index>/<size>', method=['OPTIONS', 'GET'], name="bondIsInRing")
def bondIsInRing(ctab, index, size):
    """
A predicate, checking if bond is in ring. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}bondIsInRing/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")/1/3

    """

    data = base64.urlsafe_b64decode(ctab)

    params = request.params
    params['index'] = index
    params['size'] = size

    return bondIsInRingView(data, params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/bondIsInRing', method=['OPTIONS', 'POST'], name="bondIsInRing")
def bondIsInRing():
    """
A predicate, checking if bond is in ring. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST -F "file=@rings.mol" -F "index=1" -F "size=3" ${BEAKER_ROOT_URL}bondIsInRing
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return bondIsInRingView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def numAtomRingsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_numAtomRings(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numAtomRings/<ctab>', method=['OPTIONS', 'GET'], name="numAtomRings")
def numAtomRings(ctab):
    """
Returns ring membership counts for every atom. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}numAtomRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")

    """
    data = base64.urlsafe_b64decode(ctab)
    return numAtomRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numAtomRings', method=['OPTIONS', 'POST'], name="numAtomRings")
def numAtomRings():
    """
Returns ring membership counts for every atom. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numAtomRings
    curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numAtomRings
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return numAtomRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def numBondRingsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_numBondRings(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numBondRings/<ctab>', method=['OPTIONS', 'GET'], name="numBondRings")
def numBondRings(ctab):
    """
Returns ring membership counts for every bond. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}numBondRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")

    """
    data = base64.urlsafe_b64decode(ctab)
    return numBondRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numBondRings', method=['OPTIONS', 'POST'], name="numBondRings")
def numBondRings():
    """
Returns ring membership counts for every bond. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numBondRings
    curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numBondRings
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return numBondRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def numRingsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_numRings(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numRings/<ctab>', method=['OPTIONS', 'GET'], name="numRings")
def numRings(ctab):
    """
Returns number of rings in molecule (NB SSSR set only). CTAB is either single molfile or SDF file.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}numRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
    """
    data = base64.urlsafe_b64decode(ctab)
    return numRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/numRings', method=['OPTIONS', 'POST'], name="numRings")
def numRings():
    """
Returns number of rings in molecule (NB SSSR set only). CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numRings
    curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numRings
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return numRingsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------