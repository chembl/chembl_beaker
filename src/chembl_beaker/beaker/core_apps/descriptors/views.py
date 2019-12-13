__author__ = 'mnowotka'

from beaker import app
from beaker.core_apps.descriptors.impl import _getDescriptors, _getChemblDescriptors
from beaker.utils.io import _parseFlag
from bottle import request
import base64
import json

#-----------------------------------------------------------------------------------------------------------------------

def chemblDescriptorsView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', True))
    return json.dumps(_getChemblDescriptors(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/chemblDescriptors/<ctab>', method=['OPTIONS', 'GET'], name="chemblDescriptors")
def chemblDescriptors(ctab):
    """
Returns descriptors available in ChEMBL. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}descriptors/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return chemblDescriptorsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/chemblDescriptors', method=['OPTIONS', 'POST'], name="chemblDescriptors")
def chemblDescriptors():
    """
Returns descriptors available in ChEMBL. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}descriptors
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}descriptors
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return chemblDescriptorsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def descriptorsView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', True))
    kwargs['ds'] = params.get('descrs')
    return json.dumps(_getDescriptors(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/descriptors/<ctab>', method=['OPTIONS', 'GET'], name="descriptors")
def descriptors(ctab):
    """
Returns descriptors of a compound. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:
    curl -X GET ${BEAKER_ROOT_URL}descriptors/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return descriptorsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/descriptors', method=['OPTIONS', 'POST'], name="descriptors")
def descriptors():
    """
Returns descriptors of a compound. CTAB is either single molfile or SDF file.
cURL examples:
    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}descriptors
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}descriptors
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return descriptorsView(data, request.params)