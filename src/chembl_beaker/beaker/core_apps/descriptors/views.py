__author__ = 'mnowotka'

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.descriptors.impl import _getNumAtoms,_getNumBonds, _getLogP, _getTPSA, _getMolWt, _getDescriptors
from bottle import request
import base64
import json

#-----------------------------------------------------------------------------------------------------------------------

def getNumAtomsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))

    return json.dumps(_getNumAtoms(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/getNumAtoms/<ctab>', method=['OPTIONS', 'GET'], name="getNumAtoms")
def getNumAtoms(ctab):
    """
Counts number of atoms of given compounds. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}getNumAtoms/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return getNumAtomsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/getNumAtoms', method=['OPTIONS', 'POST'], name="getNumAtoms")
def getNumAtoms():
    """
Counts number of atoms of given compounds. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}getNumAtoms
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}getNumAtoms
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return getNumAtomsView(data, request.params)
#-----------------------------------------------------------------------------------------------------------------------

def getNumBondsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))

    return json.dumps(_getNumBonds(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/getNumBonds/<ctab>', method=['OPTIONS', 'GET'], name="getNumBonds")
def getNumBonds(ctab):
    """
Counts number of bonds of given compounds. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}getNumBonds/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")

    """

    data = base64.urlsafe_b64decode(ctab)
    return getNumBondsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/getNumBonds', method=['OPTIONS', 'POST'], name="getNumBonds")
def getNumBonds():
    """
Counts number of bonds of given compounds. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}getNumBonds
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}getNumBonds
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return getNumBondsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def logPView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))

    return json.dumps(_getLogP(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/logP/<ctab>', method=['OPTIONS', 'GET'], name="logP")
def logP(ctab):
    """
Returns the logP value for a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}logP/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return logPView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/logP', method=['OPTIONS', 'POST'], name="logP")
def logP():
    """
Returns the logP value for a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}logP
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}logP
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return logPView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def tpsaView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))
    return json.dumps(_getTPSA(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/tpsa/<ctab>', method=['OPTIONS', 'GET'], name="tpsa")
def tpsa(ctab):
    """
Returns the TPSA value for a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}tpsa/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return tpsaView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/tpsa', method=['OPTIONS', 'POST'], name="tpsa")
def tpsa():
    """
Returns the TPSA value for a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}tpsa
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}tpsa
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return tpsaView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def molWtView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))
    return json.dumps(_getMolWt(data, **kwargs))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/molWt/<ctab>', method=['OPTIONS', 'GET'], name="molWt")
def molWt(ctab):
    """
Returns molecular weight of a compound. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}molWt/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return molWtView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/molWt', method=['OPTIONS', 'POST'], name="molWt")
def molWt():
    """
Returns molecular weight of a compound. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}molWt
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}molWt
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return molWtView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------

def descriptorsView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = bool(params.get('sanitize', True))
    kwargs['removeHs'] = bool(params.get('removeHs', True))
    kwargs['strictParsing'] = bool(params.get('strictParsing', True))
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

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return descriptorsView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------
