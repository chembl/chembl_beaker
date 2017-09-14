__author__ = 'mnowotka'

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.core_apps.calculations.impl import _kekulize, _sanitize, _addHs
from chembl_beaker.beaker.core_apps.calculations.impl import _removeHs, _sssr, _symmsssr
from chembl_beaker.beaker.core_apps.calculations.impl import _align
from chembl_beaker.beaker.utils.io import _parseFlag
import base64
import json

# ----------------------------------------------------------------------------------------------------------------------


def alignView(data, template, params):
    force = _parseFlag(params.get('force', False))
    ret = _align(data, template, force)
    if ret:
        return ret
    response.status = 404
    response.body = 'No template found in any of structures given.'

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/align/<template>/<ctab>', method=['OPTIONS', 'GET'], name="align")
def align(template, ctab):
    """
Generate a depiction for a molecules given as CTAB where a piece of the molecule is constrained to have the same 
coordinates as a reference molfile. TEMPLATE is urlsafe_base64 encoded string containing a single molfile.
CTAB is urlsafe_base64 encoded string containing concatenation of multiple molfiles to be aligned.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}align/$(cat pattern.mol | base64 -w 0 | tr "+/" "-_")/$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}align/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")/$(cat chembl_1999443.mol | base64 -w 0 | tr "+/" "-_")?force=true
    """

    data = base64.urlsafe_b64decode(ctab)
    t = base64.urlsafe_b64decode(template)
    return alignView(data, t, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/align', method=['OPTIONS', 'POST'], name="align")
def align():
    """
Generate a depiction for a molecules. Expects either two CTABS - one with a template fragment and another one with 
molecules to be aligned or a single CTAB with the first molecule being a template and others so be aligned.
cURL examples:

    curl -X POST -F "template=@pattern.mol" -F "ctab=@mcs.sdf" ${BEAKER_ROOT_URL}align
    curl -X POST -F "template=@aspirin.mol" -F "ctab=@chembl_1999443.mol"  -F "force=true" ${BEAKER_ROOT_URL}align
    curl -X POST -F "ctab=@pattern_and_mcs.sdf" ${BEAKER_ROOT_URL}align
    """

    number_of_files = len(request.files)

    if number_of_files == 1:
        ctab = request.files['ctab'].file.read()
        return alignView(ctab, None, request.params)

    elif number_of_files == 2:
        template = request.files['template'].file.read()
        ctab = request.files['ctab'].file.read()
        return alignView(ctab, template, request.params)

    response.status = 400
    response.body = 'This method expects one ot two ctabs.'
    print 'this should be error 400'
    return response

# ----------------------------------------------------------------------------------------------------------------------


def kekulizeView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return _kekulize(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/kekulize/<ctab>', method=['OPTIONS', 'GET'], name="kekulize")
def kekulize(ctab):
    """
Performs kekulisation on input compounds. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}kekulize/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return kekulizeView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/kekulize', method=['OPTIONS', 'POST'], name="kekulize")
def kekulize():
    """
Performs kekulisation on input compounds. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}kekulize
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}kekulize
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return kekulizeView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def sanitizeView(data, params):

    sanitizeOps = params.get('sanitizeOps')
    if sanitizeOps is None:
        return _sanitize(data)
    return _sanitize(data, int(sanitizeOps))

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sanitize/<ctab>', method=['OPTIONS', 'GET'], name="sanitize")
def sanitize(ctab):
    """
Kekulize, check valencies, set aromaticity, conjugation and hybridization. CTAB is urlsafe_base64 encoded string
containing single molfile or concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}sanitize/$(cat aromatic.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return sanitizeView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sanitize', method=['OPTIONS', 'POST'], name="sanitize")
def sanitize():
    """
Kekulize, check valencies, set aromaticity, conjugation and hybridization. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aromatic.mol ${BEAKER_ROOT_URL}sanitize
    curl -X POST -F "file=@aromatic.mol" ${BEAKER_ROOT_URL}sanitize
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return sanitizeView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def symmSSSRView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_symmsssr(data, **kwargs))

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/symmSSSR/<ctab>', method=['OPTIONS', 'GET'], name="symmSSSR")
def symmSSSR(ctab):
    """
Get a symmetrized SSSR for a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}symmSSSR/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return symmSSSRView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/symmSSSR', method=['OPTIONS', 'POST'], name="symmSSSR")
def symmSSSR():
    """
Get a symmetrized SSSR for a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}symmSSSR
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}symmSSSR
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return symmSSSRView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def sssrView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))

    return json.dumps(_sssr(data, **kwargs))

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sssr/<ctab>', method=['OPTIONS', 'GET'], name="sssr")
def sssr(ctab):
    """
Get the smallest set of simple rings for a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}sssr/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
    """

    data = base64.urlsafe_b64decode(ctab)
    return sssrView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sssr', method=['OPTIONS', 'POST'], name="sssr")
def sssr():
    """
Get the smallest set of simple rings for a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}sssr
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}sssr
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return sssrView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def addHsView(data, params):

    kwargs = dict()
    kwargs['explicitOnly'] = _parseFlag(params.get('explicitOnly', False))
    kwargs['addCoords'] = _parseFlag(params.get('addCoords', False))

    return _addHs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/addHs/<ctab>', method=['OPTIONS', 'GET'], name="addHs")
def addHs(ctab):
    """
Adds hydrogens to the graph of a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}addHs/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}addHs/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")?addCoords=1
    """

    data = base64.urlsafe_b64decode(ctab)
    return addHsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/addHs', method=['OPTIONS', 'POST'], name="addHs")
def addHs():
    """
Adds hydrogens to the graph of a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @addHs.mol ${BEAKER_ROOT_URL}addHs
    curl -X POST -F "file=@addHs.mol" -F "addCoords=1" ${BEAKER_ROOT_URL}addHs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return addHsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def removeHsView(data, params):

    kwargs = dict()
    kwargs['implicitOnly'] = _parseFlag(params.get('implicitOnly', False))
    return _removeHs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/removeHs/<ctab>', method=['OPTIONS', 'GET'], name="removeHs")
def removeHs(ctab):
    """
Removes any hydrogens from the graph of a molecule. CTAB is urlsafe_base64 encoded string containing single molfile or
concatenation of multiple molfiles.
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}removeHs/$(cat removeHs.mol | base64 -w 0 | tr "+/" "-_")
    curl -X GET ${BEAKER_ROOT_URL}removeHs/$(cat removeHs.mol | base64 -w 0 | tr "+/" "-_")?implicitOnly=1
    """

    data = base64.urlsafe_b64decode(ctab)
    return removeHsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/removeHs', method=['OPTIONS', 'POST'], name="removeHs")
def removeHs():
    """
Removes any hydrogens from the graph of a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @removeHs.mol ${BEAKER_ROOT_URL}removeHs
    curl -X POST -F "file=@removeHs.mol" -F "implicitOnly=1" ${BEAKER_ROOT_URL}removeHs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return removeHsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------
